      subroutine itrdrv (y,         ac,         banma,         
     &                   uold,      x,         
     &                   iBC,       BC,         
     &                   iper,      ilwork,     shp,       
     &                   shgl,      shpb,       shglb,
     &                   ifath,     velbar,     nsons ) 
c
c!----------------------------------------------------------------------
c
c! This iterative driver is the semi-discrete, predictor multi-corrector 
c! algorithm. It contains the Hulbert Generalized Alpha method which
c! is 2nd order accurate for Rho_inf from 0 to 1.  The method can be
c! made  first-order accurate by setting Rho_inf=-1. It uses CGP and
c! GMRES iterative solvers.
c
c! working arrays:
c!  y      (nshg,ndof)           : Y variables
c!  x      (nshg,nsd)            : node coordinates
c!  banma  (nshg,1)              : marker field 
c!  iBC    (nshg)                : BC codes
c!  BC     (nshg,ndofBC)         : BC constraint parameters
c!  iper   (nshg)                : periodicity table
c
c
c! Zdenek Johan,  Winter 1991.  (Fortran 90)
c! Alberto Figueroa, Winter 2004.  CMM-FSI
c! Irene Vignon, Fall 2004. Impedance BC
c! Jun Fang, Summer 2014. Bubble Tracking
c!----------------------------------------------------------------------
c
      USE, INTRINSIC :: ISO_C_BINDING !for calling C++ routines 
      use pvsQbi        !gives us splag (the spmass at the end of this run 
      use specialBC     !gives us itvn
      use timedata      !allows collection of time series
      use convolImpFlow !uses flow history and impedance for convolution
      use spat_var_eps  !use spatial varying eps_ls
      use cf_arrays     !access to control forces arrays
      use bub_track     !access to bubble information array
      use turbsa        !access to d2wal info
!      use readarrays    !reads in uold and acold
      use redist_freeze !access to BC arrays if freezing value of primary LS vertices ! ???
      use fncorpmod ! Arsen
      
        include "common.h"
        include "mpif.h"
        include "auxmpi.h"
!        INCLUDE "memLS.h"
        include "svLS.h"   ! Arsen, include svLS lib
c
        real*8 NewQImp(0:MAXSURF) !temporary unknown for the flow
                        !rate that needs to be added to the flow history
        
        real*8    y(nshg,ndof),              ac(nshg,ndof),           
     &	          yold(nshg,ndof),           acold(nshg,ndof),
     &            u(nshg,nsd),               uold(nshg,nsd),
     &            x(numnp,nsd),              solinc(nshg,ndof),
     &            BC(nshg,ndofBC),           tf(nshg,ndof)

c
        real*8    res(nshg,ndof)
c     
        real*8    shp(MAXTOP,maxsh,MAXQPT),  
     &            shgl(MAXTOP,nsd,maxsh,MAXQPT), 
     &            shpb(MAXTOP,maxsh,MAXQPT),
     &            shglb(MAXTOP,nsd,maxsh,MAXQPT) 
c
        integer   rowp(nshg,nnz),         colm(nshg+1),
     &            iBC(nshg),
     &            ilwork(nlwork),
     &            iper(nshg),            ifuncs(6)

        integer stopjob, iphase, recnum, numvar_rec
        character*10 cname2
        character*5  cname
        integer i_redist_counter
        real*8 redist_toler_previous
        logical iloop
c
c  stuff for dynamic model s.w.avg and wall model
c
        dimension ifath(numnp),    velbar(nfath,ndof),  nsons(nfath)

        dimension wallubar(2),walltot(2)
c     
c.... For Farzin''s Library
c
        integer eqnType, prjFlag, presPrjFlag, verbose
c
        real*8, allocatable, dimension(:,:) :: aperm,  atemp, atempS
        real*8, allocatable, dimension(:,:,:) :: apermS

        real*8, allocatable, dimension(:,:) :: lhsP, lhsK, lhsS
        real*8   almit, alfit, gamit
c
        character*1024    servername
        character*20    fname1,fmt1
        character*20    fname2,fmt2,fnamer2
        character*60    fnamepold, fvarts, fvartsb
        character*4     fieldybar
        integer         iarray(50) ! integers for headers

        real*8 rerr(nshg,numerr),ybar(nshg,5) , uhess(nshg,27),
     &         gradu( nshg, 9 )
        integer, allocatable, dimension(:) :: ivarts
        integer, allocatable, dimension(:) :: ivartsg
        real*8, allocatable, dimension(:) :: vartssoln
        real*8, allocatable, dimension(:) :: vartssolng
        real*8 elem_size(numel), elemb_size(numelb)
        real*8 elem_size_min, elem_size_mintmp
      
        real*8 xi2
c
        real*8 gradphi(nshg,3), gradphimag(nshg), maxgradphi
	integer igradphi
c
c Redistancing option of fixing phi of primary vertices
c
        real*8  primvertval(nshg,2)
        integer primvert(nshg)
        integer i_primvert,numpv,numpvset
        integer  iredist_flag, numrun

        REAL*8,  allocatable :: BCredist(:)
        integer, allocatable :: iBCredist(:)

        real*8 CFLls(nshg)
        integer ifairmod, iwaittime

!----------------------------------------------------------------------
!       The following part covers variables, arrays used in
!               1. Bubble tracking capability 
!               2. Bubble control capability
!               3. Coalescence control
!               4. Simulations steering
!               5. Mass conservation control
!               6. Bubble interface CFL number adjustment
!       These PHASTA based techniques are developed by Dr. By, ac, iBC, BC,
!       iper, ilwork,banmaolotnov group. 
!                                       Jun Fang, 2016-02-18
!----------------------------------------------------------------------
        dimension banma(nshg,1)

        real*8, allocatable ::  breakupSeederTMP(:,:)
        real*8, allocatable ::  bub2wTMP(:,:)
        real*8, allocatable ::  D_eqTMP(:)
        integer   ibreakupFlag, nBubbleLTS

        real*8    qss(nshg,ndof),       acss(nshg,ndof),
     &            uss(nshg,nsd)

        real*8    elemvol_global(nelblk,numel)   !Element volume array

        real*8 vf,              vf_tmp,         vf_initial 
        real*8 sumxcf,          sumycf,         sumzcf
!        real*8 sumxcf_o,        sumycf_o,       sumzcf_o
!        real*8 sumycf_old
        integer ierror,         numoldyhistind
!        integer ixcf,           iycf,           izcf
        integer ioldyhistst,    ioldyhisten,    indyhist
!        logical ycf_old_log

        real*8 avgxcoordf(coalest),     avgycoordf(coalest),
     &         avgzcoordf(coalest),     avgxcoordold2(coalest),
     &         avgycoordold2(coalest),  avgzcoordold2(coalest),
     &         app_time(coalest,2)
        real*8 itrtimestp
!       For simulation steering
        integer inumstart
        character*10 inpfname
!       To compute how many times level set meets convergence criteria
        integer N_Lsc  
!        real*8  oldcf_dt,       oldcf_dtlset
        character*50:: lstepc
! ---------------------------------------------------------------------------
! Arsen block
!  Setting up the svLS
      integer svLS_nFaces_sc, svLS_nFaces, gnNo, nNo, faIn, facenNo
      integer svLS_nFacesT, gnNoT, nNoT, faInT, facenNoT
!      integer, allocatable :: ltg(:)
      integer, allocatable ::  gNodes(:), gNodesT(:)
      real*8, allocatable :: sV(:,:), svT(:,:)
      
      character*128 fileName
      TYPE(svLS_commuType) communicator
      TYPE(svLS_lhsType) svLS_lhs, svLS_lhs_sc, svLS_lhsT
      TYPE(svLS_lsType) svLS_ls, svLS_sc, svLST 
      character*10 cname2nd
      real*8 sumtime
      logical ycf_old_log
      integer ixcf,iycf,izcf
      real*8 sumxcf_o(i_num_bubbles), sumycf_o(i_num_bubbles)
      real*8 oldcf_dt,oldcf_dtlset,sumycf_old
      integer svLSFlag

      real*8 :: sumzcf_o(i_num_bubbles)

!----------------------------------------------------------------------
!!     Setting up memLS
!!      INTEGER memLS_nFaces_sc, memLS_nFaces, gnNo, nNo, faIn, facenNo
!!      INTEGER memLS_nFacesT, gnNoT, nNoT, faInT, facenNoT
!!      INTEGER, ALLOCATABLE :: ltg(:), gNodes(:), gNodesT(:)
!!      REAL*8, ALLOCATABLE :: sV(:,:), sVT(:,:)

!!      CHARACTER*128 fileName
!!      TYPE(memLS_commuType) communicator
!!      TYPE(memLS_lhsType) memLS_lhs, memLS_lhs_sc, memLS_lhsT
!!      TYPE(memLS_lsType) memLS_ls, memLS_sc, memLST
!!      character*10  cname2nd
!!      REAL*8 sumtime

!----------------------------------------------------------------------
!       Initialize the current void fraction and constant used in 
!       interface adjustment 
!                                               Jun, March, 2014
!----------------------------------------------------------------------
        C_int_adjust    = 0.0d0
        vf_now          = 0.0d0
!----------------------------------------------------------------------
!       Open files for bubble analysis
!----------------------------------------------------------------------
         write (*,*) iLSet, "this is ilset !!!!!"
         flush(0)
        if(iBT.eq.1 .and. iLSet .eq. 2) then
           if(myrank .eq. master) write(*,*) 'iBT is', iBT
!       get the total number of bubbles in the entire domain (i_num_bubbles)
           call CountIniBub(banma)
           allocate( bub_cent(i_num_bubbles,5) )
		   bub_cent = 0
           ibreakupFlag         = 0
           if(iBK.eq.1) then
              allocate( breakupSeeder(i_num_bubbles,6) )
              breakupSeeder     = zero
              do i = 1, i_num_bubbles
                 breakupSeeder(i,1) = real(i)
                 breakupSeeder(i,6) = real(i)
              enddo
              allocate( D_eq(i_num_bubbles) )
              allocate( bub2w(i_num_bubbles,3) )
              D_eq  = zero
              bub2w = zero
           endif
           nBubbleLTS           = i_num_bubbles
           ts_hold              = lstep
           GhostRatio           = 0.2
           NbrhdRatio           = 0.2

           if(myrank.eq.master) then
              call OpenBubFiles(nBubbleLTS, i_num_bubbles,
     &                          C_int_adjust)
              write(*,*) 'OpenBubFiles is done!'
           endif
        endif
        if(icoalCtrl.eq.1) then
           allocate(coalCenter(1,3))
		   coalCenter = 0
           ncoalEvent = 0
        endif
!----------------------------------------------------------------------
!--------------------------------------------------------------------------------------
! Arsen block -- block added to initialise the svLS solver
      svLSFlag = 1
	  if(myrank.eq.master) write(*,*) "svLSFlag is set to ", svLSFlag
      IF (svLSFlag .EQ. 1) THEN
        !if(myrank.eq.master) write(*,*) "calling svLS_LS_CREATE"
         
        !call svLS_LS_FREE(svLS_ls)  ! MB, test
        call svLS_LS_CREATE(svLS_ls, LS_TYPE_GMRES, dimKry=Kspace,
     2   relTol=epstol(8), relTolIn=(/epstol(1), epstol(7)/),
     3   maxItr=nPrjs, maxItrIn=(/nGMRES, maxIters/))

         !if(myrank.eq.master) write(*,*) "called svLS_LS_CREATE"

         nsolt = mod(impl(1), 2)
         nsclrsol=nsolt+nsclr

         if (nsclrsol.gt.0) then
            !call svLS_LS_FREE(svLS_ls) !MB, test
            !if(myrank.eq.master) write(*,*) "calling svLS_LS_CREATE"
            call svLS_LS_CREATE(svLS_sc, LS_TYPE_GMRES, dimKry=Kspace,
     2      relTol=epstol(8), relTolIn=(/epstol(1), epstol(7)/),
     3      maxItr=nPrjs, maxItrIn=(/nGMRES, maxIters/))
            !if(myrank.eq.master) write(*,*) "called svLS_LS_CREATE"
         end if

         !if (myrank.eq.master) write(*,*) 'LS_type is ', LS_type
         !call svLS_COMMU_CREATE(communicator, MPI_COMM_WORLD)  ! MB, added to prevent MPI_ALLREDUCE problem

! Assuming the protocal to read the ltg files and set gnNo, nNO and ltg is not required
! we can simply comment all of this section out. This is what we are currently trying
!          idirstep = 512 !tmp test
!          idirtrigger = 10 ! tmp test

!          IF (numpe .GT. 1) THEN
!             if (myrank.eq.master) write(*,*) "starting to write the ltg files"
!             WRITE(fileName,*) myrank+1
!             fileName = "ltg.dat." //ADJUSTL(TRIM(fileName))

!         if (numpe.gt.idirtrigger) then
!             fileName = trim(cname2nd(int(myrank/idirstep)*idirstep))
!      1         //"-set/"//trim(fileName)
!             end if

!             open(1, FILE=fileName)
!             read(1,*) gnNo
!             read(1,*) nNo
!             allocate(ltg(nNo))
!             read(1,*) ltg
!             close(1)

!             !if (myrank.eq.master) write(*,*)"finished writing and reading the ltg files"
!          ELSE
! ! MB, I think this part of the code causes the problem with the bounds in svLS_LHS_CREATE
!             ! memLS syntax
!             gnNo = nshg
!             nNo = nshg
!             allocate(ltg(nNo))
!             do i=1, nNo
!                ltg(i) = i
!             end do
!          END IF
            ! COLORADO SYNTAX
!            nNo = nshg
!            gnNo = nshgt
	  
!!      IF (memLSFlag .EQ. 1) THEN
!!         CALL memLS_LS_CREATE(memLS_ls, LS_TYPE_GMRES, dimKry=Kspace,   ! Modified, Igor, May 2013. Was "LS_TYPE_NS"
!!     2      relTol=epstol(8), relTolIn=(/epstol(1),epstol(7)/),
!!     3      maxItr=nPrjs, maxItrIn=(/nGMRES,maxIters/))

!!      nsolt=mod(impl(1),2)      ! 1 if solving temperature
!!      nsclrsol=nsolt+nsclr      ! total number of scalars solved At
!!         if (nsclrsol.gt.0) then
!!         CALL memLS_LS_CREATE(memLS_sc, LS_TYPE_GMRES, dimKry=Kspace,
!!     2      relTol=epstol(8), relTolIn=(/epstol(1),epstol(7)/),
!!     3      maxItr=nPrjs, maxItrIn=(/nGMRES,maxIters/))
!!        end if

!!         CALL memLS_COMMU_CREATE(communicator, MPI_COMM_WORLD)

!!         IF (numpe .GT. 1) THEN
!!            WRITE(fileName,*) myrank+1
!!            fileName = "ltg.dat."//ADJUSTL(TRIM(fileName))
!!      if (numpe.gt.idirtrigger) then
!!        fileName = trim(cname2nd(int(myrank/idirstep)*idirstep))
!!     1  //"-set/"//trim(fileName)
!!      end if
!!            OPEN(1,FILE=fileName)
!!            READ(1,*) gnNo
!!            READ(1,*) nNo
!!            ALLOCATE(ltg(nNo))
!!			ltg = 0
!!            READ(1,*) ltg
!!            CLOSE(1)
!!         ELSE
!!            gnNo = nshg
!!            nNo = nshg
!!            ALLOCATE(ltg(nNo))
!!			ltg = 0
!!            DO i=1, nNo
!!               ltg(i) = i
!!            END DO
!!         END IF
!!      ELSE
!--------------------------------------------------------------------
!!!        call SolverLicenseServer(servername) ! assume never use LIBLES
      END IF

!
! only master should be verbose
!
        if(numpe.gt.0 .and. myrank.ne.master)iverbose=0  
c
c        call MPI_Barrier(MPI_COMM_WORLD)
c           if (myrank .eq. master) write(*,*) 'Starting itrdrv...' 

        inquire(file='xyzts.dat',exist=exts)
        lskeep = lstep 
        if(exts) then
           tssearch = 1
c           if (myrank .eq. master) write(*,*) 'Open the xyzts.dat...' 
           
           open(unit=626,file='xyzts.dat',status='old')
           read(626,*) ntspts, freq, tolpt, iterat, numvar, numrun 
           if (myrank .eq. master) write(*,*) 'numvar = ', numvar
           call sTD             ! sets data structures
           if (myrank .eq. master) write(*,*) 'sTD is done'
           

           do jj=1,ntspts       ! read coordinate data where solution desired
              read(626,*) ptts(jj,1),ptts(jj,2),ptts(jj,3)
           enddo
           close(626)
c           if (myrank .eq. master) write(*,*) 'ptts is read'

           statptts(:,:) = 0
           parptts(:,:) = zero
           varts(:,:) = zero

           allocate (ivarts(ntspts*numvar))
           allocate (ivartsg(ntspts*numvar))
           allocate (vartssoln(ntspts*numvar))
           allocate (vartssolng(ntspts*numvar))
		   ivarts = 0
		   ivartsg = 0
		   vartssoln = 0
		   vartssolng = 0

c          if (myrank .eq. master) write(*,*) 'allocation is done'


           if (myrank .eq. master) then
                 fvartsb='varts/varts'
                 fvartsb=trim(fvartsb)//trim(cname2(lstep))
                 fvartsb=trim(fvartsb)//'.run'
                 fvartsb=trim(fvartsb)//trim(cname2(numrun))
                 fvartsb=trim(fvartsb)//'.dat'
                 fvartsb=trim(fvartsb)
              numvar_rec = max(15, numvar)  ! Create 19 variables for LS gradient case
                 open(unit=10101, file=fvartsb, status='unknown',
     1   form='formatted', recl=2*8+3+15*numvar_rec, access='direct')   ! Structured output file
           endif ! myrank

        endif   ! exts

c        call MPI_Barrier(MPI_COMM_WORLD)
c          if (myrank .eq. master) write(*,*) 'l.176, itvn: ', itvn

	if (itvn.gt.0.and.myrank.eq.master) then
	  open(unit=1101, file='../bctinflow.dat', status='unknown')
	  close(1101)
	end if 

c
c.... open history and aerodynamic forces files
c
        if (myrank .eq. master) then
           open (unit=ihist,  file=fhist, status='unknown')
           open (unit=iforce, file=fforce, status='unknown')
           open (unit=76, file="fort.76", status='unknown')
c          write(*,*) ' l. 192; iLSet = ', iLSet
           if (iLSet .eq. 2) then
             open (unit=ivhist,  file=fvhist, status='unknown')
           endif
        endif
c
c.... initialize
c     
c         if (myrank .eq. master) write(*,*) 'l.201'
       ifuncs(:)  = 0              ! func. evaluation counter
        istep  = 0
        yold   = y
        acold  = ac

        qss = zero
        acss = zero
        uss = zero

        rerr = zero
        ybar = zero
c        call MPI_Barrier(MPI_COMM_WORLD)
c          if (myrank .eq. master) write(*,*) 'l.210'
c
!.... ---------> initialize Linear Solver parameters (memLS or LESLIB (Farzin) Library) <---------------
c
c.... assign parameter values
c     
        do i = 1, 100
           numeqns(i) = i
        enddo
        nKvecs       = Kspace
        prjFlag      = iprjFlag
        presPrjFlag  = ipresPrjFlag
        verbose      = iverbose
c
c.... determine how many scalar equations we are going to need to solve
c
      nsolt=mod(impl(1),2)      ! 1 if solving temperature
      nsclrsol=nsolt+nsclr      ! total number of scalars solved At
                                ! some point we probably want to create
                                ! a map, considering stepseq(), to find
                                ! what is actually solved and only
                                ! dimension lhs to the appropriate
                                ! size. (see 1.6.1 and earlier for a
                                ! "failed" attempt at this).


      nsolflow=mod(impl(1),100)/10  ! 1 if solving flow
c        call MPI_Barrier(MPI_COMM_WORLD)
c          if (myrank .eq. master) write(*,*) 'l. 221 itrdrv.f'

c
c.... Now, call Farzin''s lesNew routine to initialize
c     memory space
c
c        call MPI_Barrier(MPI_COMM_WORLD)
c           if (myrank .eq. master) write(*,*) 'call genadj...', icnt
      call genadj(colm, rowp, icnt, iper )  ! preprocess the adjacency list   !Igor
c          if (myrank .eq. master) write(*,*) 'l. 285 itrdrv.f'

      nnz_tot=icnt ! this is exactly the number of non-zero blocks on
                   ! this proc
      if ( numpe < 64 ) then ! Arsen mod(0,0) is wrong
         ifairmod = 1
      else
         ifairmod = (numpe/64)
      endif
      iwaittime = (mod(myrank,ifairmod) + 1)

      if (nsolflow.eq.1) then  ! CONDITION FOR FLOW SOLUTION
         lesId   = numeqns(1)
         eqnType = 1
         nDofs   = 4
!#########################################################################
! Arsen in this section we will initialse the linear solver parameters required for svLS
         ! method

         IF (svLSFlag .EQ. 1) THEN
            !write(*,*) 'myrank, gnNo =', myrank, gnNo
            !if (myrank.eq.master) write(*,*) 'calling svLS_BC_CREATE'
            !if (myrank.eq.master) write(*,*) 'gNodes is ', gNodes
            
            !if (myrank.eq.master) write(*,*) 'nNo is ', nNo
            !if (myrank .eq. master) write(*,*) 'calling svLS_LHS_CREATE'

            !-------------------------------------------------------------
            ! ------------------------------------------------------
          svLS_nFaces = 1

          CALL svLS_COMMU_CREATE(communicator, MPI_COMM_WORLD)
          nNo=nshg
          gnNo=nshgt
          
          ! Arsen
          if (.not. allocated(ltg)) then
            allocate(ltg(nNo))
          end if
          do i=1, nNo
             ltg(i) = i
          enddo
      if(myrank.eq.master) then
      write(*,*) 'about to call svLS_LHS_CREATE for field' ! Arsen
      call flush(0)
      endif
          call svLS_LHS_CREATE(svLS_lhs,communicator,gnNo,nNo,nnz_tot,
     &ltg, colm, rowp, svLS_nFaces)
            !if (myrank .eq. master) write(*,*) 'called svLS_LHS_CREATE'       
      if(myrank.eq.master) then
      write(*,*) 'done call svLS_LHS_CREATE for field' ! Arsen
      call flush(0)
      endif
            faIn = 1
            facenNo = 0

            DO i=1, nshg
               IF (IBITS(iBC(i),3,3) .NE. 0) facenNo = facenNo + 1
            END DO

            allocate(gNodes(facenNo), sV(nsd, facenNo))
            sV = 0D0
            j = 0

            DO i=1, nshg
               IF (IBITS(iBC(i),3,3) .NE. 0) THEN
                  j = j + 1
                  gNodes(j) = i
                  IF (.NOT.BTEST(iBC(i),3)) sV(1,j) = 1D0
                  IF (.NOT.BTEST(iBC(i),4)) sV(2,j) = 1D0
                  IF (.NOT.BTEST(iBC(i),5)) sV(3,j) = 1D0
               END IF
            END DO
            

            call svLS_BC_CREATE(svLS_lhs, faIn, facenNo,
     2         nsd, BC_TYPE_Dir, gNodes, sV)
            if (myrank.eq.master) write(*,*) 'called svLS_BC_CREATE'

         ELSE  ! if svLS not available

!--------------------------------------------------------------------
!     Rest of configuration of memLS is added here, where we have LHS
!     pointers
!!         IF (memLSFlag .EQ. 1) THEN
!!               memLS_nFaces = 1

!!           write(*,*) 'myrank, gnNo = ', myrank, gnNo

!!            CALL memLS_LHS_CREATE(memLS_lhs, communicator, gnNo, nNo,
!!     2         nnz_tot, ltg, colm, rowp, memLS_nFaces)

!!            faIn = 1
!!            facenNo = 0
!!            DO i=1, nshg
!!               IF (IBITS(iBC(i),3,3) .NE. 0)  facenNo = facenNo + 1
!!            END DO
!!            ALLOCATE(gNodes(facenNo), sV(nsd,facenNo))
!!			gNodes = 0
!!            sV = 0D0
!!            j = 0
!!            DO i=1, nshg
!!               IF (IBITS(iBC(i),3,3) .NE. 0) THEN
!!                  j = j + 1
!!                  gNodes(j) = i
!!                  IF (.NOT.BTEST(iBC(i),3)) sV(1,j) = 1D0
!!                  IF (.NOT.BTEST(iBC(i),4)) sV(2,j) = 1D0
!!                  IF (.NOT.BTEST(iBC(i),5)) sV(3,j) = 1D0
!!               END IF
!!            END DO
!            if(myrank.eq.master)write(*,*)'memLS_BC_CREATE first time started'
!!            CALL memLS_BC_CREATE(memLS_lhs, faIn, facenNo,
!!     2         nsd, BC_TYPE_Dir, gNodes, sV)
!            if(myrank.eq.master)write(*,*)'memLS_BC_CREATE first time done'

!!         ELSE
!--------------------------------------------------------------------
         call sleep(iwaittime)

!!!         call myfLesNew( lesId,   41994, ! assume never use LIBLES
!!!     &                 eqnType,
!!!     &                 nDofs,          minIters,       maxIters,
!!!     &                 nKvecs,         prjFlag,        nPrjs,
!!!     &                 presPrjFlag,    nPresPrjs,      epstol(1),
!!!     &                 prestol,        verbose,        statsflow,
!!!     &                 nPermDims,      nTmpDims,       servername  )

         allocate (aperm(nshg,nPermDims))
         allocate (atemp(nshg,nTmpDims))
		 aperm = 0
		 atemp = 0

        END IF  ! memLS vs. LESLIB condition

c          if (myrank .eq. master) write(*,*) 'l. 256 itrdrv.f'

         allocate (lhsP(4,nnz_tot))
         allocate (lhsK(9,nnz_tot))
		 lhsP = 0
		 lhsK = 0

!!!         call readLesRestart( lesId,  aperm, nshg, myrank, lstep, ! assume never use LIBLES
!!!     &                        nPermDims )

      else              ! must be the if (nsolflow.eq.1) alternative
         nPermDims = 0
         nTempDims = 0
      endif

!          if (myrank .eq. master) write(*,*) 'l. 262 itrdrv.f'

      if(nsclrsol.gt.0) then  ! Solving for scalars
	     svLSFlag = 1
         IF (svLSFlag .EQ. 1) THEN ! Arsen
! Possible add separate memLS_lhs definition for each scalar ? Important for
! b.c.s

        if (nsolt.gt.0) then
         if (myrank.eq.master)
     &  write(*,*) 'Setting BC for temperature'
          

! From genbc1:         where (btest(iBC,1)) BC(:,2) = BCtmp(:,2) ! temperature

               !memLS_nFacesT = 1
               svLS_nFacesT = 1 ! Arsen
!            CALL memLS_LHS_CREATE(memLS_lhsT, communicator, gnNo, nNo,
!     2         nnz_tot, ltg, colm, rowp, memLS_nFacesT)
      if(myrank.eq.master) then
      write(*,*) 'about to call svLS_LHS_CREATE for temp' ! Arsen
      call flush(0)
      endif
         CALL svLS_LHS_CREATE(svLS_lhsT, communicator, gnNo, nNo,
     2         nnz_tot, ltg, colm, rowp, svLS_nFacesT)

            faInT = 1
            facenNoT = 0
            DO i=1, nshg
               IF (btest(iBC(i),1))  facenNoT = facenNoT + 1  ! Temperature check
            END DO
            ALLOCATE(gNodesT(facenNoT), sVT(1,facenNoT))
			gNodesT = 0
            sVT = 0D0
            j = 0
            DO i=1, nshg
               IF (btest(iBC(i),1)) THEN
                  j = j + 1
                  gNodesT(j) = i
                  IF (.NOT.BTEST(iBC(i),1)) sV(1,j) = 1D0
               END IF
            END DO
!            CALL memLS_BC_CREATE(memLS_lhsT, faInT, facenNoT,
!     2         1, BC_TYPE_Dir, gNodesT, sVT)
         CALL svLS_BC_CREATE(svLS_lhsT, faInT, facenNoT, ! Arsen
     2         1, BC_TYPE_Dir, gNodesT, sVT)

        end if

        if (iLSet.gt.0) then
!               memLS_nFaces_sc = 0
               svLS_nFaces_sc = 0 ! Arsen
!            CALL memLS_LHS_CREATE(memLS_lhs_sc, communicator, gnNo, nNo,
!     2         nnz_tot, ltg, colm, rowp, memLS_nFaces_sc)
      if(myrank.eq.master) then
      write(*,*) 'about to call svLS_LHS_CREATE for scalar' ! Arsen
      call flush(0)
      endif
      CALL svLS_LHS_CREATE(svLS_lhs_sc, communicator, gnNo, nNo,
     2         nnz_tot, ltg, colm, rowp, svLS_nFaces_sc)

!            faIn = 1
!            facenNo = 0
!            sV = 0D0
!            ALLOCATE(gNodes(facenNo), sV(nsd,facenNo))
!            j = 0
!            CALL memLS_BC_CREATE(memLS_lhs_sc, faIn, facenNo,
!     2         nsd, BC_TYPE_Dir, gNodes, sV)
         if (myrank.eq.master) 
     &  write(*,*) 'Setting no-BC for level set and ls distance'
 
        end if

                ELSE   ! Acusim
       do isolsc=1,nsclrsol
         lesId       = numeqns(isolsc+1)
         eqnType     = 2
         nDofs       = 1
         presPrjflag = 0        
         nPresPrjs   = 0       
         prjFlag     = 1
         indx=isolsc+2-nsolt ! complicated to keep epstol(2) for
                             ! temperature followed by scalars
         call sleep(iwaittime)

!!!         call myfLesNew(lesId,            41994, ! assume never use LIBLES
!!!     &                 eqnType,
!!!     &                 nDofs,          minIters,       maxIters,
!!!     &                 nKvecs,         prjFlag,        nPrjs,
!!!     &                 presPrjFlag,    nPresPrjs,      epstol(indx),
!!!     &                 prestol,        verbose,        statssclr,
!!!     &                 nPermDimsS,     nTmpDimsS,      servername)
       enddo

       allocate (apermS(nshg,nPermDimsS,nsclrsol))
       allocate (atempS(nshg,nTmpDimsS))  !they can all share this
	   apermS = 0
	   atempS = 0

        END IF ! Solver choice

c
c  Assume all scalars have the same size needs
c
c          if (myrank .eq. master) write(*,*) 'l. 285 itrdrv.f'

       allocate (lhsS(nnz_tot,nsclrsol))
	   lhsS = 0
c
c actually they could even share with atemp but leave that for later
c
      else
         nPermDimsS = 0
         nTmpDimsS  = 0
      endif
c
c...  prepare lumped mass if needed
c
c      if((flmpr.ne.0).or.(flmpl.ne.0)) call genlmass(x, shp,shgl)
      call genlmass(x, shp, shgl, iBC, iper, ilwork)
c!... compute element volumes
c
      allocate(elem_local_size(numel))
	  elem_local_size = 0
      if (numelb .gt. 0) then
        allocate(elemb_local_size(numelb))
		elemb_local_size = 0
      else
	allocate(elemb_local_size(1))
	elemb_local_size = 0
      endif

      call getelsize(x,  shp,  shgl,  elem_local_size,
     &               shpb, shglb,  elemb_local_size,
     &               elemvol_global)  !Get global elem vol for Matt Thomas Control Force


c
c! Can determine psuedo time step for redistancing now that 
c! element size is known
c!
c!      if (i_dtlset_cfl .eq. 1) then
c!        elem_size_min = minval(elem_local_size)
c!        if(numpe. gt. 1) then
c!           call MPI_ALLREDUCE (elem_size_min, elem_size_mintmp, 1,
c!     &       MPI_DOUBLE_PRECISION,MPI_MIN, MPI_COMM_WORLD,ierr)
c!        else
c!           elem_size_mintmp = elem_size_min
c!        endif
c!        elem_size_min = elem_size_mintmp
c!        dtlset = dtlset_cfl * elem_size_min
c!      endif 
c
c! Initialize Level Set CFL array
c
       CFLls = zero
c
c! set flag for freezing LS Scalar 2
c
            if (iSolvLSSclr2.eq.2) then
              allocate (iBCredist(nshg))
              allocate (BCredist(nshg))
              BCredist(:) = zero
              iBCredist(:) = 0
            endif
c
c! Redistancing option of fixing phi of primary vertices
c
        i_primvert = 0
        if (i_primvert .eq. 1) then
          do inode = 1, nshg
           primvert(inode) = 0
           primvertval(inode,1) = zero 
           primvertval(inode,2) = zero
c!           write(*,*) "primvertval, inode = ", primvertval(inode,1),
c!     &                                  primvertval(inode,2), inode
          enddo
        endif
c!           if (myrank .eq. master) write(*,*) 'initialization is done'

c
c!.... -----------------> End of initialization <-----------------
!       Control force capability initialization, Jun, 2014 Fall
      if(iCForz.eq.1)then
         allocate (xcf(nstep(1)))
         allocate (ycf(nstep(1)))
         allocate (zcf(nstep(1)))
         xcf = zero
         ycf = zero
         zcf = zero

         numoldyhistind = numts_histyavg - 1
         if(myrank.eq.master) write(*,*) 'numoldyhistind =',
     &                                    numoldyhistind
         allocate (ycf_old(numoldyhistind))
		 ycf_old = 0
         CALL      CFrestar(ycf_old_log, numoldyhistind,
     &                      ixcf,        iycf,       izcf,
     &                      sumxcf_o,    sumycf_o,   sumzcf_o,
     &                      oldcf_dt,    oldcf_dtlset)

      end if      !iCForz

c!.... Initialize for Matt Talley's Coalescence Control
      if (coalcon.eq.1) then
         avgxcoordold(:) = -1.0d3 
         avgycoordold(:) = -1.0d3 
         avgzcoordold(:) = -1.0d3

         avgxcoordold2(:) = -1.0d3
         avgycoordold2(:) = -1.0d3
         avgzcoordold2(:) = -1.0d3 

         app_time(:,2) = 0.0d0
         coalcon_rem(:) = 1
         itrtimestp = 0.0d0
      endif

        N_Lsc = 0        
c.....open the necessary files to gather time series
c

      lstep0 = lstep+1
c
c.... loop through the time sequences
c

      sumtime = 0d0

      do 3000 itsq = 1, ntseq
         itseq = itsq

CAD         tcorecp1 = second(0)
CAD         tcorewc1 = second(-1)
c
c.... set up the time integration parameters
c         
         nstp   = nstep(itseq)
         nitr   = niter(itseq)
         LCtime = loctim(itseq)
         dtol(:)= deltol(itseq,:)

         call itrSetup ( y, acold )
c          if (myrank .eq. master) write(*,*) 'l. 382 itrdrv.f'

c THIS SHOULD GO IN A SUBROUTINE
c...initialize the coefficients for the impedance convolution,
c   which are functions of alphaf so need to do it after itrSetup
         if(numImpSrfs.gt.zero) then
            allocate (ConvCoef(ntimeptpT+2,2)) !same time discret. for all imp. BC
			ConvCoef = 0
            do j=1,ntimeptpT+2
                ConvCoef(j,:)=Delt(1)/2.0
            enddo
            ConvCoef(1,1)=ConvCoef(1,1)*(1.0-alfi)*(1.0-alfi)
            ConvCoef(1,2)=zero
            ConvCoef(2,2)=-ConvCoef(2,2)*(1.0-alfi)*(1.0-alfi)
            ConvCoef(ntimeptpT+1,1)=ConvCoef(ntimeptpT+1,1)*
     &                              alfi*(2.0-alfi)
            ConvCoef(ntimeptpT+2,2)=ConvCoef(ntimeptpT+2,2)*alfi*alfi
            ConvCoef(ntimeptpT+2,1)=zero  
            ConvCoef=ConvCoef/(ntimeptpT*Delt(1)) !divide by period T=N*dt
c
c...calculate the coefficients for the impedance convolution
c 
            allocate (ImpConvCoef(ntimeptpT+2,numImpSrfs))
			ImpConvCoef = 0
            do j=2,ntimeptpT+1
                ImpConvCoef(j,:) = ValueListImp(j-1,:)*ConvCoef(j,2)
     &                             + ValueListImp(j,:)*ConvCoef(j,1)  
            enddo
            ImpConvCoef(1,:) = ValueListImp(1,:)*ConvCoef(1,1)
            ImpConvCoef(ntimeptpT+2,:) = 
     &           ValueListImp(ntimeptpT+1,:)*ConvCoef(ntimeptpT+2,2)
         endif
c
c  find the last solve of the flow in the step sequence so that we will
c         know when we are at/near end of step
c
c         ilast=0
c          if (myrank .eq. master) write(*,*) 'l. 417 itrdrv.f'

         nitr=0  ! count number of flow solves in a step (# of iterations)
         do i=1,seqsize
            if(stepseq(i).eq.0) nitr=nitr+1
         enddo

         if (numpe > 1) call MPI_BARRIER(MPI_COMM_WORLD, ierr)
         if(myrank.eq.0)  then
            tcorecp1 = TMRC()
         endif

c
c.... loop through the time steps
c
         istop=0
         rmub=datmat(1,2,1)
         if(rmutarget.gt.0) then
            rmue=rmutarget
         else
            rmue=datmat(1,2,1) ! keep constant
         endif


!       do 2000 istp = 1, nstp
!        do istp = 1, nstp
        istp = 1
        do while (istp .le. nstp)
           if(iBT.eq.1 .and. myrank.eq.master)
     &        allocate( avg_info(i_num_bubbles,24) )
              avg_info = 0
           if(iBT.eq.1 .and. iBK.eq.1 .and. ibreakupFlag.eq.1) then
              nBubbleLTS = i_num_bubbles
              call breakupConfirmer(y, banma)
              if(nBubbleLTS.lt.i_num_bubbles.and.myrank.eq.master)then
                 call OpenBubFiles(nBubbleLTS, i_num_bubbles,
     &                             C_int_adjust)
                 deallocate( avg_info )
                 allocate  ( avg_info(i_num_bubbles,24) )
				 avg_info = 0
                 deallocate(D_eq)
                 deallocate(bub2w)
                 allocate  (D_eq(i_num_bubbles))
				 D_eq = 0
                 allocate  (bub2w(i_num_bubbles,3))
				 bub2w = 0

              endif
           endif

           call rerun_check(stopjob)
           if(stopjob.ne.0) goto 2001
c
c Modify time step based on CFL number
c

! Preserving the timesteps from last simulation for CF restart
           if((i_res_cf.eq.1).and.(istp.eq.1))then
            delt(itseq) = oldcf_dt
            dtlset = oldcf_dtlset
           end if
           call calc_delt(istp)


            if (iramp.gt.0) then
              xi2 = (real(time)-real(qrts0))/(real(qrts1)-real(qrts0))
               if (time.lt.qrts0) xi2 = 0.0
               if (time.gt.qrts1) xi2 = 1.0
              datmat(1,2,2) = omu*(1.0-xi2)+xi2*tmu   ! 2nd phase viscosity 
              datmat(1,1,2) = orho*(1.0-xi2)+xi2*trho   ! 2nd phase density 

                if(myrank.eq.master) then
                 if (xi2.eq.0.0) then
                   write(*,*) "Ramping properties will start at time = ", qrts0
                 else if (xi2.eq.1.0) then
c		   write(*,*) "Ramping properties is done at time = ", qrts1
		 else
                   write(*,*) "Second rho/mu",
     &    " is changed to ", datmat(1,1:2,2), " xi = ", xi2 
		 end if
                endif
            end if     
c **************************** above lines are added on 02/11/2009 by Igor Bolotnov; updated on 03/26/2009 - geometric ramp is introduced
c.... if we have time varying boundary conditions update the values of BC.
c     these will be for time step n+1 so use lstep+1
c     
c           if(itvn.gt.0) call BCint((lstep)*Delt(1), shp, shgl, 
c    &                               shpb, shglb, x, BC, iBC)

!            if(myrank.eq.master)write(*,*)'itvn, itvbc ',itvn,itvbc
            if(itvn.gt.0) then
              if (iRBCT.eq.0) then
                call BCint(time,shp, shgl, shpb, shglb, x, BC, iBC)
              endif
              if (iextsbct.eq.1) then
                call BCint(time,shp, shgl, shpb, shglb, x, BC, iBC)
              endif
            endif
            if(itvbc.gt.0) call setBC(time, x, iBC, BC)

c
c ... calculate the pressure contribution that depends on the history 
c     for the impendance BC
c
c          if (myrank .eq. master) write(*,*) 'l. 489 itrdrv.f'


            if(numImpSrfs.gt.0) call pHist(pold,QHistImp,ImpConvCoef,
     &                                          ntimeptpT,numImpSrfs)
c
c Decay of scalars
c
           if(nsclr.gt.0 .and. tdecay.ne.1) then
              yold(:,6:ndof)=y(:,6:ndof)*tdecay
              BC(:,7:6+nsclr)= BC(:,7:6+nsclr)*tdecay
           endif

           if(nosource.eq.1) BC(:,7:6+nsclr)= BC(:,7:6+nsclr)*0.8


            if(iLES.gt.0) then  !complicated stuff has moved to
                                        !routine below
               call lesmodels(yold,  acold,     shgl,      shp, 
     &                        iper,  ilwork,    rowp,      colm,
     &                        nsons, ifath,     x,   
     &                        iBC,   BC)

            
            endif

c.... set traction BCs for modeled walls
c
            if (itwmod.ne.0) then
               call asbwmod(yold,   acold,   x,      BC,     iBC,
     &                      iper,   ilwork,  ifath,  velbar)
            endif
c          if (myrank .eq. master) write(*,*) 'l. 521 itrdrv.f'

c
c.... -----------------------> predictor phase <-----------------------
c
c           if (myrank .eq. master) write(*,*) 'predictor phase'
            call itrPredict(yold, y,   acold,  ac ,  uold,  u)
            call itrBC (y,  ac, banma, iBC,  BC,  iper,ilwork)
        
c          if (myrank .eq. master) write(*,*) 'l. 529 itrdrv.f'

            if(nsolt.eq.1) then
               isclr=0
               call itrBCSclr (y, ac,  iBC, BC, iper, ilwork,banma)
!                do i=1,nshg 
!                if(banma(i,1).ne.0.0d0)write(*,*)i,banma(i,1)
!                enddo
            endif
            do isclr=1,nsclr
               call itrBCSclr (y, ac,  iBC, BC, iper, ilwork,banma)
            enddo
            iter=0
            ilss=0  ! this is a switch thrown on first solve of LS redistance
c
c set the initial tolerance for the redistance loop
c
            if (i_redist_loop_flag.eq.1) then
              redist_toler_previous = 100.0
            endif
c
c LOOP OVER SEQUENCES
c
            istepc = 1
            iloop = .true.
            i_redist_counter=0
c            do istepc=1,seqsize
             do while (iloop) 
               icode=stepseq(istepc)
c
               if(mod(icode,5).eq.0) then ! this is a solve
                  isolve=icode/10
                  if(icode.eq.0) then ! flow solve (encoded as 0)
c
                     iter   = iter+1
                     ifuncs(1)  = ifuncs(1) + 1
c     
                     Force(1) = zero
                     Force(2) = zero
                     Force(3) = zero
                     HFlux    = zero
                     entrop   = zero
                     lhs = 1 - min(1,mod(ifuncs(1)-1,LHSupd(1))) 

                     avgxcoordf(:) = -1.0d3
                     avgycoordf(:) = -1.0d3
                     avgzcoordf(:) = -1.0d3


                  call SolFlow(y,             ac,         u,
     &			       banma,          
     &                         yold,          acold,      uold,
     &                         x,             iBC,
     &                         BC,            res,
     &                         nPermDims,     nTmpDims,   aperm,
     &                         atemp,         iper,          
     &                         ilwork,        shp,        shgl,
     &                         shpb,          shglb,      rowp,     
     &                         colm,          lhsK,       lhsP,
     &                         solinc,        rerr,       sumtime,
!     &                         memLS_lhs,     memLS_ls,   memLS_nFaces,
     &                         svLS_lhs, svLS_ls, svLS_nFaces, ! Arsen
     &                         elemvol_global,
     &                         avgxcoordf,    avgycoordf, avgzcoordf)


c!....Matt Talley's Coalescence Contorl
                   if (coalcon.eq.1) then
                       if (coaltimtrak.eq.1) then

                       do k = 1, coalest
                          avgxcoordold(k) = avgxcoordf(k)
                          avgycoordold(k) = avgycoordf(k)
                          avgzcoordold(k) = avgzcoordf(k)

             if (avgxcoordold(k).gt.-1.0d3) then
             if (myrank.eq.master) write(*,*) 'Coalescence',
     &          ' Event #: ', k
             if (myrank.eq.master) write(*,*) 'x average',
     &          ' position:', avgxcoordold(k)
             if (myrank.eq.master) write(*,*) 'y average',
     &          ' position:', avgycoordold(k)
             if (myrank.eq.master) write(*,*) 'z average',
     &          ' position:', avgzcoordold(k)
             endif

                        enddo ! k

                        else

                           itrtimestp = Delt(1)

             call CoalescAppTime (avgxcoordf, avgycoordf,
     &                            avgzcoordf, avgxcoordold2,
     &                            avgycoordold2, avgzcoordold2,
     &                            app_time, itrtimestp)
                      endif ! coaltimtrak
                  endif ! coalcon

                  else          ! scalar type solve
                     if (icode.eq.5) then ! Solve for Temperature
                                ! (encoded as (nsclr+1)*10)
                        isclr=0
                        ifuncs(2)  = ifuncs(2) + 1
                        j=1
                     else       ! solve a scalar  (encoded at isclr*10)
                        isclr=isolve  
                        ifuncs(isclr+2)  = ifuncs(isclr+2) + 1
                        j=isclr+nsolt
c!  Modify psuedo time step based on CFL number for redistancing
                        if((iLSet.eq.2).and.(ilss.ge.1).and.
     &                     (i_dtlset_cfl.eq.1).and.
     &                     (isclr.eq.2)) then

                           call calc_deltau()

                           Delt(1) = dtlset ! psuedo time step for level set
                           Dtgl = one / Delt(1)

                           ilss = ilss+1
                        endif
c
                        if((iLSet.eq.2).and.(ilss.eq.0)
     &                       .and.(isclr.eq.2)) then 
                           ilss=1        ! throw switch (once per step)
                           y(:,7)=y(:,6) ! redistance field initialized
                           ac(:,7)   = zero
!----------------------------------------------------------------------
                           if (iSolvLSSclr2.eq.2)  then
                             call get_bcredist(x,y,iBCredist,BCredist,
     &                                      primvert, primvertval(:,1))
                             primvertval(:,1) = BCredist(:)
                             ib=5+isclr
                             ibb=ib+1
                             do inode = 1, nshg
                              if (iBCredist(inode).eq.1) then
                               if (btest(iBC(inode),ib)) then
                              write(*,*) "WARNING -- Bit 7 already set"
                               endif
                              endif
                             enddo
c
                             where (iBCredist.eq.1) 
                              iBC(:) = iBC(:) + 128   ! set scalar 2 (bit 7)    
                              BC(:,ibb) = BCredist(:)
                             endwhere
                             numpv = 0
                             numpvset = 0
                             do inode = 1, nshg
                               if (primvert(inode) .gt. 0) then
                               numpv = numpv + 1
                                 if (primvert(inode).eq.2) then
                                   numpvset = numpvset + 1
                                 endif
                               endif
                             enddo
                             write(*,*) lstep+1,
     &                                  " Primary Verts: set/exist = ",
     &                                  numpvset, numpv
                           endif
                           call itrBCSclr (  y,  ac,  iBC,  BC, iper,
     &                          ilwork,banma)
c     
c!....store the flow alpha, gamma parameter values and assigm them the 
c!....Backward Euler parameters to solve the second levelset scalar
c     
                           alfit=alfi
                           gamit=gami
                           almit=almi
                           Deltt=Delt(1)
                           Dtglt=Dtgl
                           alfi = 1
                           gami = 1
                           almi = 1
c!     Delt(1)= Deltt ! Give a pseudo time step
                           Delt(1) = dtlset ! psuedo time step for level set
                           Dtgl = one / Delt(1)
                        endif  ! level set eq. 2
                     endif ! deciding between temperature and scalar

                     lhs = 1 - min(1,mod(ifuncs(isclr+2)-1,
     &                                   LHSupd(isclr+2))) 

                     if((isclr.eq.1.and.iSolvLSSclr1.eq.1) .or. 
     &                  (isclr.eq.2.and.iSolvLSSclr2.eq.1) .or.
     &                  (isclr.eq.2.and.iSolvLSSclr2.eq.2)) then
                      lhs=0
c                      write(*,'(a)') 'This big block is executed'
                      call SolSclrExp(y,          ac,        yold,
     &                         acold,         x,         iBC,
     &                         BC,            nPermDimsS,nTmpDimsS,  
     &                         apermS(1,1,j), atempS,    iper,          
     &                         ilwork,        shp,       shgl,
     &                         shpb,          shglb,     rowp,     
     &                         colm,          lhsS(1,j), 
     &                         solinc(1,isclr+5), CFLls)
                     else
!                     if (myrank.eq.master) write(*,*) 'icode = ', icode
                      if (icode.eq.5) then    ! temperature
!                       write(*,*)istepc
                     call SolSclr(y,          ac,        u,
     &                         yold,          acold,     uold,
     &                         x,             iBC,
     &                         BC,            nPermDimsS,nTmpDimsS,  
     &                         apermS(1,1,j), atempS,    iper,          
     &                         ilwork,        shp,       shgl,
     &                         shpb,          shglb,     rowp,     
     &                         colm,          lhsS(1,j), 
     &                         solinc(1,isclr+5), CFLls,
     &                         svLS_lhsT, svLS_sc, svLS_nFacesT) ! Arsen
!     &               memLS_lhsT,     memLS_sc,  memLS_nFacesT)   ! Igor 10/2012
                       else    ! other scalar
                     call SolSclr(y,          ac,        u,
     &                         yold,          acold,     uold,
     &                         x,             iBC,
     &                         BC,            nPermDimsS,nTmpDimsS,
     &                         apermS(1,1,j), atempS,    iper,
     &                         ilwork,        shp,       shgl,
     &                         shpb,          shglb,     rowp,
     &                         colm,          lhsS(1,j),
     &                         solinc(1,isclr+5), CFLls,
     &               svLS_lhs_sc,     svLS_sc,  svLS_nFaces_sc)   ! Arsen
!     &               memLS_lhs_sc,     memLS_sc,  memLS_nFaces_sc)   ! Igor 10/2012
                        end if ! temp/scalar condition from icode
                     endif
                        
                  endif         ! end of scalar type solve

               else ! this is an update  (mod did not equal zero)
                  iupdate=icode/10  ! what to update
                  if(icode.eq.1) then !update flow  
                     call itrCorrect ( y,    ac,    u,   solinc)
                     call itrBC (y,  ac, banma, iBC,  BC, iper, ilwork)
                  else  ! update scalar
                     isclr=iupdate  !unless
                     if(icode.eq.6) isclr=0
                     if(iRANS.lt.0)then  ! RANS
                        call itrCorrectSclrPos(y,ac,solinc(1,isclr+5))
                     else
                        call itrCorrectSclr (y, ac, solinc(1,isclr+5))
                     endif
                     if (ilset.eq.2 .and. isclr.eq.2)  then
                        if (ivconstraint .eq. 1) then
                           call itrBCSclr (  y,  ac,  iBC,  BC, iper,
     &                          ilwork,banma)
c                    
c ... applying the volume constraint on second level set scalar
c
                           call solvecon (y,    x,      iBC,  BC, 
     &                          iper, ilwork, shp,  shgl)
c
                        endif   ! end of volume constraint calculations
                     endif      ! end of redistance calculations
c                     
                        call itrBCSclr (  y,  ac,  iBC,  BC, iper,
     &                       ilwork,banma)
c                    
c
c ... update the old value for second level set scalar
c
                     if (ilset.eq.2 .and. isclr.eq.2)  then
                         call itrUpdateDist( yold, acold, y, ac)
                     endif   

                     endif      ! end of flow or scalar update
                  endif         ! end of switch between solve or update
c
c!** Conditions for Redistancing Loop **
c! Here we test to see if the following conditions are met:
c!	no. of redistance iterations < i_redist_max_iter
c!	residual (redist_toler_curr) > redist_toler
c! If these are true then we continue in the redistance loop
c
                 if(i_redist_loop_flag.eq.1) then
                   if (icode .eq. 21) then ! only check after a redistance update
                     if((ilset.eq.2).and.(isclr.eq.2)) then !redistance condition
       if (redist_toler_curr.gt.redist_toler.or.
     &  i_redist_counter.eq.0) then !condition 1
              if (i_redist_counter.lt.i_redist_max_iter) then ! condition 2
                        i_redist_counter = i_redist_counter + 1
                        istepc = istepc - 2  ! repeat the 20 21 step
                        if(redist_toler_curr.gt.redist_toler_previous)
     &                  then
c                         if(myrank.eq.master) then
c                          write(*,*) "Warning: diverging!"
c                         endif
                        endif
                       else
                        iloop = .false. 

                        if(myrank.eq.master) then  
                         write(*,*) "Exceeded Max # of the iterations: "
     &                              , i_redist_max_iter

        write(*,*)'The redistance loop has converged ',N_Lsc,' times'

                        endif
                       endif
                       redist_toler_previous=redist_toler_curr
                      else

                       if(myrank.eq.master) then
                        write(*,*) "Redistance loop converged in ",
     &                       i_redist_counter," iterations"

                        N_Lsc = N_Lsc + 1 !this indicates that level set
                                          !met coverge criteria

       write(*,*)'The redistance loop has converged ',N_Lsc,' times'

                       endif
                       iloop = .false. 
                      endif
                     endif
                   endif !end of the redistance condition
                 endif !end of the condition for the redistance loop
c
                 if (istepc .eq. seqsize) then
                   iloop = .false.
                 endif
                 istepc = istepc + 1
c                 write(*,*)istepc
c
c**End of loop condition for Redistancing equation**
c		 		  
               end do      ! end while loop over sequence steps


c
c Check if interface has moved into region of larger interface
c
             if ((iLSet.eq.2).and.(i_check_prox.eq.1)) then 
               call check_proximity(y, stopjob)
	       if(stopjob.ne.0) then
                   lstep = lstep + 1
                   goto 2001
               endif
             endif          
c
c print out phasic volume for Level Set
c
             if (iLSet.eq.2) then
               vf =  phvol(2)/(phvol(1)+phvol(2))
!----------------------------------------------------------------------
!       Here is the interface adjustment, The C_int_adjust will be used in
!       e3source.f and cannot be too large
!                                                JUN, Feb, 2014
!       When 'Flow regime transition' is enabled, the domain volume fraction
!       will be manipulated to model the flow regime transition process, for
!       example the liquid will evaporated for the presence of volumeric heat
!       source.
!                                                JUN, Oct, 2015
!----------------------------------------------------------------------
        vf_now = vf
        IF(iFT.eq.0) THEN
        if(myrank.eq.masteri.and.lstep+1.eq.1)then
              write(*,'(1x,A,F14.7,A,2x,A,I2,A,ES14.7)')
     &        'vf_now is ',vf_now*100,'%,','iuse_vfcont_cap ',
     &        iuse_vfcont_cap, ',    C_int_cap =',C_int_cap
        end if

        C_int_adjust = C_int_adjust - (vf - vf_target)*vfcontrcoeff

        if(iuse_vfcont_cap .eq. 1)then
           if(C_int_adjust .gt. C_int_cap)then
              C_int_adjust =  C_int_cap
           else if(C_int_adjust.lt.-C_int_cap)then
              C_int_adjust = -C_int_cap
           end if
        end if
        ELSEIF(iFT.eq.1) THEN
!       We don't care about the following variables, they are kept just to
!       preseve the code structure
        i_num_bubbles = 0
        nBubbleLTS    = i_num_bubbles
        Nbubtot = 0
        Nbubdef = 0
        if(myrank.eq.master.and.istp.eq.1)then
           vf_initial = vf
           write(*,*)'Bubble volume will grow!'
           call OpenBubFiles(nBubbleLTS, i_num_bubbles,
     &                       C_int_adjust)
        end if
!       Fang: vf_tmp is used as the current expected vf in flow regime
!       transition simulations. (11/01/2015)
        vf_tmp = vf_initial +
     &          (real(istp)/real(nstp))*(vf_target-vf_initial)
        C_int_adjust = (vf_tmp-vf_now) * vfcontrcoeff

        if(myrank.eq.master) then
        write(*,*) 'The current expected alpha: ', vf_tmp
        open(unit=740,file='../bubStat/alpha.dat',status="old",
     &       action="write", position="append", iostat=ierror)
        if(ierror /= 0) STOP "Error creating file 740"
        write(740,873) lstep+1, Nbubtot, Nbubdef, C_int_adjust, vf,
     &  time
        close(740)
        endif

 873  format(I6, 1x, I4, 1x, I4, 3x, ES14.7, 1x, ES14.7, 1x, ES14.7)
        ENDIF
!----------------------------------------------------------------------

        if (myrank.eq.master) then
       write (ivhist,800) lstep+1, phvol(1), phvol(2), vf, epsilon_lsd
        if(vf_target.gt.0.0d0) then
       write (*,801) lstep+1, vf*1.0D+02, 
     &  (vf/vf_target-1.0D+00)*1.00D+02, CFLfl_max,
     &  CFLbuint_max
        else
       write (*,802) lstep+1, vf*1.0D+02,
     &  CFLfl_max,
     &  CFLbuint_max
        endif

       write(1322,*) lstep+1, real(time), vf*1.0D+02

        endif
  
      endif  !iLSet.eq.2
 800    format(1x,i6,4e15.6)
 801    format(i6,' Void Fraction: ',F15.10,
     1 '%; Mass Error: ',F15.10,'%; Domain CFL: ',
     2  ES9.3,'; Interface CFL: ',ES9.3)
 802    format(i6,' Void Fraction: ',F15.10,
     1 '%; ',' Domain CFL: ',
     2  ES9.3,'; Interface CFL: ',ES9.3)

c
c.... obtain the time average statistics
c
c	if (myrank.eq.master) write(*,*) 'GetStats is about to run' 
            if (ioform .eq. 2) then

               call stsGetStats( y,      yold,     ac,     acold,
     &                           u,      uold,     x,
     &                           shp,    shgl,     shpb,   shglb,
     &                           iBC,    BC,       iper,   ilwork,
     &                           rowp,   colm,     lhsK,   lhsP )
            endif

c     
c  Find the solution at the end of the timestep and move it to old
c
c  
c ...First to reassign the parameters for the original time integrator scheme
c
            if((iLSet.eq.2).and.(ilss.gt.0)) then 
               alfi =alfit
               gami =gamit
               almi =almit 
               Delt(1)=Deltt
               Dtgl =Dtglt
            endif          
            call itrUpdate( yold,  acold,   uold,  y,    ac,   u)
            call itrBC (yold, acold, banma, iBC,  BC,  iper,ilwork)

c        if (myrank.eq.master) write(*,*) 'itrBC is done'

!----------------------------------------------------------------------
!       Jun Fang Coalescence control based on bubble tracking information
!----------------------------------------------------------------------
        if(icoalCtrl.eq.1) then
           deallocate( coalCenter )
!       Obtain the centers of coalescence events
           if(myrank.eq.master.and.i_num_bubbles.gt.1) 
     &     call coalescenceDetection()
           if(numpe > 1) then
              call MPI_Bcast(ncoalEvent,1,MPI_INTEGER,
     &                    master,MPI_COMM_WORLD,ierr)
              if(myrank.ne.master) then
                allocate(coalCenter(ncoalEvent,3))
				coalCenter = 0
              endif
              call MPI_BARRIER(MPI_COMM_WORLD,ierr)
              call MPI_Bcast(coalCenter,ncoalEvent*3,
     &             MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
           endif
        endif !end of coalescence part
!----------------------------------------------------------------------
!       print out bubble information and some other processing at the end of
!       each time iteration
!----------------------------------------------------------------------
        if(iBT.eq.1 .and. i_num_bubbles.ne.0) then
           Nbubtot    = 0          !total number of bubbles
           Nghost     = 0          !number of ghost bubbles
           deallocate(bub_cent)
           if(myrank.eq.master) call reCountBub()
           if(numpe > 1) then
              call MPI_Bcast(Nbubtot,1,MPI_INTEGER,
     &                    master,MPI_COMM_WORLD,ierr)
              call MPI_Bcast(Nghost, 1,MPI_INTEGER,
     &                    master,MPI_COMM_WORLD,ierr)
              call MPI_Bcast(i_num_bubbles, 1,MPI_INTEGER,
     &                    master,MPI_COMM_WORLD,ierr)
           endif
           allocate(bub_cent(i_num_bubbles+Nghost,5))
           bub_cent = zero
!       breakupSeeder will be reshaped if total number of bubbles changes 
           if(iBK.eq.1 .and. nBubbleLTS .lt. i_num_bubbles) then
              allocate(breakupSeederTMP(nBubbleLTS,6))
			  breakupSeederTMP = 0
              breakupSeederTMP(1:nBubbleLTS,:) =
     &           breakupSeeder(1:nBubbleLTS,:)
              deallocate(breakupSeeder)

              allocate(D_eqTMP(nBubbleLTS))
              allocate(bub2wTMP(nBubbleLTS,3))
			  D_eqTMP = 0
			  bub2wTMP = 0
              deallocate(D_eq)
              deallocate(bub2w)

              allocate(breakupSeeder(i_num_bubbles,6))
              breakupSeeder(:,:) = 0
                 breakupSeeder(1:nBubbleLTS,:) =
     &        breakupSeederTMP(1:nBubbleLTS,:)
              deallocate(breakupSeederTMP)
  

              allocate(D_eq(i_num_bubbles))
              allocate(bub2w(i_num_bubbles,3))
              D_eq(:) = 0.0d0
              bub2w(:,:) = 0.0d0
              D_eq(1:nBubbleLTS) = D_eqTMP(1:nBubbleLTS)
              bub2w(1:nBubbleLTS,:) = bub2wTMP(1:nBubbleLTS,:)
              deallocate(D_eqTMP)
              deallocate(bub2wTMP)
              
              nBubbleLTS = i_num_bubbles
!              if(myrank.eq.master)write(*,*)D_eq,nBubbleLTS
           endif

           if(myrank.eq.master) then
!       save the bubble curvature for the detection of bubble breakup in the
!       next time step
              if(iBK.eq.1) call breakupDetector(ibreakupFlag)
!       write out bubble information, and prepare the bubble center information
!       (including real bubbles and ghost bubbles at periodic faces)
              write(*,'(1x,A,I5)') 'Nbubtot =', Nbubtot
              write(*,'(1x,A,I5)') 'Nghost  =', Nghost
              call BubPrintOut(vf)
              write(*,*) 'BubPrintOut is done!'
              deallocate(avg_info)
           endif
!       broadcast the bubbles centers and radii for the marker field updating
!       algorithm
           if(numpe > 1) then
              call MPI_Bcast(bub_cent,(i_num_bubbles+Nghost)*5,
     &             MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
              if(iBK.eq.1) call MPI_Bcast(ibreakupFlag,1,
     &             MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
              if(iBK.eq.1) call MPI_Bcast(breakupSeeder,i_num_bubbles*6,
     &             MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
              if(iBK.eq.1) call MPI_Bcast(bub2w,i_num_bubbles*3,MPI_DOUBLE_PRECISION,
     &             master,MPI_COMM_WORLD,ierr)
              if(iBK.eq.1) call MPI_Bcast(D_eq,i_num_bubbles ,MPI_DOUBLE_PRECISION,
     &             master,MPI_COMM_WORLD,ierr)
!              if(iBK.eq.1) write(*,*)D_eq

           endif
        endif     !iBT

!----------------------------------------------------------------------

            istep = istep + 1
            lstep = lstep + 1
            time  = time + delt(itseq)
c
c... Post process data and write out
c
            if ((i_gradphi.eq.1) .and. (mod(lstep, ntout).eq.0)) then
	      idflx = 0 
              if(idiff >= 1 )  idflx= (nflow-1) * nsd
              if (isurf == 1) idflx=nflow*nsd
	      call getgradphi(x, y, shp,       shgl,
     &                         shpb,          shglb, gradphi)
              gradphimag(:) = ( gradphi(:,1)**2 + 
     &                          gradphi(:,2)**2 + 
     &                          gradphi(:,3)**2 )**0.5
              maxgradphi = maxval(gradphimag)
	      write(*,1000) maxgradphi, myrank
              call MPI_BARRIER(MPI_COMM_WORLD,ierr)
	      call write_gradphi(myrank, lstep, nshg, 3, 
     &                           x, y, gradphi, gradphimag)
 1000 format ("Maximum LS gradient = ",f12.6," on proc",i6)
	    endif
c
c...Reset BC codes of primary vertices
c (need to remove dirchlet bc on primary vertices
c  by removing the 8th bit (val=128)
c
            if((iLSet.eq.2).and.(ilss.eq.1).and.(iSolvLSSclr2.eq.2))
     &      then
               where (iBCredist.eq.1)
                 iBC(:) = iBC(:) - 128   ! remove prescription on scalar 2
               endwhere
            endif

c
c... Write out Redistancing primary vertice information
c
            if ((i_primvert.eq.1) .and. (mod(lstep, ntout).eq.0)) then
                primvertval(:,2) = primvert(:) * 1.0
                call write_primvert(myrank, lstep, nshg, 2,
     &                              primvertval)
            endif

!... Processing the control forces Jun, 2014 Fall
        if(iCForz.eq.1)then

        CALL       CFpostprocs(istp,       numoldyhistind,
     &                         ixcf,       iycf,         izcf,
     &                         sumxcf_o,   sumycf_o,     sumzcf_o,
     &                         sumycf_old, ycf_old_log)

        endif !iCForz
c!
c! .. write out the solution
c!
            if ((irs .ge. 1) .and. (mod(lstep, ntout) .eq. 0)) then
!----------------------------------------------------------------------
!       Simulation Steering can be activated when PHASTA tries to write out the
!       restart files
!                                               Jun, May 2014
!----------------------------------------------------------------------
        OPEN (UNIT = 672,FILE = 'numstart.dat', STATUS = 'old')
        READ (672, *) inumstart
        CLOSE (672)
        IF ((inumstart .LT. (lstep - ntout)) .AND.
     &  (mod(inumstart, ntout) .EQ. 0)) THEN
           IF(myrank .eq. master) THEN
              WRITE(*,*)
              WRITE(*,*) 'The # in numstart.dat was changed manually!'
              WRITE(*,*) 'Restart simulation from timestep', inumstart
              WRITE(*,*)
           END IF
           istp = istp - (lstep - inumstart)
!       call input_fform.cc to get updated simulation parameters
           inpfname = 'solver.inp'
           CALL input_fform(inpfname)
!       reload the restart files at designated time step
           CALL ssread(qss, acss, uss)

           if(myrank .eq. master) write(*,*) 'istp =', istp
           if(myrank .eq. master) write(*,*) 'nstp =', nstp

           y(:,1:3)=qss(:,2:4)
           y(:,4)=qss(:,1)
           if(ndof.gt.4)  y(:,5:ndof)=qss(:,5:ndof)

           ac(:,1:3)=acss(:,2:4)
           ac(:,4)=acss(:,1)
           if(ndof.gt.4)  ac(:,5:ndof)=acss(:,5:ndof)
!        update both y and yold...
           yold = y
           acold = ac
           uold = uss
        ELSE
           IF(iBT.eq.1) call banmaCorrector(banma, ibreakupFlag, yold)
           CALL restar ('out ',  yold  ,ac)
           IF(iBT.eq.1) yold(:,7) = yold(:,6)
        END IF
!----------------------------------------------------------------------
c               if(iofieldv.ne.0) then
c                 call pstdrv (y,         ac,         x,
c     &                        iBC,       BC,
c     &                        iper,      ilwork,     shp,
c     &                        shgl,      shpb,       shglb,
c     &                        ifath,     velbar,     nsons )
c
c               endif
               if(ideformwall.eq.1) 
     &              call write_displ(myrank, lstep, nshg, 3, uold ) 
               if (ioform .eq. 2) 
     &              call stsWriteStats()
            endif
c
c SUBROUTINE
c ... update the flow rate history for the impedance convolution
c    
            if(numImpSrfs.gt.zero) then
                call GetFlowQ(NewQImp,y,nsrflistImp,numImpSrfs) !flow Q for imp BC
                do j=1, ntimeptpT
                    QHistImp(j,:)=QHistImp(j+1,:)
                enddo
                QHistImp(ntimeptpT+1,1:numImpSrfs) = 
     &          NewQImp(1:numImpSrfs)

c
c.... write out the new history of flow rates to Qhistor.dat
c      
                if (((irs .ge. 1) .and. (mod(lstep, ntout) .eq. 0)).and.
     &               (myrank .eq. zero)) then
                   open(unit=816, file='Qhistor.dat',status='replace')
                   write(816,*) ntimeptpT
                   do j=1,ntimeptpT+1
                      write(816,*) (QHistImp(j,n),n=1, numImpSrfs)
                   enddo
                   close(816)
                endif
             endif
c
c.... compute the consistent boundary flux
c
            if(abs(itwmod).ne.1) then
               call Bflux ( yold,       acold,      uold,     x,
     &                      shp,       shgl,       shpb,   
     &                      shglb,     ilwork,     iBC,
     &                      BC,        iper)
            endif

c...  dump TIME SERIES

            if (exts) then
               if (mod(lstep-1,freq).eq.0) then

                  if (numpe > 1) then
                     do jj = 1, ntspts
             vartssoln((jj-1)*numvar+1:jj*numvar)=varts(jj,1:numvar)
                        ivarts=zero
                     enddo

                     do k=1,numvar*ntspts
                        if(vartssoln(k).ne.zero) ivarts(k)=1
                     enddo
           call MPI_REDUCE(vartssoln, vartssolng, numvar*ntspts,
     &                    MPI_DOUBLE_PRECISION, MPI_SUM, master,
     &                    MPI_COMM_WORLD, ierr)

              call MPI_REDUCE(ivarts, ivartsg, numvar*ntspts,
     &                    MPI_INTEGER, MPI_SUM, master,
     &                    MPI_COMM_WORLD, ierr)

                     if (myrank.eq.zero) then
                        do jj = 1, ntspts

                           indxvarts = (jj-1)*numvar
                           do k=1, numvar
                              if(ivartsg(indxvarts+k).ne.0) then ! none of the vartssoln(parts) were non zero
                                 varts(jj,k)=vartssolng(indxvarts+k)/
     &                                ivartsg(indxvarts+k)
                              endif
                           enddo ! do k
                        enddo ! do jj
                     endif !only on master
                  endif !only if numpe > 1

c Do not search anymore:
	tssearch = 0

                  if (myrank.eq.zero) then
                     do jj = 1, ntspts
			if (numvar.ge.15) then
			   if (varts(jj,15).le.0) then
			     iphase = 1
			   else
			     iphase = 0
			   end if
			else
			 iphase = 0
			end if
			varts(jj,5) = delt(itseq)*real(freq)   ! delta t has to be recorded; multiply it by the varts stepping
			recnum = (istep-1)/freq*ntspts+jj
                        if (numvar.eq.19) then
			varts(jj,19) = one  ! Put the Heps value here (depends upon varts(jj,15) = distance field value)
             if  (abs(varts(jj,15)) .le. epsilon_lsd_tmp) then
          varts(jj,19) = 0.5*(one + varts(jj,15) / epsilon_lsd_tmp +
     &                    (sin(pi*varts(jj,15)/epsilon_lsd_tmp))/pi)
              elseif (varts(jj,15) .lt. - epsilon_lsd_tmp) then
                 varts(jj,19) = zero
              endif
                         write(10101,'(2I8, I3, 19E15.7)',REC=recnum)
     &            lstep-1,jj,iphase,(varts(jj,k), k=1, 19)   !  including the distance field (k = 15), gradient (16..18) and Heps value (19)
			else if (numvar.eq.15) then
			 write(10101,'(2I8, I3, 15E15.7)',REC=recnum) 
     &            lstep-1,jj,iphase,(varts(jj,k), k=1, 15)   !  including the distance field (k = 15)
			else
c			if (jj.eq.1) write(*,*) 'numvar = ', numvar, '; point is there'
                         write(10101,'(2I8, I3, 15E15.7)',REC=recnum)
     &            lstep-1,jj,iphase,(varts(jj,k), k=1, 14), 0.0E0  
			end if
                        ifile = 1000+jj
                        if (ntspts.gt.80) then                  ! no more than 96 files can be opened simulatneously
                           fvarts='varts/varts'
                           fvarts=trim(fvarts)//trim(cname2(jj))
                           fvarts=trim(fvarts)//trim(cname2(lskeep))
                           fvarts=trim(fvarts)//'.dat'
                           fvarts=trim(fvarts)
c                           open(unit=ifile, file=fvarts,
c     &                          position='append')
                        end if
c            write(ifile,555) lstep-1, (varts(jj,k), k=1, numvar) 
c                        call flush(ifile)
                        if (ntspts.gt.80) then
c                           close(ifile)
                        end if

                       if (((irs .ge. 1) .and. (ntspts.le.80) .and.
     &                       (mod(lstep, ntout) .eq. 0))) then
c                           close(ifile)
                           fvarts='varts/varts'
                           fvarts=trim(fvarts)//trim(cname2(jj))
                           fvarts=trim(fvarts)//trim(cname2(lskeep))
                           fvarts=trim(fvarts)//'.dat'
                           fvarts=trim(fvarts)
c                           open(unit=ifile, file=fvarts,
c     &                          position='append')
                        endif !only when dumping restart
                     enddo
                  endif !only on master

                  varts(:,:) = zero ! reset the array for next step

 555              format(i6,20(2x,E12.5e2))

               endif
            endif

c
c.... update and the aerodynamic forces
c
!	if (myrank.eq.master) write(*,*) 'calling forces, l. 974 itrdrv'
            call forces ( yold,  ilwork )
            
            if((irscale.ge.0).or.(itwmod.gt.0).or.(iDNS.ne.0)) 
     &           call getvel (yold,     ilwork, iBC,
     &                        nsons,    ifath, velbar)

            if((irscale.ge.0).and.(myrank.eq.master)) then
               call genscale(yold,       x,    iBC)   ! Corrected by Igor A. Bolotnov, 07/06/2009
c               call genscale(yold,       x,       iper, 
c     &                       iBC,     ifath,   velbar,
c     &                       nsons)
            endif
c
c....  print out results.
c
            ntoutv=max(ntout,100)   ! velb is not needed so often
            if ((irs .ge. 1) .and. (mod(lstep, ntout) .eq. 0)) then
               if( (mod(lstep, ntoutv) .eq. 0) .and.
     &              ((irscale.ge.0).or.(itwmod.gt.0) .or. 
     &              ((nsonmax.eq.1).and.(iLES.gt.0))))
     &              call rwvelb  ('out ',  velbar  ,ifail)
            endif
c
c.... end of the NSTEP and NTSEQ loops
c
c
c.... -------------------> error calculation  <-----------------
c 
            if(ierrcalc.eq.1 .or. ioybar.eq.1) then
c$$$c
c$$$c compute average
c$$$c
c$$$               tfact=one/istep
c$$$               ybar =tfact*yold + (one-tfact)*ybar

c compute average
c ybar(:,1) - ybar(:,3) is average velocity components
c ybar(:,4) is average pressure
c ybar(:,5) is average speed
c averaging procedure justified only for identical time step sizes
c istep is number of time step
c
               tfact=one/istep

c ybar to contain the averaged ((u,v,w),p)-field
c and speed average, i.e., sqrt(u^2+v^2+w^2)

               ybar(:,1) = tfact*yold(:,1) + (one-tfact)*ybar(:,1)
               ybar(:,2) = tfact*yold(:,2) + (one-tfact)*ybar(:,2)
               ybar(:,3) = tfact*yold(:,3) + (one-tfact)*ybar(:,3)
               ybar(:,4) = tfact*yold(:,4) + (one-tfact)*ybar(:,4)
               ybar(:,5) = tfact*sqrt(yold(:,1)**2+yold(:,2)**2+
     &                     yold(:,3)**2) + (one-tfact)*ybar(:,5)
c
c compute rms
c
               rerr(:, 7)=rerr(:, 7)+(yold(:,1)-ybar(:,1))**2
               rerr(:, 8)=rerr(:, 8)+(yold(:,2)-ybar(:,2))**2
               rerr(:, 9)=rerr(:, 9)+(yold(:,3)-ybar(:,3))**2
               rerr(:,10)=rerr(:,10)+(yold(:,4)-ybar(:,4))**2
c
c.....smooth the error indicators
c
               do i=1,ierrsmooth
                 call errsmooth( rerr, x, iper, ilwork, shp, shgl, iBC )
               end do
c
c.... open the output file
c
               if ((irs .ge. 1) .and. (mod(lstep, ntout) .eq. 0)) then
                 iqoldsiz=nshg*ndof*2
c                 call write_error(myrank, lstep, nshg, numerr, rerr )   ! Avoid to double write the error (see below)
               endif
            endif
            
            if(istop.eq.1000) exit ! stop when delta small (see rstatic)
!2000 continue
        istp = istp + 1
        enddo 
 2001    continue
        

CAD         tcorecp2 = second(0)
CAD         tcorewc2 = second(-1)
         
CAD         write(6,*) 'T(core) cpu-wallclock = ',tcorecp2-tcorecp1,
CAD     &                                        tcorewc2-tcorewc1

         if (numpe > 1) call MPI_BARRIER(MPI_COMM_WORLD, ierr)
         if(myrank.eq.0)  then
            tcorecp2 = TMRC()
            write(6,*) 'T(core) cpu = ',tcorecp2-tcorecp1
         endif

 3000 continue
c      open(unit=87,file="spmass.dat",status="unknown")
c      write(87,197)(splag(j),j=1,10)
c      close(87)
 197  format(10(2x,e14.7))
         if ( myrank.eq.master) print *, sumtime
c
c.... ---------------------->  Post Processing  <----------------------


         if(iCForz.eq.1)then
           if(myrank.eq.master)then
             write(*,*)'Closing lift force files' !for matt CF algorithm
             close(741)
             close(742)
             close(743)
             close(744)
             close(745)
             close(746)
             close(747)
             close(748)
           end if
         deallocate(xcf)
         deallocate(ycf)
         deallocate(zcf)

         end if

	if (myrank.eq.master.and.exts)
     1   close(10101)		! Close the binary output file

c
c=====================================================================
c.... print out the last step and save marker field
c
      if ((irs .ge. 1) .and. ((mod(lstep, ntout) .ne. 0) .or.
     &     (nstp .eq. 0) .or. (stopjob.ne.0) )) then
         if(
     &              ((irscale.ge.0).or.(itwmod.gt.0) .or.
     &              ((nsonmax.eq.1).and.iLES.gt.0)))
     &              call rwvelb  ('out ',  velbar  ,ifail)

         IF(iBT.eq.1) call banmaCorrector(banma, ibreakupFlag, yold) 
         CALL restar ('out ',  yold  ,ac)
         IF(iBT.eq.1) yold(:,7) = yold(:,6)

!       JF: write out the d2wall info in the last restart files
         if (iRANS .lt. 0. or. (iDNS.ne.0.and.abs(itwmod).eq.1)
     &      .or. iBT .eq.1) then
            call write_field(myrank, 'a', 'dwal', 4, d2wal, 'd', nshg,
     &                    4, lstep)
         endif

c=====================================================================
c         if(iofieldv.ne.0)
c     &           call pstdrv (yold,      ac,         x,
c     &                        iBC,       BC,
c     &                        iper,      ilwork,     shp,
c     &                        shgl,      shpb,       shglb,
c     &                        ifath,     velbar,     nsons )
         if (i_gradphi.eq.1) then
	    idflx = 0 
            if(idiff >= 1 )  idflx= (nflow-1) * nsd
            if (isurf == 1) idflx=nflow*nsd
	    call getgradphi(x, y, shp,       shgl,
     &                      shpb,          shglb, gradphi)
            gradphimag(:) = ( gradphi(:,1)**2 + 
     &                        gradphi(:,2)**2 + 
     &                        gradphi(:,3)**2 )**0.5
            maxgradphi = maxval(gradphimag)
	    write(*,1000) maxgradphi, myrank
            call MPI_BARRIER(MPI_COMM_WORLD,ierr)
	    call write_gradphi(myrank, lstep, nshg, 3, 
     &                           x, y, gradphi, gradphimag)
	 endif
c
c... Write out Redistancing primary vertice information
c
         if (i_primvert.eq.1) then
             primvertval(:,2) = primvert(:) * 1.0
             call write_primvert(myrank, lstep, nshg, 2,
     &                           primvertval)
         endif
         if(ideformwall.eq.1) 
     &        call write_displ(myrank, lstep, nshg, 3, u ) 
      endif


c
c Free memory for freezing LS Scalar 2
c
            if (iSolvLSSclr2.eq.2) then
              deallocate (iBCredist)
              deallocate (BCredist)
            endif

         lesId   = numeqns(1)
!      IF (memLSflag .NE. 1) THEN
      svLSflag = 1
      IF (svLSflag .NE. 1) THEN ! Arsen
!!!         call saveLesRestart( lesId,  aperm , nshg, myrank, lstep, ! ASSUME NEVER use LIBLES
!!!     &                        nPermDims )
      END IF


      if(ierrcalc.eq.1) then
c
c.....smooth the error indicators
c
        do i=1,ierrsmooth
            call errsmooth( rerr, x, iper, ilwork, shp, shgl, iBC )
        end do
c
c.... open the output file
c
         iqoldsiz=nshg*ndof*2
c           call write_error(myrank, lstep, nshg, numerr, rerr ) 
	call write_field(myrank, 'a', 'errors', 6, 
     &  rerr, 'd', nshg,numerr,lstep)
      endif

      if(ioybar.eq.1) then

         call write_field(myrank,'a','ybar',4,
     &                    ybar,'d',nshg,ndof,lstep)


c         write (fmt2,"('(''restart.'',i',i1,',1x)')") itmp
c         write (fname2,fmt2) lstep
c
c         fname2 = trim(fname2) // cname(myrank+1)
c
c.... open  files
c
c         call openfile(  fname2,  'append?', irstin )
c
c         fnamer2 = 'ybar'
c         isize = nshg*5
c         nitems = 3
c         iarray(1) = nshg
c         iarray(2) = 5
c         iarray(3) = lstep
c         call writeheader(irstin, fnamer2,iarray, nitems, isize,
c     &        'double', iotype )
c
c         nitems = nshg*5
c         call writedatablock(irstin, fnamer2,ybar, nitems,
c     &        'double', iotype)
c
c         call closefile( irstin, "append" )

      endif

      if ( ( ihessian .eq. 1 ) .and. ( numpe < 2 )  )then
          uhess = zero
          gradu = zero
          tf = zero

          do ku=1,nshg
c           tf(ku,1) = x(ku,1)**2+2*x(ku,1)*x(ku,2)
            tf(ku,1) = x(ku,1)**3
          end do

          call hessian( yold, x,     shp,  shgl,   iBC, 
     &                  shpb, shglb, iper, ilwork, uhess, gradu )

          call write_hessian( uhess, gradu, nshg )
      endif
c
c write phasta status file
c
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      if (myrank .eq. master) then
         open (unit=istat,  file=fstat,  status='unknown')
         write(istat,*) "done"
      endif
      close(istat)
c
c.... close history and aerodynamic forces files
c
      if (myrank .eq. master) then
         close (ihist)
         close (iforce)
         if (iLSet .eq. 2) then
           close (ivhist)
         endif
         if(exts) then
            do jj=1,ntspts
               close(1000+jj)
            enddo
         endif
      endif
      do isrf = 0,MAXSURF
         if ( nsrflist(isrf).ne.0 ) then
            iunit=60+isrf
            close(iunit)
         endif
      enddo
 5    format(1X,F15.10,3X,F15.10,3X,F15.10,3X,F15.10)
 444  format(6(2x,e14.7))
c
c.... end
c
      if(nsolflow.eq.1) then
         deallocate (lhsK)
         deallocate (lhsP)
       IF (svLSFlag .NE. 1) THEN ! Arsen
         deallocate (aperm)
         deallocate (atemp)
       ENDIF
      endif
      if(nsclrsol.gt.0) then
         deallocate (lhsS)
       IF (svLSFlag .NE. 1) THEN ! Arsen
         deallocate (apermS)
         deallocate (atempS)
       ENDIF
      endif
      
      if(iabc==1) deallocate(acs)
      
      if (allocated(ltg)) then
        deallocate(ltg)
      end if

      return

      ! ################################################################
      ! Arsen, add the nuemann BC
      CONTAINS

      SUBROUTINE AddNeumannBCTosvLS(srfID, faIn)

      INTEGER, INTENT(IN) :: srfID, faIn
      TYPE(svLS_lhsType) svLS_lhs
      INTEGER facenNo, i, j

      facenNo = 0
      DO i = 1, nshg
         IF (srfID .EQ. ndsurf(i)) THEN
            facenNo = facenNo + 1
         END IF
      END DO
      IF (ALLOCATED(gNodes)) DEALLOCATE(gNodes, sV)
      ALLOCATE(gNodes(facenNo), sV(nsd,facenNo))
      sV = 0D0
      j = 0
      DO i = 1, nshg
         IF (srfID .EQ. ndsurf(i)) THEN
            j = j + 1
            gNodes(j) = i
            sV(:,j) = NABI(i,1:3)
         END IF
      END DO
!      CALL memLS_BC_CREATE(memLS_lhs, faIn, facenNo,
!     2   nsd, BC_TYPE_Neu, gNodes, sV)

      RETURN
      END SUBROUTINE AddNeumannBCTosvLS

      ! ################################################################

!!      CONTAINS

!!      SUBROUTINE AddNeumannBCTomemLS(srfID, faIn)

!!      INTEGER, INTENT(IN) :: srfID, faIn
!!      TYPE(memLS_lhsType) memLS_lhs 
!!      INTEGER facenNo, i, j

!!      facenNo = 0
!!      DO i = 1, nshg
!!         IF (srfID .EQ. ndsurf(i)) THEN
!!            facenNo = facenNo + 1
!!         END IF
!!      END DO
!!      IF (ALLOCATED(gNodes)) DEALLOCATE(gNodes, sV)
!!      ALLOCATE(gNodes(facenNo), sV(nsd,facenNo))
!!	  gNodes = 0
!!      sV = 0D0
!!      j = 0
!!      DO i = 1, nshg
!!         IF (srfID .EQ. ndsurf(i)) THEN
!!            j = j + 1
!!            gNodes(j) = i
!!            sV(:,j) = NABI(i,1:3)
!!         END IF
!!      END DO
!      CALL memLS_BC_CREATE(memLS_lhs, faIn, facenNo,
!     2   nsd, BC_TYPE_Neu, gNodes, sV)

!!      RETURN
!!      END SUBROUTINE AddNeumannBCTomemLS

      end
      
      subroutine lesmodels(y,     ac,        shgl,      shp, 
     &                     iper,  ilwork,    rowp,      colm,    
     &                     nsons, ifath,     x,   
     &                     iBC,   BC)
      
      include "common.h"

      real*8    y(nshg,ndof),              ac(nshg,ndof),           
     &            x(numnp,nsd),
     &            BC(nshg,ndofBC)
      real*8    shp(MAXTOP,maxsh,MAXQPT),  
     &            shgl(MAXTOP,nsd,maxsh,MAXQPT)

c
      integer   rowp(nshg,nnz),         colm(nshg+1),
     &            iBC(nshg),
     &            ilwork(nlwork),
     &            iper(nshg)
      dimension ifath(numnp),    nsons(nfath)

      real*8, allocatable, dimension(:) :: fwr2,fwr3,fwr4
      real*8, allocatable, dimension(:) :: stabdis,cdelsq1
      real*8, allocatable, dimension(:,:) :: xavegt, xavegt2,xavegt3

      if( (iLES.gt.1) )   then ! Allocate Stuff for advanced LES models
         allocate (fwr2(nshg))
         allocate (fwr3(nshg))
         allocate (fwr4(nshg))
         allocate (xavegt(nfath,12))
         allocate (xavegt2(nfath,12))
         allocate (xavegt3(nfath,12))
         allocate (stabdis(nfath))
		 fwr2 = 0
		 fwr3 = 0
		 fwr4 = 0
		 xavegt = 0
		 xavegt2 = 0
		 xavegt3 = 0
		 stabdis = 0
      endif

c.... get dynamic model coefficient
c
      ilesmod=iLES/10  
c
c digit bit set filter rule, 10 bit set model
c
      if (ilesmod.eq.0) then    ! 0 < iLES< 10 => dyn. model calculated
                                ! at nodes based on discrete filtering


         if(isubmod.eq.2) then
            call SUPGdis(y,      ac,        shgl,      shp, 
     &                   iper,   ilwork,    
     &                   nsons,  ifath,     x,   
     &                   iBC,    BC, stabdis, xavegt3)
         endif

         if( ((isubmod.eq.0).or.(isubmod.eq.2)))then ! If no
                                                     ! sub-model
                                                     ! or SUPG
                                                     ! model wanted

            if(i2filt.eq.0)then ! If simple filter
              
               if(modlstats .eq. 0) then ! If no model stats wanted
                  call getdmc (y,       shgl,      shp, 
     &                         iper,       ilwork,    nsons,
     &                         ifath,      x)
               else             ! else get model stats 
                  call stdfdmc (y,       shgl,      shp, 
     &                          iper,       ilwork,    nsons,
     &                          ifath,      x)
               endif            ! end of stats if statement  

            else                ! else if twice filtering

               call widefdmc(y,       shgl,      shp, 
     &                       iper,       ilwork,    nsons,
     &                       ifath,      x)

               
            endif               ! end of simple filter if statement

         endif                  ! end of SUPG or no sub-model if statement


         if( (isubmod.eq.1) ) then ! If DFWR sub-model wanted
            call cdelBHsq (y,       shgl,      shp, 
     &                     iper,       ilwork,    nsons,
     &                     ifath,      x,         cdelsq1)
            call FiltRat (y,       shgl,      shp, 
     &                    iper,       ilwork,    nsons,
     &                    ifath,      x,         cdelsq1,
     &                    fwr4,       fwr3)

            
            if (i2filt.eq.0) then ! If simple filter wanted
               call DFWRsfdmc(y,       shgl,      shp, 
     &                        iper,       ilwork,    nsons,
     &                        ifath,      x,         fwr2, fwr3) 
            else                ! else if twice filtering wanted 
               call DFWRwfdmc(y,       shgl,      shp, 
     &                        iper,       ilwork,    nsons,
     &                        ifath,      x,         fwr4, fwr4) 
            endif               ! end of simple filter if statement
             
         endif                  ! end of DFWR sub-model if statement

         if( (isubmod.eq.2) )then ! If SUPG sub-model wanted
            call dmcSUPG (y,           ac,         shgl,      
     &                    shp,         iper,       ilwork,    
     &                    nsons,       ifath,      x,
     &                    iBC,    BC,  rowp,       colm,
     &                    xavegt2,    stabdis)
         endif

         if(idis.eq.1)then      ! If SUPG/Model dissipation wanted
            call ediss (y,        ac,      shgl,      
     &                  shp,      iper,       ilwork,    
     &                  nsons,    ifath,      x,
     &                  iBC,      BC,  xavegt)
         endif

      endif                     ! end of ilesmod
      
      if (ilesmod .eq. 1) then  ! 10 < iLES < 20 => dynamic-mixed
                                ! at nodes based on discrete filtering
         call bardmc (y,       shgl,      shp, 
     &                iper,    ilwork,    
     &                nsons,   ifath,     x) 
      endif
      
      if (ilesmod .eq. 2) then  ! 20 < iLES < 30 => dynamic at quad
                                ! pts based on lumped projection filt. 

         if(isubmod.eq.0)then
            call projdmc (y,       shgl,      shp, 
     &                    iper,       ilwork,    x) 
         else
            call cpjdmcnoi (y,      shgl,      shp, 
     &                      iper,   ilwork,       x,
     &                      rowp,   colm, 
     &                      iBC,    BC)
         endif

      endif

      if( (iLES.gt.1) )   then ! Deallocate Stuff for advanced LES models
         deallocate (fwr2)
         deallocate (fwr3)
         deallocate (fwr4)
         deallocate (xavegt)
         deallocate (xavegt2)
         deallocate (xavegt3)
         deallocate (stabdis)
      endif
      return
      end
