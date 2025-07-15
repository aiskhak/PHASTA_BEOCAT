c
c  Copyright (c) 2000-2007, Stanford University,
c     Rensselaer Polytechnic Institute, Kenneth E. Jansen,

      subroutine SolFlow(y,          ac,          u,
     &                   banma,      
     &                   yold,       acold,       uold,
     &                   x,          iBC,
     &                   BC,         res,             
     &                   nPermDims,  nTmpDims,    aperm,
     &                   atemp,      iper,       
     &                   ilwork,     shp,         shgl, 
     &                   shpb,       shglb,       rowp,     
     &                   colm,       lhsK,        lhsP, 
     &                   solinc,     rerr,        sumtime,
!     &                   memLS_lhs,  memLS_ls,    memLS_nFaces,
     &                   svLS_lhs, svLS_ls, svLS_nFaces, ! Arsen
     &                   elemvol_global,
     &                   avgxcoordf, avgycoordf,  avgzcoordf)

c
c!----------------------------------------------------------------------
c
c! This is the 2nd interface routine to the Farzin''s linear equation 
c! solver library that uses the CGP and GMRES methods.
c
c! input:
c!  y      (nshg,ndof)           : Y-variables at n+alpha_f
c!  ac     (nshg,ndof)           : Primvar. accel. variable n+alpha_m
c!  yold   (nshg,ndof)           : Y-variables at beginning of step
c!  acold   (nshg,ndof)          : Primvar. accel. at beginning of step
c!  x      (numnp,nsd)            : node coordinates
c!  iBC    (nshg)                : BC codes
c!  BC     (nshg,ndofBC)         : BC constraint parameters
c!  iper   (nshg)                : periodic nodal information
c
c! output:
c!  res    (nshg,nflow)           : preconditioned residual
c!  y      (nshg,ndof)           : Y-variables at n+alpha_f
c!  ac     (nshg,ndof)           : Primvar. accel. variable n+alpha_m
c
c
c! The followings are preliminary steps required to use Farzin''s
c! solver library.  New way of writing has to be used such as
c
c!          |  K     G | | du |    | Rmom  |
c!          |          | |    | =  |       |
c!          | G^t    C | | dp |    | Rcon  |
c
c!          |     E    | | dT | =  | Rtemp |
c
c!     where
c
c!      xKebe : K_ab = dRmom_a/du_b    xTe : E_ab = dRtemp_a/dT_b 
c
c!              G_ab = dRmom_a/dp_b
c!      xGoC  :
c!              C_ab = dRcon_a/dp_b       
c
c!              resf = Rmon Rcon       rest = Rtemp
c
c  
c! Zdenek Johan,  Winter 1991.  (Fortran 90)
c! Juin Kim, Summer 1998. (Incompressible flow solver interface)
c! Alberto Figueroa.  CMM-FSI
c----------------------------------------------------------------------
c
      use pointer_data
c!#ifdef AMG      
c!      use ramg_data
c!#endif     
        
      include "common.h"
      include "mpif.h"
      include "auxmpi.h"
!      INCLUDE "memLS.h"
      include "svLS.h"  ! Arsen
c     
!      TYPE(memLS_lhsType) memLS_lhs
!      TYPE(memLS_lsType) memLS_ls
      TYPE(svLS_lhsType) svLS_lhs   ! Arsen
      TYPE(svLS_lsType) svLS_ls     ! Arsen

      real*8    y(nshg,ndof),             ac(nshg,ndof),
     &          yold(nshg,ndof),          acold(nshg,ndof),
     &          u(nshg,nsd),              uold(nshg,nsd),
     &          x(numnp,nsd),             BC(nshg,ndofBC),
     &          res(nshg,nflow),          tmpres(nshg,nflow),
     &          flowDiag(nshg,4),         
     &          aperm(nshg,nPermDims),    atemp(nshg,nTmpDims),
     &          sclrDiag(nshg,1),         
     &          lhsK(9,nnz_tot),	  lhsP(4,nnz_tot)
c
      dimension banma(nshg,1)
c
      real*8    shp(MAXTOP,maxsh,MAXQPT),  
     &          shgl(MAXTOP,nsd,maxsh,MAXQPT), 
     &          shpb(MAXTOP,maxsh,MAXQPT),
     &          shglb(MAXTOP,nsd,maxsh,MAXQPT) 
c
      integer   usr(100),                 eqnType,temp,
     &          rowp(nshg*nnz),           colm(nshg+1),
     &          iBC(nshg),                ilwork(nlwork),
     &          iper(nshg)
c
      real*8    yAlpha(nshg,ndof),        acAlpha(nshg,ndof),
     &          uAlpha(nshg,nsd),
     &          lesP(nshg,4),             lesQ(nshg,4),
     &          solinc(nshg,ndof),        CGsol(nshg)

! Arsen CHANGE
      real*8     tcorecp(2)
! Arsen CHANGE END
      
      real*8    rerr(nshg,numerr),            rtmp(nshg,4),rtemp
      
      real*8    msum(4),mval(4),cpusec(10)

! BEGIN Arsen ----------------------------------------------------------------
      integer dof, svLS_nfaces, i, j, k, lesId
      integer, allocatable :: incL(:)
      real*8, allocatable :: faceRes(:), Res4(:,:), Val4(:,:)
	  integer svLSFlag
! END Arsen ------------------------------------------------------------------

      REAL*8 sumtime, timekeeper
!      INTEGER dof, memLS_nFaces, i, j, k, l
!      INTEGER, ALLOCATABLE :: incL(:)
!      REAL*8, ALLOCATABLE :: faceRes(:), Res4(:,:), Val4(:,:)

      integer sparseloc 

c!.... Matt Talley's Bubble Coal Control
        real*8, dimension (100) :: avgxcoordf, 
     &avgycoordf, avgzcoordf ! ??? hardcoded coalest

!      if (.not.allocated(avgxcoordf)) allocate (avgxcoordf(coalest))
!	  if (.not.allocated(avgycoordf)) allocate (avgycoordf(coalest))
!	  if (.not.allocated(avgzcoordf)) allocate (avgzcoordf(coalest))
c     
c.... *******************>> Element Data Formation <<******************
c
c
c.... set the parameters for flux and surface tension calculations
c
c
        temp = npro
        

        idflx = 0 
        if(idiff >= 1 )  idflx= (nflow-1) * nsd  
        if (isurf == 1)  idflx=nflow*nsd
	if (idiff >= 1)  idflx = idflx + 10     ! used in the dissipation timeseries computation  
c.... compute solution at n+alpha
c
c           if (myrank .eq. master) write(*,*) 'idflx = ', idflx

      call itrYAlpha( uold,    yold,    acold,
     &                u,       y,       ac,  
     &                uAlpha,  yAlpha,  acAlpha)
c           if (myrank .eq. master) write(*,*) 'itrYAlpha is done'

c
c.... form the LHS matrices, the residual vector (at alpha)
c
      call ElmGMR (uAlpha,    yAlpha,     acAlpha,   
     &             x,         banma,      
     &             shp,       shgl,       iBC,       
     &             BC,        shpb,       shglb,
     &             res,       iper,       ilwork,   
     &             rowp,      colm,       lhsK,      
     &             lhsP,      rerr,       elemvol_global,
     &             avgxcoordf,avgycoordf, avgzcoordf)

            tmpres(:,:) = res(:,:)
            iblk = 1

      timekeeper = cput()

! BEGIN Arsen ---------------------------------------------------------
!###################################################################
!  Calling svLs to solve
! The following block of code has been implemented just as the svLS
! is implemented in the developBoiling_IB verison for the code
! this will likely lead to mistakes if unedited with syncIO format of
! the current code -- MB, 01 Aug 2024

      svLSFlag = 1
      IF (svLSFlag.eq.1) THEN
      
      ALLOCATE(faceRes(svLS_nFaces), incL(svLS_nFaces))  ! Arsen
      faceRes(:) = 0.d0 ! Neumann will not work ???
      incL(:) = 1
      dof = 4
      if (.not.allocated(Res4)) then
            allocate(Res4(dof,nshg), Val4(dof*dof, nnz_tot))
      end if

      do i=1, nshg
            Res4(1:dof,i) = res(i,1:dof)
      end do

      do i=1, nnz_tot
            Val4(1:3,i) = lhsK(1:3,i)
            Val4(5:7,i) = lhsK(4:6,i)
            Val4(9:11,i) = lhsK(7:9,i)
            Val4(13:15,i) = lhsP(1:3,i)
            Val4(16,i) = lhsP(4,i)
      end do

      do i=1, nshg
            do j=colm(i), colm(i+1) - 1
                  k =rowp(j)
                  do l=colm(k), colm(k+1) - 1
                        if (rowp(l).eq.i) then
                              Val4(4:12:4,l) = -lhsP(1:3,j) !val4(4:12:4) are coupling terms for u and p
                              exit
                        end if
                  end do
            end do
      end do

      do j=1, 0 !nshg      ! 0 --> nshg ! loop over number of elements on proc
            i = iper(j) ! iper(j) = j if not slave, if slave it does not
                  !i will only point to nodes that are not a slave
!            if (myrank.eq.master) write(*,*) i

            do k1=colm(i), colm(i+1) - 1  ! colm(i) gives the index of the vals, that start at row i
                  if (rowp(k1).eq.i) then ! if the row index matches the column index, i,e,m the value is on a digonal
                        i1=k1             ! if the diagonal element is tsored, its index in val is stored as i1
                        exit
                  end if
            end do
            
            ! this does the same but with j1 instead
            do k1=colm(j), colm(j+1) - 1
                  if (rowp(k1).eq.j) then
                        j1=k1
                        exit
                  end if 
            end do

! in progress   
            if (i.ne.j) then        ! periodic condition on node
                                    ! if i does not equal j, then the index of the element j, does not match
                                    ! its slave status j
                  write(*,*) 'i does not equal j'     ! identified the values, j must indicate the slave (?)
                                                      ! i indicates the master (?)

                  do k1=colm(j), colm(j+1) - 1        ! colm(j) is the starting index of the
                                                      ! nonzero elements in row j in the val array
                        !write(*,*)'k1 = ', k1, ': rowp(k1)=', rowp(k1)
                        if (rowp(k1).eq.j) then
                              !write(*,*) 'unmodified val4 for node', i, 'is', val4(:,k1)
                              Val4(1:16, k1) = zero
                              Val4(1:16:5, k1) = one
!                              Res4(1:4,k1)=zero
                              !val4(1:16:5,k1) = -one
                              
                              write(*,*), 'set pos diagonal'
      write(*,*), 'modified val4 for node ',j,'is',Val4(:,k1)
                        ! above will print if the node belongs to a boundary element on the master from what i can see,
                        ! however perhaps i want the other way around?
                        else
                              Val4(1:16, k1) = zero
                        end if
                        if (rowp(k1).eq.i) then
                              !write(*,*) 'element is a slave'
                              Val4(1:16:5, k1) = -one
                              write(*,*) ' component (j,i) is modified'
                       end if
                  end do        

            end if ! end periodic node condition


      
 101   FORMAT(1x, 'VAL4(1,6,11,16) @', 2I7, 4F12.7)
 102   FORMAT(1x, 'Res4(1:4) @', I7, 4F12.7)
 103   FORMAT(1x, 'val4 row # ', I7, 104F12.7)
      
      end do            ! nshg loop
       !write(*,*) 'ncorp',ncorp
      !if (myrank.eq.master) write(*,*) 'calling svLS_SOLVE'
      !if (myrank.eq.master) write(*,*) 'dof = ', dof
      !if (myrank.eq.master) write(*,*) 'Res4 = ', Res4
      !if (myrank.eq.master) write(*,*) 'Val4 = ', Val4
      !call svLS_SOLVE(svLS_lhs,svLS_ls,dof,Res4,Val4)
      call svLS_SOLVE(svLS_lhs,svLS_ls,dof,Res4,Val4,incL,faceRes) ! Arsen
      !if(myrank.eq.master)write(*,*)'svLS_SOLVE done'

      DO i=1, nshg
            solinc(i,1:dof) = Res4(1:dof,i)
      END DO

      deallocate(faceRes, incL)
      ELSE  ! MB, use lesSolve instead
           

! ###############################################################################
! END Arsen ------------------------------------------------------------------------

!      IF (memLSFlag .EQ. 1) THEN
!####################################################################
!     Here calling memLS

!      ALLOCATE(faceRes(memLS_nFaces), incL(memLS_nFaces))
!      CALL AddElmpvsQFormemLS(faceRes, memLS_nFaces)

!!      incL = 1
!!      dof = 4
!!      IF (.NOT.ALLOCATED(Res4)) THEN
!!        ALLOCATE (Res4(dof,nshg), Val4(dof*dof,nnz_tot))
!!      END IF

!!      DO i=1, nshg
!!         Res4(1:dof,i) = res(i,1:dof)
!!      END DO

!!      DO i=1, nnz_tot
!!         Val4(1:3,i)   = lhsK(1:3,i)
!!         Val4(5:7,i)   = lhsK(4:6,i)
!!         Val4(9:11,i)  = lhsK(7:9,i)
!!         Val4(13:15,i) = lhsP(1:3,i)
!!         Val4(16,i)    = lhsP(4,i)
!!      END DO

      !Val4(4:12:4,:) = -lhsP(1:3,:)^t
!!      DO i=1, nshg
!!        Do j=colm(i), colm(i+1) - 1  
!!            k = rowp(j)
!!            DO l=colm(k), colm(k+1) - 1
!!               IF (rowp(l) .EQ. i) THEN
!!                  Val4(4:12:4,l) = -lhsP(1:3,j)
!!                  EXIT
!!               END IF
!!            END DO
!!         END DO
!!      END DO

! Check here what is the diagonal values for the j = iper(i) nodes:

!!      DO j=1, 0  !nshg   ! Loop over "global number of shape functions"   LOOP is TURNED OFF ! IGOR, DECEMBER 2012 
!!        i = iper(j)
! How to convert j from nshg to nnz_total ??? See above !!
! Define i1 and j1 here !!
!            i1 = sparseloc( rowp(colm(i)), colm(i+1)-colm(i), i )
!     &       + colm(i)-1
!            j1 = sparseloc( rowp(colm(j)), colm(j+1)-colm(j), j )
!     &       + colm(j)-1

! Above way:
!!        Do k1=colm(i), colm(i+1) - 1
!          write(*,*) 'k1 = ', k1, '; rowp(k)=', rowp(k1), j
!!             IF (rowp(k1) .EQ. i) THEN
!!                  i1 = k1
!!                  EXIT
!!            END if
!!         END DO
!!        Do k1=colm(j), colm(j+1) - 1
!!               IF (rowp(k1) .EQ. j) THEN
!!                  j1 = k1
!!                  EXIT
!!               END IF
!!         END DO

! Print the diagonal element value for both i and j:
!!        if (i.ne.j) then   ! if periodicity is on this node
!          write(*,*) 'i,j: ', i,j
!          write(*,101) rowp(i1), i1, Val4(1:16:5,i1)
!          write(*,102) i, Res4(1:dof,i)
!          write(*,101) rowp(j1), j1, Val4(1:16:5,j1)
!          write(*,102) j, Res4(1:dof,j)

! Next steps:
! a) assign the row to zero !
! b) r.h.s to zero !
! c) A(j,i) = -1.0
! d) A(j,j) = 1.0 

! Make sure this is done correctly in the sparse matrix format !

! Make the whole row "j" equal to zero (the solver should complain at this point)
!  .. and 1.0 for j column


!!        Do k1=colm(j), colm(j+1) - 1
!!          write(*,*) 'k1 = ', k1, '; rowp(k)=', rowp(k1)
!!               IF (rowp(k1) .EQ. j) THEN
!!                Val4(1:16, k1) = zero
!!                Val4(1:16:5, k1) = one
!!               else
!!                Val4(1:16, k1) = zero
!!               END IF
!!               IF (rowp(k1) .EQ. i) THEN
!!                Val4(1:16:5, k1) = -one
!!                write(*,*) ' component (j,i) is modified '
!!               END IF
!!        END DO


         
!          write(*,*) 'After modification : '
!          write(*,101) rowp(j1), j1, Val4(1:16:5,j1)
!          write(*,103) j, Val4(1,colm(j):colm(j+1) - 1)
!          write(*,102) i, Res4(1:dof,i)


!!        end if   ! Periodic node condition

!! 101   FORMAT(1x, 'VAL4(1,6,11,16) @', 2I7, 4F12.7)
!! 102   FORMAT(1x, 'Res4(1:4) @', I7, 4F12.7)
!! 103   FORMAT(1x, 'val4 row # ', I7, 104F12.7)

!!      END DO            ! nshg loop



!      CALL memLS_SOLVE(memLS_lhs, memLS_ls, dof, Res4, Val4, incL,
!     2   faceRes)

!      if(myrank.eq.master)write(*,*)'memLS_SOLVE called'
!!      CALL memLS_SOLVE(memLS_lhs, memLS_ls, dof, Res4, Val4)
!      if(myrank.eq.master)write(*,*)'memLS_SOLVE done'

!!      DO i=1, nshg
!!         solinc(i,1:dof) = Res4(1:dof,i)
!!      END DO

c####################################################################
!!      ELSE

c.... lesSolve : main matrix solver
c
      lesId   = numeqns(1)
      eqnType = 1
c
c.... setup the linear algebra solver
c
      rtmp = res(:,1:4)
!!!      call usrNew ( usr,        eqnType,          aperm, ! assume never use LIBLES
!!!     &              atemp,      rtmp,             solinc,          
!!!     &              flowDiag,   sclrDiag,         lesP,   
!!!     &              lesQ,       iBC,              BC,
!!!     &              iper,       ilwork,           numpe,
!!!     &              nshg,       nshl,             nPermDims,  
!!!     &              nTmpDims,   rowp,             colm,     
!!!     &              lhsK,       lhsP,             rdtmp,      
!!!     &              nnz_tot,    CGsol )
c
c.... solve linear system
c

!!!      call myfLesSolve ( lesId, usr ) ! assume never use LIBLES
      
      call getSol ( usr, solinc )

      if (numpe > 1) then
         call commu ( solinc, ilwork, nflow, 'out')
      endif

        END IF   ! memLS / lesLIB choice condition

! Arsen CHANGE 
      tlescp2 = TMRC()
      impistat=0
      impistat2=0

      tcorecp(1) = tcorecp(1) + telmcp2-telmcp1 ! elem. formation
      tcorecp(2) = tcorecp(2) + tlescp2-tlescp1 ! linear alg. solution
! Arsen CHANGE END
      
      call rstatic (res, y, solinc) ! output flow stats
      sumtime = sumtime + cput() - timekeeper

c     
c.... end
c     
      return
      end

      subroutine SolSclr(y,          ac,         u,
     &                   yold,       acold,      uold,
     &                   x,          iBC,
     &                   BC,         nPermDimsS,  nTmpDimsS,  
     &                   apermS,     atempS,     iper,       
     &                   ilwork,     shp,        shgl, 
     &                   shpb,       shglb,      rowp,     
     &                   colm,       lhsS,       solinc,
     &                   cfl, 
     &                   svLS_lhs_sc, svLS_sc, svLS_nFaces)
!     &                   memLS_lhs_sc,  memLS_sc,   memLS_nFaces)
c
c----------------------------------------------------------------------
c
c This is the 2nd interface routine to the Farzin''s linear equation 
c solver library.
c
c input:
c  y      (nshg,ndof)           : Y-variables at n+alpha_f
c  ac     (nshg,ndof)           : Primvar. accel. variable n+alpha_m
c  yold   (nshg,ndof)           : Y-variables at beginning of step
c  acold  (nshg,ndof)           : Primvar. accel. variable at begng step
c  x      (numnp,nsd)            : node coordinates
c  iBC    (nshg)                : BC codes
c  BC     (nshg,ndofBC)         : BC constraint parameters
c  iper   (nshg)                : periodic nodal information
c
c output:
c  y      (nshg,ndof)           : Y-variables at n+alpha_f
c  ac     (nshg,ndof)           : Primvar. accel. variable n+alpha_m
c
c
c The followings are preliminary steps required to use Farzin''s
c solver library.  New way of writing has to be used such as
c
c          |     E    | | dS | =  | RScal |
c
c----------------------------------------------------------------------
c
      use pointer_data
        
      include "common.h"
      include "mpif.h"
      include "auxmpi.h"
!      INCLUDE "memLS.h"
      INCLUDE "svLS.h"
c
!      TYPE(memLS_lhsType) memLS_lhs_sc
!      TYPE(memLS_lsType) memLS_sc
      TYPE(svLS_lhsType) svLS_lhs_sc
      TYPE(svLS_lsType) svLS_sc

c     
      real*8    y(nshg,ndof),             ac(nshg,ndof),
     &          yold(nshg,ndof),          acold(nshg,ndof),
     &          u(nshg,nsd),              uold(nshg,nsd),
     &          x(numnp,nsd),             BC(nshg,ndofBC),
     &          res(nshg,1),
     &          flowDiag(nshg,4),
     &          sclrDiag(nshg,1),         lhsS(nnz_tot),
     &          apermS(nshg,nPermDimsS),  atempS(nshg,nTmpDimsS)

c
      real*8    shp(MAXTOP,maxsh,MAXQPT),  
     &          shgl(MAXTOP,nsd,maxsh,MAXQPT), 
     &          shpb(MAXTOP,maxsh,MAXQPT),
     &          shglb(MAXTOP,nsd,maxsh,MAXQPT) 
c
      integer   usr(100),                 eqnType,
     &          rowp(nshg*nnz),           colm(nshg+1),
     &          iBC(nshg),                ilwork(nlwork),
     &          iper(nshg)
c
      real*8    yAlpha(nshg,ndof),        acAlpha(nshg,ndof),
     &          uAlpha(nshg,nsd),
     &          lesP(nshg,1),             lesQ(nshg,1),
     &          solinc(nshg,1),           cfl(nshg)
      
      REAL*8 sumtime, timekeeper
!      INTEGER dof, memLS_nFaces, i, j, k, l
      INTEGER dof, svLS_nFaces, i, j, k, l
      INTEGER, ALLOCATABLE :: incL(:)
      REAL*8, ALLOCATABLE :: faceRes(:), Res4(:,:), Val4(:,:)
      INTEGER svLSFlag

c     
c.... *******************>> Element Data Formation <<******************
c
c.... compute solution at n+alpha
c
      call itrYAlpha( uold,    yold,    acold, 
     &                u,       y,       ac,  
     &                uAlpha,  yAlpha,  acAlpha)
c
c.... form the LHS matrices, the residual vector (at alpha)
c
      call ElmGMRSclr(yAlpha,acAlpha,    x,
     &             shp,       shgl,       iBC,       
     &             BC,        shpb,       shglb,
     &             res,       iper,       ilwork,   
     &             rowp,      colm,       lhsS,
     &             cfl )

! ################################################################
! Arsen, call the svLS solver
      svLSFlag = 1
      IF (svLSFlag .EQ. 1) THEN 
      ALLOCATE(faceRes(svLS_nFaces), incL(svLS_nFaces))  ! Arsen
      faceRes(:) = 0.d0 ! Neumann will not work ???
      incL(:) = 1
            lesId   = numeqns(1+nsolt+isclr)
      !        if (myrank.eq.master) write(*,*) 'lesID = ', lesId
            dof = 1   ! Should be equal to 1 ??? (was 4 for N.-S.)
            IF (.NOT.ALLOCATED(Res4)) THEN
                  ALLOCATE (Res4(dof,nshg), Val4(dof*dof,nnz_tot))
            END IF

            DO i=1, nshg
                  Res4(1:dof,i) = res(i,1)
            END DO

            DO i=1, nnz_tot
                  Val4(1,i)   = lhsS(i)
            END DO

      !      if (lesId.eq.2) then   ! Temperature
      !         CALL svLS_SOLVE(svLS_lhsT, svLS_sc, dof, Res4, Val4)
      !      else   ! Level set 
!                  CALL svLS_SOLVE(svLS_lhs_sc, svLS_sc, dof, Res4, Val4)
      !      end if
      CALL svLS_SOLVE(svLS_lhs_sc,svLS_sc,dof,Res4,Val4,incL,faceRes) ! Arsen
            DO i=1, nshg
                  solinc(i,1:dof) = Res4(1:dof,i)
            END DO
      deallocate (faceRes, incL)
      ELSE  ! not using svLS
     

! ################################################################


!!      IF (memLSFlag .EQ. 1) THEN
!####################################################################
!     Here calling memLS

!      ALLOCATE(faceRes(memLS_nFaces), incL(memLS_nFaces))
!      CALL AddElmpvsQFormemLS(faceRes, memLS_nFaces)

!!      lesId   = numeqns(1+nsolt+isclr)
!        if (myrank.eq.master) write(*,*) 'lesID = ', lesId
!!      incL = 1
!!      dof = 1   ! Should be equal to 1 ? (was 4 for N.-S.)
!!      IF (.NOT.ALLOCATED(Res4)) THEN
!!         ALLOCATE (Res4(dof,nshg), Val4(dof*dof,nnz_tot))
!!      END IF

!!      DO i=1, nshg
!!         Res4(1:dof,i) = res(i,1)
!!      END DO

!!      DO i=1, nnz_tot
!!         Val4(1,i)   = lhsS(i)
!!      END DO

!      if (lesId.eq.2) then   ! Temperature
!         CALL memLS_SOLVE(memLS_lhsT, memLS_sc, dof, Res4, Val4)
!      else   ! Level set 
!!         CALL memLS_SOLVE(memLS_lhs_sc, memLS_sc, dof, Res4, Val4)
!      end if

!!      DO i=1, nshg
!!         solinc(i,1:dof) = Res4(1:dof,i)
!!      END DO

!####################################################################
!!      ELSE

c
c.... lesSolve : main matrix solver
c
      lesId   = numeqns(1+nsolt+isclr)
      eqnType = 2
c
c.... setup the linear algebra solver
c
!!!      call usrNew ( usr,        eqnType,          apermS, ! assume never use LIBLES
!!!     &              atempS,     res,              solinc,          
!!!     &              flowDiag,   sclrDiag,         lesP,   
!!!     &              lesQ,       iBC,              BC,
!!!     &              iper,       ilwork,           numpe,
!!!     &              nshg,       nshl,             nPermDimsS,  
!!!     &              nTmpDimsS,  rowp,             colm,     
!!!     &              rlhsK,      rlhsP,            lhsS,      
!!!     &              nnz_tot )
c
c.... solve linear system
c 
!!!      call myfLesSolve ( lesId, usr ) ! assume never use LIBLES
      call getSol ( usr, solinc )

      if (numpe > 1) then
         call commu ( solinc, ilwork, 1, 'out')
      endif
      
      END IF   ! memLS / lesLIB choice condition

      nsolsc=5+isclr
      call rstaticSclr (res, y, solinc, nsolsc) ! output scalar stats
c     
c.... end
c     
      return
      end


!####################################################################
c     This function will return the current time in sec
      FUNCTION CPUT()

      IMPLICIT NONE

      INTEGER timeArray(8)

      Real*8 CPUT
!         END DO

      CALL DATE_AND_TIME (VALUES=timeArray)
      CPUT = ((timeArray(4)*24D0 !day
     2   + timeArray(5))*6D1     !hr
     3   + timeArray(6))*6D1     !min
     4   + timeArray(7)*1D0      !Sec
     5   + timeArray(8)*1D-3     !mSec

      RETURN
      END FUNCTION CPUT


