!----------------------------------------------------------------------
!
!       This file contains several subroutines which are dealing with
!       bubble tracking capability and advanced co-processing. The
!       subroutines are called in many subroutines, such as itrdrv.f,
!       elmgmr.f and so on. 
!
!       Jun Fang,                               Summer, 2015
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!       Name              Description
!       procs_dataset   : data collected from a single processor
!       unive_dataset   : data collected from all processors
!       bub_info        : bubble information array for each element
!       avg_info        : bubble information array for each bubble
!       one_procs       : dummy of procs_dataset
!       all_procs       : dummy of unive_dataset
!       phi_min         : the min level set value in a bubble
!       phi_tmp         : dummy of phi_min
!       eq_rad          : equivalent radius for deformed bubble
!       deform          : deformability factor
!       i_mrk           : index of colored bubbles
!       elemvol_global  : array contains element volumes
!       elemvol_local   : array contains element volumes for each block
!       Spnt_min        : local shear calculation point with min y
!       coordinate
!       Spnt_max        : local shear calculation point with max y
!       coordinate
!       Spnt            : two points array used in local shear
!                         calculation
!----------------------------------------------------------------------

        subroutine CountIniBub(banma)
!----------------------------------------------------------------------
!
!       This subroutine is used to count the total number of bubbles in
!       the initial profile. 
!
!----------------------------------------------------------------------
      USE, INTRINSIC :: ISO_C_BINDING !for calling C++ routines
        include "common.h"
        include "mpif.h"
        include "auxmpi.h"

        dimension banma(nshg,1)
        real*8    local_bmmax, global_bmmax


        if(myrank.eq.master) 
     &     write(*,'(1x,A,ES14.7)') 'epsilonBT = ', epsilonBT
        local_bmmax     = 0.0d0
        global_bmmax    = 0.0d0
        do i = 1, nshg
           if(local_bmmax .lt. banma(i,1))
     &        local_bmmax = banma(i,1)
        enddo
        call MPI_ALLREDUCE (local_bmmax, global_bmmax, 1,
     &       MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)

        i_num_bubbles = INT(global_bmmax)
        if(myrank.eq.master) write(*,*) 'Number of bubbles: ',
     &  i_num_bubbles

        end     !CountIniBub ends 
c======================================================================
c======================================================================
c======================================================================
        subroutine OpenBubFiles(i_num_bubbles, C_int_adjust)
!----------------------------------------------------------------------
!
!       Open files for saving the data about bubble properties or
!       behaviors (which is part of bubble tracking algorithm)
!
!----------------------------------------------------------------------

        integer:: i_num_bubbles
        INTEGER:: file_number, il, irstart, ierror
        integer:: nl_alpha, inistep
        integer,allocatable :: alpha01(:,:)
        real*8  C_int_adjust
        real*8, allocatable :: alpha02(:,:), bubmatx(:,:)
        CHARACTER(LEN=100) :: path1, path2, path3, fname_bub
        logical:: ObjExt1, ObjExt2


        C_int_adjust = 0.0d0
        write(*,*)
        inquire(directory='../bubStat/',exist=ObjExt1)
        if(ObjExt1.eqv..False.) then
           call system('mkdir ../bubStat')
           write(*,*) 'bubStat not found, plz create one manually!'
           write(*,*)
        else
           write(*,*) 'bubStat already exists'
           write(*,*)
        end if
        write(*,*) 'In OpenBubFiles, i_num_bubbles = ', i_num_bubbles
!       open file to record void fraction during simulations and trim
!       the existing files such that new data can fit without any
!       overlap with the previous data
        inquire(file='../bubStat/alpha.dat',exist=ObjExt2)
        IF(ObjExt2.eqv..False.) THEN   !no data file found, create them

        open(unit=740,file='../bubStat/alpha.dat',status="replace",
     &       action="write",iostat=ierror)
        if(ierror /= 0) STOP "Error creating file 740"
        close(740)

        if (i_num_bubbles.ne.0) then
           path1 = "../bubStat/bubble."
           do i=1,i_num_bubbles
              write(path2,'(i4.4)') i
              fname_bub = trim(path1)//trim(path2)//".dat"
              file_number = 950 + i
           OPEN(unit=file_number, file=fname_bub, status="replace",
     &          action="write", iostat=ierror)
           IF(ierror /= 0)
     &       STOP "Error creating files:/bubStat/bubble.*.dat"
           close(file_number)
           end do
        end if
!----------------------------------------------------------------------
        ELSE                            !old data files found, trim them
                                        !if necessary
        call CountAlphaLines(nl_alpha,inistep)
        write(*,*) 'Number of timesteps recorded = ', nl_alpha
        write(*,*) 'initial timestep in alpha.dat = ',inistep

        open(unit=72,file='numstart.dat',status='old')
        read(72,*) irstart
        close(72)

        write(*,*) 'Time step when processing bubble files is', 
     &             irstart+1

!----------------------------------
        if(irstart+1 .lt. (inistep+nl_alpha) 
     &                  .and. irstart+1.ne.inistep) then 
          write(*,*) 'Bubble files are being trimed!'

          allocate (alpha01(irstart+1-inistep,3))
          allocate (alpha02(irstart+1-inistep,3))
          allocate (bubmatx(irstart+1-inistep,21))

          open(unit=740,file='../bubStat/alpha.dat',status="old",
     &       action="read",iostat=ierror)
          if(ierror /= 0) STOP "Error reading file 740"
          do il = 1, (irstart+1-inistep)
             read(740,813) alpha01(il,1:3), alpha02(il,1:3)
          enddo
          close(740)
          open(unit=740,file='../bubStat/alpha.dat',status="replace",
     &       action="write",iostat=ierror)
          if(ierror /= 0) STOP "Error creating file 740"
          do il = 1, (irstart+1-inistep)
             write(740,813) alpha01(il,1:3), alpha02(il,1:3)
!       Use the C_int_adjust from previous simulation that can better
!       preserve void fraction 
             if(il.eq.(irstart+1-inistep)) C_int_adjust=alpha02(il,1)
          enddo
          close(740)

          deallocate(alpha01)
          deallocate(alpha02)
!
          if (i_num_bubbles.ne.0) then
             path1 = "../bubStat/bubble."
             do i=1,i_num_bubbles
               bubmatx = zero
               write(path2,'(i4.4)') i
               fname_bub = trim(path1)//trim(path2)//".dat"
               file_number = 950 + i
!               write(*,*) 'fname_bub =', fname_bub
!               write(*,*) 'file_number = ', file_number
!               write(*,*) 'irstart+1-inistep =', irstart+1-inistep
!       Reading the old bubble information into temporary matrix
               OPEN(unit=file_number, file=fname_bub, 
     &          status="old", action="read", iostat=ierror)
               IF(ierror /= 0) STOP "Error reading Bubble.*.dat"
!               write(*,*) 'ierror of bubble file =', ierror
               do il = 1, (irstart+1-inistep)
                  read(file_number,812) bubmatx(il,1:21)
               enddo
!               write(*,812) bubmatx(1,1:10)
!               write(*,812) bubmatx(irstart+1-inistep,1:10)
               close(file_number) 
 
!       Writing the temporary matrix into newly created files
               OPEN(unit=file_number, file=fname_bub, status="replace",
     &         action="write", iostat=ierror)
               IF(ierror /= 0) STOP "Error creating Bubble.*.dat"
               do il = 1, (irstart+1-inistep)
                  write(file_number,812) bubmatx(il,1:21)
               enddo
               close(file_number)
             end do     !i_num_bubbles
          end if!i_num_bubbles.ne.0

          deallocate(bubmatx)
!----------------------------------
        elseif(irstart+1.eq.inistep) then
!       restart the case from the every beginning so we can totally
!       replace the previous data files        
        open(unit=740,file='../bubStat/alpha.dat',status="replace",
     &       action="write",iostat=ierror)
        if(ierror /= 0) STOP "Error creating file 740"
        close(740)

        if (i_num_bubbles.ne.0) then
           path1 = "../bubStat/bubble."
           do i=1,i_num_bubbles
              write(path2,'(i4.4)') i
              fname_bub = trim(path1)//trim(path2)//".dat"
              file_number = 950 + i
           OPEN(unit=file_number, file=fname_bub, status="replace",
     &          action="write", iostat=ierror)
           IF(ierror /= 0)
     &       STOP "Error creating files:/bubStat/bubble.*.dat"
           close(file_number)
           end do
        end if

        endif   !irstart .lt. (inistep+nl_alpha)
        ENDIF   !existence of alpha.dat

 812  format(ES14.7, 1x, ES14.7, 1x, ES14.7, 1x, ES14.7, 1x,
     &       ES14.7, 1x, ES14.7, 1x, ES14.7, 1x, ES14.7, 1x,
     &       ES14.7, 1x, ES14.7, 1x, ES14.7, 1x, ES14.7, 1x,
     &       ES14.7, 1x, F9.1) 
 813  format(I6, 1x, I4, 1x, I4, 3x, ES14.7, 1x, ES14.7, 1x, ES14.7)

        end     !OpenBubFiles ends

c======================================================================
c======================================================================
c======================================================================
        subroutine BubASSY()
!----------------------------------------------------------------------
!       Called in elmgmr.f
!       This subroutine is used to assembly bubbles' information after
!       the loop over blocks on each processor
!
!----------------------------------------------------------------------
        use bub_track   ! access to bubble information array
        include "common.h"
        include "mpif.h"
        include "auxmpi.h"

        integer i_mrk

        do i = 1, npro
           i_mrk = INT(bub_info(i,11))
           if (i_mrk .gt. 0) then
!       j=1-8: x coord, y coord, z coord, elem_vol, elem_mass (for
!       bubble), x vel, y vel, z vel (for local liquid)
              do j = 1, 8
              procs_dataset(i_mrk,j) = procs_dataset(i_mrk,j)
     &                                 + bub_info(i,j)
              enddo
!       Here the bubble velocities are brought back
              do j = 12, 14
              procs_dataset(i_mrk,j) = procs_dataset(i_mrk,j)
     &                                 + bub_info(i,j)
              enddo
!===============================================================================
!-------------------------------------------------------------------------------
!	This part is for temperature collection
!	j=18: temperature
!	j=19: temperature gradient
!	j=20: The number of element in the innner shell 
!	j=21: The number of element in the outer shell
!-------------------------------------------------------------------------------

              do j=18, 19 
              procs_dataset(i_mrk,j) = procs_dataset(i_mrk,j)
     &                                 + bub_info(i,j)
              enddo
!              if(bub_info(i,18).ne.0.0d0) then
!              write(*,*)'procs_dataset(1,18)', procs_dataset(1,18)
!              endif
              if (bub_info(i,20).ne.0.0d0) then
!              if (bub_info(i,14).eq. epsilon_ls_tmp(i)) then
              procs_dataset(i_mrk,20)=procs_dataset(i_mrk,20)
     &                                +1.0
              endif
              if (bub_info(i,21).ne.0.0d0) then
!              if (bub_info(i,15).eq. epsilon_ls_tmp(i)) then
              procs_dataset(i_mrk,21)=procs_dataset(i_mrk,21)
     &                                +1.0 
              endif
!              write(*,*)'bub_info 20,21',bub_info(i,20),bub_info(i,21) 
!-------------------------------------------------------------------------------
!===============================================================================
!       Find out the upper and lower bounds of element coordinates
!       inside each bubble
              IF(bub_info(i,4).ne.0.0d0) THEN   !Inside bubble
              if(bub_info(i,1).lt.procs_coordDn(i_mrk,1)) 
     &           procs_coordDn(i_mrk,1) = bub_info(i,1)
              if(bub_info(i,2).lt.procs_coordDn(i_mrk,2))
     &           procs_coordDn(i_mrk,2) = bub_info(i,2)
              if(bub_info(i,3).lt.procs_coordDn(i_mrk,3))
     &           procs_coordDn(i_mrk,3) = bub_info(i,3)

              if(bub_info(i,1).gt.procs_coordUp(i_mrk,1))
     &           procs_coordUp(i_mrk,1) = bub_info(i,1)
              if(bub_info(i,2).gt.procs_coordUp(i_mrk,2))
     &           procs_coordUp(i_mrk,2) = bub_info(i,2)
              if(bub_info(i,3).gt.procs_coordUp(i_mrk,3))
     &           procs_coordUp(i_mrk,3) = bub_info(i,3)
              ENDIF

!       The minimum level set value inside a bubble
              if(bub_info(i,10).lt.procs_dataset(i_mrk,9))
     &        procs_dataset(i_mrk,9) = bub_info(i,10)

!       Number of elements in each bubble (non-zero elem vol)
              if(bub_info(i,4).ne.0.0d0) then
                 procs_dataset(i_mrk,10) = procs_dataset(i_mrk,10)
     &                                     + 1.0
              endif

              if(bub_info(i,6).ne.0.0d0) then
!       Number of elements in the liquid shell (non-zero velocity)
                 procs_dataset(i_mrk,11) = procs_dataset(i_mrk,11)
     &                                     + 1.0d0
!       Find out the points with min and max d2wall and the associated
!       velocities for each bubble

                 if(bub_info(i,9).lt.Shear_NodeMin(i_mrk,1)) then
                    Shear_NodeMin(i_mrk,1) = bub_info(i,9)         !min d2w
                    Shear_NodeMin(i_mrk,2) = bub_info(i,6)         !x vel
                 endif
                 if(bub_info(i,9).gt.Shear_NodeMax(i_mrk,1)) then
                    Shear_NodeMax(i_mrk,1) = bub_info(i,9)         !max d2w
                    Shear_NodeMax(i_mrk,2) = bub_info(i,6)         !x vel
                 endif
              endif

           endif !i_mrk
        enddo !npro

        end     !BubASSY = bubble assembly

c======================================================================
c======================================================================
c======================================================================
        subroutine BubMPIprocess()
!----------------------------------------------------------------------
!       Called in elmgmr.f
!       This subroutine is dealing with the MPI processing of bubble
!       information from different processors
!
!----------------------------------------------------------------------
        use bub_track   ! access to bubble information array
        use bubboil_info  ! for boiling & condensation model
        include "common.h"
        include "mpif.h"
        include "auxmpi.h"

        real*8 one_procs,       all_procs
        real*8 phi_tmp,         phi_min
        real*8 eq_rad
        real*8 Shear_dummy01(2), Shear_dummy02(2), Shear_dummy03
        real*8 diffinY
        real*8 bubble_vol_temp,  numshell_temp1,   numshell_temp2,
     &         bubble_tempG_temp
        allocate ( unive_dataset(i_num_bubbles, 21) )
        allocate ( unive_coordDn(i_num_bubbles, 3)  )
        allocate ( unive_coordUp(i_num_bubbles, 3)  )
        allocate ( bubbl_coordDf(i_num_bubbles, 3)  )
        allocate ( Shear(i_num_bubbles,4)           )
!        allocate ( bubble_vol(i_num_bubbles)        )
!        allocate ( numshell(i_num_bubbles, 2)       )

        unive_dataset   = zero
        unive_coordDn   = zero
        unive_coordUp   = zero
        bubbl_coordDf   = zero
        Shear           = zero

        if (numpe .gt. 1) then
           do i = 1, i_num_bubbles
              do j = 1, 21
                 if (j.eq.9) then
!       Find out the minimum level set value for each bubble
                 phi_tmp = procs_dataset(i,9)
                 call MPI_ALLREDUCE (phi_tmp, phi_min, 1,
     &                MPI_DOUBLE_PRECISION, MPI_MIN, 
     &                MPI_COMM_WORLD, ierr)
                 unive_dataset(i,9) = phi_min
                 else
!       The rest includes:
!       x, y, z coord, elem_vol, elem_mass, x,y,z vel (for bubble),
!       x vel, y vel, z vel (for local liquid), curvature
!       bub_elem_count
                 one_procs = procs_dataset(i,j)
                 call MPI_ALLREDUCE (one_procs, all_procs, 1,
     &                MPI_DOUBLE_PRECISION,MPI_SUM, MPI_COMM_WORLD,ierr)
                 unive_dataset(i,j) = all_procs
                 endif
              enddo !index j
!------------------------------------------------------------------------
!       Here, we determine the x bounds, y bounds, z bounds for each
!       bubble
!------------------------------------------------------------------------
              do j = 1, 3
                 one_procs = procs_coordDn(i,j)
                 call MPI_ALLREDUCE (one_procs, all_procs, 1,
     &                MPI_DOUBLE_PRECISION, MPI_MIN,
     &                MPI_COMM_WORLD, ierr)
                 unive_coordDn(i,j) = all_procs

                 one_procs = procs_coordUp(i,j)
                 call MPI_ALLREDUCE (one_procs, all_procs, 1,
     &                MPI_DOUBLE_PRECISION, MPI_MAX,
     &                MPI_COMM_WORLD, ierr)
                 unive_coordUp(i,j) = all_procs
!       This range of elments coordinates inside a bubble will be used
!       is crossing the periodic planes
                 bubbl_coordDf(i,j) = unive_coordUp(i,j) -
     &                                unive_coordDn(i,j)
              enddo
!------------------------------------------------------------------------
!       Find the two points to calculate local shear
!------------------------------------------------------------------------
!       The one with minimum d2wall
              Shear_dummy01(1) = Shear_NodeMin(i,1)
              Shear_dummy01(2) = myrank  ! myrank is to be coerced to a real*8

              call MPI_ALLREDUCE (Shear_dummy01, Shear_dummy02, 1,
     &                MPI_2DOUBLE_PRECISION, MPI_MINLOC,
     &                MPI_COMM_WORLD, ierr)

              Shear(i,1) = Shear_dummy02(1)
              Shear_dummy03 = Shear_NodeMin(i,2)
!       Broadcast the x vel (for local liq) associated with min d2wall
              call MPI_Bcast(Shear_dummy03,1,MPI_DOUBLE_PRECISION,
     &                       INT(Shear_dummy02(2)),MPI_COMM_WORLD,ierr)
              Shear(i,2) = Shear_dummy03
!       The one with maximum d2wall
              Shear_dummy01(1) = Shear_NodeMax(i,1)
              Shear_dummy01(2) = myrank  ! myrank is coerced to a real*8

              call MPI_ALLREDUCE (Shear_dummy01, Shear_dummy02, 1,
     &                MPI_2DOUBLE_PRECISION, MPI_MAXLOC,
     &                MPI_COMM_WORLD, ierr)

              Shear(i,3) = Shear_dummy02(1)
              Shear_dummy03 = Shear_NodeMax(i,2)
!       Broadcast the x vel (for local liq) associated with max d2wall
              call MPI_Bcast(Shear_dummy03,1,MPI_DOUBLE_PRECISION,
     &                       INT(Shear_dummy02(2)),MPI_COMM_WORLD,ierr)
              Shear(i,4) = Shear_dummy03

           enddo !i_num_bubbles

!------------------------------------------------------------------------
!       Calculate the average properties
!------------------------------------------------------------------------
        avg_info = zero
        if(myrank.eq.master) then

        do i = 1, i_num_bubbles
!       The average coordinates for each bubble
           if(unive_dataset(i,10).ne.0.0) then
           avg_info(i,1:3)  = unive_dataset(i,1:3)  /unive_dataset(i,10)
           avg_info(i,11:13)= unive_dataset(i,12:14)/unive_dataset(i,10)
           endif
           avg_info(i,14)   = unive_dataset(i,10)

           eq_rad = (3.0*unive_dataset(i,4)/(4.0*pi))**(1.0/3.0)
           avg_info(i,4)    = eq_rad
!       The bubble mass
           avg_info(i,5)    = unive_dataset(i,5)
!           bubble_vol(i)    = unive_dataset(i,4)
!       The deformability factor is defined by the ratio b/w min level
!       set value with equivalent radius, which is 1.0 for spherical
!       bubbles ideally.
           if(eq_rad.ne.0.0d0)
     &     avg_info(i,6)    = abs(unive_dataset(i,9))/eq_rad

!       The local liquid velocity components near the bubble
           if(unive_dataset(i,11).ne.0.0d0)
     &     avg_info(i,7:9)  = unive_dataset(i,6:8)/unive_dataset(i,11)

!       The local liquid shear around the bubble
           diffinY = Shear(i,3) - Shear(i,1)
           if(diffinY.ne.0.0d0)
     &     avg_info(i,10)   = (Shear(i,4) - Shear(i,2))/diffinY

!       Adjust the bubble centers coordinates for those crossing
!       periodic planes (1% tolerance is allowed)
           if(abs(bubbl_coordDf(i,1)-XLEN).lt.0.01d0*XLEN) then
              if(avg_info(i,1) .ge. Xmid ) then
                 avg_info(i,1) = DomainSize(2) - eq_rad *
     &          (avg_info(i,1) - Xmid)/(0.5d0*XLEN-eq_rad)
              else
                 avg_info(i,1) = DomainSize(1) + eq_rad *
     &          (Xmid - avg_info(i,1))/(0.5d0*XLEN-eq_rad)
              endif
           endif

           if(abs(bubbl_coordDf(i,2)-YLEN).lt.0.01d0*YLEN) then
              if(avg_info(i,2) .ge. Ymid ) then
                 avg_info(i,2) = DomainSize(4) - eq_rad *
     &          (avg_info(i,2) - Ymid)/(0.5d0*YLEN-eq_rad)
              else
                 avg_info(i,2) = DomainSize(3) + eq_rad *
     &          (Ymid - avg_info(i,2))/(0.5d0*YLEN-eq_rad)
              endif
           endif

           if(abs(bubbl_coordDf(i,3)-ZLEN).lt.0.01d0*ZLEN) then
              if(avg_info(i,3) .ge. Zmid ) then
                 avg_info(i,3) = DomainSize(6) - eq_rad *
     &          (avg_info(i,3) - Zmid)/(0.5d0*ZLEN-eq_rad)
              else
                 avg_info(i,3) = DomainSize(5) + eq_rad *
     &          (Zmid - avg_info(i,3))/(0.5d0*ZLEN-eq_rad)
              endif
           endif

!===============================================================================
!-------------------------------------------------------------------------------
!	This part is for temperature gradient collection
!						-Mengnan Li
!-------------------------------------------------------------------------------
!... Store temperature and temperature gradient around the bubble
!           if(unive_dataset(i,18).ne.0.0)then
           avg_info(i,18) = unive_dataset(i,18)/unive_dataset(i,20)
           avg_info(i,19) = unive_dataset(i,19)
!... Store the number of cell outside the bubble
           avg_info(i,20) = unive_dataset(i,20)
	   avg_info(i,21) = unive_dataset(i,21)
!           endif
        if(bubboil.eq.1.0.or.bubgrow.eq.1.0) then 
           write(*,*)'# of elem(out),# of elem(in)',
     &               avg_info(i,20),avg_info(i,21)
        endif
        if(bubboil.eq.1.0) then
           write(*,*)'temp,tempG',
     &               avg_info(i,18),avg_info(i,19)
        endif
!-------------------------------------------------------------------------------
        enddo  !i_num_bubbles

        endif !myrank

        else
           write(*,*) 'Bubble tracking was not coded for 1-procs case!'
        endif
       if (numpe.gt.1) then
!       call MPI_BARRIER(MPI_COMM_WORLD,ierr)
       do i=1,i_num_bubbles
        if(myrank.eq.master)then
        bubble_vol_temp   =  unive_dataset(i,4)
        bubble_tempG_temp =  unive_dataset(i,19) 
        numshell_temp1    =  avg_info(i,20)
        numshell_temp2    =  avg_info(i,21)
        endif
!        write(*,*) 'I am here', bubble_vol(i)
          call MPI_Bcast(bubble_tempG_temp,1,MPI_DOUBLE_PRECISION,
     &                       master,MPI_COMM_WORLD,ierr)
!          for average temperature gradient
          call MPI_Bcast(numshell_temp1,1,MPI_DOUBLE_PRECISION,
     &                       master,MPI_COMM_WORLD,ierr)
!          for number of element outside the bubble
          call MPI_Bcast(numshell_temp2,1,MPI_DOUBLE_PRECISION,
     &                       master,MPI_COMM_WORLD,ierr)
!          for number of element inside the bubble
          call MPI_Bcast(bubble_vol_temp,1,MPI_DOUBLE_PRECISION,
     &                   master,MPI_COMM_WORLD,ierr)
!          for each bubble volume
          bubble_vol(i)  =  bubble_vol_temp
          bubble_tempG(i)=  bubble_tempG_temp
          numshell(i,1)  =  numshell_temp1
          numshell(i,2)  =  numshell_temp2
       enddo
!       call MPI_BARRIER(MPI_COMM_WORLD,ierr)
       endif

!       if(myrank.eq.10) write(*,*) 'bubble_vol =', bubble_vol,
!     &                                         bubble_vol_temp

        deallocate ( unive_dataset )
        deallocate ( unive_coordDn )
        deallocate ( unive_coordUp )
        deallocate ( bubbl_coordDf )
        deallocate ( Shear         )
!        deallocate ( bubble_vol    )
!        deallocate ( numshell      )

        end     !BubMPIprocess ends

c======================================================================
c======================================================================
c======================================================================
        subroutine banmaUpdate(xl, yl, banma, ien, bml)
!----------------------------------------------------------------------
!       Called in asigmr.f
!       This subroutine will update the marker field in each time
!       iteration.
!       Banma is the Chinese phonetic translation of Zebra, an African
!       wild horse with black-and-white stripes (markers) and an
!       erect mane.
!       The code will first find out the max of banma values in each
!       element and update the rest based on their levelset value on a
!       certain processor, the markers are updated from one element to
!       another, and this update is conducting samutaniously on
!       different processors
!
!----------------------------------------------------------------------
        use bub_track   ! access to bubble information array
        include "common.h"

        dimension banma(nshg,1),              bml(npro,nshl,1),
     &            yl(npro,nshl,ndofl),        xl(npro,nenl,nsd)
        real*8  bmtmp, xltmp(3)
        real*8  bmmax, Mrk_dist, Bub_radi
        real*8  deltS           !the thickness of liquid shell in the
                                !unit of epsilon
        integer:: imrk, bmflag

        deltS = 3.0d0

        IF(lstep.eq.ts_hold) THEN
!c... localization of marker field
        call local (banma,  bml,     ien,     1,      'gather  ')

        do i = 1,npro
           bmmax = 0.0d0
           do n = 1,nshl
              if(bml(i,n,1).gt.bmmax) bmmax = bml(i,n,1)
           enddo

           do n = 1,nshl
!c... if the nodal point is in the liquid
              if(yl(i,n,6).gt.0.0d0) then ! 0.0d0
                bml(i,n,1) = 0.0d0
!c... for the nodal points in the gas and near interface region...
              else
                bml(i,n,1) = bmmax
              endif
           enddo
        enddo
!
!c.... assemble the marker field
        call local (banma,    bml,   ien,    1,  'globaliz')

        ELSE
!       Do not use this part in the very first timestep because the
!       bub_cent is needed 
!       Localization of marker field
        call local (banma,  bml,     ien,     1,      'gather  ')

        do i = 1,npro
           bmmax = 0.0d0
!       bmflag = 1 means there is only one marker in the element        
           bmflag= 1
           xltmp = zero
           do n = 1,nshl
              if(bml(i,n,1).gt.bmmax) bmmax = bml(i,n,1)
!              if(bml(i,n,1).ne.bmmax) bmflag= bmflag+1
           enddo
!       Update the maker field differently in bubble elments and near
!       interface shell
           do n = 1,nshl
!       If the node is in the liquid (eg. levelset .gt. 3epsilon)
              if(yl(i,n,6).gt.deltS*epsilonBT) then
                bml(i,n,1) = 0.0d0
!       If the node is in the bubble region
              elseif(yl(i,n,6).lt.0.0d0) then
                bml(i,n,1) = bmmax
!       If the node is in the near interface liquid shell
              else
                xltmp(1:3) = xl(i,n,1:3)
                call ReColor(xltmp,bmtmp)
                bml(i,n,1) = bmtmp

              endif

           enddo        !nshl
        enddo           !npro

c.... assemble the marker field
        call local (banma,    bml,   ien,    1,  'globaliz')

        ENDIF   !ts_hold
        end     !banmaUpdate ends

c======================================================================
c======================================================================
c======================================================================
        subroutine BubCollect(u1,       u2,     u3,     Sclr, dist2w,
     &                        xx,       yl,     bml,    elemvol_local,
     &                        rho,      Tempb,  gytemp, shell_elem)
!----------------------------------------------------------------------
!       Called in e3ivar.f
!       This subroutine is dealing with the bubble information
!       collection at the very bottom level.
!       In bub_info(i_num_bubbls, 10)
!       bubble-wise: x coord, y coord, z coord, elem vol, levelset;
!       local liq  : x vel, y vel, z vel, d2wall;
!       BT field   : marker(bubble ID)
!
!----------------------------------------------------------------------
        use  bub_track
        use  spat_var_eps ! for spatially varying epsilon_ls
        include "common.h"

        dimension u1(npro),  u2(npro),  u3(npro)
        dimension dist2w(npro), rho(npro) 
        dimension yl(npro,nshl,ndof),   xx(npro,nsd),
     &            Sclr(npro)
        dimension bml(npro,nshl,1)
        real*8    elemvol_local(ibksiz)
        real*8    bmmax 
        real*8    denswght
        dimension Tempb(npro), gytemp(npro,5), shell_elem(npro,2) 
        real*8    Tempb, gytemp, shell_elem, epsilon_lst_tmp
                  ! Temperature Collection Mengnan 9/6/15
        rholiq=datmat(1,1,1)
        rhogas=datmat(1,1,2) 

        bub_info = zero

        do i = 1, npro
        epsilon_lst_tmp = epsilon_ls *
     &                       elem_local_size(lcblk(1,iblk)+i-1)
!... collect LS value & markers for bubble and liquid shell
           if(Sclr(i) .le. 3.0d0*epsilonBT) then
              if(rholiq.eq.rhogas) then
                 denswght = 1.0d0
              else
                 denswght  =(rholiq-rho(i))/(rholiq-rhogas)
              endif
              bmmax = 0.0
              do n = 1, nshl
                if(bml(i,n,1).gt.bmmax) bmmax = bml(i,n,1)
              enddo
              bub_info(i,11) = bmmax
!       Find out the minimum level set value
              do n = 1, nshl
                if(yl(i,n,6).lt.bub_info(i,10))
     &             bub_info(i,10) = yl(i,n,6)
              enddo
!       collect the local liquid velocity and the y coord of liquid
!       shell elments around the bubble
              if(Sclr(i).gt.epsilonBT) then
                bub_info(i,6)  = u1(i)
                bub_info(i,7)  = u2(i)
                bub_info(i,8)  = u3(i)
                bub_info(i,9)  = dist2w(i) 
!              elseif(Sclr(i).le.-1.0d0*epsilonBT) then  !collet info from
                                                        !smaller sphere
              elseif(Sclr(i).le.0.0d0) then
!... collect the bubble information in details
                bub_info(i,1)  = xx(i,1)
                bub_info(i,2)  = xx(i,2)
                bub_info(i,3)  = xx(i,3)
                bub_info(i,4)  = elemvol_local(i)
                bub_info(i,5)  = elemvol_local(i)*denswght*rhogas
                bub_info(i,12) = u1(i)
                bub_info(i,13) = u2(i)
                bub_info(i,14) = u3(i)
              endif ! Sclr(i)
           endif !Sclr(i) for the bubble region and liquid shell
!===============================================================================
!-------------------------------------------------------------------------------

!	Collect the temperature and its gradient for boiling
!	18: temperature around the bubble
!	19: temperature gradient around the bubble
!	20: number of elements in the shell outside the bubble
!	21: number of elements in the shell inside the bubble	
!-------------------------------------------------------------------------------
               if (bubboil.eq.1.or.bubgrow.eq.1)then
                bub_info(i,18)= Tempb(i)
                bub_info(i,19)= gytemp(i,5) ! test for tempgradient
                bub_info(i,20)= shell_elem(i,1) ! count element in the shell
                                                ! outside the bubble
                bub_info(i,21)= shell_elem(i,2) ! count element in the shell
                                                ! inside the bubble 
!                if (Tempb(i).ne.0.0d0)then
!               write(*,*)'bub_info(i,18)',bub_info(i,18)
!                endif
                endif
        enddo !npro

        end     !BubCollect ends
c======================================================================
c======================================================================
c======================================================================
        subroutine BubPrintOut(vf)
!----------------------------------------------------------------------
!       Called in itrdrv.f
!       This subroutine will print out the bubble information into
!       external files which can be further processed
!
!       bub_cent        :x coord, y coord, z coord, radius, bubble ID
!----------------------------------------------------------------------
        use bub_track   ! access to bubble information array
        include "common.h"

        INTEGER:: file_number, istp, iNghost
        integer:: Nbubdef
        CHARACTER(LEN=100) :: path1, path2, path3, fname_bub
        real*8  vf

        write(*,'(1x,A,F14.7,A)') 'Current void: ',vf*100.0d0,'%'
        path1 = "../bubStat/bubble."
        Nbubdef = 0             !number of deformed bubbles
        iNghost = 0
        
        do i = 1, i_num_bubbles
        IF(avg_info(i,4).gt.0.0d0) THEN
           bub_cent(i,1:4) = avg_info(i,1:4)
           bub_cent(i,5)   = REAL(i)

!       Find out the ghost bubbles and save them into bub_cent
           if(avg_info(i,1)-DomainSize(1).lt.GhostRatio*XLEN) then
              iNghost = iNghost + 1
              bub_cent(i_num_bubbles+iNghost,1)   = bub_cent(i,1)+XLEN
              bub_cent(i_num_bubbles+iNghost,2:3) = bub_cent(i,2:3)
              bub_cent(i_num_bubbles+iNghost,4)   = bub_cent(i,4)
              bub_cent(i_num_bubbles+iNghost,5)   = bub_cent(i,5)
           endif

           if(avg_info(i,2)-DomainSize(3).lt.GhostRatio*YLEN) then
              iNghost = iNghost + 1
              bub_cent(i_num_bubbles+iNghost,1)   = bub_cent(i,1)
              bub_cent(i_num_bubbles+iNghost,2)   = bub_cent(i,2)+YLEN
              bub_cent(i_num_bubbles+iNghost,3)   = bub_cent(i,3)
              bub_cent(i_num_bubbles+iNghost,4)   = bub_cent(i,4)
              bub_cent(i_num_bubbles+iNghost,5)   = bub_cent(i,5)
           endif

           if(avg_info(i,3)-DomainSize(5).lt.GhostRatio*ZLEN) then
              iNghost = iNghost + 1
              bub_cent(i_num_bubbles+iNghost,1:2) = bub_cent(i,1:2)
              bub_cent(i_num_bubbles+iNghost,3)   = bub_cent(i,3)+ZLEN
              bub_cent(i_num_bubbles+iNghost,4)   = bub_cent(i,4)
              bub_cent(i_num_bubbles+iNghost,5)   = bub_cent(i,5)
           endif

           if(DomainSize(2)-avg_info(i,1).lt.GhostRatio*XLEN) then
              iNghost = iNghost + 1
              bub_cent(i_num_bubbles+iNghost,1)   = bub_cent(i,1)-XLEN
              bub_cent(i_num_bubbles+iNghost,2:3) = bub_cent(i,2:3)
              bub_cent(i_num_bubbles+iNghost,4)   = bub_cent(i,4)
              bub_cent(i_num_bubbles+iNghost,5)   = bub_cent(i,5)
           endif

           if(DomainSize(4)-avg_info(i,2).lt.GhostRatio*YLEN) then
              iNghost = iNghost + 1
              bub_cent(i_num_bubbles+iNghost,1)   = bub_cent(i,1)
              bub_cent(i_num_bubbles+iNghost,2)   = bub_cent(i,2)-YLEN
              bub_cent(i_num_bubbles+iNghost,3)   = bub_cent(i,3)
              bub_cent(i_num_bubbles+iNghost,4)   = bub_cent(i,4)
              bub_cent(i_num_bubbles+iNghost,5)   = bub_cent(i,5)
           endif

           if(DomainSize(6)-avg_info(i,3).lt.GhostRatio*ZLEN) then
              iNghost = iNghost + 1
              bub_cent(i_num_bubbles+iNghost,1:2) = bub_cent(i,1:2)
              bub_cent(i_num_bubbles+iNghost,3)   = bub_cent(i,3)-ZLEN
              bub_cent(i_num_bubbles+iNghost,4)   = bub_cent(i,4)
              bub_cent(i_num_bubbles+iNghost,5)   = bub_cent(i,5)
           endif
        ENDIF !The bubble exists physically

!       the bubble is recognized to be deformed when the factor is less
!       than certain limit
           if(avg_info(i,4).gt.0.0d0 .and. avg_info(i,6).lt.0.8d0) 
     &        Nbubdef = Nbubdef + 1

!       Record the detailed bubble information into files
           write(path2,'(i4.4)') i
           fname_bub = trim(path1)//trim(path2)//".dat"
           file_number = 950 + i
           OPEN(unit=file_number, file=fname_bub, status="old",
     &          action="write", position="append", iostat=ierror)
           write(file_number,832) avg_info(i,:)
           close(file_number)

        end do

        open(unit=740,file='../bubStat/alpha.dat',status="old",
     &       action="write", position="append", iostat=ierror)
        if(ierror /= 0) STOP "Error creating file 740"
        write(740,833) lstep+1, Nbubtot, Nbubdef, C_int_adjust, vf,
     &  time
        close(740)

 832  format(ES14.7, 1x, ES14.7, 1x, ES14.7, 1x, ES14.7, 1x,
     &       ES14.7, 1x, ES14.7, 1x, ES14.7, 1x, ES14.7, 1x,
     &       ES14.7, 1x, ES14.7, 1x, ES14.7, 1x, ES14.7, 1x,
     &       ES14.7, 1x, F9.1)
 833  format(I6, 1x, I4, 1x, I4, 3x, ES14.7, 1x, ES14.7, 1x, ES14.7)

        end     !BubPrintOut ends
c======================================================================
c======================================================================
c======================================================================
        subroutine CountAlphaLines(LineCount, inistep)
!----------------------------------------------------------------------
!       This subroutine will return the intial timestep in alpha.dat and
!       total number of lines of which
!
!----------------------------------------------------------------------
        use, intrinsic :: iso_fortran_env

        integer, intent(out) :: LineCount
        integer :: Read_Code
        character (len=200) :: line
        character (len=200) :: filename
        integer :: inistep, na1, na2
        real*8  :: dummya(3)

        open(unit=740,file='../bubStat/alpha.dat',status="old",
     &       action="read",iostat=ierror)
        read(740,823) inistep, na1, na2, dummya(:)
        close(740)

        open(unit=51,file='../bubStat/alpha.dat',
     &       status="old",access='sequential',
     &       form='formatted',action='read')

        LineCount = 0

        ReadLoop: do

        read (51, '(A)', iostat=Read_Code)  line

        if(Read_Code /= 0) then
           if(Read_Code == iostat_end ) then
              exit ReadLoop    ! end of file --> line count found
           else
              write ( *, '( / "read error: ", I0 )' )  Read_Code
              stop
           endif
        endif

        LineCount = LineCount + 1

!        write (*, '( I0, ": ''", A, "''" )' )  LineCount, trim (line)
        if(len_trim(line)==0)then
           write(*,'("The above is an empty or all blank line.")')
        endif

        end do ReadLoop

!        write (*, *) "found", LineCount, " lines"
        close(51)

 823  format(I6, 1x, I4, 1x, I4, 3x, ES14.7, 1x, ES14.7, 1x, ES14.7)
        end

c======================================================================
c======================================================================
c======================================================================
        subroutine reCountBub()
!----------------------------------------------------------------------
!       Called in itrdrv.f
!       This subroutine will count how many bubble exist in current time
!       iteration and determine the number of ghost bubbles
!
!----------------------------------------------------------------------
        use bub_track   ! access to bubble information array
        include "common.h"

!       the bubble is recognized when its volume is larger than zero
        do i = 1, i_num_bubbles
           if(avg_info(i,4).gt.0.0d0) then
              Nbubtot = Nbubtot + 1
              if(avg_info(i,1)-DomainSize(1).lt.GhostRatio*XLEN)
     &              Nghost = Nghost + 1
              if(avg_info(i,2)-DomainSize(3).lt.GhostRatio*YLEN)
     &              Nghost = Nghost + 1
              if(avg_info(i,3)-DomainSize(5).lt.GhostRatio*ZLEN)
     &              Nghost = Nghost + 1

              if(DomainSize(2)-avg_info(i,1).lt.GhostRatio*XLEN)
     &              Nghost = Nghost + 1
              if(DomainSize(4)-avg_info(i,2).lt.GhostRatio*YLEN)
     &              Nghost = Nghost + 1
              if(DomainSize(6)-avg_info(i,3).lt.GhostRatio*ZLEN)
     &              Nghost = Nghost + 1
           endif !The bubble exists physically
        enddo


        end     !RecountBub ends
c======================================================================
c======================================================================
c======================================================================
        subroutine ReColor(xltmp,bmtmp)
!----------------------------------------------------------------------
!       This subroutine is used to recolor the node coordinates and
!       current bubble centers 
!
!----------------------------------------------------------------------
        use bub_track   ! access to bubble information array
        include "common.h"

        real*8  bmtmp, xltmp(3)
        real*8  Mdist, MdistMin
        integer ib

        bmtmp = 0.0d0
        Mdist = 0.0d0
        MdistMin = 1.0E6

        do ib = 1,(i_num_bubbles+Nghost) 
!       The bubble must exist physically and in the neighborhood       
           if(bub_cent(ib,4).gt.0.0d0 .and.  
     &        abs(xltmp(1)-bub_cent(ib,1)).le.NbrhdRatio*XLEN) then     

              Mdist = sqrt((xltmp(1)-bub_cent(ib,1))**2 +
     &                  (xltmp(2)-bub_cent(ib,2))**2 +
     &                  (xltmp(3)-bub_cent(ib,3))**2)-bub_cent(ib,4)

              if(MdistMin.gt.Mdist) then
                 MdistMin = Mdist
                 bmtmp = bub_cent(ib,5)         !update the marker
              endif

           endif



        enddo

        end     !ReColor ends
c======================================================================
c======================================================================
c======================================================================

