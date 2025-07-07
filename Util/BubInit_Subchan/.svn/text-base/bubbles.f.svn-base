        PROGRAM bubbles
        IMPLICIT NONE
!----------------------------------------------------------------------
!
! This routine is developed to created the required bubble information
! such as, bubble ID, bubble coordiantes and radii. And the Monta Carlo
! rejection sampling method is used to avoid the intersection between
! bubbles and the walls or between different bubbles
!
! Some variables:
! bundle(2)             The number of subchannels in y, z directions
! mixvanes(2)           The range of mixing vanes region
! coordrep              The flag for coordinate replacement
! bubdupl               The bubble array for duplication in 2x2 case
!
! Jun Fang, 2015.
!----------------------------------------------------------------------

        INTEGER  nbub,    ierror,        i,j,k
        integer::clock,   n_dimension=4, size
        integer  bundle(2)
        REAL*8   rbub,    rbox1(6),  rbox2(3), mixvanes(2)
        real*8   distance, mindist,  mindistmp, r
        real*8   ghostupplim, ghostlowlim
        real*8   random_replace(4),  temp(1,3), distw(4)
        REAL*8, ALLOCATABLE:: bubcoord(:,:), minw(:),
     &                        random1d(:),   random2d(:,:)
        LOGICAL::gate=.TRUE., check,         gothroughlist,
     &           coordrep=.FALSE.

        INTEGER::nghost,ntotal
        REAL*8:: ghost_area_ratio,outline(3),bubtmp(3)
        REAL*8, ALLOCATABLE:: bublist(:,:)
        real*8, allocatable:: bubdupl(:,:)        

        INTEGER::xres,yres,zres
        REAL*8:: xcoord=0,ycoord=0,zcoord=0,
     &           xstep,ystep,zstep,A_bubble=0
        REAL*8, ALLOCATABLE:: xvoid(:),yvoid(:),zvoid(:)

        CHARACTER (LEN=3) :: DomainType
!----------------------------------------------------------------------
!       Creating the output files 
!       Reading input parameters from files (options.dat)
!----------------------------------------------------------------------

        OPEN(UNIT=1,FILE='options.dat',IOSTAT=ierror)
        IF(ierror .NE. 0)STOP"Error opening file options.txt"

        OPEN(UNIT=2,FILE='bubbles.inp',IOSTAT=ierror)
        IF(ierror .NE. 0)STOP"Error opening file output.txt"

        OPEN(UNIT=6,FILE='Logs.out',IOSTAT=ierror)
        IF(ierror .NE. 0)STOP"Error opening file output.txt"

!        OPEN(unit=5,FILE='void_x.dat',IOSTAT=ierror)
!        IF(ierror .NE. 0)STOP"Error opening file void_x.dat"

        READ(1,*)
        READ(1,*)DomainType             !Domain type
        READ(1,*)
        READ(1,*)nbub                   !number of bubbles
        READ(1,*)
        READ(1,*)rbub                   !bubble radius
        READ(1,*)
        READ(1,*)rbox1(:)               !domain dimensions
        READ(1,*)
        READ(1,*)mixvanes(:)            !range of mixing vanes
        READ(1,*)
        READ(1,*)r                      !rod radius
        READ(1,*)
        READ(1,*)mindist                !min distance b/w bubbles
        READ(1,*)
        READ(1,*)ghost_area_ratio       !ghost region ratio
        READ(1,*)
        READ(1,*)xres                   !x resolution
        READ(1,*)
        READ(1,*)yres
        READ(1,*)
        READ(1,*)zres

        if (DomainType .eq. "Sub") then
           write(*,*) 'Domain Type: PWR subchannel'
        elseif (DomainType .eq. "2x2") then
           write(*,*) 'Domain Type: PWR 2x2 structure'
           bundle(:) = 2
           write(*,"(A, ES16.6, A, ES16.6)") 
     &                ' Mixing region is from', mixvanes(1), 
     &                ' to', mixvanes(2)
        endif

!----------------------------------------------------------------------
!       The initialization of bubble center coordinates
!----------------------------------------------------------------------

!Subtracting the radius from the box walls, gives room, lets me not worry
! about crossing at wall

        DO i=2,6,2
           rbox2(i/2)=rbox1(i)-rbox1(i-1)
        END DO
        outline(:)=rbox2(:)
        if (DomainType .eq. "2x2") then 
!...The mixing vanes region will also get shifted with the rest, and the
!array mixvanes becomes the relative distances from the channel inlet.
           mixvanes(:) = mixvanes(:) - rbox1(1)
           write(*,"(A, 3ES16.6)")' The outline of a quarter domain:',
     &          outline(:)
        endif

        ghost_area_ratio=ghost_area_ratio/100.0

!Forming the random number arrays and the bubble coordinate arrays.
        allocate( random1d(nbub*3), random2d(nbub,3) )  
        ALLOCATE( bubcoord(nbub,3), minw(nbub) )

!forming the random array
        CALL SYSTEM_CLOCK(COUNT=clock)
        CALL random_number_init(nbub*3,clock)
        CALL RANDOM_NUMBER(random1d)

        DO i=1,nbub
           random2d(i,:)=random1d((i*3-2):i*3)
        END DO

!...forming the initial bubble coordinate array

        DO i=1,nbub
           DO j=1,3
              bubcoord(i,j)=random2d(i,j)*rbox2(j) !+ rbox1(2*j-1)
           END DO
        END DO
!----------------------------------------------------------------------
!       check if the first bubble center is in the right place, replace it if
!       not
!       The radius of fuel rod (r) is increased a little to avoid that
!       bubble's initial position is too close to walls
!----------------------------------------------------------------------
        r = r + rbub

        distw(1)=sqrt((bubcoord(1,3))**2 + (bubcoord(1,2))**2)
        distw(2)=sqrt((bubcoord(1,3)-rbox2(3))**2 + (bubcoord(1,2))**2)
        distw(3)=sqrt((bubcoord(1,3))**2 + (bubcoord(1,2)-rbox2(2))**2)
        distw(4)=sqrt((bubcoord(1,3)-rbox2(3))**2 + 
     &                (bubcoord(1,2)-rbox2(2))**2)
        minw(1)=distw(1)
        DO j = 2, 4
           IF(minw(1).GE. distw(j)) minw(1)=distw(j)
        END DO

!...Check if the bubble is in mixing vanes region
        if (DomainType .eq. "2x2") then
            call checkmixing(bubcoord(1,1), coordrep, mixvanes, rbub)
        endif

        IF(minw(1).LE.(r+rbub) .or. coordrep) THEN
          L0:DO
          IF(gate)THEN
           CALL random_number_init(n_dimension,int(random1d(i)*1000000))
           CALL RANDOM_NUMBER(random_replace)
           gate=.FALSE.
          ELSE
           CALL random_number_init( n_dimension, 
     &                          int(random_replace(2)*1000000))
           CALL RANDOM_NUMBER(random_replace)
          END IF

          DO k=1,3
           bubcoord(1,k)=rbox2(k)*random_replace(k+1) !+ rbox1(2*k-1)
          END DO
!...Check if the bubble is in mixing vanes region
          if (DomainType .eq. "2x2") then
             call checkmixing(bubcoord(1,1), coordrep, mixvanes, rbub)
          endif
!          write(*,*) coordrep
!...  Check the wall again !
          distw(1)=sqrt((bubcoord(1,3))**2 + (bubcoord(1,2))**2)
          distw(2)=sqrt((bubcoord(1,3)-rbox2(3))**2 + 
     &                  (bubcoord(1,2))**2)
          distw(3)=sqrt((bubcoord(1,3))**2 + 
     &                  (bubcoord(1,2)-rbox2(2))**2)
          distw(4)=sqrt((bubcoord(1,3)-rbox2(3))**2 + 
     &                  (bubcoord(1,2)-rbox2(2))**2)
          minw(1)=distw(1)
          DO j = 2, 4
             IF(minw(1).GE. distw(j)) minw(1)=distw(j)
          END DO

          IF( minw(1).LE.(r+rbub) .or. coordrep) CYCLE L0

          EXIT L0
          END DO L0
        ENDIF
        write(*,*) 'The first bubble center coordinates: '
        write(*,"(3(1xES16.6))") bubcoord(1,:)
        if(minw(1).gt.rbub) write(*,*) 'Good coordinates!!!'
        write(*,*)

!----------------------------------------------------------------------
!       Determine the rest bubble centers one by one
!----------------------------------------------------------------------
        DO i=2,nbub
!...Check if the bubble is in mixing vanes region
           if (DomainType .eq. "2x2") then
              call checkmixing(bubcoord(i,1), coordrep, mixvanes, rbub)
           endif

!...Find the minimum distance to the walls
           distw(1)=sqrt((bubcoord(i,3))**2 + (bubcoord(i,2))**2)
           distw(2)=sqrt((bubcoord(i,3)-rbox2(3))**2 +
     &                  (bubcoord(i,2))**2)
           distw(3)=sqrt((bubcoord(i,3))**2 +
     &                  (bubcoord(i,2)-rbox2(2))**2)
           distw(4)=sqrt((bubcoord(i,3)-rbox2(3))**2 +
     &                  (bubcoord(i,2)-rbox2(2))**2)
           minw(i)=distw(1)
           DO j = 2, 4
              IF(minw(i).GE. distw(j)) minw(i)=distw(j)
           END DO
!...Find the minimum distance to other bubbles
           mindistmp = 1e9
           do j=1,i-1  !loop over the centers that have already been found
              distance=sqrt(((bubcoord(i,1)-bubcoord(j,1))**2)+
     &                      ((bubcoord(i,2)-bubcoord(j,2))**2)+
     &                      ((bubcoord(i,3)-bubcoord(j,3))**2))
              if(mindistmp.gt.distance) mindistmp = distance
           end do

           gate = .TRUE.
!...Check the overlap
        if ((mindistmp .LE. mindist) .OR. (minw(i).LE.(r+rbub)) 
     &       .or. coordrep) then
        L1: do

           IF(gate)THEN
              CALL random_number_init(n_dimension,
     &                          int(random1d(i)*1000000))
             CALL RANDOM_NUMBER(random_replace)
              gate=.FALSE.
           ELSE
              CALL random_number_init(n_dimension,
     &                          int(random_replace(2)*1000000))
              CALL RANDOM_NUMBER(random_replace)
           END IF
           DO k=1,3
              bubcoord(i,k)=rbox2(k)*random_replace(k+1) !+ rbox1(2*k-1)
           END DO
!...Check if the bubble is in mixing vanes region
           if (DomainType .eq. "2x2") then
              call checkmixing(bubcoord(i,1), coordrep, mixvanes, rbub)
           endif
!...Find the minimum distance to the walls
           distw(1)=sqrt((bubcoord(i,3))**2 + (bubcoord(i,2))**2)
           distw(2)=sqrt((bubcoord(i,3)-rbox2(3))**2 +
     &                  (bubcoord(i,2))**2)
           distw(3)=sqrt((bubcoord(i,3))**2 +
     &                  (bubcoord(i,2)-rbox2(2))**2)
           distw(4)=sqrt((bubcoord(i,3)-rbox2(3))**2 +
     &                  (bubcoord(i,2)-rbox2(2))**2)
           minw(i)=distw(1)
           DO j = 2, 4
              IF(minw(i).GE. distw(j)) minw(i)=distw(j)
           END DO
!...Find the minimum distance to other bubbles
           mindistmp = 1e9
           do j=1,i-1  !loop over the centers that have already been found
              distance=sqrt(((bubcoord(i,1)-bubcoord(j,1))**2)+
     &                      ((bubcoord(i,2)-bubcoord(j,2))**2)+
     &                      ((bubcoord(i,3)-bubcoord(j,3))**2))
              if(mindistmp.gt.distance) mindistmp = distance
           end do
           if ((mindistmp .LE. mindist).OR.(minw(i).LE.(r+rbub))
     &         .or. coordrep) then
               cycle L1
           endif

        exit L1
        end do L1
        end if

        END DO  !for i loop

!...Coordinates adjustment
!        DO i=1,3
!           rbox2(i)=rbox2(i)+3.0*rbub
!        END DO

!        DO i=1,nbub
!        DO j=1,3
!           bubcoord(i,j)=bubcoord(i,j)+1.5*rbub    
! Why do we do this ?  - to center the smaller box (which was created for the clearance !)
!        END DO
!        END DO
!----------------------------------------------------------------------
!       Void fraction calculation (Not Used)
!----------------------------------------------------------------------
!        xstep=rbox2(1)/real(xres)
!        ystep=rbox2(2)/real(yres)
!        zstep=rbox2(3)/real(zres)
!        ALLOCATE(xvoid(xres),yvoid(yres),zvoid(zres))
!
!        DO i=1,xres
!           xcoord=xcoord+xstep
!           DO j=1,nbub
!              IF((abs(xcoord-bubcoord(j,1))) .LT. rbub)THEN
!              A_bubble=A_bubble+(3.14159*(abs((rbub**2)-
!     &                          ((xcoord-bubcoord(j,1))**2))))
!              END IF
!           END DO
!           xvoid(i)=A_bubble/(rbox2(2)*rbox2(3))
!           A_bubble=0
!        END DO
!
!        DO i=1,yres
!           ycoord=ycoord+ystep
!           DO j=1,nbub
!              IF((abs(ycoord-bubcoord(j,2))) .LT. rbub)THEN
!              A_bubble=A_bubble+(3.14159*(abs((rbub**2)-
!     &                          ((ycoord-bubcoord(j,2))**2))))
!              END IF
!           END DO
!           yvoid(i)=A_bubble/(rbox2(1)*rbox2(3))
!           A_bubble=0
!        END DO
!
!        DO i=1,zres
!           zcoord=zcoord+zstep
!           DO j=1,nbub
!              IF((abs(zcoord-bubcoord(j,3))) .LT. rbub)THEN
!              A_bubble=A_bubble+(3.14159*(abs((rbub**2)-
!     &                          ((zcoord-bubcoord(j,3))**2))))
!              END IF
!           END DO
!           zvoid(i)=A_bubble/(rbox2(1)*rbox2(2))
!           A_bubble=0
!        END DO
!
!----------------------------------------------------------------------
!       Post-processing
!----------------------------------------------------------------------
        DO i=1,nbub   
           DO j=1,3
              IF(j .EQ. 1)THEN
                 bubcoord(i,j)=bubcoord(i,j)+rbox1(1)
              ELSE IF(j .EQ. 2)THEN
                 bubcoord(i,j)=bubcoord(i,j)+rbox1(3)
              ELSE IF(j .EQ. 3)THEN
                 bubcoord(i,j)=bubcoord(i,j)+rbox1(5)
              END IF
           END DO
        END DO

!...Double check the coordinates we got 
        check=.false.
        DO i=1,nbub
           distw(1)=sqrt((bubcoord(i,3)-rbox1(5))**2 + 
     &                  (bubcoord(i,2)-rbox1(3))**2)
           distw(2)=sqrt((bubcoord(i,3)-rbox1(6))**2 + 
     &                  (bubcoord(i,2)-rbox1(3))**2)
           distw(3)=sqrt((bubcoord(i,3)-rbox1(5))**2 + 
     &                  (bubcoord(i,2)-rbox1(4))**2)
           distw(4)=sqrt((bubcoord(i,3)-rbox1(6))**2 + 
     &                  (bubcoord(i,2)-rbox1(4))**2)
           minw(i)=distw(1)
           DO k=2,4
              IF(minw(i).GE. distw(k)) minw(i)=distw(k)
           ENDDO

        IF(minw(i).LE.(r+rbub))THEN
           check=.true.
           WRITE(*,*)
           WRITE(*,*) 'Bubble touches the wall!!!'
           WRITE(*,*) 'The problem point is',bubcoord(i,:)
           WRITE(*,*) 'Distance is',(minw(i)-(r+rbub))
           WRITE(*,*)
        END IF

        DO j=1,nbub
           IF(i.EQ.j)CYCLE
              distance=sqrt(((bubcoord(i,1)-bubcoord(j,1))**2)+
     &                      ((bubcoord(i,2)-bubcoord(j,2))**2)+
     &                      ((bubcoord(i,3)-bubcoord(j,3))**2))
              IF(distance .LE.mindist)THEN
                 check=.true.
                 WRITE(*,*) 'Bubble touches another bubble!!!'
              END IF
           END DO
        ENDDO

        if(check) then
           WRITE(*,*)"Problem was found!!!"
        else
           WRITE(*,*)"So far, so good!!!"
        endif
!----------------------------------------------------------------------
!       Bubble sorting (based on x coordinates) 
!----------------------------------------------------------------------
        size=nbub
        gothroughlist=.TRUE.
        DO WHILE (gothroughlist)
                gothroughlist=.FALSE.
                DO i=1,size-1
                        IF (bubcoord(i,1) .GT. bubcoord(i+1,1)) THEN
                                temp(1,:)=bubcoord(i,:)
                                bubcoord(i,:)=bubcoord(i+1,:)
                                bubcoord(i+1,:)=temp(1,:)
                                gothroughlist=.TRUE.
                        END IF
                END DO
                size=size-1
        END DO

!----------------------------------------------------------------------
!       Bubble duplication for the 2x2 geometry and ghost bubble
!       generation
!----------------------------------------------------------------------
        if (DomainType .eq. "2x2") then
           allocate ( bubdupl(nbub*4,3) )
           do i = 1, nbub
              bubdupl(i,1) = bubcoord(i,1)
              bubdupl(i,2) = bubcoord(i,2) 
              bubdupl(i,3) = bubcoord(i,3)
           enddo
           do i = nbub+1, 2*nbub
              bubdupl(i,1) = bubcoord(i-nbub,1)
              bubdupl(i,2) = bubcoord(i-nbub,2)   + outline(2)
              bubdupl(i,3) = bubcoord(i-nbub,3)
           enddo
           do i = 2*nbub+1, 3*nbub
              bubdupl(i,1) = bubcoord(i-2*nbub,1)
              bubdupl(i,2) = bubcoord(i-2*nbub,2) 
              bubdupl(i,3) = bubcoord(i-2*nbub,3) + outline(3)
           enddo
           do i = 3*nbub+1, 4*nbub
              bubdupl(i,1) = bubcoord(i-3*nbub,1)
              bubdupl(i,2) = bubcoord(i-3*nbub,2) + outline(2)
              bubdupl(i,3) = bubcoord(i-3*nbub,3) + outline(3)
           enddo
!...Enlarge the outline to the full size (we start from a quarter), but
!the x dimension is the same
           outline(2:3) = 2 * outline(2:3)

!...Find the ghost bubbles 
        WRITE(*,"(A, f5.2, A)") 'The size of ghost domain is ',
     &                  ghost_area_ratio*100,'%'
        nghost=0
        DO i=1,4*nbub
           DO j=1,3
             if (j.eq.1) then
              ghostupplim = rbox1(2*j)   -ghost_area_ratio*outline(j)
             else
              ghostupplim = rbox1(2*j-1) -ghost_area_ratio*outline(j)
     &                    + outline(j) 
             endif 

              ghostlowlim = rbox1(2*j-1) +ghost_area_ratio*outline(j)
              IF (bubdupl(i,j) .GT. ghostupplim) then
                  nghost=nghost+1
              ELSEIF (bubdupl(i,j) .LT. ghostlowlim) then
                  nghost=nghost+1
              ENDIF
           END DO
        END DO

        WRITE(*,*) nghost, 'ghost bubbles generated!!!'
        ntotal= 4*nbub + nghost
        ALLOCATE (bublist(ntotal,4))
        bublist(:,:) = 0.0

        DO i=1, 4*nbub
           bublist(i,1) = REAL(i)
           DO j=1,3
              bublist(i,j+1) = bubdupl(i,j)
           ENDDO
        ENDDO

        nghost = 1
        DO i=1,4*nbub
           DO j=1,3
             if (j.eq.1) then
              ghostupplim = rbox1(2*j)   -ghost_area_ratio*outline(j)
             else
              ghostupplim = rbox1(2*j-1) -ghost_area_ratio*outline(j)
     &                    + outline(j)
             endif
              ghostlowlim = rbox1(2*j-1) +ghost_area_ratio*outline(j)
              IF (bubdupl(i,j) .GT. ghostupplim) THEN
                 bublist(4*nbub+nghost,2)=bubdupl(i,1)
                 bublist(4*nbub+nghost,3)=bubdupl(i,2)
                 bublist(4*nbub+nghost,4)=bubdupl(i,3)
                 bublist(4*nbub+nghost,j+1)=bubdupl(i,j)-outline(j)
                 WRITE(6,*) 'Bubble',i,'has the ghost bubble'
                 nghost=nghost+1
              ELSEIF (bubdupl(i,j) .LT. ghostlowlim) THEN
                 bublist(4*nbub+nghost,2)=bubdupl(i,1)
                 bublist(4*nbub+nghost,3)=bubdupl(i,2)
                 bublist(4*nbub+nghost,4)=bubdupl(i,3)
                 bublist(4*nbub+nghost,j+1)=bubdupl(i,j)+outline(j)
                 WRITE(6,*) 'Bubble',i,'has the ghost bubble'
                 nghost=nghost+1
              ENDIF
           ENDDO
        ENDDO

        write(*,*) 'X, Y coordinates are exchanged!!!'

        DO i=1,ntotal
           WRITE(2,1011) int(bublist(i,1)), bublist(i,3), bublist(i,2),
     &                   bublist(i,4), rbub
1011    FORMAT(I8, ES16.6, ES16.6, ES16.6, ES16.6)
        END DO

        elseif (DomainType .eq. "Sub") then
!----------------------------------------------------------------------
!       Ghost bubble generation for subchannel
!----------------------------------------------------------------------
        WRITE(*,"(A, 6(1xES16.6))") 'box size: ',rbox1(:)
        WRITE(*,"(A, 3(1xES16.6))") 'outline: ',outline
        WRITE(*,"(A, f5.2, A)") 'The size of ghost domain is ',
     &                  ghost_area_ratio*100,'%'
        nghost=0
        DO i=1,nbub
           DO j=1,3
              ghostupplim = rbox1(2*j)  -ghost_area_ratio*outline(j)
              ghostlowlim = rbox1(2*j-1)+ghost_area_ratio*outline(j)
              IF (bubcoord(i,j) .GT. ghostupplim) then
                  nghost=nghost+1
              ELSEIF (bubcoord(i,j) .LT. ghostlowlim) then
                  nghost=nghost+1
              ENDIF
           END DO
        END DO

        WRITE(*,*) nghost, 'ghost bubbles generated!!!'
        ntotal=nbub+nghost
        ALLOCATE (bublist(ntotal,4))
        bublist(:,:) = 0.0

        DO i=1,nbub
           bublist(i,1) = REAL(i)
           DO j=1,3
              bublist(i,j+1)=bubcoord(i,j)
           ENDDO
        ENDDO

        nghost = 1
        DO i=1,nbub
           DO j=1,3
              ghostupplim = rbox1(2*j)  -ghost_area_ratio*outline(j)
              ghostlowlim = rbox1(2*j-1)+ghost_area_ratio*outline(j)
              IF (bubcoord(i,j) .GT. ghostupplim) THEN
                 bublist(nbub+nghost,2)=bubcoord(i,1)
                 bublist(nbub+nghost,3)=bubcoord(i,2)
                 bublist(nbub+nghost,4)=bubcoord(i,3)
                 bublist(nbub+nghost,j+1)=bubcoord(i,j)-outline(j)
                 WRITE(6,*) 'Bubble',i,'has the ghost bubble'
                 nghost=nghost+1
              ELSEIF (bubcoord(i,j) .LT. ghostlowlim) THEN
                 bublist(nbub+nghost,2)=bubcoord(i,1)
                 bublist(nbub+nghost,3)=bubcoord(i,2)
                 bublist(nbub+nghost,4)=bubcoord(i,3)
                 bublist(nbub+nghost,j+1)=bubcoord(i,j)+outline(j)
                 WRITE(6,*) 'Bubble',i,'has the ghost bubble'
                 nghost=nghost+1
              ENDIF
           ENDDO
        ENDDO

        DO i=1,ntotal
           WRITE(2,1010) int(bublist(i,1)), bublist(i,2:4), rbub
1010    FORMAT(I8, ES16.6, ES16.6, ES16.6, ES16.6)
        END DO

        endif !DomainType
!-----------------------------------------------------------------------------
!       The following part is for the void fraction calculation (Not Used)
!-----------------------------------------------------------------------------
!WRITE(5,'(//)')
!WRITE(5,*)"Void fraction in the x direction with appropriate x coordinates (x,vf(x))"
!xcoord=0
!DO i=1,xres
!	xcoord=xcoord+xstep
!	WRITE(5,*)xcoord+rbox1(1),xvoid(i)
!END DO
!WRITE(2,'(//)')
!WRITE(2,*)"Void fraction in the y direction with appropriate y coordinates (y,vf(y))"

!ycoord=0
!DO i=1,yres
!	ycoord=ycoord+ystep
!	!WRITE(2,*)ycoord+rbox1(3),yvoid(i)
!END DO

!WRITE(2,'(//)')
!WRITE(2,*)"Void fraction in the z direction with appropriate z coordinates (z,vf(z))"
!zcoord=0
!DO i=1,zres
!	zcoord=zcoord+zstep
!	WRITE(2,*)zcoord+rbox1(5),zvoid(i)
!END DO
!----------------------------------------------------------------------
!       The following part is to write the distance field expression for old
!       bubble initialization used in simmodeller, which is already obsolete and
!       is not used any more. I will just keep it as the evidence of growth of
!       this algorithm.                         Jun Fang, Jan 2015
!----------------------------------------------------------------------

!WRITE(2,'(//)')
!WRITE(2,*)"Distance field equation is below"
!DISTANCE FIELD EQUATION FORMER
!10 FORMAT('sqrt((x1-',F7.4,')^2+')
!20 FORMAT('sqrt(x1+',F7.4,')^2+')
!30 FORMAT('(x2-',F7.4,')^2+')
!40 FORMAT('(x2+',F7.4,')^2+')
!50 FORMAT('(x3-',F7.4,')^2)-',F7.4)
!60 FORMAT('(x3+',F7.4,')^2)-',F7.4)
!70 FORMAT('(x3-',F7.4,')^2)-',F7.4)
!80 FORMAT('(x3+',F7.4,')^2)-',F7.4)
!90 FORMAT(')')
!100 FORMAT(',')
!DO i=1,nbub-1
!	IF(mod(i,15) .EQ. 0 .OR. i .EQ. nbub-1)THEN
!		WRITE(3,'(A)')'min('
!	ELSE
!		WRITE(3,'(A)',ADVANCE='NO')'min('
!	END IF
!END DO

!DO i=1,nbub
!	IF(bubcoord(i,1) .GE. 0)THEN
!		WRITE(3,10,ADVANCE='NO')bubcoord(i,1)
!	ELSE IF(bubcoord(i,1) .LT. 0)THEN
!		WRITE(3,20,ADVANCE='NO')-bubcoord(i,1)
!	END IF
!	IF(bubcoord(i,2) .GE. 0)THEN
!		WRITE(3,30,ADVANCE='NO')bubcoord(i,2)
!	ELSE IF(bubcoord(i,2) .LT. 0)THEN
!		WRITE(3,40,ADVANCE='NO')-bubcoord(i,2)
!	END IF
!	IF(i .LT. nbub)THEN
!		IF(bubcoord(i,3) .GE. 0)THEN
!			WRITE(3,50,ADVANCE='no')bubcoord(i,3),rbub
!		ELSE IF(bubcoord(i,3) .LT. 0)THEN
!			WRITE(3,60,ADVANCE='no')-bubcoord(i,3),rbub
!		END IF
!	ELSE IF(i .EQ. nbub)THEN
!		IF(bubcoord(i,3) .GE. 0)THEN
!			WRITE(3,70,ADVANCE='no')bubcoord(i,3),rbub
!		ELSE IF(bubcoord(i,3) .LT. 0)THEN
!			WRITE(3,80,ADVANCE='no')-bubcoord(i,3),rbub
!		END IF
!	END IF
!	IF(i .EQ. 1)THEN
!		WRITE(3,100,ADVANCE='NO')
!	ELSE IF(mod(i,2) .EQ. 0 .AND. i .LT. nbub)THEN
!		WRITE(3,90,ADVANCE='NO')
!		WRITE(3,100)
!	ELSE IF(mod(i,2) .EQ. 1 .AND. i .LT. nbub)THEN
!		WRITE(3,90,ADVANCE='NO')
!		WRITE(3,100,ADVANCE='NO')
!	ELSE IF(i .EQ. nbub)THEN
!		WRITE(3,90)
!	END IF
!END DO

      END PROGRAM


!----------------------------------------------------------------------
!       The subroutine for random seed 
!----------------------------------------------------------------------

      SUBROUTINE random_number_init(n,input)
        IMPLICIT NONE
        INTEGER::i,n,input
        INTEGER, ALLOCATABLE::seed(:)

        CALL RANDOM_SEED(size=n)
        ALLOCATE(seed(n))

        seed=input+37*(/(i-1,i=1,n)/)
        CALL RANDOM_SEED(PUT=seed)

        DEALLOCATE(seed)
      END SUBROUTINE


!----------------------------------------------------------------------
!       The subroutine to check if a bubble is in mixing vanes or not
!----------------------------------------------------------------------
      subroutine checkmixing(xpostn, coordrep, mixvanes, rbub)
        implicit none
        logical :: coordrep
        real*8     xpostn, mixvanes(2)
        real*8     rbub

        if(xpostn.ge.(mixvanes(1)-rbub) .and.
     &     xpostn.le.(mixvanes(2)+rbub)) then
           coordrep = .True.
        else
!           write(*,*) xpostn, mixvanes
           coordrep = .False.
        endif
        return
      end subroutine




