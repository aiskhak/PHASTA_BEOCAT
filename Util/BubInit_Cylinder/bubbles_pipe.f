        PROGRAM bubbles
        IMPLICIT NONE
c!----------------------------------------------------------------------
c
c! This routine is developed to created the required bubble information
c! such as, bubble ID, bubble coordiantes and radii. And the Monta Carlo
c! rejection sampling method is used to avoid the intersection between
c! bubbles and the walls or between different bubbles 
c!
c! key variables:
c! nbub         number of bubbles
c! rbub         radius of bubbles (constant for now)
c! rbox         dimensions of box
c!
c!
c! Jun Fang, 2014. 
c!----------------------------------------------------------------------
        INTEGER::nbub,ierror,i,j,k,clock,n_dimension=4,size
        REAL*8  rbub,rbox1(6),rbox2(3),distance,mindist,
     &          random_replace(4),temp(1,3),distw(4),l,w,h,r
        real*8, ALLOCATABLE:: random1d(:),   random2d(:,:)
        real*8, allocatable:: bubcoord(:,:), minw(:)
        LOGICAL::gate=.TRUE.,check,gothroughlist

        INTEGER::nghost,ntotal
        REAL*8  ghost_area_ratio,outline(3),bubtmp(3)
        REAL*8, ALLOCATABLE:: bublist(:,:)
        real*8 Pipe_R

        OPEN(UNIT=1,FILE='options.dat',IOSTAT=ierror)
        IF(ierror .NE. 0)STOP"Error opening file options.txt"

        OPEN(UNIT=2,FILE='bubbles.inp',IOSTAT=ierror)
        IF(ierror .NE. 0)STOP"Error opening file output.txt"

        OPEN(UNIT=3,FILE='rawbubbles.inp',IOSTAT=ierror)
        IF(ierror .NE. 0)STOP"Error opening file output.txt"

!----------------------------------------------------------------------
!       Read from options.txt
!       rbox1(:) is in the following format:
!       (x_min,x_max,y_min,y_max,z_min,z_max)
!----------------------------------------------------------------------
        READ(1,*)
        READ(1,*)nbub
        READ(1,*)
        READ(1,*)rbub
        READ(1,*)
        READ(1,*)rbox1(:)
        READ(1,*)
        READ(1,*)r
        READ(1,*)
        READ(1,*)mindist
        READ(1,*)
        READ(1,*)ghost_area_ratio
!----------------------------------------------------------------------
!       Initialization
!----------------------------------------------------------------------
        DO i=2,6,2
           rbox2(i/2) = rbox1(i) - rbox1(i-1)
        END DO
        ghost_area_ratio = ghost_area_ratio/100.0
        outline(:)=rbox2(:)

        ALLOCATE( random1d(nbub*3), random2d(nbub,3) )
        allocate( bubcoord(nbub,3), minw(nbub) )
!...Forming the random array...
        CALL SYSTEM_CLOCK(COUNT=clock)
        CALL random_number_init(nbub*3,clock)
        CALL RANDOM_NUMBER(random1d)
!Putting the random array into a two dimensional array, making it 
!easier to use with the bubble coordinate array
        DO i=1,nbub
           random2d(i,:)=random1d((i*3-2):i*3)
        END DO
!forming the bubble coordinate array
        write(*,*)
        DO i=1,nbub
           do j=1,3
              bubcoord(i,j)=random2d(i,j)*rbox2(j)+rbox1(2*j-1)   
           enddo
        END DO

!... Calculate the distance b/w first bubble center to the pipe wall...
        Pipe_R = 0.1
        minw(1)= Pipe_R - 
     &            sqrt((bubcoord(1,3))**2 + (bubcoord(1,2))**2)

!... check if the first bubble is in the right position 
        IF( minw(1) .LE. rbub ) THEN
            L0:DO
            IF(gate)THEN

        CALL random_number_init(n_dimension,
     &          int(random1d(1)*1.0E05))
        CALL RANDOM_NUMBER(random_replace)
               gate=.FALSE.
            ELSE

        CALL random_number_init( n_dimension, 
     &          int(random_replace(2)*1.0E05))
        CALL RANDOM_NUMBER(random_replace)

            END IF

            DO k=1,3
               bubcoord(1,k)=rbox2(k)*random_replace(k+1)+rbox1(2*k-1)
            END DO

            minw(1)= Pipe_R -
     &            sqrt((bubcoord(1,3))**2 + (bubcoord(1,2))**2)
            IF( minw(1).LE.rbub )THEN
                CYCLE L0
            END IF
            EXIT L0
            END DO L0
        ENDIF
        write(*,*) 'The first center: ', bubcoord(1,:)
        if(minw(1).gt.rbub) write(*,*) 'First center is good!'
        write(*,*) 

        DO i=2,nbub
           gate = .TRUE.
           DO j=1,i-1  !loop over the centers that have been found
              minw(i)= Pipe_R -
     &            sqrt((bubcoord(i,3))**2 + (bubcoord(i,2))**2)
              distance = sqrt(((bubcoord(i,1)-bubcoord(j,1))**2)
     &                  +((bubcoord(i,2)-bubcoord(j,2))**2) 
     &                  +((bubcoord(i,3)-bubcoord(j,3))**2))
              IF((distance.LE.mindist).OR.
     &                  (minw(i).LE.rbub)) THEN


                L1:DO
                IF(gate)THEN
                  CALL random_number_init(n_dimension,
     &                  int(random1d(i)*1000000))
                  CALL RANDOM_NUMBER(random_replace)
                  gate=.FALSE.
                ELSE
                  CALL random_number_init(n_dimension,
     &                  int(random_replace(2)*1000000))
                  CALL RANDOM_NUMBER(random_replace)
                END IF

                DO k=1,3
                   bubcoord(i,k)=rbox2(k)*random_replace(k+1)
     &                  +rbox1(2*k-1)
                END DO
! Check the wall here as well !
                minw(i)= Pipe_R -
     &            sqrt((bubcoord(i,3))**2 + (bubcoord(i,2))**2)
                DO k=1,i-1
                   distance=sqrt(((bubcoord(i,1)-bubcoord(k,1))**2)
     &                  +((bubcoord(i,2)-bubcoord(k,2))**2)
     &                  +((bubcoord(i,3)-bubcoord(k,3))**2))
                   IF((distance.LE.mindist).OR.
     &                   (minw(i).LE.rbub))THEN
                       CYCLE L1


                   END IF
                END DO
                EXIT L1
                END DO L1
             END IF !for the first check

           END DO !for j loop
           if(minw(i).le.rbub) write(*,*) 'The point is bad!'
        END DO  !for i loop

!----------------------------------------------------------------------
!      Double-check the bubblen centers  
!----------------------------------------------------------------------
        check=.false.
        DO i=1,nbub
           minw(i)= Pipe_R -
     &            sqrt((bubcoord(i,3))**2 + (bubcoord(i,2))**2)
           IF(minw(i).LE.rbub)THEN
                check=.true.
                WRITE(*,*)
                WRITE(*,*) 'Problem:  check min dist. to walls'
                WRITE(*,*) 'The problem point is',bubcoord(i,:)
                WRITE(*,*) 'Distance is',minw(i)
                WRITE(*,*)
           END IF

           DO j=1,nbub
              IF(i.EQ.j)CYCLE
              distance = sqrt(((bubcoord(i,1)-bubcoord(j,1))**2)+
     &                  ((bubcoord(i,2)-bubcoord(j,2))**2)+
     &                  ((bubcoord(i,3)-bubcoord(j,3))**2))
              IF(distance .LE.mindist)THEN
                 check=.true.
                 WRITE(*,*) 'Problem: check min dist b/w bubbles'
              END IF
           END DO
        ENDDO

        if(check)THEN
        WRITE(*,*)"There is a problem"
        else
        WRITE(*,*)"No problem found w/ random generator or replac."
        END IF

!-----------------------------------------------------------------------------
!       Find Ghost bubbles
!-----------------------------------------------------------------------------
        WRITE(*,*) 'box size',rbox1(:)
        WRITE(*,*) 'outline is',outline
        WRITE(*,*) 'Ghost domain Size: ',ghost_area_ratio*100,'%'
        nghost=0
        DO i=1,nbub
           DO j=1,1
! Only streamwise ghosts this time !
              IF (bubcoord(i,j).GT.(rbox1(2*j)
     &                  - ghost_area_ratio*outline(j))) then
                  nghost=nghost+1
        !      WRITE(*,*) 'nghost is ', nghost
              ELSEIF (bubcoord(i,j).LT.(rbox1(2*j-1)
     &                  + ghost_area_ratio*outline(j))) then
                  nghost=nghost+1
        !      WRITE(*,*) 'nghost is ', nghost
              ENDIF
           END DO
        END DO

        WRITE(*,*) nghost, 'ghost bubbles will be generated !'
        ntotal=nbub+nghost
        ALLOCATE (bublist(ntotal,4))
        bublist(:,:) = 0.0
!-----------------------------------------------------------------------------
!Sorting according to x
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
!-----------------------------------------------------------------------------

        do i = 1,nbub
           write(3,1011) bubcoord(i,:)
        enddo
1011    FORMAT(ES16.6, ES16.6, ES16.6)

        DO i=1,nbub
           bublist(i,1) = REAL(i)
           bublist(i,2:4) = bubcoord(i,:)
        ENDDO

        nghost = 1
        DO i=1,nbub
           DO j=1,1
              IF (bubcoord(i,j) .GT. (rbox1(2*j) 
     &                  - ghost_area_ratio*outline(j))) THEN
                bublist(nbub+nghost,2)=bubcoord(i,1)
                 bublist(nbub+nghost,3)=bubcoord(i,2)
                 bublist(nbub+nghost,4)=bubcoord(i,3)
                 bublist(nbub+nghost,j+1)=bubcoord(i,j)-outline(j)
!                 WRITE(*,*) 'Bubble',i,'has the ghost bubble'
                 nghost=nghost+1
              ELSEIF (bubcoord(i,j) .LT. (rbox1(2*j-1)
     &                  + ghost_area_ratio*outline(j))) THEN
                 bublist(nbub+nghost,2)=bubcoord(i,1)
                 bublist(nbub+nghost,3)=bubcoord(i,2)
                 bublist(nbub+nghost,4)=bubcoord(i,3)
                 bublist(nbub+nghost,j+1)=bubcoord(i,j)+outline(j)
!                 WRITE(*,*) 'Bubble',i,'has the ghost bubble'
                 nghost=nghost+1
             ENDIF
           ENDDO
        ENDDO

        DO i=1,ntotal
                WRITE(2,1010) int(bublist(i,1)), bublist(i,2:4), rbub
        END DO
1010    FORMAT(I8, ES16.6, ES16.6, ES16.6, ES16.6)
!-----------------------------------------------------------------------------

        END PROGRAM

        SUBROUTINE random_number_init(n,input)
        IMPLICIT NONE
        INTEGER*8::i,n,input
        INTEGER*8, ALLOCATABLE::seed(:)

        CALL RANDOM_SEED(size=n)
        ALLOCATE(seed(n))

        seed=input+13*(/(i-1,i=1,n)/)
        CALL RANDOM_SEED(PUT=seed)

        DEALLOCATE(seed)
        END SUBROUTINE
