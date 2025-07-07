	program Convert_Inflow 

	implicit none

	character*80 ipath

	real tol, varts(1:20)

	integer, parameter :: maxp = 100000
	real*8, parameter  :: eps = 1.0D-010
        real*8, parameter  :: epsd = 1.0D-30   ! Distance SQUARED !
	real*8, parameter :: pi = 3.1415926535	

	logical iflag

	integer i, j, k, myrank, stat
	integer iostatus, ntv, ntime
        
        character*50 :: filenum
        character*50 :: filenum1
        character*50 :: buffer
	real*8 x, y, z, w, t, r
	real*8 time1, time2
	real*8 x1(1:maxp), y1(1:maxp), z1(1:maxp)
        
       ! open(1,file='myrank.dat')
       ! NR1=0
       ! do while(.true.)
       !    read(1,fmt="(A79)",iostat=stat) buffer
        !   if (stat==0) then 
        !       NR1=NR1+1
         !    else if(stat>0) then
        !       write(*,*) 'something is wrong with Reading operation'
        !     else
        !       exit
        !     endif
       ! enddo
        !close(1)
       ! write(*,*) NR1
! Read the input data:
        open(1,file='myrank.txt')
        do while(.true.)
           read(1,*,iostat=stat) myrank
           if(stat==0) then
                 if( myrank.ge.0 .and. myrank.le.9)then
                    write(filenum1,'(i1.1)')myrank
                 else if(myrank.ge.10 .and. myrank.le.99)then
                    write(filenum1,'(i2.2)')myrank
                 else if(myrank.ge.100 .and. myrank.le.999)then
                    write(filenum1,'(i3.3)')myrank
                 else if(myrank.ge.1000 .and. myrank.le.9999)then
                    write(filenum1,'(i4.4)')myrank
                 else if(myrank.ge.10000 .and. myrank.le.99999)then
                    write(filenum1,'(i5.5)')myrank
                 else if(myrank.ge.100000 .and. myrank.le.999999)then
                    write(filenum1,'(i6.6)')myrank
                 else if(myrank.ge.1000000 .and. myrank.le.9999999)then
                    write(filenum1,'(i7.7)')myrank
                 else if(myrank.ge.10000000 .and. myrank.le.99999999)then
                    write(filenum1,'(i8.8)')myrank
                 else if(myrank.ge.100000000 .and. myrank.le.999999999)then
                    write(filenum1,'(i9.9)')myrank
                 else if(myrank.ge.1000000000 .and. myrank.le.9999999999)then
                    write(filenum1,'(i10.10)')myrank
                 end if
        
	open(1+myrank, file = 'bct.dat.'//trim(filenum1)//'')
	i = 0

! Read the coordinates from the PHASTA-produced file:
	do 
	 read(1+myrank, *, IOSTAT = iostatus) myrank, x, y, z
	 if (iostatus.ne.0) then
	   goto 23
	 else
! Here we will first check the new point has been already recorded:
	   iflag = 1
	   do j = 1, i
!	    if (abs(x1(j)-x).lt.eps.and.abs(y1(j)-y).lt.eps.and.abs(z1(j)-z).lt.eps) then 
            if (((x1(j)-x)**2.0D0 + (y1(j)-y)**2.0D0 + (z1(j)-z)**2.0D0).lt.epsd) then
		iflag = 0
		goto 22
	    end if ! abs
	   end do

	   if (iflag) then 
             i = i + 1
	     x1(i) = x
	     y1(i) = y
	     z1(i) = z
	   end if ! iflag
   22   continue
	 end if ! iostatus

	end do
   23   continue

	write(*,*) i, 'unique points loaded from the original PHASTA file'
       
         if( myrank.ge.0 .and. myrank.le.9)then
                    write(filenum,'(i1.1)')myrank
                 else if(myrank.ge.10 .and. myrank.le.99)then
                    write(filenum,'(i2.2)')myrank
                 else if(myrank.ge.100 .and. myrank.le.999)then
                    write(filenum,'(i3.3)')myrank
                 else if(myrank.ge.1000 .and. myrank.le.9999)then
                    write(filenum,'(i4.4)')myrank
                 else if(myrank.ge.10000 .and. myrank.le.99999)then
                    write(filenum,'(i5.5)')myrank
                 else if(myrank.ge.100000 .and. myrank.le.999999)then
                    write(filenum,'(i6.6)')myrank
                 else if(myrank.ge.1000000 .and. myrank.le.9999999)then
                    write(filenum,'(i7.7)')myrank
                 else if(myrank.ge.10000000 .and. myrank.le.99999999)then
                    write(filenum,'(i8.8)')myrank
                 else if(myrank.ge.100000000 .and. myrank.le.999999999)then
                    write(filenum,'(i9.9)')myrank
                 else if(myrank.ge.1000000000 .and. myrank.le.9999999999)then
                    write(filenum,'(i10.10)')myrank
                 end if


        open(200+myrank, file = 'bctin.dat.'//trim(filenum)//'',status='new')
         write(200+myrank, *) i
          do j = 1, i
!          write(2, 10) x1(j), y1(j), z1(j)
           write(200+myrank, 10) 0.0, y1(j), z1(j)
          end do
        close(200+myrank)
10      format(10E16.7)


	close(1+myrank)
        else if(stat>0) then
          write(*,*) 'something is wrong with reading operation'
        else 
          exit
        endif
        enddo
        close(1)

! Create a short bctin.dat file with the unique points:
!	open(2+filenum, file = 'bctin.dat')	
!	 write(2+filenum, *) i
!	  do j = 1, i
!	   write(2, 10) x1(j), y1(j), z1(j)
!           write(2, 10) 0.0, y1(j), z1(j)
!	  end do
!	close(2)	
!10	format(10E16.7)

	open(3, file = 'bct.dat')

	ntv = i

! The next part of the code is generating the time dependent sin wave as an inflow BC for PHASTA
! (to test the procedure)
	ntime = 31   ! Number of time steps 
	time1 = 0.0D0 ! Begin time	
	time2 = 0.3D0 ! End time
	
	write(3,*) ntv, ntime
	do i = 1, ntv
	  write(3,11) x1(i), y1(i), z1(i), ntime
	  do j = 1, ntime
! Compute t and w here:
	   t = time1 + real(j-1)/real(ntime-1)*(time2-time1)
	   r = sqrt(x1(i)**2.0+y1(i)**2.0)
	   w = max((1.0**2.0-(r/4.999E-04)**2.0)*0.5*(1.1+cos(2.0*pi*t/time2)), 0.0)
	   write(3,12) 0.0, 0.0, w, t
	  end do
	end do  ! i, ntv
11 	format (3E16.7, I7)
12	format (4x, 4E16.7)
	close(3)

	end program  


	character*(1) function MYCHAR1(i)
	integer i
	 MYCHAR1(1:1) = ACHAR(ICHAR('0')+i)
	end function

	character*(2) function MYCHAR2(i)
	integer i
	 MYCHAR2(1:1) = ACHAR(ICHAR('0')+INT(i/10))
	 MYCHAR2(2:2) = ACHAR(ICHAR('0')+i-10*INT(i/10))
	end function

	character*(3) function MYCHAR3(i)
	integer i
	 MYCHAR3(1:1) = ACHAR(ICHAR('0')+INT(i/100))
	 MYCHAR3(2:2) = ACHAR(ICHAR('0')+INT(i/10)-10*INT(i/100))
	 MYCHAR3(3:3) = ACHAR(ICHAR('0')+i-10*INT(i/10))
	end function

	character*(4) function MYCHAR4(i)
	integer i
	 MYCHAR4(1:1) = ACHAR(ICHAR('0')+INT(i/1000))
	 MYCHAR4(2:2) = ACHAR(ICHAR('0')+INT(i/100)-10*INT(i/1000))
	 MYCHAR4(3:3) = ACHAR(ICHAR('0')+INT(i/10)-10*INT(i/100))
	 MYCHAR4(4:4) = ACHAR(ICHAR('0')+i-10*INT(i/10))
	end function

	character*(5) function MYCHAR5(i)
	integer i
	 MYCHAR5(1:1) = ACHAR(ICHAR('0')+INT(i/10000))
	 MYCHAR5(2:2) = ACHAR(ICHAR('0')+INT(i/1000)-10*INT(i/10000))
	 MYCHAR5(3:3) = ACHAR(ICHAR('0')+INT(i/100)-10*INT(i/1000))
	 MYCHAR5(4:4) = ACHAR(ICHAR('0')+INT(i/10)-10*INT(i/100))
	 MYCHAR5(5:5) = ACHAR(ICHAR('0')+i-10*INT(i/10))
	end function

	character*(6) function MYCHAR6(i)
	integer i
	 MYCHAR6(1:1) = ACHAR(ICHAR('0')+INT(i/100000))
	 MYCHAR6(2:2) = ACHAR(ICHAR('0')+INT(i/10000)-10*INT(i/100000))
	 MYCHAR6(3:3) = ACHAR(ICHAR('0')+INT(i/1000)-10*INT(i/10000))
	 MYCHAR6(4:4) = ACHAR(ICHAR('0')+INT(i/100)-10*INT(i/1000))
	 MYCHAR6(5:5) = ACHAR(ICHAR('0')+INT(i/10)-10*INT(i/100))
	 MYCHAR6(6:6) = ACHAR(ICHAR('0')+i-10*INT(i/10))
	end function
