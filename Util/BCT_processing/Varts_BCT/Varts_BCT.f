! Convertion code developed for generating the turbulent inflow conditions using PHASTA tools
! Igor A. Bolotnov, April 2010
	program Varts_to_BCT

	implicit none

	character*1 MyChar1
	external MyChar1
        character*2 MyChar2
        external MyChar2
        character*3 MyChar3
        external MyChar3
        character*4 MyChar4
        external MyChar4
        character*5 MyChar5
        external MyChar5
        character*6 MyChar6
        external MyChar6

	character*80 ipath
        character*80 filenum

	integer, parameter :: mp = 20000

	real tol, varts(1:20), ctime, tstart, tstop
	real xx(1:mp), yy(1:mp), zz(1:mp)
	real deltat, phdt

	integer i1, i2, i, j, k, maxts(1:100), ts1(1:100), ts2(1:100)
	integer Nrun, Nfiles, np, nskip, nd1, nd2 
	integer reclength, lstep, jj, icount, iphase
	integer istep(1:100), istart, istop
	integer fout, fin, isn, iread
	integer nbct, mbct, ioverlap
	
! read the current case path:
	open(101, file = 'path.dat')
	 read(101, 50) ipath
   50    format(A80)
	close(101)

	write(*,*) 'Processing the case located in ', trim(ipath)

! Read the input data:
	open(1, file = trim(ipath)//'/vbct.inp')

	do i = 1, 6		! skip the header
	  read(1,*)
	end do

        read(1,*) tstart		! Range of time in the result
	read(1,*) tstop
	read(1,*)
	read(1,*) deltat
        read(1,*)
	read(1,*) phdt
	read(1,*)
        read(1,*) istart                ! Range of timesteps in the source varts file
        read(1,*) istop
        read(1,*)
        read(1,*) ioverlap
	
	close(1)
	
	j = 1
	istep(j) = istart

! Read the point data:
	open(2, file = trim(ipath)//'/xyzts.dat')
	
	read(2, *) np, nskip, tol, nd1, nd2, NRun
! Read the point coordinates:
	 do i = 1, np
	  read(2, *) xx(i), yy(i), zz(i)
	 end do
	close(2)

	write(*,*) 'Number of points per step = ', np

	write(*,*) 'Processing run #', Nrun

	reclength = 1*8+4*15		! record length

	Nbct = int(1.0E0*(istop - istart)/(deltat/phdt))+1
	write(*,*) 'Anticipated number of bct.dat files: ', Nbct

! loop over files:
	icount = 0
	 write(*,*) 'Processing the run number: ', Nrun 
          open(3, file = trim(ipath)//'/varts_run'//MyChar2(Nrun)//'.dat'
     1     , status='unknown', form='formatted', recl=reclength, access='direct')

! read/write the information
         iread = 0
	 ts1(1) = 1
	 do i2 = 1, nbct
	 if (i2.lt.nbct) then
		maxts(i2) = nint(deltat/phdt)
	  else
		maxts(i2) = mod(istop - istart, maxts(i2-1))
	 end if
	 ts2(i2) = ts1(i2) + maxts(i2) - 1
! Expand the timestep range to provide slight overlap:
	 if (i2.eq.1) then
	   ts2(i2) = ts2(i2) + ioverlap
	 else if (i2.eq.nbct) then
	   ts1(i2) = ts1(i2) - ioverlap
	 else
           ts1(i2) = ts1(i2) - ioverlap
           ts2(i2) = ts2(i2) + ioverlap
	 end if
	 write(*,*) 'File ', i2, ' range is timesteps: ', ts1(i2), ts2(i2)
	 ts1(i2+1) = ts1(i2) + maxts(i2)
	 end do
        ctime = 0.0
! Loop over bct.dat files:
	do mbct = 1, nbct
	 write(*,*) 'Processing file #', mbct, ' of ', nbct
                 if( mbct.ge.0 .and. mbct.le.9)then
                    write(filenum,'(i1.1)')mbct
                 else if(mbct.ge.10 .and. mbct.le.99)then
                    write(filenum,'(i2.2)')mbct
                 else if(mbct.ge.100 .and. mbct.le.999)then
                    write(filenum,'(i3.3)')mbct
                 else if(mbct.ge.1000 .and. mbct.le.9999)then
                    write(filenum,'(i4.4)')mbct
                 end if
         open(20+mbct, file = trim(ipath)//'/bct.'//trim(filenum)//'.dat'
     1     , status='unknown')
	 write(20+mbct, 13) np, ts2(mbct) - ts1(mbct) + 1   !, maxts(mbct) 
! Loop over points:
     	do i1 = 1, np
         ctime = real(ts1(mbct)-1)*phdt
	!if (mod(i, 50).eq.0)write(*,*) 'Point ', i1, ' out of ', np
	 do i = ts1(mbct), ts2(mbct)   ! Loop over number of steps   
	  read(3, 11, REC=np*(i-1)+i1)
     1            lstep,(varts(k), k=2, 5)
! Record the point coordinates and the number of time steps
	  if (i.eq.ts1(mbct)) then
          write(20+mbct, 12) xx(i1), yy(i1), zz(i1), ts2(mbct) - ts1(mbct) + 1   !, maxts(mbct)
	  end if
! Accumulate time:
	  if (i.gt.1) ctime = ctime + varts(5)
	  fout = 0
	  if (ctime.ge.tstart.and.ctime.lt.tstop) fout = 1
	  if (fout) then
	  icount = icount + 1 
	  write(20+mbct,10) varts(2:4), ctime  !, i1, i, lstep, jj, iphase 
	  end if

	 end do   ! i, isteps
         end do ! i1, np
	 close(20+mbct)
	end do ! mbct, nbct

	 close(3)

	close(20)

10	format(1x, 4E15.7, 10I10)
11	format(1x, I8, 4E15.7)
12	format(1x, 3E15.7, I10)
13      format(1x, 2I10)
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
