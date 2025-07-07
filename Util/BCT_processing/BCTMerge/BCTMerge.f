! This is the BCT version of the merge code which prints out ONLY the velocity components to save space and time for the Varts_BCT code
! Created 05/05/2010
	program Merge1 

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

	real tol, varts(1:20)

	integer i1, i, j, k
	integer Nrun, Nfiles, np, nskip, nd1, nd2 
	integer reclength, lstep, jj, icount, iphase
	integer istep(1:100), istart, istop
	integer fout, fin, isn, iread, recl2

! read the current case path:
	open(101, file = 'path.dat')
	 read(101, 50) ipath
   50    format(A80)
	close(101)

	write(*,*) 'Processing the case located in ', trim(ipath)

! Read the input data:
	open(1, file = trim(ipath)//'/merge.inp')

	do i = 1, 6		! skip the header
	  read(1,*)
	end do

	read(1,*) Nrun
	read(1,*)
        read(1,*) istart		! Range of timesamples in the result
	read(1,*) istop
        read(1,*)
        read(1,*) Nfiles
	read(1,*)
	do i = 1, Nfiles+1
	 read(1,*) istep(i)
	end do

	close(1)

! Read the point data:
	open(2, file = trim(ipath)//'/xyzts.dat')
	
	read(2, *) np, nskip, tol, nd1, nd2
	close(2)

	write(*,*) 'Number of points per step = ', np

	write(*,*) 'Processing run #', Nrun

	reclength = 2*8+3+15*15		! record length
	recl2 = 1*8 + 4*15  ! merged file rec length
          open(20, file = trim(ipath)//'/varts_run'//MyChar2(Nrun)//'.dat'
     1     , status='unknown', form='formatted', recl=recl2, access='direct')

! loop over files:
	icount = 0
	do j = 1, Nfiles
	 write(*,*) 'Processing the step number: ', istep(j) 
	if (Nrun.gt.9) then
        if (istep(j).lt.1000) then
          open(3, file = trim(ipath)//'/varts.'//MyChar3(istep(j))//'.run.'//MyChar2(Nrun)//'.dat'
     1     , status='unknown', form='formatted', recl=reclength, access='direct')
	 else if (istep(j).lt.10000) then
          open(3, file = trim(ipath)//'/varts.'//MyChar4(istep(j))//'.run.'//MyChar2(Nrun)//'.dat'
     1     , status='unknown', form='formatted', recl=reclength, access='direct')
	 else if (istep(j).lt.100000) then
          open(3, file = trim(ipath)//'/varts.'//MyChar5(istep(j))//'.run.'//MyChar2(Nrun)//'.dat'
     1     , status='unknown', form='formatted', recl=reclength, access='direct')
	 else
          open(3, file = trim(ipath)//'/varts.'//MyChar6(istep(j))//'.run.'//MyChar2(Nrun)//'.dat'
     1     , status='unknown', form='formatted', recl=reclength, access='direct')
	 end if
	else
        if (istep(j).lt.1000) then
          open(3, file = trim(ipath)//'/varts.'//MyChar3(istep(j))//'.run.'//MyChar1(Nrun)//'.dat'
     1     , status='unknown', form='formatted', recl=reclength, access='direct')
         else if (istep(j).lt.10000) then
          open(3, file = trim(ipath)//'/varts.'//MyChar4(istep(j))//'.run.'//MyChar1(Nrun)//'.dat'
     1     , status='unknown', form='formatted', recl=reclength, access='direct')
         else if (istep(j).lt.100000) then
          open(3, file = trim(ipath)//'/varts.'//MyChar5(istep(j))//'.run.'//MyChar1(Nrun)//'.dat'
     1     , status='unknown', form='formatted', recl=reclength, access='direct')
         else
          open(3, file = trim(ipath)//'/varts.'//MyChar6(istep(j))//'.run.'//MyChar1(Nrun)//'.dat'
     1     , status='unknown', form='formatted', recl=reclength, access='direct')
         end if
	end if

! read/write the information
         iread = 0
	 do i = 1, istep(j+1)-istep(j)   ! Loop over number of steps   !/nskip
          isn = istep(j) + i - 1
! We check if the step is dividable by Nskip:
	  fin = 0
	  if (mod(isn,NSkip).eq.0) fin = 1	  
	  if (fin) then
	  iread = iread + 1
!	 write(*,*) 'start reading time step number : ', isn
	   do i1 = 1, np
	  read(3, '(2I8, I3, 15E15.7)', REC=np*(iread-1)+i1)
     1            lstep,jj,iphase,(varts(k), k=1, 15)
!	  if (i1.eq.np) write(*,*) i, lstep, isn, jj, iphase
	  fout = 0
	  if (lstep.ge.istart.and.lstep.lt.istop) fout = 1
!	  if (i.eq.100) write(*,*) varts(1:15)
	  if (fout) then
	  icount = icount + 1 
	  write(20,'(1I8, 4E15.7)', REC=icount)
     1            lstep,(varts(k), k=2, 5) 
	  end if

	  end do ! i1, np
	  end if ! fin
	 end do   ! i, isteps
	 close(3)
	end do  ! j, Nfiles

	close(20)

	end program Merge1 


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
