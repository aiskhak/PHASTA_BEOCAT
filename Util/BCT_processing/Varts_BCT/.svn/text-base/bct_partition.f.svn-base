      program bct_partition

      implicit none
     
      character*80 filenum1
      character*80 filenum
      character*80 filenum2
      character*80 ipath

      integer np, nts, myrank, stat, np_bct, i, j, k, m, l, p
      integer, parameter :: maxp=10000000
      integer, parameter :: bctnum=10 !number of bct files for different time window

      real*8 u(1:maxp), v(1:maxp), w(1:maxp), t(1:maxp)
      real*8 x(1:maxp), y(1:maxp), z(1:maxp)
      real*8 x_p, y_p, z_p
      real*8, parameter :: epsd=1E-7

      open(30000,file='path.dat')
         read(30000,50) ipath
50       format(A80)
      close(30000)

      do p=1, bctnum
                 if( p.ge.0 .and. p.le.9)then
                    write(filenum2,'(i1.1)')p
                 else if(p.ge.10 .and.p .le.99)then
                    write(filenum2,'(i2.2)')p
                 else if(p.ge.100 .and.p.le.999)then
                    write(filenum2,'(i3.3)')p
                 else if(p.ge.1000 .and.p.le.9999)then
                    write(filenum2,'(i4.4)')p
                 end if
      open(p, file=trim(ipath)//'/bct.'//trim(filenum2)//'.dat')
      read(p,*) np, nts
      do m=1, np*(nts+1) !lopp over the points
         read(p,*) u(m), v(m), w(m), t(m)
      enddo
      close(p)

      open(20000,file=trim(ipath)//'/myrank.txt')
      do while(.true.)
         read(20000,*,iostat=stat) myrank
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
                 else if(myrank.ge.10000000 .and.myrank.le.99999999)then
                    write(filenum1,'(i8.8)')myrank
                 else if(myrank.ge.100000000 .and.myrank.le.999999999)then
                    write(filenum1,'(i9.9)')myrank
                 else if(myrank.ge.1000000000 .and.myrank.le.9999999999)then
                    write(filenum1,'(i10.10)')myrank
                 end if

        open(bctnum+2+myrank, file = trim(ipath)//'/bctin.dat.'//trim(filenum1)//'')
        read(bctnum+2+myrank, *) np_bct
        do l=1, np_bct
           read(bctnum+2+myrank, *) x(l), y(l), z(l)
        enddo
        close(bctnum+2+myrank)
         
        open(1000+myrank,file=trim(ipath)//'/bct_input/bct.dat.'//trim(filenum2)//'.'
     & //trim(filenum1)//'', status='unknown')
        write(*,*) 'Processing processor', myrank,' of time series', p
           write(1000+myrank,13) np_bct, nts
        do i=1, np_bct
           x_p=x(i)
           y_p=y(i)
           z_p=z(i)
           write(1000+myrank,12) x_p, y_p, z_p, nts
           do j=1, np
              if((abs(y_p-v((j-1)*(nts+1)+1)).lt.epsd).and.(abs(z_p-w((
     & j-1)*(nts+1)+1)).lt.epsd)) then
                do k=1, nts
                write(1000+myrank,11) u((j-1)*(nts+1)+1+k), v((j-1)*(nts
     & +1)+1+k), w((j-1)*(nts+1)+1+k), t((j-1)*(nts+1)+1+k)
                enddo
              else
                continue
              end if
            enddo !end of search identical points
         enddo !end of loop over points on processor
         close(1000+myrank)
         else if(stat>0) then
            write(*,*) 'something is wrong with reading operation'
         else
            exit
         endif
         enddo
         close(20000)
      !   enddo !end of current processor
       !  close(bctnum+1)
        ! close(p) !end of current bct series
         enddo
13       format(1x, 2I10)
12       format(1x, 3E15.7, I10)
11       format(1x, 4E15.7, 10I10)
         end program
