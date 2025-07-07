      program procstats

      integer, allocatable :: imap(:), invmap(:), mapb(:)
      real*8,  allocatable :: qstat(:,:), qstr(:,:), qsta(:,:),
     &                        xyz(:,:), qread(:,:), qstat_new(:,:),
     &			      qsta1(:,:)  
c      real*8 wght,wghta
      integer      numnp, nshg, i, j, k, itmp, itmpo, nsteps,
     &             fkmap, iorig, ifirstf, ilastf, nstepst,
     &             Nx, Ny, Nz, iout, ifirst, ilast,istatq,ind,
     &             isize, nitems, irstin, ind2

      character*30 fname1,fmt1, fnamer
      character*1  dash
      integer      nnp, nsd,ii, iout2, iout3, iout4, iout5
      integer iarray(50)
c
c.... read coordinates from geom.old
c
    
      open (unit=33,file='geom.old',status='unknown')
      read (33,*) nnp, nsd
      allocate (xyz(nnp,3))
      do i = 1, nnp
         read (33,*) j, j, xyz(i,1), xyz(i,2), xyz(i,3)
      enddo
      close (33)
      
c
c.... read the time averaged data
c
      write(*,*) 'enter nx, ny, and nz '
      read(*,*) Nx, Ny, Nz
      numnp = Nx * Ny * Nz
      allocate (imap(numnp), invmap(numnp), mapb(numnp) )
      allocate (qstat(numnp,19))
      allocate (qstr(numnp,19))
      allocate (qsta(Nx*Ny,19))
      allocate (qsta1(Nx*Ny,19))
      allocate (qread(numnp,19))
      allocate (qstat_new(numnp,19))
  
      qstr=0.0
c      wghta=0.0
      ifirstf=1000000000
      ilastf=0

 169  write(*,*) 'enter first & last time steps '
      read(*,*) ifirst,ilast
      
      dash='-'
      if (ifirst .gt. 0) itmpo = int(log10(float(ifirst)))+1
      if (ilast .gt. 0)  itmp  = int(log10(float(ilast)))+1
c     
      write (fmt1,
     &     "('(''stats.'',i',i1,',a1,i',i1,',1x)')")
     &     itmpo,itmp
      write (fname1,fmt1) ifirst,dash,ilast
      fname1 = trim(fname1) // '.out'
      call openfile(  fname1,  'read?', irstin )

      fnamer = 'statistics'
      nitems = 5
      
      call readheader (irstin, fnamer, iarray, nitems,
     &     'integer', iotype)

      nshg = iarray(1)     
      isize = nshg*19


      call readdatablock (irstin, fnamer, qread, isize, 'double',iotype)

      qstat(1:nshg,:)=qread(1:nshg,:)

      qstr=qstr+qstat
      ifirstf=min(ifirstf,ifirst)
      ilastf=max(ilastf,ilast)

       ifirst=ifirstf
      ilast=ilastf
c      wghta=1.0/wghta
      qstat=qstr
c      write(*,*) "wghta=",wghta
      
c
c.... read the mapping file to get from mesh-database to (i,j,k)
c     numbering
c
      fkmap = 18
      open(fkmap,file='geom.kmap',status='unknown')
      do i=1,numnp
         read(fkmap,*) iorig
         imap(i)=iorig+1  ! imap takes arguement of original and returns new
         invmap(iorig+1)=i
      enddo 
      close(fkmap)
c
c.... convert to original structured numbering (i,j,k) (and time average)
c
      do i = 1, numnp
         qstr(i,:) = qstat(invmap(i),:)
      enddo
c
c.... perform the homogeneous averaging       
c
      qsta = 0.0
      do j = 1, Ny
         do i = 1, Nx
            do k = 1, Nz
               ind       = (i-1)*Ny*Nz + (j-1)*Nz + k
	       ind2      = (i-1)*Ny + j
               qsta(ind2,:) = qsta(ind2,:) + qstr(ind,:)
            enddo
         enddo
      enddo
      do ind2 = 1, Nx*Ny
         qsta(ind2,:) = qsta(ind2,:)/real(Nz,8)
      enddo
      
      qsta1 = qsta
      do ind2 = 1, Nx*Ny
         qsta1(ind2,10) = qsta(ind2,11)
         qsta1(ind2,11) = qsta(ind2,12)
         qsta1(ind2,12) = qsta(ind2,13)
      enddo

c
c.... output the results
c
      iout = 39
      iout2= 40
      iout3= 41
      iout4= 42
      iout5= 43

      if (ifirst .gt. 0) itmpo = int(log10(float(ifirst)))+1
      if (ilast .gt. 0)  itmp  = int(log10(float(ilast)))+1
c     

      write (fmt1,
     &     "('(''sts_ui.'',i',i1,',a1,i',i1,',1x)')")
     &     itmpo,itmp
      write (fname1,fmt1) ifirst,dash,ilast
      fname1 = trim(fname1)
      open (unit=iout, file=fname1, status='unknown' ) 
c
      write (fmt1,
     &     "('(''sts_uiuj.'',i',i1,',a1,i',i1,',1x)')")
     &     itmpo,itmp
      write (fname1,fmt1) ifirst,dash,ilast
      fname1 = trim(fname1)
      open (unit=iout2, file=fname1, status='unknown') 
c
      write (fmt1,
     &     "('(''sts_stress.'',i',i1,',a1,i',i1,',1x)')")
     &     itmpo,itmp
      write (fname1,fmt1) ifirst,dash,ilast
      fname1 = trim(fname1)
      open (unit=iout3, file=fname1, status='unknown')
c
      write (fmt1,
     &     "('(''sts_p.'',i',i1,',a1,i',i1,',1x)')")
     &     itmpo,itmp
      write (fname1,fmt1) ifirst,dash,ilast
      fname1 = trim(fname1)
      open (unit=iout4, file=fname1, status='unknown') 
c
      write (fmt1,
     &     "('(''sts_psq.'',i',i1,',a1,i',i1,',1x)')")
     &     itmpo,itmp
      write (fname1,fmt1) ifirst,dash,ilast
      fname1 = trim(fname1)
      open (unit=iout5, file=fname1, status='unknown') 

      do i = 1, Nx
         do j = 1, Ny        
            ind  = (i-1)*Ny*Nz + (j-1)*Nz + 1  
	    ind2 =  (i-1)*Ny + j
            write (iout, 111) 
     &		(xyz(ind,1),xyz(ind,2),(qsta(ind2,k), k=2,4))
            write (iout2,112) 
     &		(xyz(ind,1),xyz(ind,2),(qsta1(ind2,k), k=7,12))
c            write (iout3,112) 
c     &		(xyz(ind,1),xyz(ind,2),(qsta(ind2,k), k=14,19))
            write (iout4,110) 
     &		(xyz(ind,1),xyz(ind,2),qsta(ind2,1))
            write (iout5,110) 
     &		(xyz(ind,1),xyz(ind,2),qsta(ind2,5))
         enddo
      enddo

   
c      do j = 1, Ny 
c         if (j.eq.1) then
c           ind = 1
c         else if (j.eq.2) then
c           ind = 10 
c         else 
c           ind = 3
c         endif
c         print *, 'ind', ind
c         write (iout,111) (xyz(ind,2),(qsta(j,i), i=2,4))
c         write (iout2,112) (xyz(ind,2),(qsta(j,i), i=7,12))
c         write (iout3,112) (xyz(ind,2),(qsta(j,i), i=14,19))
c         write (iout4,110) (xyz(ind,2),qsta(j,1))
c         write (iout5,110) (xyz(ind,2),qsta(j,6))
c      enddo


      close (iout)
      close (iout2)
      close (iout3)
      close (iout4)
      close (iout5)

 110  format(1p,3e12.4)
 111  format(1p,5e12.4)
 112  format(1p,8e12.4)
 113  format(1p,7e24.16)

      close(iout)
      
      stop
c

c.... ---------------------->  Error Handling  <-----------------------
c
c.... Error handling
c
ci
c.... file error handling
c
995	write(*,*) 'input file error'
996	write(*,*) 'output file error'
c
c.... end
c
	end
