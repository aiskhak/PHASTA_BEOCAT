      subroutine openf(fname,unt)
c---------------------------------------------------------------------
c
c open the Fortran file
c      
c---------------------------------------------------------------------
      integer unt
      character*(*) fname
      
      open(unit=unt, file=fname, status='unknown', form='unformatted')
      
      return
      end
      
      subroutine closef(unt)
c---------------------------------------------------------------------
c
c close the Fortran file
c      
c---------------------------------------------------------------------
      integer unt
      
      close(unt)
      
      return
      end
     
      subroutine readfdbl(unt,ar,n1,n2,nth)
c---------------------------------------------------------------------
c
c Read the Fortran file (real*8)
c      
c---------------------------------------------------------------------
      real*8 ar(n1*n2)
      integer n1,n2,nth,unt
c
c.... skip to the desired record
c      
      do i=1,nth
         read (unt) 
      enddo

      read (unt) ((ar((i-1)*n2+j),i=1,n1),j=1,n2)
      
      return
      end
      
     
      subroutine readfhd(unt,ar)
c---------------------------------------------------------------------
c
c Read the Fortran file (restart header)
c      
c---------------------------------------------------------------------
      integer ar(2),unt
      character*8 machine

      read (unt) machine, ar(1),ar(2)
      
      return
      end
      
      subroutine readfint(unt,ar,n1,n2,nth)
c---------------------------------------------------------------------
c
c Read the Fortran file (integer)
c      
c---------------------------------------------------------------------
      integer ar(n1*n2)
      integer n1,n2,nth,unt
c
c.... skip to the desired record
c      
      do i=1,nth
         read (unt) 
      enddo

      read (unt) ((ar((i-1)*n2+j),i=1,n1),j=1,n2)
      
      return
      end

      subroutine writefdbl(unt,ar,n1,n2)
c---------------------------------------------------------------------
c
c Write the Fortran file (real*8)
c      
c---------------------------------------------------------------------
      real*8 ar(n1*n2)
      integer n1,n2,unt

      write (unt) ((ar((i-1)*n2+j),i=1,n1),j=1,n2)
      
      return
      end
      
     
      subroutine writefhd(unt,ar)
c---------------------------------------------------------------------
c
c Write the Fortran file (restart header)
c      
c---------------------------------------------------------------------
      integer ar(2),unt
      character*8 machine

      machine = 'noname  '
      
      write (unt) machine, ar(1),ar(2)
      
      return
      end
      
      subroutine writefint(unt,ar,n1,n2)
c---------------------------------------------------------------------
c
c Write the Fortran file (integer)
c      
c---------------------------------------------------------------------
      integer ar(n1*n2)
      integer n1,n2,unt

      write (unt) ((ar((i-1)*n2+j),i=1,n1),j=1,n2)
      
      return
      end
      
     
     
