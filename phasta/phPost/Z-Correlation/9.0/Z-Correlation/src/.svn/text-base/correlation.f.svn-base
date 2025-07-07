      subroutine Correlation(gvelavg,ur,vr,wr,Lx,Lz,
     &                        nx,ny,nz,
     &                        utausq,Ruu,ynum,t)

      integer nz,ny,nx,i,j,k,x,y,z,nn,n,e,itmp

      integer kcount, t, ynum

      real*8 ur(nx+1,ny+1,nz+1)   ! x-velocity
      real*8 vr(nx+1,ny+1,nz+1)   ! y-velocity
      real*8 wr(nx+1,ny+1,nz+1)   ! z-velocity
      
      real*8 urE(nx+1,ny+1,2*nz)
      real*8 vrE(nx+1,ny+1,2*nz)
      real*8 wrE(nx+1,ny+1,2*nz)

      real*8  Ruup(3,nx,ny+1,nz)
      real*8 Ruup2(3,nx,ny+1,nz)
      real*8 Ruu( 3,ny+1,nz+1)  ! Average of correlation 
      real*8 Ruu2( 3,ny+1,nz+1)
      real*8 Ruu3( 3,ny+1,nz+1)  
      real*8 avgvel(3)
      real*8 avgvel2(nx+1,ny+1,3), avgvel3(3,ny+1)    
      real*8 gvelavg(3,ny+1)

      real*8 pi, Lx, Lz, Area, fwall, utausq

      character*10 cname
      character*30 fname,fmt

      pi = 4.0*atan(1.0)

      Area = 2.d0*Lx*Lz
c      utausq = 0.004219353d0 !   fwall/Area ! = tau_wall 

c---------------------------Correlations----------------------------

c... First obtain periodic extensions of ur,vr,wr to facilitate correlations

         urE = 0.0d0
         vrE = 0.0d0
         wrE = 0.0d0

      do k = 1, nz
         do j = 1, ny+1
            do i = 1, nx+1
               urE(i,j,k) = ur(i,j,k)
               vrE(i,j,k) = vr(i,j,k)
               wrE(i,j,k) = wr(i,j,k)
            enddo
         enddo
      enddo

      kcount = 0
      do k = nz+1, 2*nz
         kcount = kcount + 1
         do j = 1, ny+1
            do i = 1, nx+1
               urE(i,j,k) = urE(i,j,kcount)
               vrE(i,j,k) = vrE(i,j,kcount)
               wrE(i,j,k) = wrE(i,j,kcount)
               
            enddo
         enddo
      enddo

      avgvel2 = 0.0d0

      do i = 1, nx+1
         do j = 1, ny+1            
            do z = 1, nz

               avgvel2(i,j,1) = avgvel2(i,j,1) + ur(i,j,z)
               avgvel2(i,j,2) = avgvel2(i,j,2) + vr(i,j,z)
               avgvel2(i,j,3) = avgvel2(i,j,3) + wr(i,j,z) 

            enddo            
         enddo
      enddo

      avgvel2 = avgvel2/dfloat(nz)
      avgvel3 = 0.0d0

      do j = 1, ny+1
         do k = 1, nz
            do i = 1, nx
               
               avgvel3(1,j) = avgvel3(1,j) + ur(i,j,k)
               avgvel3(2,j) = avgvel3(2,j) + vr(i,j,k)
               avgvel3(3,j) = avgvel3(3,j) + wr(i,j,k)
               
            enddo
         enddo
      enddo
      
      gvelavg =  gvelavg + avgvel3 ! average of velocity

      do k = 1, 2*nz            ! Substract average velocity to obtain fluctuation
         do j = 1, ny+1
            do i = 1, nx+1
               
               urE(i,j,k) = urE(i,j,k)  
     &              - gvelavg(1,j)/( dfloat( nz*nx ) )
               vrE(i,j,k) = vrE(i,j,k) 
     &              - gvelavg(2,j)/( dfloat( nz*nx ) )
               wrE(i,j,k) = wrE(i,j,k) 
     &              - gvelavg(3,j)/( dfloat( nz*nx ) )
               
            enddo
         enddo
      enddo

c... Get correlations

      Ruup = 0.0d0
      Ruup2 = 0.0d0

      do j = 1, ny+1            ! Fix y-plane
         do i = 1, nx           ! Fix a j
            do n = 1, nz        ! Loop over dx              
               do k = 1, nz     ! Loop over i (loop over products)                 
                  if (n .eq. 1) then

                     Ruup(1,i,j,n) = Ruup(1,i,j,n) + 
     &                    urE(i,j,k)*urE(i,j,k)
                     Ruup(2,i,j,n) = Ruup(2,i,j,n) + 
     &                    vrE(i,j,k)*vrE(i,j,k)
                     Ruup(3,i,j,n) = Ruup(3,i,j,n) + 
     &                    wrE(i,j,k)*wrE(i,j,k)                         
                     
                  else

                     nn = k + (n-1)
                     
                     Ruup(1,i,j,n) = Ruup(1,i,j,n) + 
     &                    urE(i,j,k)*urE(i,j,nn)
                     Ruup(2,i,j,n) = Ruup(2,i,j,n) + 
     &                    vrE(i,j,k)*vrE(i,j,nn)
                     Ruup(3,i,j,n) = Ruup(3,i,j,n) + 
     &                    wrE(i,j,k)*wrE(i,j,nn)             

                  endif                  
               enddo               
            enddo
         enddo
      enddo


      do j = 1, ny+1            ! Fix y-plane
         do i = 1, nx           ! Fix a j
            do n = 1, nz        ! Loop over dx               
               do k = 2*nz, nz+1, -1 ! Loop over i (loop over products)                  
                  if (n .eq. 1) then
                     
                     Ruup2(1,i,j,n) = Ruup2(1,i,j,n) + 
     &                    urE(i,j,k)*urE(i,j,k)
                     Ruup2(2,i,j,n) = Ruup2(2,i,j,n) + 
     &                    vrE(i,j,k)*vrE(i,j,k)
                     Ruup2(3,i,j,n) = Ruup2(3,i,j,n) + 
     &                    wrE(i,j,k)*wrE(i,j,k)                         
                     
                  else
                     
                     nn = k - (n-1)
                     
                     Ruup2(1,i,j,n) = Ruup2(1,i,j,n) + 
     &                    urE(i,j,k)*urE(i,j,nn)
                     Ruup2(2,i,j,n) = Ruup2(2,i,j,n) + 
     &                    vrE(i,j,k)*vrE(i,j,nn)
                     Ruup2(3,i,j,n) = Ruup2(3,i,j,n) + 
     &                    wrE(i,j,k)*wrE(i,j,nn)             
                     
                  endif                  
               enddo               
            enddo
         enddo
      enddo
      
      Ruup =   Ruup/utausq
      Ruup2 = Ruup2/utausq


      do k = 1, 3
         do n = 1, nz
            do j = 1, ny+1               
               do i = 1, nx

                  Ruu(k,j,n) = Ruu(k,j,n) +
     &       ( Ruup(k,i,j,n)+Ruup2(k,i,j,n) )/2.d0

               enddo               
            enddo
         enddo
      enddo


c***************************************************************
         WRITE(*,*)" Writing Z Correlation file "
         WRITE(*,*)" "

      itmp = 1

      if (t.gt. 0) itmp = int(log10(float(t)))+1

      write (fmt,
     &     "('(''Zcorrelation.'',i',i1,',1x)')") itmp
      write (fname,fmt) t
      fname = trim(fname)
      open (unit=800, file=fname, status='unknown' )

      do k = 1, nz+1
            write(800,999)float(k-1)*(Lz/dfloat(nz)),
     &        Ruu(1,ynum,k)/
     &        ( dfloat( nz*nx ) ) ,
     &        Ruu(2,ynum,k)/
     &        ( dfloat( nz*nx ) ) ,
     &        Ruu(3,ynum,k)/
     &        ( dfloat( nz*nx ) )
         enddo

      close(800)

 999     format(1x,4e18.8)

c********************************************************
         return
         end
        

