      subroutine Correlation(gvelavg,ur,vr,wr,Lx,Lz,
     &                        nx,ny,nz,
     &                        utausq,Ruu,ynum,t)

      integer nz,ny,nx,i,j,k,x,y,z,nn,n,e,itmp

      integer icount, t, ynum

      real*8 ur(nx+1,ny+1,nz+1)   ! x-velocity
      real*8 vr(nx+1,ny+1,nz+1)   ! y-velocity
      real*8 wr(nx+1,ny+1,nz+1)   ! z-velocity
      
      real*8 urE(2*nx,ny+1,nz+1)
      real*8 vrE(2*nx,ny+1,nz+1)
      real*8 wrE(2*nx,ny+1,nz+1)

      real*8  Ruup(3,nx,ny+1,nz)
      real*8 Ruup2(3,nx,ny+1,nz)
      real*8 Ruu( 3,nx+1,ny+1)  ! Average of correlation 
      real*8 Ruu2( 3,nx+1,ny+1)
      real*8 Ruu3( 3,nx+1,ny+1)  
      real*8 avgvel(3)
      real*8 avgvel2(3,ny+1,nz+1), avgvel3(3,ny+1)    
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

      do i = 1, nx
         do j = 1, ny+1
            do k = 1, nz+1
               urE(i,j,k) = ur(i,j,k)
               vrE(i,j,k) = vr(i,j,k)
               wrE(i,j,k) = wr(i,j,k)
            enddo
         enddo
      enddo

      icount = 0
      do i = nx+1, 2*nx
         icount = icount + 1
         do j = 1, ny+1
            do k = 1, nz+1
               urE(i,j,k) = urE(icount,j,k)
               vrE(i,j,k) = vrE(icount,j,k)
               wrE(i,j,k) = wrE(icount,j,k)
               
            enddo
         enddo
      enddo

      avgvel2 = 0.0d0

      do k = 1, nz+1
         do j = 1, ny+1            
            do x = 1, nx

               avgvel2(1,j,k) = avgvel2(1,j,k) + ur(x,j,k)
               avgvel2(2,j,k) = avgvel2(2,j,k) + vr(x,j,k)
               avgvel2(3,j,k) = avgvel2(3,j,k) + wr(x,j,k) 

            enddo            
         enddo
      enddo

      avgvel2 = avgvel2/dfloat(nx)
      avgvel3 = 0.0d0

      do j = 1, ny+1
         do i = 1, nx
            do k = 1, nz
               
               avgvel3(1,j) = avgvel3(1,j) + ur(i,j,k)
               avgvel3(2,j) = avgvel3(2,j) + vr(i,j,k)
               avgvel3(3,j) = avgvel3(3,j) + wr(i,j,k)
               
            enddo
         enddo
      enddo
      
      gvelavg =  gvelavg + avgvel3 ! average of velocity

      do i = 1, 2*nx            ! Substract average velocity to obtain fluctuation
         do j = 1, ny+1
            do k = 1, nz+1
               
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
c      Ruu = 0.0d0

      do j = 1, ny+1            ! Fix z-plane
         do k = 1, nz           ! Fix a j
            do n = 1, nx        ! Loop over dx              
               do i = 1, nx     ! Loop over i (loop over products)                 
                  if (n .eq. 1) then

                     Ruup(1,n,j,k) = Ruup(1,n,j,k) + 
     &                    urE(i,j,k)*urE(i,j,k)
                     Ruup(2,n,j,k) = Ruup(2,n,j,k) + 
     &                    vrE(i,j,k)*vrE(i,j,k)
                     Ruup(3,n,j,k) = Ruup(3,n,j,k) + 
     &                    wrE(i,j,k)*wrE(i,j,k)                         
                     
                  else

                     nn = i + (n-1)
                     
                     Ruup(1,n,j,k) = Ruup(1,n,j,k) + 
     &                    urE(i,j,k)*urE(nn,j,k)
                     Ruup(2,n,j,k) = Ruup(2,n,j,k) + 
     &                    vrE(i,j,k)*vrE(nn,j,k)
                     Ruup(3,n,j,k) = Ruup(3,n,j,k) + 
     &                    wrE(i,j,k)*wrE(nn,j,k)             

                  endif                  
               enddo               
            enddo
         enddo
      enddo


      do j = 1, ny+1            ! Fix z-plane
         do k = 1, nz           ! Fix a j
            do n = 1, nx        ! Loop over dx               
               do i = 2*nx, nx+1, -1 ! Loop over i (loop over products)                  
                  if (n .eq. 1) then
                     
                     Ruup2(1,n,j,k) = Ruup2(1,n,j,k) + 
     &                    urE(i,j,k)*urE(i,j,k)
                     Ruup2(2,n,j,k) = Ruup2(2,n,j,k) + 
     &                    vrE(i,j,k)*vrE(i,j,k)
                     Ruup2(3,n,j,k) = Ruup2(3,n,j,k) + 
     &                    wrE(i,j,k)*wrE(i,j,k)                         
                     
                  else
                     
                     nn = i - (n-1)
                     
                     Ruup2(1,n,j,k) = Ruup2(1,n,j,k) + 
     &                    urE(i,j,k)*urE(nn,j,k)
                     Ruup2(2,n,j,k) = Ruup2(2,n,j,k) + 
     &                    vrE(i,j,k)*vrE(nn,j,k)
                     Ruup2(3,n,j,k) = Ruup2(3,n,j,k) + 
     &                    wrE(i,j,k)*wrE(nn,j,k)             
                     
                  endif                  
               enddo               
            enddo
         enddo
      enddo
      
      Ruup =   Ruup/utausq
      Ruup2 = Ruup2/utausq


      do i = 1, 3
         do n = 1, nx
            do j = 1, ny+1               
               do k = 1, nz

                  Ruu(i,n,j) = Ruu(i,n,j) +
     &       ( Ruup(i,n,j,k)+Ruup2(i,n,j,k) )/2.d0

               enddo               
            enddo
         enddo
      enddo

c***************************************************************
         WRITE(*,*)" Writing X Correlation file "
         WRITE(*,*)" "

      itmp = 1

      if (t.gt. 0) itmp = int(log10(float(t)))+1

      write (fmt,
     &     "('(''Xcorrelation.'',i',i1,',1x)')") itmp
      write (fname,fmt) t
      fname = trim(fname)
      open (unit=800, file=fname, status='unknown' )

      do i = 1, nx+1
            write(800,999)float(i-1)*(Lx/dfloat(nx)),
     &        Ruu(1,i,ynum)/
     &        ( dfloat( nz*nx ) ) ,
     &        Ruu(2,i,ynum)/
     &        ( dfloat( nz*nx ) ) ,
     &        Ruu(3,i,ynum)/
     &        ( dfloat( nz*nx ) )
         enddo

      close(800)

 999     format(1x,4e18.8)

c********************************************************
         return
         end
        

