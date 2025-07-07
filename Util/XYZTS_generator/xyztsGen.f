        Program XYZTS_generator
!
!----------------------------------------------------------------------
!
!       This routine is developed to generate probes which can be used to
!       collect the data during PHASTA simulations for post- statistical
!       analysis. 
!       variables:
!       npro            total number of probes
!       nproCont        total number of probes along a contour
!       nproBL          total number of probes along initial BL contour
!
!       Jun Fang, Igor Bolotnov                               Fall, 2014
!----------------------------------------------------------------------
        use xyztsInp
        implicit none

        real*8  r, dr, dx
        real*8  Lx, xc, yc
        real*8  x, y, z, pi
        real*8  theta, theta1, theta2, theta3, theta4, t
        real*8  dtheta, dtheta1, dtheta2, dtheta3, dtheta4

        integer i, j, k 
        integer npro, nproCont, nproBL, nprobelayer  
        integer number_probes

        integer Nar, Nrad, Nlayer1      ! For pipe domain
        integer Nregion(1:20), i_max(20)
        
        real*8  Angles(1:10, 1:256)     ! Angles for different regions
        real*8  Radii(1:256)
        logical istrue

c!... Initialization...
        pi = 4.0*atan(1.0)

        call xyztsReader() 
        nproBL = 30              
        t = 0.0002                  !increment in theta to fit desired number of points
        dtheta = 0.0

        open(13,file="xyzts.dat")

        if (idomain .eq. 1) then
        nprobelayer    = 0
        nproCont= nproBL
        r       = rpin
        dr      = yplus
        npro    = 0
        do while (r.lt.rcentre)         !rpin+bl_thick)!rcentre)!start radial distance loop

           do i =1,S                    !start loop over subchannels, find centre of the subchannel
              xc=(Cx(i,1)+Cx(i,2)+Cx(i,3)+Cx(i,4))/4.0
              yc=(Cy(i,1)+Cy(i,2)+Cy(i,3)+Cy(i,4))/4.0

              do j = 1,4                !loop over centres of a subchannel
                 theta1 = 0.0
                 theta2 = pi/2.0 
                 theta3 = pi
                 theta4 = 3.0*pi/2.0
                 theta  = pi/2.0        ! angle total available
                 dtheta = pi/2.0/real(nproBL - 1) ! initial guess for dtheta
                 dtheta1= dtheta
                 dtheta2= dtheta
                 dtheta3= dtheta
                 dtheta4= dtheta
                             
!       first find out dtheta for this sector
                 istrue=.false.                        
                    
                 do while(.not.istrue)
                            
                    if(j.eq.1) x=Cx(i,j)+r*cos(theta1) !test parameter
                    if(j.eq.1) y=Cy(i,j)+r*sin(theta1)

                    if(j.eq.2) x=Cx(i,j)+r*cos(theta2) !testparameter
                    if(j.eq.2) y=Cy(i,j)+r*sin(theta2) !testparameter

                    if(j.eq.3) x=Cx(i,j)+r*cos(theta3) !test parameter
                    if(j.eq.3) y=Cy(i,j)+r*sin(theta3) !testparameter

                    if(j.eq.4) x=Cx(i,j)+r*cos(theta4) !test parameter
                    if(j.eq.4) y=Cy(i,j)+r*sin(theta4) !test arameter
               

                    if(j.eq.1.and.((x.gt.xc).or.(y.gt.yc))) then
                       theta1 = theta1+t
                       theta  = theta-2.0*t ! decreasing available angle              
                       nproCont = int(theta/(pi/2)*real(nproBL))
                       if(nproCont.lt.3) nproCont = 3
                       dtheta1= theta/real(nproCont - 1)
                    endif

                    if(j.eq.2.and.((x.lt.xc).or.(y.gt.yc))) then
                       theta2=theta2+t
                       theta=theta-2.0*t
                       nproCont = int(theta/(pi/2)*real(nproBL))
                       if(nproCont.lt.3) nproCont = 3
                       dtheta2= theta/real(nproCont - 1)
                    endif   

                    if(j.eq.3.and.((x.lt.xc).or.(y.lt.yc))) then
                       theta3=theta3+t
                       theta=theta-2.0*t
                       nproCont = int(theta/(pi/2)*real(nproBL))
                       if(nproCont.lt.3) nproCont = 3
                       dtheta3= theta/real(nproCont - 1)
                    endif

                    if(j.eq.4.and.((x.gt.xc).or.(y.lt.yc))) then
                       theta4=theta4+t
                       theta=theta-2.0*t
                       nproCont = int(theta/(pi/2)*real(nproBL))
                       if(nproCont.lt.3) nproCont = 3
                       dtheta4= theta/real(nproCont - 1)
                    endif

                    if((j.eq.1.and.(x.lt.xc.and.y.lt.yc)).or.
     &                 (j.eq.2.and.(x.gt.xc.and.y.lt.yc)).or.
     &                 (j.eq.3.and.(x.gt.xc.and.y.gt.yc)).or.
     &                 (j.eq.4.and.(x.lt.xc.and.y.gt.yc))) then
                      
                        istrue=.true. ! point within limits
                             
                    endif

                       
                    enddo ! dtheta has been found out for this sector

                    do k = 1, nproCont  ! start along circle seont =

                       if(j.eq.1) then  ! centre type 1
                               
                          x=Cx(i,j)+r*Cos(theta1)
                          y=Cy(i,j)+r*Sin(theta1)         
                          theta1=theta1+dtheta1
                                 
                       elseif(j.eq.2) then
                                        
                          x=Cx(i,j)+r*Cos(theta2)
                          y=Cy(i,j)+r*Sin(theta2)   
                          theta2=theta2+dtheta2

                       elseif(j.eq.3) then     
                                 
                          x=Cx(i,j)+r*Cos(theta3)
                          y=Cy(i,j)+r*Sin(theta3)          
                          theta3=theta3+dtheta3
                                        
                       elseif(j.eq.4) then

                          x=Cx(i,j)+r*Cos(theta4)
                          y=Cy(i,j)+r*Sin(theta4)            
                          theta4=theta4+dtheta4
                       endif                                   
                    
                       write(13,*) axial, x,y
                       npro = npro + 1
                    enddo ! end loop along sectors

                 enddo ! end loop over centres of subchannel

              enddo ! end loop over subchannels

              if (r .le. (rpin + bl_thick)) then
                 dr = dr * ratio
              else
                 dr = 10.0 * yplus
              endif
              r = r + dr
              nprobelayer = nprobelayer + 1
!                Print *, r,rmax
        enddo ! end radial distance loop

        Print *, "radial layers", nprobelayer

        write(*,*) 'Total number of probes: ', npro
        
        Print * , "N,s",N,s, rcentre-rpin

        end if

        if (idomain.eq.2) then  ! Pipe option here

! Create an array for radial locations
        
! Number of regions:
         Nar = 3
           
! Create an array for theta angles for all regions:
         do j = 1, Nar
          i_max(j) = int(real(nhom)/2.0**real(j-1))
          do i = 1, i_max(j)
            Angles(j,i) = real(i-1)/real(i_max(j))*2.0*pi           
!           write(*,*) j, i, Angles(j,i)
          end do
         end do

! Radial locations: (ignore BL for now!)
! For uniform in R spacing and preserving total homogeneous directions
! we will need to have a geometric progression in terms of number of
! layers in each "angle" region
         
! First let's compute number of the radius coordinates
        Nrad = 0
! Number of layers for the dense near wall region (will be spaced as BL
! in the future):
        Nlayer1 = 4
        do j = 1, Nar
         Nregion(j) = Nlayer1*int(2.0**(j-1))
         Nrad = Nrad + Nregion(j)
        end do
        write(*,*) 'Nrad = ', Nrad

! Split and fill with coordinates:
        do j = 1, Nrad
!         Radii(j) = rmax*(1.0 - (real(j)-0.5)/real(Nrad))
         Radii(j) = rmax*cos((real(j)-0.5)/real(Nrad)*0.5*pi)
        end do

! See if we can print out the final result:        
! (check the sequence later)
        Nrad = 0
        do j = 1, Nar   ! Loop over the regions
         do k = 1, Nregion(j)  ! Loop over the radial positions within the region
          do i = 1, i_max(j)   ! Loop over the tangential angles
           write(13,*) axial, Radii(Nrad+k)*cos(Angles(j,i))
     &                      , Radii(Nrad+k)*sin(Angles(j,i))
          end do
         end do
         Nrad = Nrad + Nregion(j)
        end do

        end if        
     
        close(13)

        end program XYZTS_generator   
