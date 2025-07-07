!-------------------------------------------------------------------------------
!
!	This file contains several subrountines which are related to bubble boi-
!	ling phenomenon including bubble evaporation capability and condensation
!	capability.
!	                                                                        
!	Mengnan Li,                                                Spring, 2016
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!	Name                    Description
!       dVolume   :             Volume change due to boiling or condensation
!       numshell  :             numbers of element in bubble info collection   
!                               shell
!                               1. shell outside the bubble
!                               2. shell inside  the bubble
!       bubble_vol:             volume of each bubble
!       B_factor  :             Factor controls bubble growth rate from
!                               analytical solution
!       Tempb     :             temperature in each element
!       gytemp    :             temperature gradient in each element
!
!
!
!-------------------------------------------------------------------------------

        subroutine Bubvolgenera(dVolume, shell_num,bml,sclr_ls,elemvol_local) 
!-------------------------------------------------------------------------------
!
!	This subroutine is calculated the volume generation term during boiling
!	and condensation process.(caleld by e3res.f)
!
!-------------------------------------------------------------------------------
        use spat_var_eps ! for spatially varying epsilon_ls	
!        use bub_track
!        use bubboil_info       
        include "common.h"

        dimension dVolume(npro), elem_shell_num(100,2), bubvol(100), R(100)
        dimension R_1(npro), B_factor(100), shell_num(npro),sclr_ls(npro)
        dimension bml(npro,nshl,1), bubdVolume(npro,100),elemvol_local(ibksiz)
        dimension R_0(npro), shell_num_old(npro), shell_num_new(npro)

        integer i, j
	real*8 Rho_l,     Rho_v,     cp_l,       k_l,            elem_shell_num
	real*8 R,         a_l,       B_factor,   dVolume,        epsilon_ls_tmp
	real*8 R_1,       bubvol,    sclr_ls,    elemvol_local,  R_0

        Rho_l = datmat(1,1,1)  !958
        Rho_v = datmat(1,1,2) !0.579
        cp_l = datmat(1,3,1)!1.22
        k_l = datmat(1,4,1)!0.679
        dVolume(:) = 1.0e-10
        shell_num(:)=1.0e-10

!        if(myrank.eq.master)write(*,*)bubboil
        do i=1, i_num_bubbles
!        if (lstep .eq. int(1.0))then
         if (lstep.le.irstart)then
          elem_shell_num(i,1) = 2000
          elem_shell_num(i,2) = 2000
          bubvol(i)=(4.0/3.0)*pi*((2E-4)**3.0)
          bubble_tempG(i)=1.0E-12
!          avg_info(1,19)=1.0E-15
        else
          elem_shell_num(i,1) = numshell_in(i)
          elem_shell_num(i,2) = numshell_out2(i)
          bubvol(i)=bubble_vol(i)
        endif 

 
        R(i) =(((3.0E0/4.0E0)*bubvol(i))/pi)**(1.0E0/3.0E0)
     
        !if(elem_shell_num(i,1).le.2.0d0)then
         ! bubboil = 0.0
        !else
         ! bubboil = 1.0
        !endif
              
!         write(*,*)'R', R
!        endif
        if(bubboil.eq.0.and.bubgrow.eq.1.0d0)then
!         if(bubgrow.eq.1)then  ! for test only
          a_l = k_l/(Rho_l*cp_l) !1.0e3
          B_factor(i) = ((12.0E0*a_l/pi)**(0.5E0))* 
     &         ((delt_T(i)*cp_l*Rho_l)/(h_fg*Rho_v))
!         if(myrank.eq.master)write(*,*)h_fg,T_sat
!        if (myrank.eq.master)write(*,*)'delt_T', delt_T(1,i)
        do j = 1, npro
          bubdVolume(j,i) = 2.0E0*pi*R(i)*(B_factor(i)**(2.0E0))
        enddo
        endif

        if(bubboil.eq.1.0d0)then  ! loop over elements (at the local level)
        do j = 1, npro 
!          epsilon_lsd_tmp = epsilon_ls*
!     &      elem_local_size(lcblk(1,iblk)+j-1)
!          epsilon_lsd_tmp = epsilon_lsd*
!     &         (elemvol_local(j)**(1.0/3.0))
!          epsilon_lsd_tmp = 
!     &         (elemvol_local(j)**(1.0/3.0))
!          R_1(i)=R(i)+epsilonBT*1.0   ! sqrt(3.0) for structure mesh
                                                    ! sqrt(2.0) for parasolid mesh
!        write(*,*)'elemvol_local',elemvol_local(j)   
        bubdVolume(j,i) = bubble_tempG(i)*k_l*(1.0/Rho_v)
     &                    /h_fg
        bubdVolume(j,i) = bubdVolume(j,i)*2.0d0*pi*(R(i)**(2.0))
        enddo
                
        endif

        enddo ! i_num_bubbles

        bubboil = 1.0   ! reset boiling flag
!       Assembly the source matrix
        do i = 1, i_num_bubbles
           do j = 1, npro

          epsilon_ls_tmp = epsilon_lst*
     &      elem_local_size(lcblk(1,iblk)+j-1)

      
!          if ((sclr_ls(j).GT.-2.0E0*epsilon_ls_tmp).and.
!     &    (sclr_ls(j).LT.-1.0E0*epsilon_ls_tmp)) then
!              do n = 1, nshl
!                if(INT(bml(j,n,1)).eq.i) then
!                if((INT(bml(j,n,1)).eq.i).or.(INT(bml(j,n,1)).eq.(3*i_num_bubbles+i)))then
!               specially for single nucleation site
!                 dVolume(j) = bubdVolume(j,i)
!                 shell_num(j) = elem_shell_num(i,1)
!                endif
!              enddo          


          if ((sclr_ls(j).GT.1.0E0*epsilonBT).and.
     &    (sclr_ls(j).LT.2.0E0*epsilonBT)) then
              do n = 1, nshl
!                if(INT(bml(j,n,1)).eq.i) then
                if((INT(bml(j,n,1)).eq.i).or.(INT(bml(j,n,1)).eq.(3*i_num_bubbles+i)))then
!               specially for single nucleation site
                 dVolume(j) = bubdVolume(j,i)*Rho_v/Rho_l
                 shell_num(j) = elem_shell_num(i,2)
                endif
              enddo
          elseif ((sclr_ls(j).GT.-2.0E0*epsilon_ls_tmp).and.
     &    (sclr_ls(j).LT.-1.0E0*epsilon_ls_tmp)) then
              do n = 1, nshl
!                if(INT(bml(j,n,1)).eq.i) then
                if((INT(bml(j,n,1)).eq.i).or.(INT(bml(j,n,1)).eq.(3*i_num_bubbles+i)))then
!               specially for single nucleation site
                 dVolume(j) = bubdVolume(j,i)
                 shell_num(j) = elem_shell_num(i,1)
                endif
              enddo
          else

                 dVolume(j) = 0.0d0
                 shell_num(j) = 0.0d0
           endif

!        if(shell_num(j).ne.0.0)write(*,*)shell_num(j)

           if(isnan(dVolume(j))) dVolume(j) = 0.0d0
         
!        if(lstep.eq.2)then
!        if(shell_num(j).ne.0.0d0.and.shell_num(j).lt.1.0)then
!         write(*,*)shell_num(j)
!        write(*,*) bubble_tempG(i)
!        endif
!        endif

          enddo
        
        enddo
        return
        end
!-------------------------------------------------------------------------------
c===============================================================================
c===============================================================================

        subroutine Bubheatflux(yl, shpfun, shg, elemvol_local,Tempb, gytemp,bml)
!-------------------------------------------------------------------------------
!
!	This subroutine is used to calculated total heat flux flowing into bubb-
!	-le due to temperature difference(called in e3ivar.f)
!
!-------------------------------------------------------------------------------
        use  spat_var_eps ! for spatially varying epsilon_ls
!        use  bubboil_info!        
        use  bub_track
        include "common.h"
        integer    i,  n
        dimension  yl(npro,nshl,ndof),       shpfun(npro,nshl),
     &             shg(npro,nshl,nsd),       Temp(npro),
     &             elemvol_local(ibksiz),    
     &             gyti(npro,nsd),
     &             Sclr(npro),               Tempb(npro),
     &             gytemp(npro,5),           R_0(npro),
     &             R_1(npro),                bubvol(100)
        dimension  bml(npro,nshl,1)
	real*8     elemvol_local,            elemvol            
        real*8     epsilon_lsd_tmp
	real*8     epsilon_ls_tmp,           Tempb, 
     &             gytemp,                     
     &             R_0,                      R_1,
     &             num,
     &             local_volume,             bml,
     &             bubvol
        Sclr = zero
        Tempb = zero
        gyti = zero
        gytemp(:,:) = zero
        R_0 = zero
        
        do i = 1, npro
          do n = 1, nshl
            Sclr(i) = Sclr(i) + shpfun(i,n) * yl(i,n,6) !scalar
c     
c!     .... compute the global gradient of Scalar variable
c     
            gyti(i,1) = gyti(i,1) + shg(i,n,1) * yl(i,n,6) 
            gyti(i,2) = gyti(i,2) + shg(i,n,2) * yl(i,n,6)
            gyti(i,3) = gyti(i,3) + shg(i,n,3) * yl(i,n,6)
c     
           enddo
         enddo
c
c!  .... compute the global gradient of Temperature outside bubble 2 epsilon
c

!        if(myrank.eq.master)write(*,*)'lstep', lstep

       
        do i=1, npro
         do k=1, i_num_bubbles

             if (lstep.le.irstart) then
                 bubvol(k)=(4.0/3.0)*pi*((2E-4)**3.0)
             else
                 bubvol(k)=bubble_vol(k)
        !           write(*,*) "R_0", R_0(i)
             endif
             


             do n = 1, nshl
!                if(INT(bml(i,n,1)).eq.k) then
                if((INT(bml(i,n,1)).eq.k).or.(INT(bml(i,n,1)).eq.(3*i_num_bubbles+k)))then
!               specially for single nucleation site
                R_0(i)=(((3.0E0/4.0E0)*bubvol(k))/pi)
                R_0(i)=(R_0(i))**(1.0d0/3.0d0)
!                R_1(i)=R_0(i)+epsilon_ls_tmp*2.5/sqrt(3.0)  
!                R_1(i)=R_0(i)+epsilonBT*2.5

                endif
              enddo

            enddo
          
         enddo 
!         elemvol = 0.0d0
!         do i=1, npro
!           if (Sclr(i).le. 1.0E0*epsilonBT
!     &       .and. Sclr(i).ge. 0.0E0*epsilonBT) then
!           elemvol=elemvol+elemvol_local(i)
!           endif
!         enddo
!           R_1=0.0
!           R_0=(((3.0E0/4.0E0)*bubble_vol(k))/pi)**(1.0E0/3.0E0)
!           R_1=R_0+epsilon_lsd_tmp*2.5/sqrt(2.0)

!          write(*,*) 'bubble_vol,R_0',bubble_vol(1),R_0
!      .... R_0 is zero in the begining of the run due to the sequence issue
!        if (lstep+1 .gt. int(10+inistep))then
        do i = 1, npro

!          epsilon_ls_tmp = epsilon_lst*
!     &      elem_local_size(lcblk(1,iblk)+i-1)

              do n = 1, nshl
!                   Tempb(i) = Tempb(i) + shpfun(i,n) * yl(i,n,5)
!                   R_1(i)=R_0(i)+ shpfun(i,n)*yl(i,n,6)
!                    R_1(i)=R_0(i)+ 2.0d0*epsilon_ls_tmp
                    R_1(i) = R_0(i) + 2.0d0*epsilonBT    
              enddo
!           if (Sclr(i).le. 2.0E0*epsilon_ls_tmp
!     &       .and. Sclr(i).ge. 0.0E0*epsilon_ls_tmp) then
           if (Sclr(i).le. 2.0E0*epsilonBT
     &       .and. Sclr(i).ge. 0.0E0*epsilonBT.and.R_0(i).gt.0.0d0) then

              do n = 1, nshl
                   Tempb(i) = Tempb(i) + shpfun(i,n) * yl(i,n,5) !Temperature
                   gytemp(i,1) = gytemp(i,1) + shg(i,n,1) * yl(i,n,5)
                   gytemp(i,2) = gytemp(i,2) + shg(i,n,2) * yl(i,n,5)
                   gytemp(i,3) = gytemp(i,3) + shg(i,n,3) * yl(i,n,5)
             enddo ! for nshl
!           gytemp(i,4) = (gytemp(i,1)*gyti(i,1)+gytemp(i,2)*gyti(i,2)
!     &   +gytemp(i,3)*gyti(i,3))/sqrt(gyti(i,1)**2+gyti(i,2)**2+gyti(i,3)**2)
           gytemp(i,4) = (gytemp(i,1)*gyti(i,1)+gytemp(i,2)*gyti(i,2)
     &   +gytemp(i,3)*gyti(i,3))
!           gytemp(i,5)=gytemp(i,4)*((R_1(i)/R_0(i))**2.0)*elemvol_local(i)/
!     &       (pi*R_0(i)*R_0(i)*epsilonBT)
           gytemp(i,5)=gytemp(i,4)*((R_1(i)/R_0(i))**2.0)*elemvol_local(i)/
     &       (pi*R_0(i)*R_0(i)*epsilonBT)
            

           else
            gytemp(i,5)=0.0d0
            Tempb(i)=0.0d0
           endif
        
!        if(R_0(i).ne.0.0d0)then
!            write(*,*) R_0(i)
!        endif

        if(isnan(gytemp(i,5)))  gytemp(i,5) = 0.0d0
 
        enddo

        return
        end

!-------------------------------------------------------------------------------
c===============================================================================
c===============================================================================
