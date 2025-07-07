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
!
!
!
!
!
!
!
!
!
!
!
!
!-------------------------------------------------------------------------------

        subroutine Bubvolgenera(dVolume, elem_shell_num) 
!-------------------------------------------------------------------------------
!
!	This subroutine is calculated the volume generation term during boiling
!	and condensation process.(caleld by e3res.f)
!
!-------------------------------------------------------------------------------
        use spat_var_eps ! for spatially varying epsilon_ls	
!        use bub_track
        use bubboil_info       
        include "common.h"

        dimension dVolume(npro)
        integer i
	real*8 Rho_l,     Rho_v,     cp_l,       k_l,            elem_shell_num
	real*8 R,         a_l,       B_factor,   dVolume,        epsilon_lsd_tmp
	real*8 R_1

        Rho_l = datmat(1,1,1)  !958
        Rho_v = datmat(1,1,2) !0.579
        cp_l = datmat(1,3,1)!4.22
        k_l = datmat(1,4,1)!0.679
        if (lstep+1 .eq. int(1+inistep))then
          elem_shell_num = 2000
          bubble_vol(1)=(4/3)*pi*((2E-4)**3)
!          avg_info(1,19)=1.0E-15
        else
          elem_shell_num = numshell(1,2)
        endif
        R =(((3.0E0/4.0E0)*bubble_vol(1))/pi)**(1.0E0/3.0E0)
!        write(*,*)'R,numshell',R, elem_shell_num
        if(bubboil.eq.0.and.bubgrow.eq.1)then
          a_l = k_l/(Rho_l*cp_l*1.0E3)
          B_factor = ((12.0E0*a_l/pi)**(0.5E0))* 
     &               ((delt_T*cp_l*Rho_l)/(h_fg*Rho_v))
        do i = 1, npro
          dVolume(i) = 2.0E0*pi*R*(B_factor**(2.0E0))
        enddo
!        write(*,*)'dVolume',dVolume
        endif
        if(bubboil.eq.1)then  ! loop over elements (at the local level)
        do i = 1, npro 
          epsilon_lsd_tmp = epsilon_lsd*elem_local_size(lcblk(1,iblk)+i-1) 
          R_1=R+epsilon_lsd_tmp*2.5/sqrt(2.0)
!          dVolume(i) = k_l/(h_fg*Rho_v)/(epsilon_lsd_tmp)
          dVolume(i) = bubble_tempG(1)*k_l/(h_fg*Rho_v)/(epsilon_lsd_tmp)
        enddo
        endif
        return
        end
!-------------------------------------------------------------------------------
c===============================================================================
c===============================================================================

        subroutine Bubheatflux(yl, shpfun, shg, Tempb,  gytemp,  shell_elem)
!-------------------------------------------------------------------------------
!
!	This subroutine is used to calculated total heat flux flowing into bubb-
!	-le due to temperature difference(called in e3ivar.f)
!
!-------------------------------------------------------------------------------
        use  spat_var_eps ! for spatially varying epsilon_ls
        use  bubboil_info
!        use  bub_track
        include "common.h"
        integer    i,  n
        dimension  yl(npro,nshl,ndof),       shpfun(npro,nshl),
     &             shg(npro,nshl,nsd),       Temp(npro),
     &             elemvol_local(ibksiz),    
     &             gyti(npro,nsd),
     &             Sclr(npro),               Tempb(npro),
     &             gytemp(npro,5),           shell_elem(npro,2)
	real*8     elemvol_local            
        real*8     epsilon_lsd_tmp
	real*8     epsilon_ls_tmp,           Tempb, 
     &             gytemp,                     
     &             R_0,                      R_1,
     &             shell_elem,               num,
     &             local_volume
        Sclr = zero
        Tempb = zero
        gyti = zero
        shell_elem= zero
!         write(*,*) 'bubble_vol',bubble_vol(1)
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
        do i = 1, npro
            epsilon_lsd_tmp = epsilon_ls *
     &                       elem_local_size(lcblk(1,iblk)+i-1)
!        count the number of shell's element outside the bubble
           if (Sclr(i).le. 3.0E0*epsilon_lsd_tmp
     &       .and. Sclr(i).ge. 2.0E0*epsilon_lsd_tmp) then
           shell_elem(i,1)= epsilon_lsd_tmp
           endif
!        count the number of shell's element inside the bubble
           if (Sclr(i).le. -1.0E0*epsilon_lsd_tmp
     &       .and. Sclr(i).ge. -2.0E0*epsilon_lsd_tmp) then
           shell_elem(i,2)= epsilon_lsd_tmp
           endif 
!           if (shell_elem(i,1).ne. 0.0) then
!            write(*,*)'i,shell_elem(i,1)',i,shell_elem(i,1)
!           endif 
       enddo
c
c!  .... compute the global gradient of Temperature outside bubble 2 epsilon
c
        if (lstep+1 .eq. int(1+inistep))then
            bubble_vol(1)=(4/3)*pi*((2E-4)**3)
        endif

        do i = 1, npro
           epsilon_lsd_tmp = epsilon_ls *
     &                       elem_local_size(lcblk(1,iblk)+i-1)
           gytemp(i,:) = zero
           R_1=0.0
           R_0=(((3.0E0/4.0E0)*bubble_vol(1))/pi)**(1.0E0/3.0E0)
           R_1=R_0+epsilon_lsd_tmp*2.5/sqrt(2.0)

!          write(*,*) 'bubble_vol,R_0',bubble_vol(1),R_0
c!      .... R_0 is zero in the begining of the run due to the sequence issue

           if (Sclr(i).le. 3.0E0*epsilon_lsd_tmp
     &       .and. Sclr(i).ge. 2.0E0*epsilon_lsd_tmp) then
              do n = 1, nshl
                   Tempb(i) = Tempb(i) + shpfun(i,n) * yl(i,n,5) !Temperature
                   gytemp(i,1) = gytemp(i,1) + shg(i,n,1) * yl(i,n,5)
                   gytemp(i,2) = gytemp(i,2) + shg(i,n,2) * yl(i,n,5)
                   gytemp(i,3) = gytemp(i,3) + shg(i,n,3) * yl(i,n,5)
             enddo ! for nshl
           gytemp(i,4) = (gytemp(i,1)*gyti(i,1)+gytemp(i,2)*gyti(i,2)
     &   +gytemp(i,3)*gyti(i,3))/sqrt(gyti(i,1)**2+gyti(i,2)**2+gyti(i,3)**2)
           gytemp(i,5)=gytemp(i,4)*elemvol_local(i)
           else
            gytemp(i,:)=zero
            Tempb(i)=0.0
           endif
!           if (Tempb(i).ne.0.0d0) write(*,*)'Tempb(i)', Tempb(i)
        enddo
        return
        end

!-------------------------------------------------------------------------------
c===============================================================================
c===============================================================================

	subroutine Bubibc(y)
!-------------------------------------------------------------------------------
!
!	This subroutine is used to set up the boudary condition inside the bubb-
!	le(called in itrbc.f)
!
!-------------------------------------------------------------------------------
        include "common.h"
	integer i
	dimension y(nshg,ndof)
	real*8    y
	if(isclr.eq.0) then
	 do i=1, nshg
	   if (y(i,6).LE.0.0) then
		y(i,id) = 373.15 ! Kelvin
	   endif
	 enddo
	endif

	return
	end

