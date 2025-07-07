      subroutine e3ivar (yl,          acl,       shpfun,
     &                   shgl,        xl,        bml,
     &                   aci,         g1yi,      g2yi,    
     &                   g3yi,        shg,       dxidx,   
     &                   WdetJ,       rho,       pres, 
     &                   u1,          u2,        u3,              
     &                   ql,          rLui,      src,
     &                   rerrl,       rlsl,      rlsli,
     &                   dwl,         xwl,       ywl,
     &                   zwl,         elemvol_local, 
     &                   Sclr,        epsilon_ls_tmp, 
     &                   xarray,      yarray,    zarray,
     &                   bubradius,   bubradius2,
     &                   coordtag)
c
c!----------------------------------------------------------------------
c
c!  This routine computes the variables at integration point.
c
c! input:
c!  yl     (npro,nshl,ndof)      : primitive variables
c!  acl    (npro,nshl,ndof)      : prim.var. accel. 
c!  shp    (nen)                 : element shape-functions
c!  shgl   (nsd,nen)             : element local-grad-shape-functions
c!  xl     (npro,nenl,nsd)       : nodal coordinates at current step
c!  ql     (npro,nshl,nsd*nsd) : diffusive flux vector
c!  rlsl   (npro,nshl,6)       : resolved Leonard stresses
c
c! output:
c!  aci    (npro,3)              : primvar accel. variables 
c!  g1yi   (npro,ndof)           : grad-y in direction 1
c!  g2yi   (npro,ndof)           : grad-y in direction 2
c!  g3yi   (npro,ndof)           : grad-y in direction 3
c!  shg    (npro,nshl,nsd)       : element global grad-shape-functions
c!  dxidx  (npro,nsd,nsd)        : inverse of deformation gradient
c!  WdetJ  (npro)                : weighted Jacobian
c!  rho    (npro)                : density
c!  pres   (npro)                : pressure
c!  u1     (npro)                : x1-velocity component
c!  u2     (npro)                : x2-velocity component
c!  u3     (npro)                : x3-velocity component
c!  rLui   (npro,nsd)            : xi-momentum residual
c!  src    (npro,nsd)            : body force term (not density weighted)
c!  rlsli  (npro,6)              : resolved Leonard stresses at quad pt

c! locally calculated and used
c!  divqi  (npro,nsd+isurf)      : divergence of reconstructed quantity
c
c! Zdenek Johan, Summer 1990. (Modified from e2ivar.f)
c! Zdenek Johan, Winter 1991. (Fortran 90)
c! Kenneth Jansen, Winter 1997. Primitive Variables
c! Christian Whiting, Winter 1999. (uBar formulation)
c! Jun Fang, Spring 2015. (Bubble Tracking)
c
c! This routine identifies coalescense locations and implements the force
c! at the center of the event
c
c! output:
c!  xarray   (ibksiz)            : x coordinates of coalescence event
c!  yarray   (ibksiz)            : y coordinates of coalescence event
c!  zarray   (ibksiz)            : z coordinates of coalescence event
c!  coordtag (ibksiz)            : flag of a set of coordinates
c
c!----------------------------------------------------------------------
c
      use  spat_var_eps ! for spatially varying epsilon_ls
c      use  bub_track

      include "common.h"
c
c!  passed arrays
c
      dimension   yl(npro,nshl,ndof),        dwl(npro,nenl),       
     &            acl(npro,nshl,ndof),       shpfun(npro,nshl),
     &            shgl(npro,nsd,nshl),       xl(npro,nenl,nsd),
     &            aci(npro,nsd),             g1yi(npro,ndof),
     &            g2yi(npro,ndof),           g3yi(npro,ndof),
     &            shg(npro,nshl,nsd),        dxidx(npro,nsd,nsd),
     &            WdetJ(npro),               
     &            rho(npro),                 pres(npro),
     &            u1(npro),                  u2(npro),
     &            u3(npro),                  divqi(npro,nflow-1+isurf),
     &            ql(npro,nshl,idflx),       rLui(npro,nsd),
     &            src(npro,nsd), Temp(npro), xx(npro,nsd),
     &            xwl(npro,nenl), ywl(npro,nenl), zwl(npro,nenl)
c
        dimension tmp(npro), dkei(npro), dist2w(npro),
     &            x2w(npro), y2w(npro), z2w(npro)   ! Igor: vector towards a wall
c
        dimension rlsl(npro,nshl,6),         rlsli(npro,6)
c
        real*8    rerrl(npro,nshl,6), omega(3), divu(npro)

        dimension bml(npro,nshl,1)  ! Bubble tracking, Jun, Jul 2014 
        dimension bml_dummy(npro,nshl,1)
        real*8    elemvol_local(ibksiz)

        dimension gyti(npro,nsd),            gradh(npro,nsd),
     &            sforce(npro,3),            weber(npro),
     &            Sclr(npro)
        real*8    epsilon_ls_tmp, bml, bml_dummy
        real*8    alfaT, SurfDPind(npro,3),          !Farhad
     &            Xdis(npro), Ydis(npro), Zdis(npro), Rdis(npro) !Farhad
        real*8    VelMag, BubForce   ! Igor
        real*8    MaxCurv, Ccurv, MinCurv   ! Surface tension force limiter, Igor, April 2010   
        real*8    CurvInfo(npro)
	real*8    A11, B11  ! Wall repellant force now depends upon the Weber number, Igor, April 2010

        real*8    xarray(ibksiz), yarray(ibksiz), zarray(ibksiz) ! Matt Talley: Coalesc. Control
        real*8    bubradius, bubradius2, PtsPerDiam, Curvlim, Curvsel,
     &            appvolume(ibksiz,nenl,coalest) ! Matt Talley: Coalesc. Control
        integer   coordtag(ibksiz) ! Matt Talley: Coalesc. Control

        real*8    xxtmp(nsd), Wetmp !Jun Fang: Coal Ctrl

        real*8    fit
        dimension Tempb(npro), gytemp(npro,5) ! Mengnan: temperature
        real*8    Tempb, gytemp
        real*8    num
!c ...   Contact Angle Mengnan
        real*8    cos_contact, contact_angle,gcntang !From Anand: goal contact
!angle
        real*8    F_app,avg_CA,interval, mid_taninv,delta_contact !From Anand:
!Applied Force 040712
        real*8, dimension(npro, 3):: cntfdir ! From Anand: contact force
!direction
        logical, dimension(npro)::CA_achieved1,CA_achieved2 !From Anand
        integer s,app_points
        real*8    Fapp_thick,Fapp_height,constant,min_delta,stretch,offset
        integer  istp_dummy, istp !From Anand       
        real*8, dimension(ibksiz):: CA_vector, CAnode !From Anand
        real*8, dimension(ibksiz):: contact_force_vector
        integer, dimension(ibksiz)::tag
        real*8 measure_thick, measure_height, eps_tmp, mod_surften !surface
!tension
        real*8 mod_cntfdir, mod_wnrm,eps_correction
        REAL*8, DIMENSION(3) :: F_applied
        integer flag
!        integer contactangle_flag !LBW
        real*8 delta_high, delta_low, mid_high, mid_low, CA_constant, CA_temp
        real*8 delta_const_min, delta_const_max, offset_high, offset_low
        real*8 CLS_max, CLS, proj_length !contact line speed
        real*8 avgspeed !average bubble speed
        integer bubnode !nodes within the bubble
        integer velodel ! velocity theta dependence model
        real*8  eps_ave, eps_sum, upperlimit, CA_alpha, lowerlimit
        integer cnt
        real*8 cntfdir2(npro,3)
        real*8 wnrm_new(3)
        real*8, dimension(npro):: site_loct
!        real*8, dimension(npro):: site_1, site_2 !Anna 
 


c... Contact Angle Mengnan
! -----------------------------------------------------------------------
!------------------------------------------------------------------------
!       if (istp.eq.1.and.iblk.eq.1)then                  !From Anand

!            open(12,file='contact_angle_parameter.dat',status='old') !From
!Anand

!            read(12,*) gcntang,interval,min_delta,mid_taninv,Fapp_thick,
!           &Fapp_height,constant

!           contactangle_flag = 0 !LBW
!           contactangle_flag = 0
!           if(myrank.eq.master)write(*,*) Forcecont,stretch
           if (CA_flag.eq.1.0d0) then !LBW Summer 2014
           gcntang = 0.0 !60.0 !110.0 !60.0

           interval=5.0 !15.0 !20.0 !10.0
           min_delta=25.0 !25.0
           mid_taninv=15.0 !15.0
!           Fapp_thick= 2.0 !2.0      !should be 2.0 !1.0/2.5 worked well
!           Fapp_height= 2.0 !2.0         !should be 4.0 !0.5!1.0
           CAnode = zero
                                   !8E11/10**4/10 !4E9/10**3!4E9 !4.5E9 
!           stretch= 100 !1.0
           flag=1
!           theta_adv=40.0  !90.0
!           theta_rec=40.0 !90.0
           velmodel=1 !decide which vel model 1: Sine interpolation, 2: Linear
!3: Step
!###   Implementaiton of dynamic contact angle

                delta_high=-0.5   !5!20.0 !5 worked well !
                mid_high=20
                delta_const_max=40

                delta_low=-0.5 !5!20.0 !-5 worked well
                     ! signed minimum deviation from static contact angle = 40
                mid_low=-20
                delta_const_min=-40
                CA_alpha = 0.2  ! Relaxation factor

               offset_low=(1/2-1/3.141*atan((delta_low-mid_low)/stretch))
!offset along the y axis to ensure continuity of the contact fforce funciton

               offset_high=(1/2+1/3.141*atan((delta_high-mid_high)/stretch))

              ! Cavities location
              !if (lstep.eq.4500) then
              !if (myrank.eq.master) write(*,*) 'Assign sites locations '
              site_loct(:) = 0.d0
              do i = 1, npro
                   if ((xx(i,2).ge.2.85d-3).or.(xx(i,2).le.0.15d-3)) then
                      site_loct(i) = 1.d0
                   endif
              !    ! #1 and #3
              !    if (xx(i,2).ge.2.7d-3) then
              !        ! #1
              !        !if ((xx(i,1).ge.-5.7d-3).and.(xx(i,1).le.-5.0d-3)) then
              !        if ((xx(i,1).ge.-5.95d-3).and.(xx(i,1).le.-4.55d-3)) then
              !            site_loct(i) = 1.d0
              !        endif
              !        ! #3
              !        !if ((xx(i,1).ge.-2.5d-3).and.(xx(i,1).le.-1.8d-3)) then
              !        if ((xx(i,1).ge.-2.95d-3).and.(xx(i,1).le.-1.55d-3)) then
              !            site_loct(i) = 1.d0
              !        endif
              !     ! #2 and #4
              !     elseif (xx(i,2).le.0.3d-3) then
              !        ! #2
              !        !if ((xx(i,1).ge.-3.9d-3).and.(xx(i,1).le.-3.2d-3)) then
              !        if ((xx(i,1).ge.-4.45d-3).and.(xx(i,1).le.-3.05d-3)) then
              !            site_loct(i) = 1.d0
              !        endif
              !        ! #4
              !        !if ((xx(i,1).ge.-1.1d-3).and.(xx(i,1).le.-0.4d-3)) then
              !        if ((xx(i,1).ge.-1.45d-3).and.(xx(i,1).le.-0.05d-3)) then
              !            site_loct(i) = 1.d0
              !        endif
              !    endif
              enddo
              !!endif

!the elusive offset of the tan inverse function along the y axis
            endif
!Ramping constact angle control force
           if (CAramp.eq.1.0d0)then
              if (CA_temp.lt.Forcecont)CA_temp = Forcecont
              CA_temp = CA_temp + CAF_interval*(lstep-CA_startT)
              if(CA_temp .gt. CAF_upper)then
                CA_constant = CAF_upper
              else
                CA_constant = CA_temp
              endif
           else
              CA_constant = Forcecont
           endif
          
           
c
c!.... ------------->  Primitive variables at int. point  <--------------
c
c!.... compute primitive variables
c
       SurfDPind        = zero  !Farhad
       CurvInfo         = zero
       xx               = zero
       pres             = zero
       u1               = zero
       u2               = zero
       u3               = zero
c
       do n = 1, nshl 
          pres = pres + shpfun(:,n) * yl(:,n,1)
          u1   = u1   + shpfun(:,n) * yl(:,n,2)
          u2   = u2   + shpfun(:,n) * yl(:,n,3)
          u3   = u3   + shpfun(:,n) * yl(:,n,4)
!       xx can be used in  user-specified body force or coriolis force specified
!       , coalescence control, lift/drag control and also bubble tracking
          xx(:,1) = xx(:,1)  + shpfun(:,n) * xl(:,n,1)
          xx(:,2) = xx(:,2)  + shpfun(:,n) * xl(:,n,2)
          xx(:,3) = xx(:,3)  + shpfun(:,n) * xl(:,n,3)
       enddo
       if(matflg(5,1).eq.2) then ! boussinesq body force
          Temp = zero
          do n = 1, nshl
             Temp = Temp + shpfun(:,n) * yl(:,n,5)
          enddo
       endif
c
       if(iRANS.eq.-2.or.(iDNS.ne.0.and.abs(itwmod).eq.1)
     &    .or.iBT.eq.1) then ! kay-epsilon & Slip-velocity DNS
          dist2w = zero
          do n = 1, nenl
             dist2w = dist2w + shpfun(:,n) * dwl(:,n)
          enddo
       endif

c! Caculate distance to the wall for Contact Angle:
       if(CA_flag.eq.1.0d0) then
          dist2w = zero
          do n = 1, nenl
             dist2w = dist2w + shpfun(:,n) * dwl(:,n)
          enddo
       endif


c! Slip velocity wall DNS distance to the wall/wall vector:
       if (iDNS.ne.0.and.abs(itwmod).eq.1) then
          x2w = zero
          y2w = zero
          z2w = zero
          do n = 1, nenl
             x2w = x2w + shpfun(:,n) * xwl(:,n)
             y2w = y2w + shpfun(:,n) * ywl(:,n)
             z2w = z2w + shpfun(:,n) * zwl(:,n)
          enddo
       end if

        
!          if (lstep.lt.1) then
!            do i = 1, npro
!                   open(unit=1991, file='d2wall.txt')
!                   write(*,118) dist2w(i), x2w(i), y2w(i),z2w(i)
!            enddo
!          endif
       

 118  format(ES14.7, 1x, ES14.7, 1x, ES14.7, 1x, ES14.7)


c! Caculate distance to the wall for Contact Angle:
       if(CA_flag.eq.1.0d0) then
          x2w = zero
          y2w = zero
          z2w = zero
          do n = 1, nenl
             x2w = x2w + shpfun(:,n) * xwl(:,n)
             y2w = y2w + shpfun(:,n) * ywl(:,n)
             z2w = z2w + shpfun(:,n) * zwl(:,n)
          enddo
       end if


c
       if((iLES.gt.10).and.(iLES.lt.20))  then  ! bardina
          rlsli = zero
          do n = 1, nshl 

             rlsli(:,1) = rlsli(:,1) + shpfun(:,n) * rlsl(:,n,1)
             rlsli(:,2) = rlsli(:,2) + shpfun(:,n) * rlsl(:,n,2)
             rlsli(:,3) = rlsli(:,3) + shpfun(:,n) * rlsl(:,n,3)
             rlsli(:,4) = rlsli(:,4) + shpfun(:,n) * rlsl(:,n,4)
             rlsli(:,5) = rlsli(:,5) + shpfun(:,n) * rlsl(:,n,5)
             rlsli(:,6) = rlsli(:,6) + shpfun(:,n) * rlsl(:,n,6)

          enddo
       else
          rlsli = zero
       endif
c
c!.... ----------------------->  accel. at int. point  <----------------------
c
       aci = zero
       do n = 1, nshl
          aci(:,1) = aci(:,1) + shpfun(:,n) * acl(:,n,2)
          aci(:,2) = aci(:,2) + shpfun(:,n) * acl(:,n,3)
          aci(:,3) = aci(:,3) + shpfun(:,n) * acl(:,n,4)
       enddo
c
c!.... --------------------->  Element Metrics  <-----------------------
c
       call e3metric( xl,         shgl,       dxidx,  
     &                shg,        WdetJ)
c
c!.... compute the global gradient of u and P
c
c
       g1yi = zero
       g2yi = zero
       g3yi = zero
       do n = 1, nshl
          g1yi(:,1) = g1yi(:,1) + shg(:,n,1) * yl(:,n,1)
          g1yi(:,2) = g1yi(:,2) + shg(:,n,1) * yl(:,n,2)
          g1yi(:,3) = g1yi(:,3) + shg(:,n,1) * yl(:,n,3)
          g1yi(:,4) = g1yi(:,4) + shg(:,n,1) * yl(:,n,4)
c
          g2yi(:,1) = g2yi(:,1) + shg(:,n,2) * yl(:,n,1)
          g2yi(:,2) = g2yi(:,2) + shg(:,n,2) * yl(:,n,2)
          g2yi(:,3) = g2yi(:,3) + shg(:,n,2) * yl(:,n,3)
          g2yi(:,4) = g2yi(:,4) + shg(:,n,2) * yl(:,n,4)
c
          g3yi(:,1) = g3yi(:,1) + shg(:,n,3) * yl(:,n,1)
          g3yi(:,2) = g3yi(:,2) + shg(:,n,3) * yl(:,n,2)
          g3yi(:,3) = g3yi(:,3) + shg(:,n,3) * yl(:,n,3)
          g3yi(:,4) = g3yi(:,4) + shg(:,n,3) * yl(:,n,4)
       enddo

       divqi = zero
       idflow = 3
       if ( idiff >= 1 .or. isurf==1 ) then
c     
c!.... compute divergence of diffusive flux vector, qi,i
c     
          if(idiff >= 1) then
             do n=1, nshl
                divqi(:,1) = divqi(:,1) + shg(:,n,1)*ql(:,n,1 ) 
     &                                  + shg(:,n,2)*ql(:,n,4 )
     &                                  + shg(:,n,3)*ql(:,n,7 )

                divqi(:,2) = divqi(:,2) + shg(:,n,1)*ql(:,n,2 ) 
     &                                  + shg(:,n,2)*ql(:,n,5 )
     &                                  + shg(:,n,3)*ql(:,n,8)

                divqi(:,3) = divqi(:,3) + shg(:,n,1)*ql(:,n,3 ) 
     &                                  + shg(:,n,2)*ql(:,n,6 )
     &                                  + shg(:,n,3)*ql(:,n,9 )

          enddo

          endif                 !end of idiff
c     
          if (isurf .eq. 1) then   
c!     .... divergence of normal calculation (curvature)
c!  note ql(:,:,idflx-2 through idflx) contains the normal vector
c!  e3q computed grad(phi) and qpbc divided by the magnitude of grad(phi).
             do n=1, nshl
                divqi(:,idflow+1) = divqi(:,idflow+1) 
     &               + shg(:,n,1)*ql(:,n,idflx-2)
     &               + shg(:,n,2)*ql(:,n,idflx-1)
     &               + shg(:,n,3)*ql(:,n,idflx)
             enddo 
c!     .... initialization of some variables
             Sclr = zero
             gradh= zero
             gyti = zero
             sforce=zero
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
c
c!..Now we have to use spatially varied epsilon_ls
c
                epsilon_ls_tmp = epsilon_ls * 
     &               elem_local_size(lcblk(1,iblk)+i-1)

  
                if (abs (sclr(i)) .le. epsilon_ls_tmp) then
                   gradh(i,1) = 0.5/epsilon_ls_tmp * (1.0 
     &                  + cos(pi*Sclr(i)/epsilon_ls_tmp)) * gyti(i,1)
                   gradh(i,2) = 0.5/epsilon_ls_tmp * (1.0 
     &                  + cos(pi*Sclr(i)/epsilon_ls_tmp)) * gyti(i,2)
                   gradh(i,3) = 0.5/epsilon_ls_tmp * (1.0 
     &                  + cos(pi*Sclr(i)/epsilon_ls_tmp)) * gyti(i,3)
                endif
             enddo              !end of the loop over npro

!        do n = 1,nshl
!        if (bml(i,n,1).ne.0.0d0)then
!        write(*,*)bml(i,n,1)
!        endif
!        enddo

c     
c! .. surface tension force calculation
c! .. divide by density now as it gets multiplied in e3res.f, as surface
c!    tension force is already in the form of force per unti volume
c    
c!....Temperature gradient calculation --Mengnan 9/4/2015......................!
       if (bubboil.eq.1.0d0.or.bubgrow.eq.1.0d0)then
!       write(*,*)"I'm here"
        call Bubheatflux(yl, shpfun,  shg, elemvol_local,Tempb, gytemp,bml)

!        do i = 1,npro
!        do n = 1,nshl
!        if (bml(i,n,1).ne.0.0d0)then
!        write(*,*)bml(i,n,1)
!        endif
!        enddo
!        enddo

        else
        Tempb=zero
        gytemp(:,:)=zero
       endif
c!.....end of bubble boiling...................................................!
c!.... Matt Talley's Bubble Coalescence Control
      if (coalcon.eq.1) then
         bubradius = coalbubrad + 6.0d0*(epsilon_ls_tmp)
         bubradius2 = coalbubrad + epsilon_ls_tmp
         PtsPerDiam = (2.0d0*coalbubrad) / (epsilon_ls_tmp/1.8d0)
         Curvlim = (-1.8083d2 * PtsPerDiam + 1.2873d3) 
      endif

c!	MaxCurv = 30.0   ! Max curvature allowed for survace tension force usage, Igor, April 2010
        do i = 1, npro
           weber(i) = Bo
           Ccurv = divqi(i,idflow+1)
c!.... Matt Talley's Bubble Coalescence Contorl
        if (coalcon.eq.1) then
              if (((Sclr(i).ge.(1.0d0*epsilon_ls_tmp)).and.
     &        (Sclr(i).le.(6.0d0*epsilon_ls_tmp)))
     &        .and.((Ccurv.ge.(-1.45d0*Curvlim)).or.
     &        (Ccurv.le.(1.45d0*Curvlim)))) then

                 xarray(i) = xx(i,1)
                 yarray(i) = xx(i,2)
                 zarray(i) = xx(i,3)
                 coordtag(i) = 1

              endif

         if (coaltimtrak.eq.1) then
             do k = 1, coalest
                do n = 1,nenl
                  appvolume(i,n,k) = 
     &            sqrt((xl(i,n,1)-avgxcoordold(k))**2 +
     &                 (xl(i,n,2)-avgycoordold(k))**2 +
     &                 (xl(i,n,3)-avgzcoordold(k))**2)
                  if (appvolume(i,n,k).le.5.0d0*epsilon_ls_tmp) then
                       weber(i) = CoalInvSigma !Tension
                  endif
                enddo
              enddo
         else
            do k = 1, coalest
               if(coalcon_rem(k).eq.0) then
                  do n = 1,nenl
                  appvolume(i,n,k) = 
     &            sqrt((xl(i,n,1)-avgxcoordold(k))**2 +
     &                 (xl(i,n,2)-avgycoordold(k))**2 +
     &                 (xl(i,n,3)-avgzcoordold(k))**2)
                  if(appvolume(i,n,k).le.5.0d0*epsilon_ls_tmp) then
                     weber(i) = CoalInvSigma !Tension
                  endif
                  enddo
               endif
            enddo
         endif ! coaltimtrak
        endif ! coalcon

c!      MaxCurv = 700.0   ! Max curvature allowed for survace tension force usage, Igor, April 2010
c        do i = 1, npro
c           weber(i) = Bo
c           Ccurv = divqi(i,idflow+1)
c!.... Matt Talley's Bubble Coalescence Contorl
        if (coalcon.eq.1) then
              if (((Sclr(i).ge.(1.0d0*epsilon_ls_tmp)).and.
     &        (Sclr(i).le.(6.0d0*epsilon_ls_tmp)))
     &        .and.((Ccurv.ge.(-1.45d0*Curvlim)).or.
     &        (Ccurv.le.(1.45d0*Curvlim)))) then

                 xarray(i) = xx(i,1)
                 yarray(i) = xx(i,2)
                 zarray(i) = xx(i,3)
                 coordtag(i) = 1

              endif

         if (coaltimtrak.eq.1) then
             do k = 1, coalest
                do n = 1,nenl
                  appvolume(i,n,k) =
     &            sqrt((xl(i,n,1)-avgxcoordold(k))**2 +
     &                 (xl(i,n,2)-avgycoordold(k))**2 +
     &                 (xl(i,n,3)-avgzcoordold(k))**2)
                  if (appvolume(i,n,k).le.5.0d0*epsilon_ls_tmp) then
                       weber(i) = CoalInvSigma !Tension
                  endif
                enddo
              enddo
         else
            do k = 1, coalest
               if(coalcon_rem(k).eq.0) then
                  do n = 1,nenl
                  appvolume(i,n,k) =
     &            sqrt((xl(i,n,1)-avgxcoordold(k))**2 +
     &                 (xl(i,n,2)-avgycoordold(k))**2 +
     &                 (xl(i,n,3)-avgzcoordold(k))**2)
                  if(appvolume(i,n,k).le.5.0d0*epsilon_ls_tmp) then
                     weber(i) = CoalInvSigma !Tension
                  endif
                  enddo
               endif
            enddo
         endif ! coaltimtrak
        endif ! coalcon
!       Jun Fang Coalescence Control based on bubble tracking
        if(icoalCtrl.eq.1.and.ncoalEvent.ge.1) then
!        if(myrank.eq.master)write(*,*)'I am in e3ivar'
           xxtmp(:)     = xx(i,:)
           Wetmp        = weber(i)
           call coalCtrlApp(xxtmp, Wetmp)
           weber(i)     = Wetmp
        endif



!        if(abs(sclr(i)).le.eps_thick.and.(dist2w(i).le.eps_height))then
!         MaxCurv = 10000.0d0
!         MinCurv = -10000.0d0
!         if (Ccurv.gt.MaxCurv) then
!            Ccurv = MaxCurv
!         elseif(Ccurv.lt.MinCurv)then
!            Ccurv = MinCurv
!         end if
!        endif
! New version to apply the 2x surface tension ONLY inside the gas part of the
! interface:
!          if (Sclr(i).le.0.0E0) then
!             sforce(i,1) = -(2.0/weber(i)) * Ccurv !x-direction
!     &            *gradh(i,1) /rho(i)
!             sforce(i,2) = -(2.0/weber(i)) * Ccurv !y-direction
!     &            *gradh(i,2) /rho(i)
!             sforce(i,3) = -(2.0/weber(i)) * Ccurv !z-direction
!     &            *gradh(i,3) /rho(i)
!          end if
c... Contact Angle Mengnan
           if (CA_flag.eq.1.0d0) then
           node = 0
           avgspeed = 0.0d0
           bubnode = 0
           endif

! Original version for the surface tension model:
             sforce(i,1) = (-1.0/weber(i)) * Ccurv !x-direction
     &            *gradh(i,1) /rho(i)
             sforce(i,2) = (-1.0/weber(i)) * Ccurv !y-direction
     &            *gradh(i,2) /rho(i)
             sforce(i,3) = (-1.0/weber(i)) * Ccurv !z-direction
     &            *gradh(i,3) /rho(i)
!        end do   ! i, npro
c ...   Contact Angle Mengnan  --------------------------------------------
        if(CA_flag.eq.1.0d0) then
        
        if (dist2w(i).ne.0.d0) then
        wnrm_new(1) = (-x2w(i))/(dist2w(i))
        wnrm_new(2) = (-y2w(i))/(dist2w(i))
        wnrm_new(3) = (-z2w(i))/(dist2w(i))
        endif

!          if (lstep.lt.1) then
!                   open(unit=1991, file='d2wall.txt')
!                   write(1991,118) dist2w(i), x2w(i), y2w(i),z2w(i)
!          endif

! 118  format(ES14.7, 1x, ES14.7, 1x, ES14.7, 1x, ES14.7)


!        if(wnrm_new(2).ne.0.0d0)then
!        if(myrank.eq.master) write(*,*)'wnrm_new(2)', wnrm_new(2)
!        if(myrank.eq.master) write(*,*)'y2w(i)', y2w(i)
!        endif
        proj_length=wnrm_new(1)*gyti(i,1)+wnrm_new(2)*gyti(i,2)+
     &              wnrm_new(3)*gyti(i,3) 
        cntfdir(i,:)=gyti(i,:)-wnrm_new(:)*proj_length
! opoints towards outward direction of bubble
        cntfdir(i,:)=cntfdir(i,:)
     &  /sqrt(cntfdir(i,1)**2+cntfdir(i,2)**2+cntfdir(i,3)**2)
!dynamic contact angle implentation contact line speed
       CLS=u1(i)*cntfdir(i,1)+u2(i)*cntfdir(i,2)+u3(i)*cntfdir(i,3)

!model
       CLS_max=1.0E-2
       if(velmodel.eq.1) then ! sine model

         if(CLS.ge.CLS_max) gcntang=theta_adv
         if(CLS.le.-1.0*CLS_max) gcntang=theta_rec
        if(CLS.gt.-1.0*CLS_max.and.CLS.lt.CLS_max) then
       gcntang=theta_adv/2.0*(sin(3.141*CLS/(2.0*CLS_max))+1)
     &         -theta_rec/2.0*(sin(3.141*CLS/(2.0*CLS_max))-1)
        endif
!      if(myrank.eq.master) Write(*,*)  "Sine model",velmodel

       elseif(velmodel.eq.2) then ! linear model

         if(CLS.ge.CLS_max) gcntang=theta_adv
         if(CLS.le.-1.0*CLS_max) gcntang=theta_rec
        if(CLS.gt.-1.0*CLS_max.and.CLS.lt.CLS_max) then
       gcntang=(theta_adv-theta_rec)/(2*CLS_max)*CLS
     &         +(theta_rec +theta_adv)/2.0
        endif
       elseif(velmodel.eq.3) then ! step model

!      if(myrank.eq.master) Write(*,*)  "step model",velmodel

         if(CLS.ge.0.0) gcntang=theta_adv
         if(CLS.le.0.0) gcntang=theta_rec

       endif !end selection of models
!      Print *, "Contact Angle Section Line 515"
!      if(gcntang.gt.theta_rec) Print *, gcntang

! initially before looping over the nodes in a block avgspeed was zero now it
! will add up within the bubble

       if(sclr(i).le.0.0) then
          avgspeed=avgspeed+sqrt(u1(i)*u1(i)+u2(i)*u2(i)+u3(i)*u3(i))
          bubnode=bubnode+1
       endif
!       bubnode=bubnode+1


! find contact angle
       if (dist2w(i).ne.0.d0) then
       cos_contact=(gyti(i,1)*x2w(i)+gyti(i,2)*y2w(i)+gyti(i,3)*z2w(i))
     & /((gyti(i,1))**2.0D0+(gyti(i,2))**2.0D0+(gyti(i,3))**2.0D0)**0.5
     & /(dist2w(i))
        endif

        !if(isnan(cos_contact)) cos_contact = 0.0d0

!       if(myrank.eq.master)write(*,*)'cos_contact',cos_contact

       if  (cos_contact.lt.0.d0) then

           contact_angle=180.0-acos(-cos_contact)*180.0/3.141
       else
           contact_angle=acos(cos_contact)*180.0/3.141
       endif

       !upperlimit = max(bubble_avgCA, gcntang)
       !lowerlimit = min(gcntang, bubble_avgCA)
!       if (bubble_avgCA.ge.gcntang)then
! 
         !if((contact_angle.gt.upperlimit).and.(upperlimit.ne.0.0))then
         !contact_angle = (1-CA_alpha)*bubble_avgCA + CA_alpha*contact_angle
         ! elseif(contact_angle.lt.lowerlimit)then
         !contact_angle = (1-CA_alpha)*bubble_avgCA + CA_alpha*contact_angle         
         !endif
         
!        else

!         if((contact_angle.gt.upperlimit).and.(upperlimit.ne.0.0))then
!         contact_angle = CA_alpha*bubble_avgCA + (1-CA_alpha)*contact_angle
!         elseif(contact_angle.lt.gcntang)then
!         contact_angle = CA_alpha*bubble_avgCA + (1-CA_alpha)*contact_angle
!         endif

!       endif
  
!       if((contact_angle.gt.upperlimit).and.(upperlimit.ne.0.0))then
!         contact_angle = upperlimit
!        elseif(contact_angle.le.lowerlimit)then
!         contact_angle = lowerlimit
!       elseif((contact_angle.gt.(bubble_avgCA*0.88))
!     &  .and.(contact_angle.lt.bubble_avgCA))then
!         contact_angle = lowerlimit
!        else
!         contact_angle = bubble_avgCA
!       endif 

!        if(myrank.eq.master)write(*,*)"bubble_avgCA",bubble_avgCA

       delta_contact=(gcntang-contact_angle)

       if (i_spat_var_eps_flag.eq.0) then
            eps_thick = epsilon_ls*Fapp_thick
            eps_height = epsilon_ls*Fapp_heigh
        else
            epsilon_ls_tmp = epsilon_ls*
     &      elem_local_size(lcblk(1,iblk)+i-1)
            eps_thick = epsilon_ls_tmp*Fapp_thick
            eps_height = epsilon_ls_tmp*Fapp_heigh
        endif
       
 
            do n = 1, nshl
                !if ((bml(i,n,1).gt.0.d0).and.() then
                if ((bml(i,n,1).gt.0.d0).and.(site_loct(i).eq.1.d0)) then
            !if ((xx(i,2).lt.0.d0).or.(xx(i,2).gt.2.8d-3)) then
                eps_thick = epsilon_ls_tmp*1.4d0
                eps_height = epsilon_ls_tmp*1.4d0

                goto 333
          ! elseif (bml(i,n,1).eq.3) then
                !eps_thick = epsilon_ls_tmp*2.5d0
                !eps_height = epsilon_ls_tmp*3.5d0
                !goto 334
            else
                delta_contact = 0.d0
            endif
            enddo

 333    continue
! 334    continue

         !if (xx(i,1).ge.4.0d-4) then
         !    delta_contact = 0.d0
         !endif
        
        !if ((xx(i,1).ge.-3.0d-3).and.(xx(i,2).le.1.5d-3)) then
         !    delta_contact = 0.d0
         !endif
 
!-----------------------initial guess ends------------------------------------
!        Print *, "Contact Angle Section Line 543"

! find if current contact angle is what we want
        if((delta_contact.gt.delta_low)
     &     .and.(delta_contact.lt.delta_high))then
                          CA_achieved2(i)=.true.
        elseif (delta_contact.eq.0.d0) then
                          CA_achieved2(i)=.true.
        else
                          CA_achieved2(i)=.false.
        endif
! if we are not in the force application region then we don't consider contact
! angle and node. so putting them to zero:
                    F_app=0.d0
                    CA_vector(i)=0.d0
                    contact_force_vector(i)=0.d0

                    tag(i)=0  ! tag 0 is not in advancing region, 1 advancing, 2
                              ! receding
   !                     Print *, "Contact Angle Section Line 560"
c------------------------------------------------------------------------------
c_______________________________________________________________________________
!             Print *, "Contact Angle Section Line 567"
! FORCE APPLICATION REGION
!          if (1.eq.0) then ! Initial force formulation to disable force
!          if ((sclr(i)).le.Fapp_thick*eps.and.

!           if(sclr(i).le.Fapp_thick*eps.and.sclr(i).ge.0.0.and.
           if(abs(sclr(i)).le.eps_thick.and.(dist2w(i).le.eps_height)
     &         .and.(.not.CA_achieved2(i)))then

              if(delta_contact.le.delta_const_min) then
                 delta_contact=delta_const_min
              elseif(delta_contact.ge.delta_const_max) then
                 delta_contact=delta_const_max
              endif

!............INITIAL FPRCE FORMULATION.........................................


              if (delta_contact.le.delta_low) then ! push interface in
                 F_app=(1.d0/(eps_thick*eps_height**3))*CA_constant*
     &                 (-1.0/2.0+1.0/3.141*atan((delta_contact-mid_low)/stretch)
     &                 +offset_low)*cos(sclr(i)*3.141/2.0/(eps_thick))
     &                 *(eps_height-dist2w(i))**2.0*rho(i)
              elseif(delta_contact.ge.delta_high) then ! pull interface out
                 F_app=(1.d0/(eps_thick*eps_height**3))*CA_constant*
     &                 (1.0/2.0+1.0/3.141*atan((delta_contact-mid_high)/stretch)
     &                  -offset_high)*cos(sclr(i)*3.141/2.0/(eps_thick))
     &                  *(eps_height-dist2w(i))**2.0*rho(i)
               !write(*,*) 'TH^3',
!     &                     eps_thick*eps_height**3 !Anna

!   initial static formulation:
!      F_app=constant*(1.0/2.0+1.0/3.141*atan((delta_contact
!     &-mid_taninv)/stretch)-offset)
!     &*cos(sclr(i)*3.141/2/(2*Fapp_thick*eps))
!     &*(Fapp_height*eps-dist2w(i))**2.0*rho(i)

               else

                       F_app=0.0
               endif
!        if(myrank.eq.master)write(*,*)'F_app', F_app
 !         Print *, "Contact Angle Section Line 607"
!        Print *, "gcntang",gcntang
! Transferring force to sforce so that its incorporated

!          if (lstep.eq.1) then
!                   open(unit=3470, file='raw_sforce.txt')
!                   write(3470,*) sforce(i,1),sforce(i,2), sforce(i,3)
!          endif
         ! end of the Contact angle force formulation, out of region now
              !if ((xx(i,1).ge.3.8d-4).or.(xx(i,2).lt.-1.d-8)) then
              !    F_app=0.d0
              !endif
            !do n = 1, nshl
            !if (bml(i,n,1).gt.1) then
            !    F_app = 0.d0                
            !endif
            !enddo
  
              sforce(i,:) = sforce(i,:) + F_app*cntfdir(i,:)


         if(abs(sclr(i)).le.eps_thick.and.
     &     (dist2w(i).le.eps_height).and.
     &     (.not.CA_achieved2(i))) then
         if (lstep.eq.5) then
!              open(unit=3460, file="FAR_Coord.txt")
!              open(unit=3461, file="CA_Force.txt")
!              open(unit=3462, file="wnrm.txt")
!              open(unit=3463, file="wnrm_new.txt")
!              open(unit=3464, file="TotalForce.txt")
!              open(unit=3465, file="Eps.txt")
!              open(unit=3466, file='Contact_angle.txt')
!              open(unit=3467, file='CAF_direction.txt')
!     
!              write(3460,*) xl(i,1,1),xl(i,1,2),xl(i,1,3)
!              write(3461,*) F_app
!              write(3462,*) wnrm(i,1), wnrm(i,2), wnrm(i,3)
!              write(3463,*) wnrm_new(1), wnrm_new(2), wnrm_new(3)
!              write(3464,*) sforce(i,1), sforce(i,2), sforce(i,3)
!              write(3465,*) epsilon_ls_tmp, eps_thick, eps_height
!              write(3466,*) contact_angle, delta_contact
!              write(3467,*) cntfdir(i,1),cntfdir(i,2), cntfdir(i,3)
         endif
!        if(lstep.eq.10)then
!              open(unit=3467, file='Contact_angle10.txt')
!              write(3467,*) contact_angle, delta_contact
!        endif

!        if(lstep.eq.200)then
!              open(unit=3468, file='Contact_angle200.txt')
!              write(3468,*) contact_angle, delta_contact
!        endif

        
!        if(lstep.eq.500)then
!              open(unit=3469, file='Contact_angle500.txt')
!              write(3469,*) contact_angle, delta_contact
!        endif


!        if(lstep.eq.1500)then
!              open(unit=3470, file='Contact_angle1500.txt')
!              write(3470,*) contact_angle, delta_contact
!        endif


         ! out of Force application region, end of writing
         endif
! taking care of force in the direction normal to the wall
!             sforce(i,:)=sforce(i,:)
!     &   -(sforce(i,1)*wnrm(i,1)+sforce(i,2)*wnrm(i,2)
!     &   +sforce(i,3)*wnrm(i,3))/(dist2w(i))**2*wnrm(i,:)
!          Print *, "Contact Angle Section Line 649"

!      ENDIF ! 1eq 0OMMENTING OUT THE INITIAL FORCE FORMULATION
!-----------------INIITIAL FORCE
!FORMULATION-------------------------------------------


       endif ! ENDING FORCE APPLICATION REGION IF STATEMENT
!                Print *, F_app,rho(i),cntang,gcntang
!                 Print *, "gcntang", gcntang !LBW
  
!        if(sclr(i).le.Fapp_thick*eps.and.sclr(i).ge.0.0.and.
!     &          (dist2w(i).le.Fapp_height*eps)) then
         if(abs(sclr(i)).le.eps_thick.and.
     &     (dist2w(i).le.eps_height).and.
     &     (.not.CA_achieved2(i))) then 
                CAnode(i)=1.0d0
                CA_vector(i)=contact_angle
                contact_force_vector(i)=F_app
         endif
         
         !do n = 1, nshl
         !   if (bml(i,n,1).gt.1) then
         !       CAnode(i)=0.d0
         !       CA_vector(i)=0.d0
         !       contact_force_vector(i)=0.d0
         !    endif
         !enddo         

         !if ((xx(i,1).ge.3.8d-4).or.(xx(i,2).lt.-1.d-8)) then
             !write (*,*)'y=', xx(i,2)
         !    CAnode(i)=0.d0
         !    CA_vector(i)=0.d0
         !    contact_force_vector(i)=0.d0
         !endif

        endif !contactangle_flag
        
!       Jun Fang Coalescence Control based on bubble tracking
        if(icoalCtrl.eq.1.and.ncoalEvent.ge.1) then
!        if(myrank.eq.master)write(*,*)'I am in e3ivar'
           xxtmp(:)     = xx(i,:)
           Wetmp        = weber(i) 
           call coalCtrlApp(xxtmp, Wetmp)
           weber(i)     = Wetmp
        endif
!	  if (Ccurv.gt.MaxCurv) then
!	     Ccurv = MaxCurv
!	  end if
! New version to apply the 2x surface tension ONLY inside the gas part of the interface:
!          if (Sclr(i).le.0.0E0) then
!             sforce(i,1) = -(2.0/weber(i)) * Ccurv !x-direction
!     &            *gradh(i,1) /rho(i)
!             sforce(i,2) = -(2.0/weber(i)) * Ccurv !y-direction
!     &            *gradh(i,2) /rho(i)
!             sforce(i,3) = -(2.0/weber(i)) * Ccurv !z-direction
!     &            *gradh(i,3) /rho(i)
!          end if
! Original version for the surface tension model:
!             sforce(i,1) = (-1.0/weber(i)) * Ccurv !x-direction
!     &            *gradh(i,1) /rho(i)
!             sforce(i,2) = (-1.0/weber(i)) * Ccurv !y-direction
!     &            *gradh(i,2) /rho(i)
!             sforce(i,3) = (-1.0/weber(i)) * Ccurv !z-direction
!     &            *gradh(i,3) /rho(i)          

        end do   ! i, npro

!----------------------------------------------------------------------
!       collect information for lift/drag control
!----------------------------------------------------------------------
        if(iCForz.eq.1)then
           call CFCollect(xx,        u1,     u2,     u3,
     &                    Sclr,      epsilon_ls_tmp, elemvol_local,
     &                    rho,       sforce)
        end if !iCForz
!----------------------------------------------------------------------
!       collect information for advanced analysis. 
!----------------------------------------------------------------------
        if(iBT .eq. 1) then
           CurvInfo(:) = divqi(:,idflow+1)
           call BubCollect(u1,       u2,     u3,     Sclr, dist2w,
     &                     xx,       yl,     bml,    elemvol_local,
     &                     rho,      Tempb,  gytemp, CurvInfo,
     &                     CAnode,   CA_vector)
        bml_dummy = bml
        endif !iBT 

!----------------------------------------------------------------------
!....Temperature gradient calculation --Mengnan 9/4/2015......................!
!       if (bubboil.eq.1.0d0.or.bubgrow.eq.1.0d0)then
!        write(*,*)"I'm here"
!        call Bubheatflux(yl, shpfun,  shg, elemvol_local,Tempb, gytemp,bml_dummy)
!
!        else
!        Tempb=zero
!        gytemp(:,:)=zero
!       endif
!-----------------------------------------------------------------------------

c
c!            Xdis=zero                                               !Farhad
c!            Ydis=zero                                               !Farhad
c!            Zdis=zero                                               !Farhad
c!             do n = 1, nshl                                         !Farhad
c!                Xdis(:)= Xdis(:) +  shpfun(:,n)*xl(:,n,1)           !Farhad
c!                Ydis(:)= Ydis(:) +  shpfun(:,n)*xl(:,n,2)           !Farhad
c!                Zdis(:)= Zdis(:) +  shpfun(:,n)*xl(:,n,3)           !Farhad
c!             end do                                                 !Farhad
c!             do i = 1, npro                                         !Farhad
c!                if ((Xdis(i).ge.0.37).and.(Xdis(i).le.0.395)) then  !Farhad
c!                   alfaT = (Xdis(i)-0.37)/(0.395-0.37)              !Farhad
c!                   alfaT = 2.0*(alfaT **3) - 3.0*(alfaT**2) + 1.0   !Farhad
c!                   sforce(i,1) = alfaT * sforce(i,1)                !Farhad
c!                   sforce(i,2) = alfaT * sforce(i,2)                !Farhad
c!                   sforce(i,3) = alfaT * sforce(i,3)                !Farhad
c!                elseif (Xdis(i).gt.0.395) then                      !Farhad
c!                   sforce(i,1) = 0.0                                !Farhad
c!                   sforce(i,2) = 0.0                                !Farhad
c!                   sforce(i,3) = 0.0                                !Farhad
c!                endif                                               !Farhad
c!             enddo                                                  !Farhad


          endif        ! end of the surface tension force calculation
       endif           ! diffusive flux computation
c
c! Calculate strong form of pde as well as the source term
c     
c!      SurfDPind= zero         
c!          if (Xdis(i)>0.0077) then                                    !Farhad
c!             if ((Sclr(i)<0).and.(abs(Sclr(i))>epsilon_ls_tmp)) then !Farhad
c!                SurfDPind(i,1)= g1yi(i,1)                      !Farhad
cc!                SurfDPind(i,2)= g2yi(i,1)                      !Farhahd
cc!                SurfDPind(i,3)= g3yi(i,1)                      !Farhad
c!              elseif (abs (sclr(i)) .le. epsilon_ls_tmp) then         !Farhad
c!                SurfDPind(i,1)=(g1yi(i,1)-(rho(i)*sforce(i,1))) !Farhad
c!                SurfDPind(i,2)=(g2yi(i,1)-(rho(i)*sforce(i,2))) !Farhad
c!                SurfDPind(i,3)=(g3yi(i,1)-(rho(i)*sforce(i,3))) !Farhad
c!             endif                                                   !Farhad
c!          endif                                                      !Farhad
c!       enddo                                                         !Farhad


    ! Coordinates of nucleation sites

!        if (lstep.eq.0) then
!        if (myrank.eq.master) then
!       site_1(:) = zero
!          site_2(:) = zero
!           do i=1, npro
!           if ((xx(i,1).ge.-4.0d-3).and.
!     &          (xx(i,1).le.-3.42d-3)) then
!             site_1(i) = xx(i,1)
!            endif
!
!            if ((xx(i,1).ge.-2.6d-3).and.
!     &          (xx(i,1).le.-1.89d-3)) then
!             site_2(i) = xx(i,1)
!            endif
!            enddo
!
!           call MPI_Bcast(site_1,npro,MPI_DOUBLE_PRECISION,
!     &                   master,MPI_COMM_WORLD,ierr)
!           call MPI_Bcast(site_2,npro,MPI_DOUBLE_PRECISION,
!     &                   master,MPI_COMM_WORLD,ierr) 
!    
!          endif
!          endif


!   This is a subgrid wall force which prevents bubbles from attaching to the wall (Igor A. Bolotnov, Summer 2009) 
      SurfDPind= zero   
      if (iDNS.ne.0.and.abs(itwmod).eq.1) then ! works in case of DNS and slip velocity wall modeling
!       if (myrank.eq.master) write(*,*) 'e3ivar version 07/22/2009'
       do i=1, npro
            if (site_loct(i).eq.0.d0) then
          ! if ((xx(i,1).le.4.6d-3).or.(xx(i,1).ge.5.4d-3)) then
           !if ((xx(i,2).gt.0.0d0).and.(xx(i,2).lt.3.0d-3)) then
!        write(*,*) 'sclr, d2w, y2w, eps = ', sclr(i), dist2w(i), y2w(i), epsilon_ls_tmp       
!              if (abs (sclr(i)).le.epsilon_ls_tmp.and.sclr(i).gt.0.0.and.dist2w(i).lt.4.0*epsilon_ls_tmp) then  ! liquid part of the interface
        if (abs (sclr(i)).le.epsilon_ls_tmp.and.dist2w(i)
     1   .lt.8.0*epsilon_ls_tmp) then  ! whole interface, about 10 elements from the wall 
!              if (abs (sclr(i)).le.epsilon_ls_tmp.and.dist2w(i).lt.16.0*epsilon_ls_tmp) then  ! whole interface, about 10 elements from the wall
!              if (abs (sclr(i)).le.epsilon_ls_tmp.and.dist2w(i).lt.4.0*epsilon_ls_tmp) then  ! whole interface
!              if (sclr(i).le.0.0.and.dist2w(i).lt.4.0*epsilon_ls_tmp) then    ! The force is applied to the internal volume of the bubbles
                VelMag = 1.0  ! Temp fix (!)  sqrt(u1(i)**2.0+u2(i)**2.0+u3(i)**2.0) ! Should be parallel to the wall, however that may cause trouble in the jet case
                 BubForce = datmat(1,2,1)*VelMag*
     1 BubRad*(550.0/dist2w(i)+35.0/(dist2w(i)**2.0))   ! 1/d2 dependence is added here
! Currently disabled (04/29/2011) to make sure Re = 180 runs are done for BA70 cases (Bubble/Air at 70 atm)
! A11 / B11 introduce Weber number dependence in the wall force computation:, Igor 04/26/2010
!              A11 = 7200.0 / sqrt(Bo)  !15000 10000  9200  fixed rad/diam 3600.0 / sqrt(Bo)   ! Bo is the Weber number from the input
!              B11 = 1200.0 / sqrt(Bo)    ! 2000.0 1200;    600, fixed radius/diameter issue Increased from 200.0
!                BubForce = datmat(1,2,1)*VelMag*BubRad*(A11/dist2w(i)+B11/(dist2w(i)**2.0))   ! 1/d2 dependence is added here
!              write(*,*) 'myrank, i, dist2w, BubForce = ', myrank, i, dist2w(i), BubForce, BubRad
                if (dist2w(i).ne.0.d0) then
                SurfDPind(i,1) = -1.0E0*BubForce*x2w(i)/dist2w(i) *rho(i)   !   A negative since is needed since the force
                SurfDPind(i,2) = -1.0E0*BubForce*y2w(i)/dist2w(i) *rho(i)   !   is acting away from the wall whilst
                SurfDPind(i,3) = -1.0E0*BubForce*z2w(i)/dist2w(i) *rho(i)   !   the normal vector is directed towards the wall
                endif
!	  if (myrank.eq.16.and.i.gt.800) write(*,*) 'mr, i, y2w, SurfDPind, sforce', myrank, i, y2w(i), SurfDPind(i,2), sforce(i,2)*rho(i)
! Increase surface tension in this region as well:
!		sforce(i,1:3) = 2.0*sforce(i,1:3)
             endif
             endif
       enddo                                                  
      end if  ! iDNS condition

       call e3resStrongPDE(
     &      aci,  u1,   u2,   u3,   Temp, rho,  xx, SurfDPind,
     &            g1yi, g2yi, g3yi,
     &      rLui, src, divqi, sforce)
      
c
c!.... -------------------> error calculation  <-----------------
c     
       if((ierrcalc.eq.1).and.(nitr.eq.iter)) then
          do ia=1,nshl
             tmp=shpfun(:,ia)*WdetJ(:)
             rerrl(:,ia,1) = rerrl(:,ia,1) +
     &                       tmp(:)*rLui(:,1)*rLui(:,1)
             rerrl(:,ia,2) = rerrl(:,ia,2) +
     &                       tmp(:)*rLui(:,2)*rLui(:,2)
             rerrl(:,ia,3) = rerrl(:,ia,3) +
     &                       tmp(:)*rLui(:,3)*rLui(:,3)

             rerrl(:,ia,4) = rerrl(:,ia,4) +
     &                       tmp(:)*divqi(:,1)*divqi(:,1)
             rerrl(:,ia,5) = rerrl(:,ia,5) +
     &                       tmp(:)*divqi(:,2)*divqi(:,2)
             rerrl(:,ia,6) = rerrl(:,ia,6) +
     &                       tmp(:)*divqi(:,3)*divqi(:,3)
          enddo
       endif
       distcalc=0  ! return to 1 if you want to compute T-S instability
       if(distcalc.eq.1) then
c
c!.... ----------------------->  dist. kin energy at int. point  <--------------
c
       
       if (ires .ne. 2 .and. iter.eq.1)  then  !only do at beginning of step
c
c! calc exact velocity for a channel at quadrature points.
c
       dkei=0.0
c
       do n = 1, nenl 
          dkei = dkei + shpfun(:,n) * (1.0-xl(:,n,2)**2) !u_ex^~ (in FEM space)
       enddo
          dkei = (u1(:)-dkei)**2 +u2(:)**2  ! u'^2+v'^2
          dkei = dkei*WdetJ  ! mult function*W*det of jacobian to
c!                              get this quadrature point contribution
          dke  = dke+sum(dkei) ! we move the sum over elements inside of the
c!                              sum over quadrature to save memory (we want
c!                              a scalar only)
       endif
       endif
c     
c!.... return
c
       return
       end

c!-----------------------------------------------------------------------
c 
c!     Calculate the variables for the scalar advection-diffusion
c!     equation.
c
c!-----------------------------------------------------------------------
      subroutine e3ivarSclr (yl,          acl,       shpfun,
     &                      shgl,        xl,        xmudmi,
     &                      Sclr,        Sdot,      gradS,  
     &                      shg,         dxidx,     WdetJ,
     &                      u1,          u2,        u3,              
     &                      ql,          rLS ,       SrcR,
     &                      SrcL,        uMod,      dwl,
     &                      diffus,      srcRat,    evl,
     &                      cfll )
c
      use spat_var_eps   ! use spatially-varying epl_ls
c
      include "common.h"
c
c!  passed arrays
c
      dimension yl(npro,nshl,ndof),        acl(npro,nshl,ndof), 
     &          Sclr(npro),                Sdot(npro),
     &          gradS(npro,nsd),           shpfun(npro,nshl),
     &          shgl(npro,nsd,nshl),       xl(npro,nenl,nsd),
     &          shg(npro,nshl,nsd),        dxidx(npro,nsd,nsd),
     &          WdetJ(npro),              
     &          u1(npro),                  u2(npro),
     &          u3(npro),                  divS(npro),
     &          ql(npro,nshl,nsd),         rLS(npro),
     &          SrcR(npro),                 SrcL(npro),
     &          dwl(npro,nshl),            diffus(npro),
     &          umod(npro,nsd),            Temp(npro),
     &          xx(npro,nsd),              evl(npro,nenl),
     &          divqi(npro,nsd),
     &          cfll(npro,nshl) 
c
      dimension tmp(npro), srcRat(npro)
      real*8 rLui(npro,nsd),     aci(npro,nsd),
     &       g1yi(npro,nflow),   g2yi(npro,nflow),
     &       g3yi(npro,nflow),
     &       src(npro,nsd),      rho(npro),
     &       rmu(npro),          cp(npro),
     &       k_T(npro)
      real*8 uBar(npro,nsd), xmudmi(npro,ngauss)
      real*8 epsilon_ls_tmp, SurfDPind(npro,3), Xdis(npro) !Farhad
c
c!.... ------------->  Primitive variables at int. point  <--------------
c
c!.... compute primitive variables
c
      u1   = zero
      u2   = zero
      u3   = zero
      Sclr = zero
      SurfDPind = zero !Farhad
      Xdis = zero      !Farhad
c
      id=isclr+5
      do n = 1, nshl 
         u1   = u1   + shpfun(:,n) * yl(:,n,2)
         u2   = u2   + shpfun(:,n) * yl(:,n,3)
         u3   = u3   + shpfun(:,n) * yl(:,n,4)
         Sclr = Sclr + shpfun(:,n) * yl(:,n,id)
      enddo
c
c
c!.... ----------------------->  dS/dt at int. point  <----------------------
c
      Sdot = zero
      do n = 1, nshl
         Sdot = Sdot + shpfun(:,n) * acl(:,n,id)
      enddo
c
c!.... --------------------->  Element Metrics  <-----------------------
c

      call e3metric( xl,         shgl,        dxidx,  
     &               shg,        WdetJ)

c
c!.... compute the global gradient of u and P
c

       gradS = zero
       do n = 1, nshl
          gradS(:,1) = gradS(:,1) + shg(:,n,1) * yl(:,n,id)
          gradS(:,2) = gradS(:,2) + shg(:,n,2) * yl(:,n,id)
          gradS(:,3) = gradS(:,3) + shg(:,n,3) * yl(:,n,id)
       enddo

       divS = zero
       if ( idiff >= 1 ) then
c
c!.... compute divergence of diffusive flux vector, qi,i
c
          do n=1, nshl
             divS(:) = divS(:) + shg(:,n,1)*ql(:,n,1 ) 
     &                         + shg(:,n,2)*ql(:,n,2 ) 
     &                         + shg(:,n,3)*ql(:,n,3 ) 
          enddo
       endif                    ! diffusive flux computation

       if(consrv_sclr_conv_vel.eq.1) then
c!         Calculate uBar = u - TauM*L, where TauM is the momentum
c!         stabilization factor and L is the momentum residual

          if(matflg(5,1).eq.2) then ! boussinesq body force
             Temp = zero
             do n = 1, nshl
                Temp = Temp + shpfun(:,n) * yl(:,n,5)
             enddo
          endif
          if(matflg(5,1).eq.3.or.matflg(6,1).eq.1) then
c!     user-specified body force or coriolis force specified
             xx = zero
             do n  = 1,nenl
                xx(:,1) = xx(:,1)  + shpfun(:,n) * xl(:,n,1)
                xx(:,2) = xx(:,2)  + shpfun(:,n) * xl(:,n,2)
                xx(:,3) = xx(:,3)  + shpfun(:,n) * xl(:,n,3)
             enddo
          endif
          aci = zero
          do n = 1, nshl
             aci(:,1) = aci(:,1) + shpfun(:,n) * acl(:,n,2)
             aci(:,2) = aci(:,2) + shpfun(:,n) * acl(:,n,3)
             aci(:,3) = aci(:,3) + shpfun(:,n) * acl(:,n,4)
          enddo
          g1yi = zero
          g2yi = zero
          g3yi = zero
          do n = 1, nshl
             g1yi(:,1) = g1yi(:,1) + shg(:,n,1) * yl(:,n,1)
             g1yi(:,2) = g1yi(:,2) + shg(:,n,1) * yl(:,n,2)
     
        g1yi(:,3) = g1yi(:,3) + shg(:,n,1) * yl(:,n,3)
             g1yi(:,4) = g1yi(:,4) + shg(:,n,1) * yl(:,n,4)
c     
             g2yi(:,1) = g2yi(:,1) + shg(:,n,2) * yl(:,n,1)
             g2yi(:,2) = g2yi(:,2) + shg(:,n,2) * yl(:,n,2)
             g2yi(:,3) = g2yi(:,3) + shg(:,n,2) * yl(:,n,3)
             g2yi(:,4) = g2yi(:,4) + shg(:,n,2) * yl(:,n,4)
c     
             g3yi(:,1) = g3yi(:,1) + shg(:,n,3) * yl(:,n,1)
             g3yi(:,2) = g3yi(:,2) + shg(:,n,3) * yl(:,n,2)
             g3yi(:,3) = g3yi(:,3) + shg(:,n,3) * yl(:,n,3)
             g3yi(:,4) = g3yi(:,4) + shg(:,n,3) * yl(:,n,4)
          enddo
c          
          if (iLSet .eq. 0)then
             rho  = datmat(1,1,1)
             rmu = datmat(1,2,1)
          else
             write(*,*) 'Not sure if we can handle level set with K-E'
             write(*,*) '(different uMods? correct value of rho?)'
          endif
          divqi=zero  ! until we reconstruct q_flow for scalar solve

          call e3resStrongPDE(
     &         aci,  u1,   u2,   u3,   Temp, rho,  x,  SurfDPind,     !Farhad
     &               g1yi, g2yi, g3yi,
     &         rLui, src, divqi)
          src(:,1)=u1           !
          src(:,2)=u2           ! store u in src memory
          src(:,3)=u3           !
c!         e3uBar calculates Tau_M and assembles uBar
          call getdiff(dwl, yl, shpfun, xmudmi, xl, rmu, rho,
     &               cp,  k_T,
     &               elem_local_size(lcblk(1,iblk):lcblk(1,iblk+1)-1),
     &               evl)
          call e3uBar(rho, src, dxidx, rLui, rmu, uBar)
          u1=ubar(:,1)          ! the entire scalar residual
          u2=ubar(:,2)          ! is based on the modified
          u3=ubar(:,3)          ! velocity for conservation
       endif
c
c!.... Initialize uMod, the modified velocity uMod
c!      We initialize it to u_i and then calculate
c!      the correction in e3sourcesclr
c

       umod(:,1) = u1
       umod(:,2) = u2
       umod(:,3) = u3
c     
c!.... compute  source terms
c
cad
cad    if we are solving the redistancing equation, the umod(:,:) are 
CAD    modified in e3sourceSclr.  
CAD
CAD  if we are redistancing levelset variable we want to use a use the  
CAD  convective term from the equation.  


       if(nosource.ne.1) then
        call e3sourceSclr ( Sclr,         Sdot,      gradS,  dwl,
     &                      shpfun,       shg,       yl,     dxidx,
     &                      diffus,       u1,        u2,     u3,
     &                      cfll,         Temp,      srcR,   srcL,
     &                      uMod,         srcRat )
       else
        srcRat = zero
        srcR   = zero
        srcL   = zero
       endif
c
c!.... -------------------> Scalar residual  <-----------------
c

        rLS(:) = ( Sdot(:) +  (u1*gradS(:,1) + 
     &                          u2*gradS(:,2) +
     &                           u3*gradS(:,3)) )
     &                                          - divS(:)           

c
c!.... return
c
       return
       end

