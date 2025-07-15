! file: levlset_mod.f90
      module levlset_mod
      use iso_c_binding, only: c_double, c_int
      implicit none

!---------------------------------------------------------------------
! This module defines the levlset common block as a C-bindable type.
! It must match exactly the layout in common_c.h's struct levlset.
!---------------------------------------------------------------------

      type, bind(C) :: levlset_type
! Fields in the same order as common_c.h
      real(c_double) epsilon_ls
      real(c_double) epsilon_lsd
      real(c_double) dtlset
      real(c_double) dtlset_cfl
      real(c_double) redist_toler
      real(c_double) redist_toler_curr
      real(c_double) r_int_buffer
      real(c_double) r_int_elem_size
      real(c_double) phvol(2)
      real(c_double) AdjRedistVelCFL
      real(c_double) BubRad
      real(c_double) vf_target
      real(c_double) C_int_adjust
      real(c_double) vf_now
      real(c_double) vfcontrcoeff
      real(c_double) C_int_cap
      real(c_double) epsilonBT
      real(c_double) xdistancesum
      real(c_double) ydistancesum
      real(c_double) zdistancesum
      real(c_double) totalxdist
      real(c_double) totalydist
      real(c_double) totalzdist
      real(c_double) avgxdistance
      real(c_double) avgydistance
      real(c_double) avgzdistance
      real(c_double) avgxdistold
      real(c_double) avgydistold
      real(c_double) avgzdistold
      real(c_double) xdistideal
      real(c_double) ydistideal
      real(c_double) zdistideal
      real(c_double) dx_new
      real(c_double) dy_new
      real(c_double) dz_new
      real(c_double) ddxvel
      real(c_double) ddyvel
      real(c_double) ddzvel
      real(c_double) xvelsum
      real(c_double) yvelsum
      real(c_double) zvelsum
      real(c_double) totalxvel
      real(c_double) totalyvel
      real(c_double) totalzvel
      real(c_double) avgxvel
      real(c_double) avgyvel
      real(c_double) avgzvel
      real(c_double) avgxvelold
      real(c_double) avgyvelold
      real(c_double) avgzvelold
      real(c_double) velwghtsum
      real(c_double) totalvelwght
      real(c_double) bubvolsum
      real(c_double) totbubvol
      real(c_double) denssum
      real(c_double) totbubdens
      real(c_double) xcforcesum
      real(c_double) ycforcesum
      real(c_double) zcforcesum
      real(c_double) totalxcforce
      real(c_double) totalycforce
      real(c_double) totalzcforce
      real(c_double) totalxcforceold
      real(c_double) totalycforceold
      real(c_double) totalzcforceold
      real(c_double) avgxcforce
      real(c_double) avgycforce
      real(c_double) avgzcforce
      real(c_double) avgxcforceold
      real(c_double) avgycforceold
      real(c_double) avgzcforceold
      real(c_double) avgxcf
      real(c_double) avgycf
      real(c_double) avgzcf
      real(c_double) x_c_f
      real(c_double) y_c_f
      real(c_double) z_c_f
      real(c_double) xforcenewtsum
      real(c_double) yforcenewtsum
      real(c_double) zforcenewtsum
      real(c_double) totxfnewtsum
      real(c_double) totyfnewtsum
      real(c_double) totzfnewtsum
      real(c_double) xcfnewtons
      real(c_double) ycfnewtons
      real(c_double) zcfnewtons
      real(c_double) rholiq
      real(c_double) rhogas
      real(c_double) coalbubrad
      real(c_double) avgxcoordold(100)
      real(c_double) avgycoordold(100)
      real(c_double) avgzcoordold(100)
      integer(c_int) i_res_cf
      integer(c_int) nzinBsum
      integer(c_int) ntotnzinB
      integer(c_int) iLSet
      integer(c_int) iuse_vfcont_cap
      integer(c_int) i_num_bubbles
      integer(c_int) ivconstraint
      integer(c_int) iSolvLSSclr1
      integer(c_int) iSolvLSSclr2
      integer(c_int) i_redist_loop_flag
      integer(c_int) i_redist_max_iter
      integer(c_int) i_spat_var_eps_flag
      integer(c_int) i_dtlset_cfl
      integer(c_int) i_check_prox
      integer(c_int) i_gradphi
      integer(c_int) i_focusredist
      integer(c_int) i_AdjRedistVel
      integer(c_int) buintcfl
      integer(c_int) iBT
      integer(c_int) iFT
      integer(c_int) icoalCtrl
      integer(c_int) coalcon
      integer(c_int) update_coalcon
      integer(c_int) coaltimtrak
      integer(c_int) coalest
      integer(c_int) coalcon_rem(100)
      end type levlset_type

! Bind the Fortran variable to the C symbol "levlset_"
      type(levlset_type), bind(C, name="levlset_") :: levlset

      end module levlset_mod