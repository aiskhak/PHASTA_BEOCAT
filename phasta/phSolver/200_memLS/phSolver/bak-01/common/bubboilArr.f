        module bubboil_info
!-------------------------------------------------------------------
!
!       This module defines some arrays and variables, which are co-
!       llected by bubble tracking and used in e3res.f, e3ivarf,etc 
!       to calculate bubble volume and number of shell element.
!       Mengnan Li,                                      Spring,2016
!
!-------------------------------------------------------------------
        real*8  bubble_vol(2) ! could contain 2 bubble at most
        real*8  numshell(2,2)
        real*8  bubble_tempG(2)
        end module

!-------------------------------------------------------------------
