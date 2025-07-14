        !COMPILER-GENERATED INTERFACE MODULE: Fri Jul 11 21:31:20 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE PHIST__genmod
          INTERFACE 
            SUBROUTINE PHIST(PRESSHIST,QHIST,BETAS,NTIMEPOINT,NSRFS)
              INTEGER(KIND=4) :: NSRFS
              INTEGER(KIND=4) :: NTIMEPOINT
              REAL(KIND=8) :: PRESSHIST(0:20)
              REAL(KIND=8) :: QHIST(NTIMEPOINT+1,NSRFS)
              REAL(KIND=8) :: BETAS(NTIMEPOINT+2,NSRFS)
            END SUBROUTINE PHIST
          END INTERFACE 
        END MODULE PHIST__genmod
