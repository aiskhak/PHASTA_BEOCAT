        !COMPILER-GENERATED INTERFACE MODULE: Fri Jul 11 21:31:03 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE CTYPES__genmod
          INTERFACE 
            SUBROUTINE CTYPES(ILWORK)
              COMMON/FRONTS/ MAXFRONT,NLWORK,IDIRSTEP,IDIRTRIGGER
                INTEGER(KIND=4) :: MAXFRONT
                INTEGER(KIND=4) :: NLWORK
                INTEGER(KIND=4) :: IDIRSTEP
                INTEGER(KIND=4) :: IDIRTRIGGER
              INTEGER(KIND=4) :: ILWORK(NLWORK)
            END SUBROUTINE CTYPES
          END INTERFACE 
        END MODULE CTYPES__genmod
