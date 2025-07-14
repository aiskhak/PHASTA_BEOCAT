        !COMPILER-GENERATED INTERFACE MODULE: Fri Jul 11 21:31:18 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE GENSIDCOUNT__genmod
          INTERFACE 
            SUBROUTINE GENSIDCOUNT(NSIDG)
              COMMON/WORKFC/ MASTER,NUMPE,MYRANK
                INTEGER(KIND=4) :: MASTER
                INTEGER(KIND=4) :: NUMPE
                INTEGER(KIND=4) :: MYRANK
              INTEGER(KIND=4) :: NSIDG
            END SUBROUTINE GENSIDCOUNT
          END INTERFACE 
        END MODULE GENSIDCOUNT__genmod
