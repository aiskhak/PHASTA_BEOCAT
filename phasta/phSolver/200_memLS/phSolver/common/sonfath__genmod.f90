        !COMPILER-GENERATED INTERFACE MODULE: Fri Jul 11 21:31:13 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SONFATH__genmod
          INTERFACE 
            SUBROUTINE SONFATH(NCORP,VNUMNPL,IFATH,NSONS,TFATH,NUMPE,   &
     &NUMNP,MAXNODE)
              INTEGER(KIND=4) :: MAXNODE
              INTEGER(KIND=4) :: NUMNP
              INTEGER(KIND=4) :: NUMPE
              INTEGER(KIND=4) :: TFATH
              INTEGER(KIND=4) :: NCORP(NUMPE,MAXNODE)
              INTEGER(KIND=4) :: VNUMNPL(NUMPE)
              INTEGER(KIND=4) :: IFATH(NUMPE,MAXNODE)
              INTEGER(KIND=4) :: NSONS(TFATH)
            END SUBROUTINE SONFATH
          END INTERFACE 
        END MODULE SONFATH__genmod
