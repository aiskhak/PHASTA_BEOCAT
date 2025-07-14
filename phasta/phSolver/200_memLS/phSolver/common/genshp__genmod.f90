        !COMPILER-GENERATED INTERFACE MODULE: Fri Jul 11 21:31:06 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE GENSHP__genmod
          INTERFACE 
            SUBROUTINE GENSHP(SHP,SHGL,NSHP,NBLK)
              REAL(KIND=8) :: SHP(6,32,125)
              REAL(KIND=8) :: SHGL(6,3,32,125)
              INTEGER(KIND=4) :: NSHP
              INTEGER(KIND=4) :: NBLK
            END SUBROUTINE GENSHP
          END INTERFACE 
        END MODULE GENSHP__genmod
