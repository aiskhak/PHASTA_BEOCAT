        !COMPILER-GENERATED INTERFACE MODULE: Fri Jul 11 21:31:08 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE READTABLE_1__genmod
          INTERFACE 
            SUBROUTINE READTABLE_1(ISLOT,TABLE,NUMDATA,DX,MAXDATA,      &
     &MAXSLOTS)
              INTEGER(KIND=4) :: MAXSLOTS
              INTEGER(KIND=4) :: MAXDATA
              INTEGER(KIND=4) :: ISLOT
              REAL(KIND=8) :: TABLE(0:MAXDATA,2,MAXSLOTS)
              INTEGER(KIND=4) :: NUMDATA
              REAL(KIND=8) :: DX
            END SUBROUTINE READTABLE_1
          END INTERFACE 
        END MODULE READTABLE_1__genmod
