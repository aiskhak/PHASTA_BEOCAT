        !COMPILER-GENERATED INTERFACE MODULE: Fri Jul 11 21:31:07 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE READTABLE_2__genmod
          INTERFACE 
            SUBROUTINE READTABLE_2(ISLOT,TABLE,NUMDATA,DX,MAXDATA,      &
     &MAXSLOTS)
              INTEGER(KIND=4) :: MAXSLOTS
              INTEGER(KIND=4) :: MAXDATA
              INTEGER(KIND=4) :: ISLOT
              REAL(KIND=8) :: TABLE(4,0:MAXDATA,0:MAXDATA,MAXSLOTS)
              INTEGER(KIND=4) :: NUMDATA(2,MAXSLOTS)
              REAL(KIND=8) :: DX(2,MAXSLOTS)
            END SUBROUTINE READTABLE_2
          END INTERFACE 
        END MODULE READTABLE_2__genmod
