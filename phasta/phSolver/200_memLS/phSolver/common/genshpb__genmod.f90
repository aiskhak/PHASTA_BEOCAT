        !COMPILER-GENERATED INTERFACE MODULE: Fri Jul 11 21:31:06 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE GENSHPB__genmod
          INTERFACE 
            SUBROUTINE GENSHPB(SHPB,SHGLB,NSHPB,NBLK)
              REAL(KIND=8) :: SHPB(6,32,125)
              REAL(KIND=8) :: SHGLB(6,3,32,125)
              INTEGER(KIND=4) :: NSHPB
              INTEGER(KIND=4) :: NBLK
            END SUBROUTINE GENSHPB
          END INTERFACE 
        END MODULE GENSHPB__genmod
