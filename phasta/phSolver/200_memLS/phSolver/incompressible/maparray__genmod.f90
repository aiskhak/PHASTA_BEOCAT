        !COMPILER-GENERATED INTERFACE MODULE: Fri Jul 11 21:31:30 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE MAPARRAY__genmod
          INTERFACE 
            SUBROUTINE MAPARRAY(X,X2,MAP,NSHL,NPROLD)
              INTEGER(KIND=4) :: NPROLD
              INTEGER(KIND=4) :: NSHL
              REAL(KIND=8) :: X(NPROLD,NSHL)
              REAL(KIND=8) :: X2(NPROLD,NSHL)
              INTEGER(KIND=4) :: MAP(NPROLD)
            END SUBROUTINE MAPARRAY
          END INTERFACE 
        END MODULE MAPARRAY__genmod
