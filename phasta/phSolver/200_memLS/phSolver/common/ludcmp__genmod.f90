        !COMPILER-GENERATED INTERFACE MODULE: Fri Jul 11 21:31:14 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE LUDCMP__genmod
          INTERFACE 
            SUBROUTINE LUDCMP(AA,N,NP,INDX,D)
              INTEGER(KIND=4) :: NP
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: AA(NP,NP)
              INTEGER(KIND=4) :: INDX(N)
              REAL(KIND=8) :: D
            END SUBROUTINE LUDCMP
          END INTERFACE 
        END MODULE LUDCMP__genmod
