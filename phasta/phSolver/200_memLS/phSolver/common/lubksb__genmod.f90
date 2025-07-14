        !COMPILER-GENERATED INTERFACE MODULE: Fri Jul 11 21:31:12 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE LUBKSB__genmod
          INTERFACE 
            SUBROUTINE LUBKSB(AA,N,NP,INDX,BB)
              INTEGER(KIND=4) :: NP
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: AA(NP,NP)
              INTEGER(KIND=4) :: INDX(N)
              REAL(KIND=8) :: BB(N)
            END SUBROUTINE LUBKSB
          END INTERFACE 
        END MODULE LUBKSB__genmod
