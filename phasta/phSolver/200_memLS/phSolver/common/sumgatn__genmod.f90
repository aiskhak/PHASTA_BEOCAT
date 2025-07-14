        !COMPILER-GENERATED INTERFACE MODULE: Fri Jul 11 21:31:19 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SUMGATN__genmod
          INTERFACE 
            SUBROUTINE SUMGATN(U,N,SUMMED,NNP)
              COMMON/FRONTS/ MAXFRONT,NLWORK,IDIRSTEP,IDIRTRIGGER
                INTEGER(KIND=4) :: MAXFRONT
                INTEGER(KIND=4) :: NLWORK
                INTEGER(KIND=4) :: IDIRSTEP
                INTEGER(KIND=4) :: IDIRTRIGGER
              INTEGER(KIND=4) :: NNP
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: U(NNP,N)
              REAL(KIND=8) :: SUMMED
            END SUBROUTINE SUMGATN
          END INTERFACE 
        END MODULE SUMGATN__genmod
