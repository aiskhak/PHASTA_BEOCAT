        !COMPILER-GENERATED INTERFACE MODULE: Fri Jul 11 21:31:21 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE USRBC__genmod
          INTERFACE 
            SUBROUTINE USRBC(CFLAG,RX,TIME1,XVAL,NSIZE)
              INTEGER(KIND=4) :: NSIZE
              CHARACTER(LEN=8) :: CFLAG
              REAL(KIND=8) :: RX(3)
              REAL(KIND=8) :: TIME1
              REAL(KIND=8) :: XVAL(NSIZE)
            END SUBROUTINE USRBC
          END INTERFACE 
        END MODULE USRBC__genmod
