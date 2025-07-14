        !COMPILER-GENERATED INTERFACE MODULE: Fri Jul 11 21:31:27 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COALESCAPPTIME__genmod
          INTERFACE 
            SUBROUTINE COALESCAPPTIME(AVGXCOORDF,AVGYCOORDF,AVGZCOORDF, &
     &AVGXCOORDOLD2,AVGYCOORDOLD2,AVGZCOORDOLD2,APP_TIME,ITRTIMESTP)
              REAL(KIND=8) :: AVGXCOORDF(100)
              REAL(KIND=8) :: AVGYCOORDF(100)
              REAL(KIND=8) :: AVGZCOORDF(100)
              REAL(KIND=8) :: AVGXCOORDOLD2(100)
              REAL(KIND=8) :: AVGYCOORDOLD2(100)
              REAL(KIND=8) :: AVGZCOORDOLD2(100)
              REAL(KIND=8) :: APP_TIME(100,2)
              REAL(KIND=8) :: ITRTIMESTP
            END SUBROUTINE COALESCAPPTIME
          END INTERFACE 
        END MODULE COALESCAPPTIME__genmod
