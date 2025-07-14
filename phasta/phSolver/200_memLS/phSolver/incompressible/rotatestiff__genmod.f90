        !COMPILER-GENERATED INTERFACE MODULE: Fri Jul 11 21:31:30 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE ROTATESTIFF__genmod
          INTERFACE 
            SUBROUTINE ROTATESTIFF(RKLOCAL,ROTATION,RKGLOBAL)
              COMMON/PROPAR/ NPRO
                INTEGER(KIND=4) :: NPRO
              REAL(KIND=8) :: RKLOCAL(NPRO,3,3)
              REAL(KIND=8) :: ROTATION(NPRO,3,3)
              REAL(KIND=8) :: RKGLOBAL(NPRO,3,3)
            END SUBROUTINE ROTATESTIFF
          END INTERFACE 
        END MODULE ROTATESTIFF__genmod
