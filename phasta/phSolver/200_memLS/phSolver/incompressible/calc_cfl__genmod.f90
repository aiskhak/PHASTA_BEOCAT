        !COMPILER-GENERATED INTERFACE MODULE: Fri Jul 11 21:31:27 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE CALC_CFL__genmod
          INTERFACE 
            SUBROUTINE CALC_CFL(RHO,U1,U2,U3,DXIDX,RMU,CFL_LOC)
              COMMON/PROPAR/ NPRO
                INTEGER(KIND=4) :: NPRO
              REAL(KIND=8) :: RHO(NPRO)
              REAL(KIND=8) :: U1(NPRO)
              REAL(KIND=8) :: U2(NPRO)
              REAL(KIND=8) :: U3(NPRO)
              REAL(KIND=8) :: DXIDX(NPRO,3,3)
              REAL(KIND=8) :: RMU(NPRO)
              REAL(KIND=8) :: CFL_LOC(NPRO)
            END SUBROUTINE CALC_CFL
          END INTERFACE 
        END MODULE CALC_CFL__genmod
