        !COMPILER-GENERATED INTERFACE MODULE: Fri Jul 11 21:31:27 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE E3DCSCLR__genmod
          INTERFACE 
            SUBROUTINE E3DCSCLR(GRADS,GIJU,GGRADS,RLS,TAUS,SRCR,DCFCT)
              COMMON/PROPAR/ NPRO
                INTEGER(KIND=4) :: NPRO
              REAL(KIND=8) :: GRADS(NPRO,3)
              REAL(KIND=8) :: GIJU(NPRO,6)
              REAL(KIND=8) :: GGRADS(NPRO,3)
              REAL(KIND=8) :: RLS(NPRO)
              REAL(KIND=8) :: TAUS(NPRO)
              REAL(KIND=8) :: SRCR(NPRO)
              REAL(KIND=8) :: DCFCT(NPRO)
            END SUBROUTINE E3DCSCLR
          END INTERFACE 
        END MODULE E3DCSCLR__genmod
