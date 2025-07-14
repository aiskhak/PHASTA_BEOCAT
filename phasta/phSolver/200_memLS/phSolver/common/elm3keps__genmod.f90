        !COMPILER-GENERATED INTERFACE MODULE: Fri Jul 11 21:31:00 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE ELM3KEPS__genmod
          INTERFACE 
            SUBROUTINE ELM3KEPS(KAY,EPSILON,DWI,GRADV,SRCRAT1,SRC1,     &
     &SRCJAC)
              COMMON/PROPAR/ NPRO
                INTEGER(KIND=4) :: NPRO
              REAL(KIND=8) :: KAY(NPRO)
              REAL(KIND=8) :: EPSILON(NPRO)
              REAL(KIND=8) :: DWI(NPRO)
              REAL(KIND=8) :: GRADV(NPRO,3,3)
              REAL(KIND=8) :: SRCRAT1(NPRO)
              REAL(KIND=8) :: SRC1(NPRO)
              REAL(KIND=8) :: SRCJAC(NPRO,4)
            END SUBROUTINE ELM3KEPS
          END INTERFACE 
        END MODULE ELM3KEPS__genmod
