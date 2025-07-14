        !COMPILER-GENERATED INTERFACE MODULE: Fri Jul 11 21:31:06 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE GET_COEFF_TET__genmod
          INTERFACE 
            SUBROUTINE GET_COEFF_TET(XC,ANOT)
              USE SPEBC
              COMMON/ELMPAR/ LELCAT,LCSYST,IORDER,NENB,NELBLK,NELBLB,   &
     &NDOFL,NSYMDL,NENL,NFACEL,NENBL,INTIND,MATTYP
                INTEGER(KIND=4) :: LELCAT
                INTEGER(KIND=4) :: LCSYST
                INTEGER(KIND=4) :: IORDER
                INTEGER(KIND=4) :: NENB
                INTEGER(KIND=4) :: NELBLK
                INTEGER(KIND=4) :: NELBLB
                INTEGER(KIND=4) :: NDOFL
                INTEGER(KIND=4) :: NSYMDL
                INTEGER(KIND=4) :: NENL
                INTEGER(KIND=4) :: NFACEL
                INTEGER(KIND=4) :: NENBL
                INTEGER(KIND=4) :: INTIND
                INTEGER(KIND=4) :: MATTYP
              REAL(KIND=8) :: XC(NELINT,NENL,3)
              REAL(KIND=8) :: ANOT(NELINT,NENL,3)
            END SUBROUTINE GET_COEFF_TET
          END INTERFACE 
        END MODULE GET_COEFF_TET__genmod
