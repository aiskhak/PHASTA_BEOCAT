        !COMPILER-GENERATED INTERFACE MODULE: Fri Jul 11 21:31:03 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE CMASS__genmod
          INTERFACE 
            SUBROUTINE CMASS(SHP,SHGL,XL,EM)
              COMMON/SHPDAT/ NSHAPE,NSHAPEB,MAXSHB,NSHL,NSHLB,NFATH,    &
     &NTOPSH,NSONMAX
                INTEGER(KIND=4) :: NSHAPE
                INTEGER(KIND=4) :: NSHAPEB
                INTEGER(KIND=4) :: MAXSHB
                INTEGER(KIND=4) :: NSHL
                INTEGER(KIND=4) :: NSHLB
                INTEGER(KIND=4) :: NFATH
                INTEGER(KIND=4) :: NTOPSH
                INTEGER(KIND=4) :: NSONMAX
              COMMON/PROPAR/ NPRO
                INTEGER(KIND=4) :: NPRO
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
              REAL(KIND=8) :: SHP(NSHL,125)
              REAL(KIND=8) :: SHGL(3,NSHL,125)
              REAL(KIND=8) :: XL(NPRO,NENL,3)
              REAL(KIND=8) :: EM(NPRO,NSHL,NSHL)
            END SUBROUTINE CMASS
          END INTERFACE 
        END MODULE CMASS__genmod
