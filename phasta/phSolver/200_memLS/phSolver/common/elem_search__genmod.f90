        !COMPILER-GENERATED INTERFACE MODULE: Fri Jul 11 21:31:02 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE ELEM_SEARCH__genmod
          INTERFACE 
            SUBROUTINE ELEM_SEARCH(XL,XTS1,XTS2,XTS3,XSIC,ELMT,IDFILE)
              USE SPEBC
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
              REAL(KIND=8) :: XL(NELINT,NENL,3)
              REAL(KIND=8) :: XTS1
              REAL(KIND=8) :: XTS2
              REAL(KIND=8) :: XTS3
              REAL(KIND=8) :: XSIC(3)
              INTEGER(KIND=4) :: ELMT
              INTEGER(KIND=4) :: IDFILE
            END SUBROUTINE ELEM_SEARCH
          END INTERFACE 
        END MODULE ELEM_SEARCH__genmod
