        !COMPILER-GENERATED INTERFACE MODULE: Fri Jul 11 21:31:00 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE E3METRIC__genmod
          INTERFACE 
            SUBROUTINE E3METRIC(XL,SHGL,DXIDX,SHG,WDETJ)
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
              REAL(KIND=8) :: XL(NPRO,NENL,3)
              REAL(KIND=8) :: SHGL(NPRO,3,NSHL)
              REAL(KIND=8) :: DXIDX(NPRO,3,3)
              REAL(KIND=8) :: SHG(NPRO,NSHL,3)
              REAL(KIND=8) :: WDETJ(NPRO)
            END SUBROUTINE E3METRIC
          END INTERFACE 
        END MODULE E3METRIC__genmod
