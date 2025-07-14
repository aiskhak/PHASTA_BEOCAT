        !COMPILER-GENERATED INTERFACE MODULE: Fri Jul 11 21:31:14 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE E3METRICF__genmod
          INTERFACE 
            SUBROUTINE E3METRICF(XL,SHGL,DXIDX,SHG,WDETJ,QWTF)
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
              COMMON/TURBVAR/ ELES,YLIMIT,RMUTARGET,PZERO,WTAVEI,DTAVEI,&
     &DKE,FWR1,FLUMP,IERRCALC,IHESSIAN,ITWMOD,NGAUSSF,IDIM,NLIST,NINTF
                REAL(KIND=8) :: ELES
                REAL(KIND=8) :: YLIMIT(3,9)
                REAL(KIND=8) :: RMUTARGET
                REAL(KIND=8) :: PZERO
                REAL(KIND=8) :: WTAVEI
                REAL(KIND=8) :: DTAVEI
                REAL(KIND=8) :: DKE
                REAL(KIND=8) :: FWR1
                REAL(KIND=8) :: FLUMP
                INTEGER(KIND=4) :: IERRCALC
                INTEGER(KIND=4) :: IHESSIAN
                INTEGER(KIND=4) :: ITWMOD
                INTEGER(KIND=4) :: NGAUSSF
                INTEGER(KIND=4) :: IDIM
                INTEGER(KIND=4) :: NLIST
                INTEGER(KIND=4) :: NINTF(6)
              REAL(KIND=8) :: XL(NPRO,NENL,3)
              REAL(KIND=8) :: SHGL(NPRO,3,NSHL)
              REAL(KIND=8) :: DXIDX(NPRO,3,3)
              REAL(KIND=8) :: SHG(NPRO,NSHL,3)
              REAL(KIND=8) :: WDETJ(NPRO)
              REAL(KIND=8) :: QWTF(NGAUSSF)
            END SUBROUTINE E3METRICF
          END INTERFACE 
        END MODULE E3METRICF__genmod
