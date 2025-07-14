        !COMPILER-GENERATED INTERFACE MODULE: Fri Jul 11 21:31:30 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE ASIGMRSCLR__genmod
          INTERFACE 
            SUBROUTINE ASIGMRSCLR(Y,AC,X,SHP,SHGL,IEN,RES,QRES,XSEBE,   &
     &XMUDMI,CFL,ICFLHITS,CFLOLD)
              COMMON/CONPAR/ NUMNP,NUMEL,NUMELB,NUMPBC,NEN,NFACES,NUMFLX&
     &,NDOF,IALE,ICOORD,NAVIER,IBLK,IRS,IEXEC,NECHO,ICHEM,IRK,NEDOF,NSHG&
     &,NNZ,ISTOP,NFLOW,NNZ_TOT,IDTN
                INTEGER(KIND=4) :: NUMNP
                INTEGER(KIND=4) :: NUMEL
                INTEGER(KIND=4) :: NUMELB
                INTEGER(KIND=4) :: NUMPBC
                INTEGER(KIND=4) :: NEN
                INTEGER(KIND=4) :: NFACES
                INTEGER(KIND=4) :: NUMFLX
                INTEGER(KIND=4) :: NDOF
                INTEGER(KIND=4) :: IALE
                INTEGER(KIND=4) :: ICOORD
                INTEGER(KIND=4) :: NAVIER
                INTEGER(KIND=4) :: IBLK
                INTEGER(KIND=4) :: IRS
                INTEGER(KIND=4) :: IEXEC
                INTEGER(KIND=4) :: NECHO
                INTEGER(KIND=4) :: ICHEM
                INTEGER(KIND=4) :: IRK
                INTEGER(KIND=4) :: NEDOF
                INTEGER(KIND=4) :: NSHG
                INTEGER(KIND=4) :: NNZ
                INTEGER(KIND=4) :: ISTOP
                INTEGER(KIND=4) :: NFLOW
                INTEGER(KIND=4) :: NNZ_TOT
                INTEGER(KIND=4) :: IDTN
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
              COMMON/INTPT/ QPT,QWT,QPTB,QWTB,NINT,NINTB,NGAUSS,NGAUSSB,&
     &INTP,MAXNINT
                REAL(KIND=8) :: QPT(6,4,125)
                REAL(KIND=8) :: QWT(6,125)
                REAL(KIND=8) :: QPTB(6,4,125)
                REAL(KIND=8) :: QWTB(6,125)
                INTEGER(KIND=4) :: NINT(6)
                INTEGER(KIND=4) :: NINTB(6)
                INTEGER(KIND=4) :: NGAUSS
                INTEGER(KIND=4) :: NGAUSSB
                INTEGER(KIND=4) :: INTP
                INTEGER(KIND=4) :: MAXNINT
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
              REAL(KIND=8) :: Y(NSHG,NDOFL)
              REAL(KIND=8) :: AC(NSHG,NDOFL)
              REAL(KIND=8) :: X(NUMNP,3)
              REAL(KIND=8) :: SHP(NSHL,NGAUSS)
              REAL(KIND=8) :: SHGL(3,NSHL,NGAUSS)
              INTEGER(KIND=4) :: IEN(NPRO,NSHL)
              REAL(KIND=8) :: RES(NSHG)
              REAL(KIND=8) :: QRES(NSHG,3)
              REAL(KIND=8) :: XSEBE(NPRO,NSHL,NSHL)
              REAL(KIND=8) :: XMUDMI(NPRO,NGAUSS)
              REAL(KIND=8) :: CFL(NSHG)
              INTEGER(KIND=4) :: ICFLHITS(NSHG)
              REAL(KIND=8) :: CFLOLD(NSHG)
            END SUBROUTINE ASIGMRSCLR
          END INTERFACE 
        END MODULE ASIGMRSCLR__genmod
