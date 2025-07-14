        !COMPILER-GENERATED INTERFACE MODULE: Fri Jul 11 21:31:33 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE E3IVARSCLR__genmod
          INTERFACE 
            SUBROUTINE E3IVARSCLR(YL,ACL,SHPFUN,SHGL,XL,XMUDMI,SCLR,SDOT&
     &,GRADS,SHG,DXIDX,WDETJ,U1,U2,U3,QL,RLS,SRCR,SRCL,UMOD,DWL,DIFFUS, &
     &SRCRAT,EVL,CFLL)
              COMMON/PROPAR/ NPRO
                INTEGER(KIND=4) :: NPRO
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
              REAL(KIND=8) :: YL(NPRO,NSHL,NDOF)
              REAL(KIND=8) :: ACL(NPRO,NSHL,NDOF)
              REAL(KIND=8) :: SHPFUN(NPRO,NSHL)
              REAL(KIND=8) :: SHGL(NPRO,3,NSHL)
              REAL(KIND=8) :: XL(NPRO,NENL,3)
              REAL(KIND=8) :: XMUDMI(NPRO,NGAUSS)
              REAL(KIND=8) :: SCLR(NPRO)
              REAL(KIND=8) :: SDOT(NPRO)
              REAL(KIND=8) :: GRADS(NPRO,3)
              REAL(KIND=8) :: SHG(NPRO,NSHL,3)
              REAL(KIND=8) :: DXIDX(NPRO,3,3)
              REAL(KIND=8) :: WDETJ(NPRO)
              REAL(KIND=8) :: U1(NPRO)
              REAL(KIND=8) :: U2(NPRO)
              REAL(KIND=8) :: U3(NPRO)
              REAL(KIND=8) :: QL(NPRO,NSHL,3)
              REAL(KIND=8) :: RLS(NPRO)
              REAL(KIND=8) :: SRCR(NPRO)
              REAL(KIND=8) :: SRCL(NPRO)
              REAL(KIND=8) :: UMOD(NPRO,3)
              REAL(KIND=8) :: DWL(NPRO,NSHL)
              REAL(KIND=8) :: DIFFUS(NPRO)
              REAL(KIND=8) :: SRCRAT(NPRO)
              REAL(KIND=8) :: EVL(NPRO,NENL)
              REAL(KIND=8) :: CFLL(NPRO,NSHL)
            END SUBROUTINE E3IVARSCLR
          END INTERFACE 
        END MODULE E3IVARSCLR__genmod
