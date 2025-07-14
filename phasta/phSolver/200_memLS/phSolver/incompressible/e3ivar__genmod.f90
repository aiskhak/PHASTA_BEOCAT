        !COMPILER-GENERATED INTERFACE MODULE: Fri Jul 11 21:31:30 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE E3IVAR__genmod
          INTERFACE 
            SUBROUTINE E3IVAR(YL,ACL,SHPFUN,SHGL,XL,BML,ACI,G1YI,G2YI,  &
     &G3YI,SHG,DXIDX,WDETJ,RHO,PRES,U1,U2,U3,QL,RLUI,SRC,RERRL,RLSL,    &
     &RLSLI,DWL,XWL,YWL,ZWL,ELEMVOL_LOCAL,SCLR,EPSILON_LS_TMP,XARRAY,   &
     &YARRAY,ZARRAY,BUBRADIUS,BUBRADIUS2,COORDTAG)
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
              COMMON/GENPAR/ E3NSD,I3NSD,NSYMDF,NDOFBC,NDIBCB,NDBCB,    &
     &JACTYP,JUMP,IRES,IPREC,IPREV,IBOUND,IDIFF,LHS,ITAU,IPORD,IPRED,   &
     &LSTRES,IEPSTM,DTSFCT,DTSFCTSCLR,TAUCFCT,IBKSIZ,IABC,ISURF,IDFLX,BO&
     &,COALINVSIGMA,PRESAVG,ENTROPYPRESSURE
                REAL(KIND=8) :: E3NSD
                INTEGER(KIND=4) :: I3NSD
                INTEGER(KIND=4) :: NSYMDF
                INTEGER(KIND=4) :: NDOFBC
                INTEGER(KIND=4) :: NDIBCB
                INTEGER(KIND=4) :: NDBCB
                INTEGER(KIND=4) :: JACTYP
                INTEGER(KIND=4) :: JUMP
                INTEGER(KIND=4) :: IRES
                INTEGER(KIND=4) :: IPREC
                INTEGER(KIND=4) :: IPREV
                INTEGER(KIND=4) :: IBOUND
                INTEGER(KIND=4) :: IDIFF
                INTEGER(KIND=4) :: LHS
                INTEGER(KIND=4) :: ITAU
                INTEGER(KIND=4) :: IPORD
                INTEGER(KIND=4) :: IPRED
                INTEGER(KIND=4) :: LSTRES
                INTEGER(KIND=4) :: IEPSTM
                REAL(KIND=8) :: DTSFCT
                REAL(KIND=8) :: DTSFCTSCLR
                REAL(KIND=8) :: TAUCFCT
                INTEGER(KIND=4) :: IBKSIZ
                INTEGER(KIND=4) :: IABC
                INTEGER(KIND=4) :: ISURF
                INTEGER(KIND=4) :: IDFLX
                REAL(KIND=8) :: BO
                REAL(KIND=8) :: COALINVSIGMA
                REAL(KIND=8) :: PRESAVG
                INTEGER(KIND=4) :: ENTROPYPRESSURE
              REAL(KIND=8) :: YL(NPRO,NSHL,NDOF)
              REAL(KIND=8) :: ACL(NPRO,NSHL,NDOF)
              REAL(KIND=8) :: SHPFUN(NPRO,NSHL)
              REAL(KIND=8) :: SHGL(NPRO,3,NSHL)
              REAL(KIND=8) :: XL(NPRO,NENL,3)
              REAL(KIND=8) :: BML(NPRO,NSHL,1)
              REAL(KIND=8) :: ACI(NPRO,3)
              REAL(KIND=8) :: G1YI(NPRO,NDOF)
              REAL(KIND=8) :: G2YI(NPRO,NDOF)
              REAL(KIND=8) :: G3YI(NPRO,NDOF)
              REAL(KIND=8) :: SHG(NPRO,NSHL,3)
              REAL(KIND=8) :: DXIDX(NPRO,3,3)
              REAL(KIND=8) :: WDETJ(NPRO)
              REAL(KIND=8) :: RHO(NPRO)
              REAL(KIND=8) :: PRES(NPRO)
              REAL(KIND=8) :: U1(NPRO)
              REAL(KIND=8) :: U2(NPRO)
              REAL(KIND=8) :: U3(NPRO)
              REAL(KIND=8) :: QL(NPRO,NSHL,IDFLX)
              REAL(KIND=8) :: RLUI(NPRO,3)
              REAL(KIND=8) :: SRC(NPRO,3)
              REAL(KIND=8) :: RERRL(NPRO,NSHL,6)
              REAL(KIND=8) :: RLSL(NPRO,NSHL,6)
              REAL(KIND=8) :: RLSLI(NPRO,6)
              REAL(KIND=8) :: DWL(NPRO,NENL)
              REAL(KIND=8) :: XWL(NPRO,NENL)
              REAL(KIND=8) :: YWL(NPRO,NENL)
              REAL(KIND=8) :: ZWL(NPRO,NENL)
              REAL(KIND=8) :: ELEMVOL_LOCAL(IBKSIZ)
              REAL(KIND=8) :: SCLR(NPRO)
              REAL(KIND=8) :: EPSILON_LS_TMP
              REAL(KIND=8) :: XARRAY(IBKSIZ)
              REAL(KIND=8) :: YARRAY(IBKSIZ)
              REAL(KIND=8) :: ZARRAY(IBKSIZ)
              REAL(KIND=8) :: BUBRADIUS
              REAL(KIND=8) :: BUBRADIUS2
              INTEGER(KIND=4) :: COORDTAG(IBKSIZ)
            END SUBROUTINE E3IVAR
          END INTERFACE 
        END MODULE E3IVAR__genmod
