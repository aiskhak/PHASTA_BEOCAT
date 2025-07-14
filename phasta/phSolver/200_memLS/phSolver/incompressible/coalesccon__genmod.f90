        !COMPILER-GENERATED INTERFACE MODULE: Fri Jul 11 21:31:28 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COALESCCON__genmod
          INTERFACE 
            SUBROUTINE COALESCCON(XARRAY,YARRAY,ZARRAY,COORDTAG,        &
     &BUBRADIUS,BUBRADIUS2,AVGXCOORDF,AVGYCOORDF,AVGZCOORDF)
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
              REAL(KIND=8) :: XARRAY(IBKSIZ)
              REAL(KIND=8) :: YARRAY(IBKSIZ)
              REAL(KIND=8) :: ZARRAY(IBKSIZ)
              INTEGER(KIND=4) :: COORDTAG(IBKSIZ)
              REAL(KIND=8) :: BUBRADIUS
              REAL(KIND=8) :: BUBRADIUS2
              REAL(KIND=8) :: AVGXCOORDF(100)
              REAL(KIND=8) :: AVGYCOORDF(100)
              REAL(KIND=8) :: AVGZCOORDF(100)
            END SUBROUTINE COALESCCON
          END INTERFACE 
        END MODULE COALESCCON__genmod
