        !COMPILER-GENERATED INTERFACE MODULE: Fri Jul 11 21:31:13 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE QPBC__genmod
          INTERFACE 
            SUBROUTINE QPBC(RMASS,QRES,IBC,IPER,ILWORK)
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
              REAL(KIND=8) :: RMASS(NSHG)
              REAL(KIND=8) :: QRES(NSHG,IDFLX)
              INTEGER(KIND=4) :: IBC(NSHG)
              INTEGER(KIND=4) :: IPER(NSHG)
              INTEGER(KIND=4) :: ILWORK
            END SUBROUTINE QPBC
          END INTERFACE 
        END MODULE QPBC__genmod
