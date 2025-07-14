        !COMPILER-GENERATED INTERFACE MODULE: Fri Jul 11 21:31:34 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE CFCOLLECT__genmod
          INTERFACE 
            SUBROUTINE CFCOLLECT(XX,U1,U2,U3,SCLR,EPSILON_LS_TMP,       &
     &ELEMVOL_LOCAL,RHO,SFORCE)
              COMMON/PROPAR/ NPRO
                INTEGER(KIND=4) :: NPRO
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
              REAL(KIND=8) :: XX(NPRO,3)
              REAL(KIND=8) :: U1(NPRO)
              REAL(KIND=8) :: U2(NPRO)
              REAL(KIND=8) :: U3(NPRO)
              REAL(KIND=8) :: SCLR(NPRO)
              REAL(KIND=8) :: EPSILON_LS_TMP
              REAL(KIND=8) :: ELEMVOL_LOCAL(IBKSIZ)
              REAL(KIND=8) :: RHO(NPRO)
              REAL(KIND=8) :: SFORCE(NPRO,3)
            END SUBROUTINE CFCOLLECT
          END INTERFACE 
        END MODULE CFCOLLECT__genmod
