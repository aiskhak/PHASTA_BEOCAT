        !COMPILER-GENERATED INTERFACE MODULE: Fri Jul 11 21:31:34 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE FILTRAT__genmod
          INTERFACE 
            SUBROUTINE FILTRAT(Y,SHGL,SHP,IPER,ILWORK,NSONS,IFATH,X,    &
     &CDELSQ1,FWR4,FWR3)
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
              COMMON/FRONTS/ MAXFRONT,NLWORK,IDIRSTEP,IDIRTRIGGER
                INTEGER(KIND=4) :: MAXFRONT
                INTEGER(KIND=4) :: NLWORK
                INTEGER(KIND=4) :: IDIRSTEP
                INTEGER(KIND=4) :: IDIRTRIGGER
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
              REAL(KIND=8) :: Y(NSHG,5)
              REAL(KIND=8) :: SHGL(6,3,32,125)
              REAL(KIND=8) :: SHP(6,32,125)
              INTEGER(KIND=4) :: IPER(NSHG)
              INTEGER(KIND=4) :: ILWORK(NLWORK)
              INTEGER(KIND=4) :: NSONS(NSHG)
              INTEGER(KIND=4) :: IFATH(NSHG)
              REAL(KIND=8) :: X(NUMNP,3)
              REAL(KIND=8) :: CDELSQ1(NFATH)
              REAL(KIND=8) :: FWR4(NSHG)
              REAL(KIND=8) :: FWR3(NSHG)
            END SUBROUTINE FILTRAT
          END INTERFACE 
        END MODULE FILTRAT__genmod
