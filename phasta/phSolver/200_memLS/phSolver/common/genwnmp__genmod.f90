        !COMPILER-GENERATED INTERFACE MODULE: Fri Jul 11 21:31:17 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE GENWNMP__genmod
          INTERFACE 
            SUBROUTINE GENWNMP(IBC,BCTMP,X,ILWORK,IPER,WNRMP)
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
              INTEGER(KIND=4) :: IBC(NSHG)
              REAL(KIND=8) :: BCTMP(NSHG,NDOF+7)
              REAL(KIND=8) :: X(NUMNP,3)
              INTEGER(KIND=4) :: ILWORK(NLWORK)
              INTEGER(KIND=4) :: IPER(NSHG)
              REAL(KIND=8) :: WNRMP(NSHG,3)
            END SUBROUTINE GENWNMP
          END INTERFACE 
        END MODULE GENWNMP__genmod
