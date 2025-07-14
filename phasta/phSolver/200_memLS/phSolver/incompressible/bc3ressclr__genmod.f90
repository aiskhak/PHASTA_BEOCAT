        !COMPILER-GENERATED INTERFACE MODULE: Fri Jul 11 21:31:30 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE BC3RESSCLR__genmod
          INTERFACE 
            SUBROUTINE BC3RESSCLR(IBC,RES,IPER,ILWORK)
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
              COMMON/FRONTS/ MAXFRONT,NLWORK,IDIRSTEP,IDIRTRIGGER
                INTEGER(KIND=4) :: MAXFRONT
                INTEGER(KIND=4) :: NLWORK
                INTEGER(KIND=4) :: IDIRSTEP
                INTEGER(KIND=4) :: IDIRTRIGGER
              INTEGER(KIND=4) :: IBC(NSHG)
              REAL(KIND=8) :: RES(NSHG)
              INTEGER(KIND=4) :: IPER(NSHG)
              INTEGER(KIND=4) :: ILWORK(NLWORK)
            END SUBROUTINE BC3RESSCLR
          END INTERFACE 
        END MODULE BC3RESSCLR__genmod
