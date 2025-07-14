        !COMPILER-GENERATED INTERFACE MODULE: Fri Jul 11 21:31:06 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE GETSTRL__genmod
          INTERFACE 
            SUBROUTINE GETSTRL(Y,X,IEN,STRNRM,SHGL,SHP)
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
              REAL(KIND=8) :: Y(NSHG,5)
              REAL(KIND=8) :: X(NUMNP,3)
              INTEGER(KIND=4) :: IEN(NPRO,NSHL)
              REAL(KIND=8) :: STRNRM(NPRO,32)
              REAL(KIND=8) :: SHGL(3,NSHL,32)
              REAL(KIND=8) :: SHP(NSHL,32)
            END SUBROUTINE GETSTRL
          END INTERFACE 
        END MODULE GETSTRL__genmod
