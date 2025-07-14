        !COMPILER-GENERATED INTERFACE MODULE: Fri Jul 11 21:31:20 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SMAGCOEFF__genmod
          INTERFACE 
            SUBROUTINE SMAGCOEFF(FRES,XDEN,XNUM,IEN,SHP,LFATHE)
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
              COMMON/TURBVAR/ ELES,YLIMIT,RMUTARGET,PZERO,WTAVEI,DTAVEI,&
     &DKE,FWR1,FLUMP,IERRCALC,IHESSIAN,ITWMOD,NGAUSSF,IDIM,NLIST,NINTF
                REAL(KIND=8) :: ELES
                REAL(KIND=8) :: YLIMIT(3,9)
                REAL(KIND=8) :: RMUTARGET
                REAL(KIND=8) :: PZERO
                REAL(KIND=8) :: WTAVEI
                REAL(KIND=8) :: DTAVEI
                REAL(KIND=8) :: DKE
                REAL(KIND=8) :: FWR1
                REAL(KIND=8) :: FLUMP
                INTEGER(KIND=4) :: IERRCALC
                INTEGER(KIND=4) :: IHESSIAN
                INTEGER(KIND=4) :: ITWMOD
                INTEGER(KIND=4) :: NGAUSSF
                INTEGER(KIND=4) :: IDIM
                INTEGER(KIND=4) :: NLIST
                INTEGER(KIND=4) :: NINTF(6)
              REAL(KIND=8) :: FRES(NSHG,22)
              REAL(KIND=8) :: XDEN(IDIM)
              REAL(KIND=8) :: XNUM(IDIM)
              INTEGER(KIND=4) :: IEN(NPRO,NSHL)
              REAL(KIND=8) :: SHP(NSHL,32)
              INTEGER(KIND=4) :: LFATHE(NPRO,32)
            END SUBROUTINE SMAGCOEFF
          END INTERFACE 
        END MODULE SMAGCOEFF__genmod
