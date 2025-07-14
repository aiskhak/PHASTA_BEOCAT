        !COMPILER-GENERATED INTERFACE MODULE: Fri Jul 11 21:31:18 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE GETNU__genmod
          INTERFACE 
            SUBROUTINE GETNU(IEN,STRL,XMUDMI,CDELSQ,LFATHE)
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
              INTEGER(KIND=4) :: IEN(NPRO,NSHL)
              REAL(KIND=8) :: STRL(NPRO,NGAUSS)
              REAL(KIND=8) :: XMUDMI(NPRO,NGAUSS)
              REAL(KIND=8) :: CDELSQ(NLIST)
              INTEGER(KIND=4) :: LFATHE(NPRO,NGAUSS)
            END SUBROUTINE GETNU
          END INTERFACE 
        END MODULE GETNU__genmod
