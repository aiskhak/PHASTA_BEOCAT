        !COMPILER-GENERATED INTERFACE MODULE: Fri Jul 11 21:31:15 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE GETSHPB__genmod
          INTERFACE 
            SUBROUTINE GETSHPB(SHP,SHGL,SGN,SHAPE,SHDRV)
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
              COMMON/PROPAR/ NPRO
                INTEGER(KIND=4) :: NPRO
              REAL(KIND=8) :: SHP(NSHL,NGAUSSB)
              REAL(KIND=8) :: SHGL(3,NSHL,NGAUSSB)
              REAL(KIND=8) :: SGN(NPRO,NSHL)
              REAL(KIND=8) :: SHAPE(NPRO,NSHL)
              REAL(KIND=8) :: SHDRV(NPRO,3,NSHL)
            END SUBROUTINE GETSHPB
          END INTERFACE 
        END MODULE GETSHPB__genmod
