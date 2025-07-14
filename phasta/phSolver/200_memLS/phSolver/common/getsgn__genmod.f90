        !COMPILER-GENERATED INTERFACE MODULE: Fri Jul 11 21:31:11 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE GETSGN__genmod
          INTERFACE 
            SUBROUTINE GETSGN(IEN,SGN)
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
              INTEGER(KIND=4) :: IEN(NPRO,NSHL)
              REAL(KIND=8) :: SGN(NPRO,NSHL)
            END SUBROUTINE GETSGN
          END INTERFACE 
        END MODULE GETSGN__genmod
