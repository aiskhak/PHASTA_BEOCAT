        !COMPILER-GENERATED INTERFACE MODULE: Fri Jul 11 21:31:06 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE GENSAV__genmod
          INTERFACE 
            SUBROUTINE GENSAV(IENTMP,MATTMP,IEN,MATER)
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
              INTEGER(KIND=4) :: IENTMP(NPRO,NSHL)
              INTEGER(KIND=4) :: MATTMP(NPRO)
              INTEGER(KIND=4) :: IEN(NPRO,NSHL)
              INTEGER(KIND=4) :: MATER(NPRO)
            END SUBROUTINE GENSAV
          END INTERFACE 
        END MODULE GENSAV__genmod
