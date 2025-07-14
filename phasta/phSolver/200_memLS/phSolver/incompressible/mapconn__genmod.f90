        !COMPILER-GENERATED INTERFACE MODULE: Fri Jul 11 21:31:30 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE MAPCONN__genmod
          INTERFACE 
            SUBROUTINE MAPCONN(IEN,IEN2,MASK,MAP,NSHL,NPRO,NPRO2,NSHG)
              INTEGER(KIND=4) :: NSHG
              INTEGER(KIND=4) :: NPRO
              INTEGER(KIND=4) :: NSHL
              INTEGER(KIND=4) :: IEN(NPRO,NSHL)
              INTEGER(KIND=4) :: IEN2(NPRO,NSHL)
              INTEGER(KIND=4) :: MASK(NSHG)
              INTEGER(KIND=4) :: MAP(NPRO)
              INTEGER(KIND=4) :: NPRO2
            END SUBROUTINE MAPCONN
          END INTERFACE 
        END MODULE MAPCONN__genmod
