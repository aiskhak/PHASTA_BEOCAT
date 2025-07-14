        !COMPILER-GENERATED INTERFACE MODULE: Fri Jul 11 21:31:28 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE CFRESTAR__genmod
          INTERFACE 
            SUBROUTINE CFRESTAR(YCF_OLD_LOG,NUMOLDYHISTIND,IXCF,IYCF,   &
     &IZCF,SUMXCF_O,SUMYCF_O,SUMZCF_O,OLDCF_DT,OLDCF_DTLSET)
              LOGICAL(KIND=4) :: YCF_OLD_LOG
              INTEGER(KIND=4) :: NUMOLDYHISTIND
              INTEGER(KIND=4) :: IXCF
              INTEGER(KIND=4) :: IYCF
              INTEGER(KIND=4) :: IZCF
              REAL(KIND=8) :: SUMXCF_O
              REAL(KIND=8) :: SUMYCF_O
              REAL(KIND=8) :: SUMZCF_O
              REAL(KIND=8) :: OLDCF_DT
              REAL(KIND=8) :: OLDCF_DTLSET
            END SUBROUTINE CFRESTAR
          END INTERFACE 
        END MODULE CFRESTAR__genmod
