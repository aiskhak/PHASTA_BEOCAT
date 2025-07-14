        !COMPILER-GENERATED INTERFACE MODULE: Fri Jul 11 21:31:34 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE CFPOSTPROCS__genmod
          INTERFACE 
            SUBROUTINE CFPOSTPROCS(ISTP,NUMOLDYHISTIND,IXCF,IYCF,IZCF,  &
     &SUMXCF_O,SUMYCF_O,SUMZCF_O,SUMYCF_OLD,YCF_OLD_LOG)
              INTEGER(KIND=4) :: ISTP
              INTEGER(KIND=4) :: NUMOLDYHISTIND
              INTEGER(KIND=4) :: IXCF
              INTEGER(KIND=4) :: IYCF
              INTEGER(KIND=4) :: IZCF
              REAL(KIND=8) :: SUMXCF_O
              REAL(KIND=8) :: SUMYCF_O
              REAL(KIND=8) :: SUMZCF_O
              REAL(KIND=8) :: SUMYCF_OLD
              LOGICAL(KIND=4) :: YCF_OLD_LOG
            END SUBROUTINE CFPOSTPROCS
          END INTERFACE 
        END MODULE CFPOSTPROCS__genmod
