   
! Copyright ©2011 UCSD
! Created by Mahdi Esmaily Moghadam
! contact memt63@gmail.com for reporting the bugs.

      USE ISO_FORTRAN_ENV
      
      IMPLICIT NONE
      
      INCLUDE "mpif.h"
      INCLUDE "STRUCT.h"

!     Communication parameters
      INTEGER, PARAMETER :: stdout = OUTPUT_UNIT
      INTEGER(KIND=MPI_ADDRESS_KIND), PARAMETER :: 
     2   mplog  = MPI_LOGICAL, mpint  = MPI_INTEGER,
     3   mpreal = MPI_DOUBLE_PRECISION, mpchar = MPI_CHARACTER
