MODULE TPMALD_DISTR

! Module for distributed memory environment.

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

TYPE ALDDISTR_TYPE

INTEGER(KIND=JPIM) ,POINTER :: NESM0(:)  ! Address in a spectral array of (m, n=m)
INTEGER(KIND=JPIM) ,POINTER :: NCPL2M(:) ! Number of complex Laplace coefficient for m given
INTEGER(KIND=JPIM) ,POINTER :: NPME(:)   ! Address for the Laplace operator and its inverse

END TYPE ALDDISTR_TYPE

TYPE(ALDDISTR_TYPE),ALLOCATABLE,TARGET :: ALDDISTR_RESOL(:)
TYPE(ALDDISTR_TYPE),POINTER     :: DALD

INTEGER(KIND=JPIM), ALLOCATABLE :: DALD_NESM0  (:)
INTEGER(KIND=JPIM), ALLOCATABLE :: DALD_NCPL2M (:)
INTEGER(KIND=JPIM), ALLOCATABLE :: DALD_NPME   (:)

END MODULE TPMALD_DISTR

