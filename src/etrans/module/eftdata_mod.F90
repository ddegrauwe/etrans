MODULE EFTDATA_MOD

USE PARKIND1, ONLY : JPRB

IMPLICIT NONE

REAL(KIND=JPRB), ALLOCATABLE, SAVE, TARGET :: ZGTF_PERM (:,:)
END MODULE
