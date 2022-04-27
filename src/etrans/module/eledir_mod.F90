MODULE ELEDIR_MOD
CONTAINS
SUBROUTINE ELEDIR(KFC,KLED2,PFFT)

!**** *ELEDIR* - Direct meridional transform.

!     Purpose.
!     --------
!        Direct meridional tranform of state variables.

!**   Interface.
!     ----------
!        CALL ELEDIR(...)

!        Explicit arguments :  KM - zonal wavenumber
!        --------------------  KFC - number of field to transform
!                              PAIA - antisymmetric part of Fourier
!                              fields for zonal wavenumber KM
!                              PSIA - symmetric part of Fourier
!                              fields for zonal wavenumber KM
!                              POA1 -  spectral
!                              fields for zonal wavenumber KM
!                              PLEPO - Legendre polonomials

!        Implicit arguments :  None.
!        --------------------

!     Method.
!     -------

!     Reference.
!     ----------

!     Author.
!     -------

!     Modifications.
!     --------------
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM, JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE TPM_DISTR       ,ONLY : D, D_NUMP
USE TPM_DIM         ,ONLY : R
USE TPMALD_FFT      ,ONLY : TALD
USE TPMALD_DIM      ,ONLY : RALD
USE ABORT_TRANS_MOD ,ONLY : ABORT_TRANS
USE TPM_FFTR        ,ONLY : CREATE_PLAN_FFTR, EXECUTE_PLAN_FFTR_INPLACE
!USE CUDA_DEVICE_MOD
!
USE, INTRINSIC :: ISO_C_BINDING

IMPLICIT NONE

INTEGER(KIND=JPIM), INTENT(IN)  :: KFC,KLED2
REAL(KIND=JPRB) ,   INTENT(INOUT)  :: PFFT(:,:,:)

INTEGER(KIND=JPIM) :: IRLEN, ICLEN
INTEGER(KIND=JPIM) :: IPLAN_R2C
TYPE(C_PTR) :: PLAN_R2C_PTR
INTEGER(KIND=JPIM) :: JM, JF, JJ
REAL (KIND=JPRB)   :: ZSCAL

integer :: istat
character(len=1024) :: frmt

!     ------------------------------------------------------------------

!*       1.       PERFORM FOURIER TRANFORM.
!                 --------------------------

IRLEN=R%NDGL+R%NNOEXTZG
ICLEN=RALD%NDGLSUR+R%NNOEXTZG

CALL CREATE_PLAN_FFTR (PLAN_R2C_PTR, -1, KN=IRLEN, KLOT=UBOUND (PFFT,2)*UBOUND (PFFT, 3), &
                    & KISTRIDE=1, KIDIST=ICLEN, KOSTRIDE=1, KODIST=ICLEN/2, LDINPLACE=.TRUE.)

! !$acc update host (PFFT)
! write (0,*) ''
! write (0,*) '=== before transform ==='
! write (0,*) '  shape(PFFT) = ',shape(PFFT)
! write (0,*) '  PFFT = '
! write (frmt,*) '(2X,',ICLEN,'F8.2)'
! write (0,frmt) PFFT
! write (0,*) ''
! call flush(0)

CALL EXECUTE_PLAN_FFTR_INPLACE(PLAN_R2C_PTR, PFFT (1, 1, 1))

! !$acc update host (PFFT)
! write (0,*) ''
! write (0,*) '=== after transform ==='
! write (0,*) '  shape(PFFT) = ',shape(PFFT)
! write (0,*) '  PFFT = '
! write (frmt,*) '(2X,',ICLEN,'F8.2)'
! write (0,frmt) PFFT
! write (0,*) ''
! call flush(0)

ZSCAL = 1._JPRB / REAL (IRLEN, JPRB)

!$acc parallel loop collapse (3) copyin (D_NUMP, KFC, ICLEN, ZSCAL) present (PFFT)
DO JF = 1, KFC
  DO JM = 1, D_NUMP
    DO JJ = 1, ICLEN
      PFFT (JJ, JM, JF) = PFFT (JJ, JM, JF) * ZSCAL
    ENDDO
  ENDDO
ENDDO
!$acc end parallel loop


!     ------------------------------------------------------------------

END SUBROUTINE ELEDIR
END MODULE ELEDIR_MOD
