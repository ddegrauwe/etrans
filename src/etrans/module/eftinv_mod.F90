MODULE EFTINV_MOD
CONTAINS
SUBROUTINE EFTINV(PREEL,KFIELDS)

!**** *FTINV - Inverse Fourier transform

!     Purpose. Routine for Fourier to Grid-point transform
!     --------

!**   Interface.
!     ----------
!        CALL FTINV(..)

!        Explicit arguments :  PREEL   - Fourier/grid-point array
!        --------------------  KFIELDS - number of fields

!     Method.
!     -------

!     Externals.  FFT992 - FFT routine
!     ----------
!

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-03-03
!        G. Radnoti 01-04-24 : 2D model (NLOEN=1)
!        D. Degrauwe  (Feb 2012): Alternative extension zone (E')
!        G. Mozdzynski (Oct 2014): support for FFTW transforms
!        G. Mozdzynski (Jun 2015): Support alternative FFTs to FFTW
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM, JPRB
USE PARKIND_ECTRANS, ONLY : JPRBT
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE TPM_DISTR       ,ONLY : D, MYSETW,  MYPROC, NPROC
USE TPM_GEOMETRY    ,ONLY : G
use tpm_gen, only: nout
USE TPM_FFT         ,ONLY : T !, TB
!USE BLUESTEIN_MOD   ,ONLY : BLUESTEIN_FFT
#ifdef WITH_FFTW
USE TPM_FFTW        ,ONLY : TW, EXEC_FFTW
#endif
USE TPM_FFTR        ,ONLY : CREATE_PLAN_FFTR, EXECUTE_PLAN_FFTR_INPLACE
USE TPM_DIM         ,ONLY : R
!USE CUDA_DEVICE_MOD
USE, INTRINSIC :: ISO_C_BINDING

IMPLICIT NONE

INTEGER (KIND=JPIM), INTENT(IN)    :: KFIELDS
REAL (KIND=JPRBT),   INTENT(INOUT) :: PREEL(:,:)

INTEGER(KIND=JPIM) :: IRLEN, ICLEN
INTEGER(KIND=JPIM) :: IPLAN_C2R
TYPE(C_PTR) :: PLAN_C2R_PTR
integer :: istat
character(len=1024) :: frmt
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('EFTINV_MOD:EFTINV',0,ZHOOK_HANDLE)

IRLEN=R%NDLON+R%NNOEXTZG
ICLEN=D%NLENGTF/D%NDGL_FS

CALL CREATE_PLAN_FFTR (PLAN_C2R_PTR, +1, KN=IRLEN, KLOT=KFIELDS*D%NDGL_FS, &
                    & KISTRIDE=1, KIDIST=ICLEN/2, KOSTRIDE=1, KODIST=ICLEN, LDINPLACE=.TRUE.)					
					
! !$acc update host (PREEL)
! write (0,*) ''
! write (0,*) '=== before transform ==='
! write (0,*) '  shape(PREEL) = ',shape(PREEL)
! write (0,*) '  PREEL = '
! write (frmt,*) '(2X,',ICLEN,'F8.2)'
! write (0,frmt) PREEL
! write (0,*) ''
! call flush(0)

CALL EXECUTE_PLAN_FFTR_INPLACE (PLAN_C2R_PTR, PREEL (1, 1))

! !$acc update host (PREEL)
! write (0,*) ''
! write (0,*) '=== after transform ==='
! write (0,*) '  shape(PREEL) = ',shape(PREEL)
! write (0,*) '  PREEL = '
! write (frmt,*) '(2X,',ICLEN,'F8.2)'
! write (0,frmt) PREEL
! write (0,*) ''
! call flush(0)

IF (LHOOK) CALL DR_HOOK('EFTINV_MOD:EFTINV',1,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

END SUBROUTINE EFTINV
END MODULE EFTINV_MOD
