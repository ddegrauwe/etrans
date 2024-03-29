MODULE ELEDIRAD_MOD
CONTAINS
SUBROUTINE ELEDIRAD(KM,KFC,KLED2,PFFT)

!**** *ELEDIRAD* - Direct Legendre transform.

!     Purpose.
!     --------
!        Direct Legendre tranform of state variables.

!**   Interface.
!     ----------
!        CALL ELEDIRAD(...)

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

!     Externals.   MXMAOP - matrix multiply
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Mats Hamrud and Philippe Courtier  *ECMWF*

!     Modifications.
!     --------------
!        Original : 88-01-28
!        Modified : 91-07-01 Philippe Courtier/Mats Hamrud - Rewrite
!                            for uv formulation
!        Modified : 93-03-19 D. Giard - NTMAX instead of NSMAX
!        Modified : 04/06/99 D.Salmond : change order of AIA and SIA
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        D. Degrauwe  (Feb 2012): Alternative extension zone (E')

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE TPM_DIM         ,ONLY : R
!USE TPM_GEOMETRY
!USE TPM_TRANS

USE TPMALD_FFT      ,ONLY : TALD
USE TPMALD_DIM      ,ONLY : RALD
USE ABORT_TRANS_MOD,ONLY: ABORT_TRANS

IMPLICIT NONE

INTEGER(KIND=JPIM), INTENT(IN)  :: KM
INTEGER(KIND=JPIM), INTENT(IN)  :: KFC
INTEGER(KIND=JPIM), INTENT(IN)  :: KLED2

REAL(KIND=JPRB),   INTENT(INOUT)  :: PFFT(:,:)

INTEGER(KIND=JPIM) :: JF, JJ
REAL(KIND=JPRB) :: ZNORM
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

call abort_trans("Error: TL/AD needs FFT992, removed from trans")

IF (LHOOK) CALL DR_HOOK('ELEDIRAD_MOD:ELEDIRAD',0,ZHOOK_HANDLE)
DO JJ=1,1
  DO JF=1,KFC
    PFFT(JJ,JF) = 2.0_JPRB * PFFT(JJ,JF)
  ENDDO
ENDDO

!CALL FFT992(PFFT,TALD%TRIGSE,TALD%NFAXE,1,RALD%NDGLSUR+R%NNOEXTZG,R%NDGL+R%NNOEXTZG,KFC,1)

ZNORM=1.0_JPRB/(2.0_JPRB*REAL(R%NDGL+R%NNOEXTZG,JPRB))
DO JJ=1,R%NDGL+R%NNOEXTZG
  DO JF=1,KFC
    PFFT(JJ,JF) = ZNORM * PFFT(JJ,JF)
  ENDDO
ENDDO
IF (LHOOK) CALL DR_HOOK('ELEDIRAD_MOD:ELEDIRAD',1,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

END SUBROUTINE ELEDIRAD
END MODULE ELEDIRAD_MOD
