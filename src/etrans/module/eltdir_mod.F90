MODULE ELTDIR_MOD
CONTAINS
SUBROUTINE ELTDIR(KF_FS,KF_UV,KF_SCALARS,KLED2,&
 & PSPVOR,PSPDIV,PSPSCALAR,PSPSC3A,PSPSC3B,PSPSC2,&
 & KFLDPTRUV,KFLDPTRSC,PSPMEANU,PSPMEANV)

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE TPM_DIM         ,ONLY : R
USE TPM_DISTR       ,ONLY : D
USE TPMALD_DIM      ,ONLY : RALD

USE EPRFI2_MOD      ,ONLY : EPRFI2
USE ELEDIR_MOD      ,ONLY : ELEDIR
USE EUVTVD_MOD
USE EUPDSP_MOD      ,ONLY : EUPDSP
USE EUVTVD_COMM_MOD 
USE EXTPER_MOD      ,ONLY : EXTPER

!
!**** *ELTDIR* - Control of Direct Legendre transform step

!     Purpose.
!     --------
!        Tranform from Fourier space to spectral space, compute
!        vorticity and divergence.

!**   Interface.
!     ----------
!        *CALL* *ELTDIR(...)*

!        Explicit arguments :
!        --------------------  IM     - zonal wavenumber
!                              JM  - local zonal wavenumber

!        Implicit arguments :  None
!        --------------------

!     Method.
!     -------

!     Externals.
!     ----------
!         EPRFI2      - prepares the Fourier work arrays for model variables
!         ELEDIR      - direct Legendre transform
!         EUVTVD      -
!         EUPDSP      - updating of spectral arrays (fields)
!         EUVTVD_COMM -
!         EXTPER      -


!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Mats Hamrud and Philippe Courtier  *ECMWF*

!     Modifications.
!     --------------
!        Original : 87-11-24
!        Modified : 91-07-01 Philippe Courtier/Mats Hamrud - Rewrite
!                            for uv formulation
!        Modified 93-03-19 D. Giard - CDCONF='T' for tendencies
!        Modified 93-11-18 M. Hamrud - use only one Fourier buffer
!        Modified 94-04-06 R. El khatib Full-POS implementation
!        M.Hamrud  : 94-11-01 New conf 'G' - vor,div->vor,div
!                             instead of u,v->vor,div
!        MPP Group : 95-10-01 Support for Distributed Memory version
!        K. YESSAD (AUGUST 1996):
!               - Legendre transforms for transmission coefficients.
!        Modified : 04/06/99 D.Salmond : change order of AIA and SIA
!            01-03-14 G. Radnoti aladin version
!     01-08-27 : R. El Khatib Fix for NPROMATR /= 0
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        F. Vana + NEC 28-Apr-2009 MPI-OpenMP fix
!        D. Degrauwe  (Feb 2012): Alternative extension zone (E')
!        N. Lopes & R. El Khatib 15-Jun-2012 Scalability enhancement
!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)   :: KF_FS,KF_UV,KF_SCALARS,KLED2

REAL(KIND=JPRB)  ,OPTIONAL, INTENT(OUT) :: PSPVOR(:,:)
REAL(KIND=JPRB)  ,OPTIONAL, INTENT(OUT) :: PSPDIV(:,:)
REAL(KIND=JPRB)  ,OPTIONAL, INTENT(OUT) :: PSPSCALAR(:,:)
REAL(KIND=JPRB)   ,OPTIONAL,INTENT(OUT) :: PSPSC2(:,:)
REAL(KIND=JPRB)   ,OPTIONAL,INTENT(OUT) :: PSPSC3A(:,:,:)
REAL(KIND=JPRB)   ,OPTIONAL,INTENT(OUT) :: PSPSC3B(:,:,:)
INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN)  :: KFLDPTRUV(:)
INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN)  :: KFLDPTRSC(:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(OUT) :: PSPMEANU(:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(OUT) :: PSPMEANV(:)

INTEGER(KIND=JPIM) :: IFC, IINDEX(2*KF_FS), JF, JDIM
INTEGER(KIND=JPIM) :: IM
INTEGER(KIND=JPIM) :: JM
INTEGER(KIND=JPIM) :: IUS,IUE,IVS,IVE,IVORS,IVORE,IDIVS,IDIVE 

REAL(KIND=JPRB) :: ZFFT(RALD%NDGLSUR+R%NNOEXTZG,KLED2,D%NUMP)
REAL(KIND=JPRB) :: ZVODI(RALD%NDGLSUR+R%NNOEXTZG,MAX(4*KF_UV,1),D%NUMP)
REAL(KIND=JPRB) :: ZFFT2(KLED2,RALD%NDGLSUR+R%NNOEXTZG)
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('ELTDIR_MOD:ELTDIR',0,ZHOOK_HANDLE)

IF( KF_UV > 0 ) THEN
  IUS = 1
  IUE = 2*KF_UV
  IVS = 2*KF_UV+1
  IVE = 4*KF_UV
  IVORS = 1
  IVORE = 2*KF_UV
  IDIVS = 2*KF_UV+1
  IDIVE = 4*KF_UV
ENDIF
IFC = 2*KF_FS

!$OMP PARALLEL DO PRIVATE(JM,IM,JF,JDIM,IINDEX,ZFFT2)
DO JM=1,D%NUMP
  IM = D%MYMS(JM)

!*     1.    PREPARE WORK ARRAYS.
!            --------------------

  ZFFT(:,:,JM)=0.0_JPRB
  ZVODI(:,:,JM)=0.0_JPRB
  CALL EPRFI2(IM,JM,KF_FS,ZFFT(:,:,JM))


!*     2.    PERIODICIZATION IN Y DIRECTION
!            ------------------------------

  IF(R%NNOEXTZG>0) THEN
    DO JF = 1,IFC
      DO JDIM = 1,R%NDGL
        ZFFT2(JF,JDIM)=ZFFT(JDIM,JF,JM)
      ENDDO
    ENDDO
    IINDEX(1)=0
    CALL EXTPER(ZFFT2(:,:),R%NDGL+R%NNOEXTZG,1,R%NDGL,IFC,1,IINDEX,0)
    DO JF = 1,IFC
      DO JDIM = 1,R%NDGL+R%NNOEXTZG
        ZFFT(JDIM,JF,JM) = ZFFT2(JF,JDIM)
      ENDDO
    ENDDO
  ENDIF


!*     3.    DIRECT LEGENDRE TRANSFORM.
!            --------------------------

  CALL ELEDIR(IM,IFC,KLED2,ZFFT(:,:,JM))


!*     4.    COMPUTE VORTICITY AND DIVERGENCE.
!            ---------------------------------

  IF( KF_UV > 0 ) THEN
    CALL EUVTVD(IM,JM,KF_UV,KFLDPTRUV,ZFFT(:,IUS:IUE,JM),ZFFT(:,IVS:IVE,JM),&
     & ZVODI(:,IVORS:IVORE,JM),ZVODI(:,IDIVS:IDIVE,JM),PSPMEANU,PSPMEANV)
  ENDIF

ENDDO
!$OMP END PARALLEL DO


!*     5.    COMMUNICATION OF MEAN WIND
!            --------------------------

IF( KF_UV > 0 ) THEN

  DO JM=1,D%NUMP
    IM = D%MYMS(JM)

    CALL EUVTVD_COMM(IM,JM,KF_UV,KFLDPTRUV,ZFFT(:,IUS:IUE,JM), &
     & ZFFT(:,IVS:IVE,JM),ZVODI(:,IVORS:IVORE,JM),ZVODI(:,IDIVS:IDIVE,JM), &
     & PSPMEANU,PSPMEANV)

  ENDDO

ENDIF


!*       6.    UPDATE SPECTRAL ARRAYS.
!              -----------------------

!$OMP PARALLEL DO PRIVATE(JM,IM)
DO JM=1,D%NUMP
  IM = D%MYMS(JM)

  CALL EUPDSP(IM,KF_UV,KF_SCALARS,ZFFT(:,:,JM),ZVODI(:,:,JM), &
   & PSPVOR,PSPDIV,PSPSCALAR,PSPSC3A,PSPSC3B,PSPSC2,KFLDPTRUV,KFLDPTRSC)

ENDDO
!$OMP END PARALLEL DO

IF (LHOOK) CALL DR_HOOK('ELTDIR_MOD:ELTDIR',1,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

END SUBROUTINE ELTDIR
END MODULE ELTDIR_MOD
