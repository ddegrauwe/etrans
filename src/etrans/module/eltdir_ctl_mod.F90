MODULE ELTDIR_CTL_MOD
CONTAINS
SUBROUTINE ELTDIR_CTL(KF_FS,KF_UV,KF_SCALARS, &
 & PSPVOR,PSPDIV,PSPSCALAR, &
 & PSPSC3A,PSPSC3B,PSPSC2, &
 & KFLDPTRUV,KFLDPTRSC,PSPMEANU,PSPMEANV,AUX_PROC)

!**** *ELTDIR_CTL* - Control routine for direct Legendre transform

!     Purpose.
!     --------
!        Direct Legendre transform

!**   Interface.
!     ----------
!     CALL ELTDIR_CTL(...)

!     Explicit arguments :
!     --------------------
!     KF_FS      - number of fields in Fourier space
!     KF_UV      - local number of spectral u-v fields
!     KF_SCALARS - local number of scalar spectral fields
!     PSPVOR(:,:) - spectral vorticity (output)
!     PSPDIV(:,:) - spectral divergence (output)
!     PSPSCALAR(:,:) - spectral scalarvalued fields (output)
!     KFLDPTRUV(:) - field pointer for vorticity and divergence (input)
!     KFLDPTRSC(:) - field pointer for scalarvalued fields (input)
!     PSPMEANU(:),PSPMEANV(:) - mean winds
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE TPM_GEN         ,ONLY : LALLOPERM
USE TPM_TRANS       ,ONLY : FOUBUF, FOUBUF_IN
USE TPM_DISTR       ,ONLY : D, MYPROC

USE ELTDIR_MOD      ,ONLY : ELTDIR
USE TRLTOM_MOD      ,ONLY : TRLTOM, TRLTOM_CUDAAWARE
!

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN) :: KF_FS,KF_UV,KF_SCALARS
REAL(KIND=JPRB)   ,OPTIONAL,INTENT(OUT) :: PSPVOR(:,:)
REAL(KIND=JPRB)   ,OPTIONAL,INTENT(OUT) :: PSPDIV(:,:)
REAL(KIND=JPRB)   ,OPTIONAL,INTENT(OUT) :: PSPSCALAR(:,:)
REAL(KIND=JPRB)   ,OPTIONAL,INTENT(OUT) :: PSPSC3A(:,:,:)
REAL(KIND=JPRB)   ,OPTIONAL,INTENT(OUT) :: PSPSC3B(:,:,:)
REAL(KIND=JPRB)   ,OPTIONAL,INTENT(OUT) :: PSPSC2(:,:)
INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN)  :: KFLDPTRUV(:)
INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN)  :: KFLDPTRSC(:)
REAL(KIND=JPRB)   ,OPTIONAL,INTENT(OUT) :: PSPMEANU(:)
REAL(KIND=JPRB)   ,OPTIONAL,INTENT(OUT) :: PSPMEANV(:)
EXTERNAL AUX_PROC
OPTIONAL AUX_PROC

INTEGER(KIND=JPIM) :: JM,IM,IBLEN,ILED2,INUL
REAL(KIND=JPRB) :: ZDUM
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!     ------------------------------------------------------------------

! Transposition from Fourier space distribution to spectral space distribution


IF (LHOOK) CALL DR_HOOK('ELTDIR_CTL_MOD:ELTDIR_CTL',0,ZHOOK_HANDLE)

IBLEN = D%NLENGT0B*2*KF_FS
#ifdef gnarls
IF (ALLOCATED(FOUBUF)) THEN
  write (77+MYPROC,*)  __FILE__, __LINE__,': foubuf already allocated'; call flush(77+MYPROC)
  IF (MAX(1,IBLEN) > SIZE(FOUBUF)) THEN
!$acc exit data delete (FOUBUF)
    DEALLOCATE(FOUBUF)
	write (77+MYPROC,*)  __FILE__, __LINE__,': reallocating foubuf'; call flush(77+MYPROC)
    ALLOCATE (FOUBUF (MAX (1,IBLEN)))
!$acc enter data create (FOUBUF)
  ENDIF
ELSE
  write (77+MYPROC,*)  __FILE__, __LINE__,': allocating foubuf'; call flush(77+MYPROC)
  ALLOCATE (FOUBUF (MAX (1,IBLEN)))
!$acc enter data create (FOUBUF)
ENDIF
#else
IF (ALLOCATED(FOUBUF)) THEN
  !$acc exit data delete (FOUBUF)
  DEALLOCATE(FOUBUF)
ENDIF

ALLOCATE (FOUBUF (MAX (1,IBLEN)))
!$acc data create(foubuf)

#endif


CALL GSTATS(153,0)
#ifdef USE_CUDA_AWARE_MPI_ELTDIR
write (77+MYPROC,*)  __FILE__, __LINE__,': calling trltom'; call flush(77+MYPROC)
CALL TRLTOM_CUDAAWARE(FOUBUF_IN,FOUBUF,2*KF_FS)
write (77+MYPROC,*)  __FILE__, __LINE__,': trltom call complete'; call flush(77+MYPROC)
#else
CALL TRLTOM(FOUBUF_IN,FOUBUF,2*KF_FS)
#endif

CALL GSTATS(153,1)

! TODO: don't deallocate if LALLOPERM
IF ( ALLOCATED(FOUBUF_IN) ) THEN
write (77+MYPROC,*)  __FILE__, __LINE__,': deallocating foubuf_in'; call flush(77+MYPROC)
!$acc exit data delete (FOUBUF_IN)
  DEALLOCATE (FOUBUF_IN)
ENDIF

! Periodization of auxiliary fields in y direction

IF (PRESENT(AUX_PROC)) THEN
  CALL AUX_PROC(ZDUM,FOUBUF,2*KF_FS,1,IBLEN,0,D%NUMP,.FALSE.,&
   & INUL,D%NPROCL,D%NSTAGT0B,D%NPNTGTB1)
ENDIF

! Direct Legendre transform

ILED2 = 2*KF_FS
CALL GSTATS(1645,0)
IF(KF_FS>0) THEN
  write (77+MYPROC,*)  __FILE__, __LINE__,': calling eltdir'; call flush(77+MYPROC)
    CALL ELTDIR(KF_FS,KF_UV,KF_SCALARS,ILED2, &
     & PSPVOR,PSPDIV,PSPSCALAR,&
     & PSPSC3A,PSPSC3B,PSPSC2 , &
     & KFLDPTRUV,KFLDPTRSC,PSPMEANU,PSPMEANV)
  write (77+MYPROC,*)  __FILE__, __LINE__,': eltdir call complete'; call flush(77+MYPROC)
ENDIF
CALL GSTATS(1645,1)

! TODO: don't deallocate if LALLOPERM
#ifdef gnarls
IF ( ALLOCATED(FOUBUF) ) THEN
write (77+MYPROC,*)  __FILE__, __LINE__,': deallocating foubuf'; call flush(77+MYPROC)
!$acc exit data delete (FOUBUF)
  DEALLOCATE (FOUBUF)
write (77+MYPROC,*)  __FILE__, __LINE__,': deallocation of foubuf done'; call flush(77+MYPROC)
ENDIF
#else

write (77+MYPROC,*)  __FILE__, __LINE__,': launching acc kernels on foubuf'; call flush(77+MYPROC)
!$acc kernels
do iled2=1,size(foubuf,1)
  foubuf(iled2)=1.
enddo
!$acc end kernels

write (77+MYPROC,*)  __FILE__, __LINE__,': leaving acc data region for foubuf'; call flush(77+MYPROC)
!$acc end data

write (77+MYPROC,*)  __FILE__, __LINE__,': deallocating foubuf'; call flush(77+MYPROC)
DEALLOCATE(FOUBUF)
write (77+MYPROC,*)  __FILE__, __LINE__,': deallocation of foubuf done'; call flush(77+MYPROC)

#endif

IF (LHOOK) CALL DR_HOOK('ELTDIR_CTL_MOD:ELTDIR_CTL',1,ZHOOK_HANDLE)

write (77+MYPROC,*)  __FILE__, __LINE__,': leaving eltdir_ctl'; call flush(77+MYPROC)

!     -----------------------------------------------------------------

END SUBROUTINE ELTDIR_CTL
END MODULE ELTDIR_CTL_MOD


