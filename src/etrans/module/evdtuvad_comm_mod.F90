MODULE EVDTUVAD_COMM_MOD
CONTAINS
SUBROUTINE EVDTUVAD_COMM(KM,KMLOC,KFIELD,KFLDPTR,PSPMEANU,PSPMEANV)

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE TPM_DIM
USE TPM_FIELDS
USE TPM_DISTR

USE TPMALD_FIELDS
USE TPMALD_GEO
USE TPMALD_DISTR

USE MPL_MODULE
USE ABORT_TRANS_MOD
USE SET2PE_MOD


!**** *EVDTUVAD_COMM* - Compute U,V in  spectral space

!     Purpose.
!     --------
!        In Laplace space communicate the mean winds
!        from vorticity and divergence.

!**   Interface.
!     ----------
!        CALL EVDTUVAD_COMM(...)

!        Explicit arguments :  KM -zonal wavenumber (input-c)
!        --------------------  KFIELD - number of fields (input-c)
!                              KFLDPTR - fields pointers
!                              PEPSNM - REPSNM for wavenumber KM (input-c)
!        Organisation within NLEI1:
!        NLEI1 = NSMAX+4+mod(NSMAX+4+1,2)
!                        overdimensioning
!        1        : n=NSMAX+2
!        2        : n=NSMAX+1
!        3        : n=NSMAX
!        .        :
!        .        :
!        NSMAX+3  : n=0
!        NSMAX+4  : n=-1

!        Implicit arguments :  Eigenvalues of inverse Laplace operator
!        --------------------  from YOMLAP

!     Method.
!     -------

!     Externals.   None.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS
!        Temperton, 1991, MWR 119 p1303

!     Author.
!     -------
!        Mats Hamrud and Philippe Courtier  *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-02-01 From VDTUVAD in IFS CY22R1
!        01-08-27 : R. El Khatib Fix for NPROMATR /= 0
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        01-Dec-2004   A. Deckmyn    Fix mean wind for NPRTRW > 1
!        N. Lopes & R. El Khatib 15-Jun-2012 Scalability enhancement +
!        thread-safety
!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM), INTENT(IN)    :: KM, KFIELD, KMLOC

INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN)  :: KFLDPTR(:)
REAL(KIND=JPRB),    OPTIONAL, INTENT(OUT) :: PSPMEANU(:),PSPMEANV(:)

INTEGER(KIND=JPIM) :: II, IJ, IR, J, JN, IFLD

INTEGER(KIND=JPIM) :: IN
INTEGER(KIND=JPIM) :: ISND, JA, ITAG, ILEN

INTEGER(KIND=JPIM) :: ISENDREQ(NPRTRW)

REAL(KIND=JPRB) :: ZSPU(2*KFIELD)
REAL(KIND=JPRB) :: ZKM
REAL(KIND=JPRB) :: ZIN
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('EVDTUVAD_COMM_MOD:EVDTUVAD_COMM',0,ZHOOK_HANDLE)

IF (NPRTRW > 1 .AND. KFIELD > 0) THEN
  IF (KM == 0) THEN
    IF (PRESENT(KFLDPTR)) THEN
      DO J=1,KFIELD
        IFLD=KFLDPTR(J)
        ZSPU(J)=PSPMEANU(IFLD)
        ZSPU(KFIELD+J)=PSPMEANV(IFLD)
      ENDDO
    ELSE
      DO J=1,KFIELD
        ZSPU(J)=PSPMEANU(J)
        ZSPU(KFIELD+J)=PSPMEANV(J)
      ENDDO
    ENDIF 
    DO JA=1,NPRTRW
      IF (JA /= MYSETW) THEN
        CALL SET2PE(ISND,0,0,JA,MYSETV)
        ISND=NPRCIDS(ISND)
        ITAG=300000+KFIELD*NPROC+ISND
        CALL MPL_SEND(ZSPU(1:2*KFIELD),KDEST=ISND,KTAG=ITAG, &
         & KMP_TYPE=JP_NON_BLOCKING_STANDARD,KREQUEST=ISENDREQ(JA), &
         & CDSTRING='EVDTUVAD_COMM:')
      ENDIF
    ENDDO
  ELSE
    IF (KMLOC == 1) THEN
      IF (D%NPROCM(0) /= MYSETW) THEN
        CALL SET2PE(ISND,0,0,D%NPROCM(0),MYSETV)
        ISND=NPRCIDS(ISND)
        ITAG=300000+KFIELD*NPROC+MYPROC
        CALL MPL_RECV(ZSPU(1:2*KFIELD),KSOURCE=ISND,KTAG=ITAG,KOUNT=ILEN,CDSTRING='EVDTUVAD_COMM:')
        IF (ILEN /= 2*KFIELD) THEN
          CALL ABORT_TRANS('EVDTUVAD_COMM: RECV INVALID RECEIVE MESSAGE LENGTH')
        ENDIF
        IF (PRESENT(KFLDPTR)) THEN
          DO J=1,KFIELD
            IFLD=KFLDPTR(J)
            PSPMEANU(IFLD)=ZSPU(J)
            PSPMEANV(IFLD)=ZSPU(KFIELD+J)
          ENDDO
        ELSE
          DO J=1,KFIELD
            PSPMEANU(J)=ZSPU(J)
            PSPMEANV(J)=ZSPU(KFIELD+J)
          ENDDO
        ENDIF
      ENDIF
    ENDIF
  ENDIF
ENDIF

IF (LHOOK) CALL DR_HOOK('EVDTUVAD_COMM_MOD:EVDTUVAD_COMM',1,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

END SUBROUTINE EVDTUVAD_COMM
END MODULE EVDTUVAD_COMM_MOD
