MODULE EFOURIER_IN_MOD
CONTAINS
SUBROUTINE EFOURIER_IN(PREEL,KFIELDS)

!**** *FOURIER_IN* - Copy fourier data from buffer to local array

!     Purpose.
!     --------
!        Routine for copying fourier data from buffer to local array

!**   Interface.
!     ----------
!     CALL FOURIER_IN(...)

!     Explicit arguments :  PREEL - local fourier/GP array
!     --------------------  KFIELDS - number of fields
!
!     Externals.  None.
!     ----------

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 2000-04-01

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM
USE PARKIND_ECTRANS, ONLY : JPRBT

USE TPM_DISTR       ,ONLY : D, MYSETW, MYPROC, NPROC, &
             & D_NSTAGTF, D_MSTABF, D_NSTAGT0B, D_NPNTGTB0, D_NPROCM, D_NPTRLS
USE TPM_TRANS       ,ONLY : FOUBUF
USE TPM_GEOMETRY    ,ONLY : G, G_NMEN,G_NMEN_MAX
!

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN) :: KFIELDS
REAL(KIND=JPRBT), INTENT(OUT) :: PREEL(:,:)

INTEGER(KIND=JPIM) :: JM, JF, IGLG, IPROC, IR, II, ISTA
INTEGER(KIND=JPIM) :: IOFF, JGL

!     ------------------------------------------------------------------

!$acc data &
!$acc& copyin(D_NPTRLS,D_NSTAGTF,D_MSTABF,D_NSTAGT0B,D_NPNTGTB0,G_NMEN,G_NMEN_MAX,D_NPROCM) &
!$acc& present(PREEL,FOUBUF)

!!$acc update host(FOUBUF)
!write (77,*) __FILE__, __LINE__
!write (77,*) 'foubuf = '; write(77,'(999F8.2)') foubuf

!$acc parallel loop collapse(2)
do jf=1,KFIELDS
  do jgl=1,size(preel,1)
    PREEL(jgl,jf)=0._JPRBT
  enddo
enddo


!$acc parallel loop collapse(3) private(IGLG,IPROC,ISTA,IOFF)
DO JGL = 1, D%NDGL_FS
   DO JM=0,G_NMEN_MAX      
      DO JF=1,KFIELDS     
         IGLG = D_NPTRLS(MYSETW)+JGL-1
         IF ( JM .LE. G_NMEN(IGLG)) THEN
            IPROC = D_NPROCM(JM)
            ISTA  = (D_NSTAGT0B(D_MSTABF(IPROC))+D_NPNTGTB0(JM,JGL))*2*KFIELDS
            IOFF  = 1+D_NSTAGTF(JGL)
            PREEL(IOFF+2*JM+0, JF) = FOUBUF(ISTA+2*JF-1) 
            PREEL(IOFF+2*JM+1, JF) = FOUBUF(ISTA+2*JF  ) 
         END IF
      ENDDO
   ENDDO
END DO

!!$acc update host(PREEL)
!write (77,*) __FILE__, __LINE__
!write (77,*) 'preel = '; write(77,'(999F8.2)') preel


!$acc end data

END SUBROUTINE EFOURIER_IN
END MODULE EFOURIER_IN_MOD

