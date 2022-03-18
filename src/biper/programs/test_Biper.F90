PROGRAM TEST_BIPER
!   purpose  :
!   --------
!    To test biperiodicization.

!    method  :
!   ---------
!    Reads horizontal input field on C U I U E, changes layout for use in ETIBIHI
!    or for use in FPBIPER (takes into account that for FPIBIPER is possible to have
!    input field on C U I U E or C U I. Makes biperiodicization. 

!   interface  :
!   ---------
!
!   externals :
!   ----------
!   ETIBIHI - Doubly-periodicisation

!   
!   references :
!    ----------

!    author :
!    -----
!    23-May-2008   Antonio Stanesic
!    ----------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

INTEGER(KIND=JPIM)               :: NNUBI       ! number of levels to biperiodicise
INTEGER(KIND=JPIM)               :: NSTART      ! first dimension in x direction of g-p array
INTEGER(KIND=JPIM)               :: NDLSM       ! second dimension in x direction of g-p array

INTEGER(KIND=JPIM)               :: NDLON       ! upper bound for the x (or longitude) dimension
                                                ! of the gridpoint array on C U I U E
INTEGER(KIND=JPIM)               :: NDGL        ! upper bound for the y (or latitude) dimension
                                                ! of the gridpoint array on C U I U E
INTEGER(KIND=JPIM)               :: NDLUX       ! upper bound for the x (or longitude) dimension
                                                ! of  C U I.
INTEGER(KIND=JPIM)               :: NDGUX       ! upper bound for the y (or latitude) dimension
                                                ! of  C U I.
REAL(KIND=JPRB),ALLOCATABLE      :: GPBI(:,:,:) ! random horizontal field for test ETIBIHI
REAL(KIND=JPRB),ALLOCATABLE      :: GGPBI(:,:)  ! random horizontal field for test FPBIPER
REAL(KIND=JPRB),ALLOCATABLE      :: HFIELD(:,:) ! input horizontal field 
LOGICAL                          :: LBIPX       ! .TRUE. biperiodicisation in x
LOGICAL                          :: LBIPY       ! .TRUE. biperiodicisation in y
LOGICAL                          :: LNZON       ! .TRUE. if input grid on C U I U E (.FALSE. if C U I)
INTEGER(KIND=JPIM)               :: JX,JY,JLEV,ISTAE,ND1,IENDX,IENDY
INTEGER(KIND=JPIM)               :: IADD        ! 1 if test of spline

#include "etibihie.h"
#include "fpbipere.h"
#include "horiz_field.h"

! ------------------------------------------------------------------

!* 1. Initialization.

!* 1.1 Grindpoint common
! -----------------------------

NNUBI=2
NSTART=1
NDLSM=288
NDLON=288
NDGL=288
NDLUX=266  
NDGUX=266 
LBIPX=.TRUE. 
!LBIPX=.FALSE. 
LBIPY=.TRUE. 
!LBIPY=.FALSE. 
LNZON=.FALSE.
!LNZON=.TRUE.
IADD=1 ! test spline




!* 1.1 Read input field.
! -----------------------------
!Field is read for C+I+E
ALLOCATE(HFIELD(NDLON,NDGL))

CALL HORIZ_FIELD(NDLON,NDGL,HFIELD)

!* 1.2 Field for ETIBIHI.
! ----------------------------- 

ALLOCATE(GPBI(NDLON+IADD,NNUBI,NDGL+IADD)) 

!change layout
DO JLEV=1,NNUBI 
  DO JY=1,NDGL
    DO JX=1,NDLON
      GPBI(JX,JLEV,JY)=HFIELD(JX,JY)
    ENDDO
  ENDDO
ENDDO


!* 1.3 Field for FPBIPER.
! -----------------------------

!It's same field horizontal fild in different layout
!depending on assumption that we have field on C U I U E or 
!C U I. It is done for testing purposes because FBIPER takes
!both kinds of inputs.
ND1=(NDGL+IADD)*(NDLON+IADD)
IF(LNZON) THEN
  IENDX=NDLON
  IENDY=NDGL
ELSE
  IENDX=NDLUX
  IENDY=NDGUX
ENDIF


ALLOCATE(GGPBI(ND1,NNUBI)) !first dimension of GGPBI must always be NDGL*NDLON,
!because output from FPBIPER is always on C U I U E
!change of layout
DO JLEV=1,NNUBI
  ISTAE=0
  DO JY=1,IENDY
    DO JX=1,IENDX
      GGPBI(ISTAE+JX,JLEV)=HFIELD(JX,JY)
    ENDDO
    ISTAE=ISTAE+IENDX
  ENDDO
ENDDO


! ------------------------------------------------------------------

!* 2. ETIBIHIE test

! -----------------------------

WRITE(*,*) " "
WRITE(*,*) "Test of external program ETIBIHIE ... "

CALL ETIBIHIE(NDLON,NDGL,NNUBI,NDLUX,NDGUX,NSTART,NDLSM,GPBI,LBIPX,LBIPY,IADD) 

OPEN(UNIT=10,FILE="etibihi_test.txt",STATUS='REPLACE',POSITION='APPEND')

DO JLEV=1,NNUBI
WRITE(10,*) "******* Left(NDLUN) - right(NDLON+1) border  *******"
 DO JY=1,NDGL
   WRITE(10,FMT='(F12.4)') GPBI(1,JLEV,JY)-GPBI(NDLON+IADD,JLEV,JY)  
 ENDDO

 WRITE(10,*) "******* Bottom(NDGUN) - top(NDGL+1) border  *******"
 DO JX=1,NDLON
   WRITE(10,FMT='(F12.4)') GPBI(JX,JLEV,1)-GPBI(JX,JLEV,NDGL+IADD)  
 ENDDO
ENDDO
WRITE(*,*) "Output from ETIBIHIE test written in file etibihi_test.txt."

! ------------------------------------------------------------------

!* 3. FBIPERE test

! -----------------------------

WRITE(*,*) " "
WRITE(*,*) "Test of external program FPBIPERE ... "

CALL FPBIPERE(NDLUX,NDGUX,NDLON,NDGL,NNUBI,ND1,GGPBI,IADD,LNZON)

OPEN(UNIT=11,FILE="fpbiper_test.txt",STATUS='REPLACE',POSITION='APPEND')

DO JLEV=1,NNUBI
WRITE(11,*) "******* Left(NDLUN) - right(NDLON+1) border  *******"
ISTAE=1
 DO JY=1,NDGL
  WRITE(11,FMT='(F12.4)') GGPBI(ISTAE,JLEV)-GGPBI(ISTAE+NDLON,JLEV)
  ISTAE=ISTAE+NDLON+IADD
 ENDDO

 WRITE(11,*) "******* Bottom(NDGUN) - top(NDGL+1) border  *******"
 DO JX=1,NDLON
  WRITE(11,FMT='(F12.4)') GGPBI(JX,JLEV)-GGPBI(ND1-NDLON-IADD+JX,JLEV)
 ENDDO
ENDDO
WRITE(*,*) "Output from FPBIPER written in file fpbiper_test.txt."

END PROGRAM TEST_BIPER
