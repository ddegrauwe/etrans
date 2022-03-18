PROGRAM TEST_EZONES
!   purpose  :
!   --------
!    To test SPECTRAL TRANSFORM.

!    method  :
!   ---------
!    Input spectral transform with two horizontal field's after biperiodicization
!    (analitical temperature filed in this test), both having same
!    C+I zone but second having larger E zone and larger truncation number.
!    Do direct and inverse transform. Compute RMSE after inverse transform
!    on C+I zone for both domains.

!   interface  :
!   ---------
!
!   externals :
!   ----------
!   SETUP_TRANS0 - General setup routine for transform package
!   ESETUP_TRANS - Setup transform package for specific resolution
!   ETRANS_INQ   - Extract information from the transform package
!   EINV_TRANS   - Inverse spectral transform (from spectral to grid-point)
!   EDIR_TRANS   - Direct spectral transform (from grid-point to spectral)
!   HORIZ_FIELD  - Calculate test analitical temperature field
!   FPBIPERE     - Interface routine for biperiodicization

!   references :
!    ----------

!    author :
!    -----
!    15-04-2008   Antonio Stanesic
!    ----------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE
INTEGER(KIND=JPIM)              :: NDLON1,NDGL1,NDGUX1,NDLUX1,NGPTOT1,NSMAX1,NMSMAX1,NSPEC21,NRESOL1
INTEGER(KIND=JPIM)              :: NDLON2,NDGL2,NDGUX2,NDLUX2,NGPTOT2,NSMAX2,NMSMAX2,NSPEC22,NRESOL2
REAL(KIND=JPRB),ALLOCATABLE     :: GGPBI1(:,:,:),GGPBIR1(:,:,:),SPECTEMP1(:,:)
REAL(KIND=JPRB),ALLOCATABLE     :: GGPBI2(:,:,:),GGPBIR2(:,:,:),SPECTEMP2(:,:)
INTEGER(KIND=JPIM),ALLOCATABLE  :: NLOEN1(:)
INTEGER(KIND=JPIM),ALLOCATABLE  :: NLOEN2(:)
INTEGER(KIND=JPIM)              :: JX,JY,ISTAE,COUNT
REAL(KIND=JPRB)                 :: SUM2,RMS
REAL(KIND=JPRB),ALLOCATABLE     :: HFIELD(:,:)

#include "setup_trans0.h"
#include "esetup_trans.h"
#include "etrans_inq.h"
#include "einv_trans.h"
#include "edir_trans.h"
#include "horiz_field.h"
#include "fpbipere.h"


!-----------------------------------------------------------
!1. GRID SETUP
!-----------------------------------------------------------

WRITE(*,*) '*****************************************************************'
WRITE(*,*) '*                                                               *'
WRITE(*,*) '*                 TEST OF SPECTRAL TRANSFORMS                   *'
WRITE(*,*) '*                                                               *'
WRITE(*,*) '*****************************************************************'
WRITE(*,*) ' '
WRITE(*,*) '************************** START SETUP **************************'
WRITE(*,*) ' '

!1.1. Setup for grid 1
!------------------------

NDLON1=100
NDGL1=100
NDLUX1=90
NDGUX1=90
ALLOCATE(NLOEN1(NDGL1))
NLOEN1(:)=NDLON1
NSMAX1=49
NMSMAX1=49
WRITE(*,*) '************ GRIDPOINT SETUP 1 **************'
WRITE(UNIT=*,FMT='(4(A,I4,2X))') 'NDLON:',NDLON1,'NDGL:',NDGL1,'NDLUX:',NDLUX1,'NDGUX:',NDGUX1
WRITE(UNIT=*,FMT='(4(A,I4,2X))') 'NSMAX:',NSMAX1,'NMSMAX:',NMSMAX1

!1.1.1 Read and biper input field
!--------------------------------
!allocate two horizontal field's
!one for spectral transformations and one for reference
ALLOCATE(GGPBI1(NDLON1*NDGL1,1,1),GGPBIR1(NDLON1*NDGL1,1,1))

!read input field
WRITE(*,*) 'READ FIELD 1'
ALLOCATE(HFIELD(NDLUX1,NDGUX1))

CALL HORIZ_FIELD(NDLUX1,NDGUX1,HFIELD)

!change layout
ISTAE=0
DO JY=1,NDGUX1
  DO JX=1,NDLUX1
   GGPBI1(ISTAE+JX,1,1)=HFIELD(JX,JY)
  ENDDO
ISTAE=ISTAE+NDLUX1
ENDDO

!biper input field
WRITE(*,*) 'BIPER FIELD 1'

CALL FPBIPERE(NDLUX1,NDGUX1,NDLON1,NDGL1,1,NDLON1*NDGL1,GGPBI1,0,.FALSE.)


!1.2. Setup for grid 2
!------------------------

NDLON2=120
NDGL2=120
NDLUX2=90
NDGUX2=90
ALLOCATE(NLOEN2(NDGL2))
NLOEN2(:)=NDLON2
NSMAX2=59
NMSMAX2=59

WRITE(*,*) '************ GRIDPOINT SETUP 2 **************'
WRITE(UNIT=*,FMT='(4(A,I4,2X))') 'NDLON:',NDLON2,'NDGL:',NDGL2,'NDLUX:',NDLUX2,'NDGUX:',NDGUX2
WRITE(UNIT=*,FMT='(4(A,I4,2X))') 'NSMAX:',NSMAX2,'NMSMAX:',NMSMAX2

!1.2.1 Read and biper input field
!----------------------------
!allocate two horizontal field's
!one for spectral transformations and one for reference
ALLOCATE(GGPBI2(NDLON2*NDGL2,1,1),GGPBIR2(NDLON2*NDGL2,1,1))

!read input field
WRITE(*,*) 'READ FIELD 2'
ALLOCATE(HFIELD(NDLUX2,NDGUX2))

CALL HORIZ_FIELD(NDLUX2,NDGUX2,HFIELD)

!change layout
ISTAE=0
DO JY=1,NDGUX2
  DO JX=1,NDLUX2
   GGPBI2(ISTAE+JX,1,1)=HFIELD(JX,JY)
  ENDDO
ISTAE=ISTAE+NDLUX2
ENDDO

!biper input field
WRITE(*,*) 'BIPER FIELD 2'

CALL FPBIPERE(NDLUX2,NDGUX2,NDLON2,NDGL2,1,NDLON2*NDGL2,GGPBI2,0,.FALSE.)


!-----------------------------------------------------------
!2. SPECTRAL SETUP
!-----------------------------------------------------------

!2.1. Resolution independent part
!---------------------------------


CALL SETUP_TRANS0(KMAX_RESOL=2,LDMPOFF=.TRUE.)


!2.2. Setup for grid 1
!-------------------------------------

WRITE(*,*) '************ SPECTRAL SETUP 1 **************'

CALL ESETUP_TRANS(NMSMAX1,NSMAX1,NDGL1,NDGUX1,NLOEN1,KRESOL=NRESOL1)

CALL ETRANS_INQ(KRESOL=NRESOL1,KSPEC2=NSPEC21,KGPTOT=NGPTOT1)

ALLOCATE(SPECTEMP1(1,NSPEC21))

WRITE(UNIT=*,FMT='(3(A,I7,2X))') 'NSPEC21:',NSPEC21,'NGPTOT1:',NGPTOT1,'NRESOL1:',NRESOL1
WRITE(*,*) 'END SETUP 1'


!2.3. Setup for grid 2
!-------------------------------------

WRITE(*,*) '************ SPECTRAL SETUP 2 **************'

CALL ESETUP_TRANS(NMSMAX2,NSMAX2,NDGL2,NDGUX2,NLOEN2,KRESOL=NRESOL2)

CALL ETRANS_INQ(KRESOL=NRESOL2,KSPEC2=NSPEC22,KGPTOT=NGPTOT2)

ALLOCATE(SPECTEMP2(1,NSPEC22))

WRITE(UNIT=*,FMT='(3(A,I7,2X))') 'NSPEC22:',NSPEC22,'NGPTOT2:',NGPTOT2,'NRESOL2:',NRESOL2
WRITE(UNIT=*,FMT=*) 'END SETUP 2'

WRITE(UNIT=*,FMT=*) '*************************** END SETUP ***************************'
WRITE(*,*) ' '


!-----------------------------------------------------------
!3. DIRECT AND INVERSE TRANSFORM + RMSE COMPUTATIONS
!-----------------------------------------------------------

!3.1. Computations for grid 1
!----------------------------

WRITE(*,*) '**********************************************'
WRITE(*,*) '*    DATA ON C+I AND SMALLER E ZONE          *'
WRITE(*,*) '**********************************************'
WRITE(*,*) ' '


!3.1.1. Go to spectral
!----------------------------

WRITE(*,*) 'DIRECT TRANSFORM OF FIELD 1'

CALL EDIR_TRANS(PSPSCALAR=SPECTEMP1,PGP=GGPBI1,KPROMA=NGPTOT1,KRESOL=NRESOL1)


!3.1.2. Go back to gridpoint
!----------------------------

WRITE(UNIT=*,FMT='(A,I2)') 'INVERSE TRANSFORM OF FIELD 1'

CALL EINV_TRANS(PSPSCALAR=SPECTEMP1,PGP=GGPBIR1,KPROMA=NGPTOT1,KRESOL=NRESOL1)


!3.1.3. Compute RMS on C U I
!----------------------------

SUM2=0
ISTAE=0
COUNT=0
DO JY=1,NDGUX1
 DO JX=1,NDLUX1
  SUM2=SUM2+(GGPBI1(JX+ISTAE,1,1)-GGPBIR1(JX+ISTAE,1,1))**2
  COUNT=COUNT+1
 ENDDO
ISTAE=ISTAE+NDLON1
ENDDO
RMS=SQRT(SUM2/COUNT)


WRITE(UNIT=*,FMT='(A,I2)') 'CALCULATION OF RMS ERROR ON C+I ZONE FOR FIELD 1'
WRITE(UNIT=*,FMT='(A,I5)') 'NUMBER OF POINTS IN CALCULATION 1:',COUNT
WRITE(*,*) 'RMS ERROR 1:',RMS


!3.2. Computations for grid 2
!----------------------------

WRITE(*,*) ' '
WRITE(*,*) '**********************************************'
WRITE(*,*) '*    DATA ON C+I AND BIGGER E ZONE           *'
WRITE(*,*) '**********************************************'
WRITE(*,*) ' '

!3.2.1. Go to spectral
!----------------------------

WRITE(*,*) 'DIRECT TRANSFORM OF FIELD 2'

CALL EDIR_TRANS(PSPSCALAR=SPECTEMP2,PGP=GGPBI2,KPROMA=NGPTOT2,KRESOL=NRESOL2)


!3.2.2. Go back to gridpoint
!----------------------------
WRITE(UNIT=*,FMT='(A,I2)') 'INVERSE TRANSFORM OF FIELD 2'

CALL EINV_TRANS(PSPSCALAR=SPECTEMP2,PGP=GGPBIR2,KPROMA=NGPTOT2,KRESOL=NRESOL2)


!3.2.3. Compute RMS on C U I
!-----------------------------

SUM2=0
ISTAE=0
COUNT=0
DO JY=1,NDGUX2
 DO JX=1,NDLUX2
  SUM2=SUM2+(GGPBI2(JX+ISTAE,1,1)-GGPBIR2(JX+ISTAE,1,1))**2
  COUNT=COUNT+1
 ENDDO
ISTAE=ISTAE+NDLON2
ENDDO

RMS=SQRT(SUM2/COUNT)

WRITE(UNIT=*,FMT='(A,I2)') 'CALCULATION OF RMS ERROR ON C+I ZONE FOR FIELD 2'
WRITE(UNIT=*,FMT='(A,I5)') 'NUMBER OF POINTS IN CALCULATION 2:',COUNT
WRITE(*,*) 'RMS ERROR 2:',RMS


!-----------------------------------------------------------
!4. END SPECTRAL TRANSFORM
!-----------------------------------------------------------

CALL TRANS_END


END PROGRAM TEST_EZONES
