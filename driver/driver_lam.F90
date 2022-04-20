PROGRAM DRIVER_LAM

!
! Spectral transform test
!
! This test performs spectral to real and real to spectral transforms repeated in
! timed loop.
!
! 
! Author : George Mozdzynski
! Modifications: Daan Degrauwe (adapted global to LAM)

USE PARKIND1, ONLY : JPIM, JPRB, JPRD

USE MPL_MODULE
!USE GRIB_API
USE YOMGSTATS

USE PARKIND_ECTRANS, ONLY : JPRBT

use cudafor

IMPLICIT NONE

INTEGER(KIND=JPIM) :: RETURN_CODE

INTEGER(KIND=JPIM),PARAMETER :: JPMLAT=8000

INTEGER(KIND=JPIM) :: ISTACK, GETSTACKUSAGE
REAL(KIND=JPRB)    :: ZERR,ZMAXERR

INTEGER(KIND=JPIM) :: NERR,NULNAM,NPROC,NLIN,NSMAX,NQ
INTEGER(KIND=JPIM) :: NOUT,NSPEC2,NGPTOT,NGPTOTG,IFLD
INTEGER(KIND=JPIM) :: NPROMA, NBLK

INTEGER(KIND=JPIM) :: ITAG,NSPEC2G,I
INTEGER(KIND=JPIM) ,ALLOCATABLE :: NLOEN(:),ITO(:),NPRCIDS(:)
INTEGER(KIND=JPIM) :: MYPROC,JJ
INTEGER   :: JSTEP
REAL(KIND=JPRB)    :: ZTINIT,ZTLOOP,ZTSTEP,TIMEF
REAL(KIND=JPRB),ALLOCATABLE :: ZSPEC(:,:),ZNORM(:),ZNORM1(:)
REAL(KIND=JPRB),ALLOCATABLE :: ZG(:,:), Z(:,:,:)
REAL(KIND=JPRD) :: ZAVEAVE(0:JPMAXSTAT)

LOGICAL :: LUSERPNM, LKEEPRPNM, LUSEFLT
LOGICAL :: LLTRACE_STATS, LLSTATS_OMP, LLSTATS_COMMS
LOGICAL :: LNORMS=.FALSE.
INTEGER(KIND=JPIM) :: IITER, NITERS
REAL(KIND=JPRB) :: ZMAXERR_CHECK=0.0_JPRB

INTEGER(KIND=JPIM) :: NLON, NLAT, NMSMAX, NDGUX, NLATE
REAL(KIND=JPRB) :: ZDX, ZDY
INTEGER(KIND=JPIM) :: NPRGPNS, NPRGPEW, NPRTRV, NPRTRW, MYSETV
INTEGER(KIND=JPIM) :: JST, JNUM, JPRV
LOGICAL :: LSPLIT
CHARACTER(LEN=256) :: CDATAFILE
INTEGER(KIND=JPIM) :: IFILE, NFLD, NGRIB, JFLD, NFLDSP
INTEGER(KIND=JPIM) :: IGRIB
INTEGER(KIND=JPIM) :: JLON, JLAT, KX, KY
REAL(KIND=JPRB), PARAMETER :: RPI=ACOS(-1._JPRB)
INTEGER(KIND=JPIM), ALLOCATABLE :: IVSET(:)

! Define namelists
NAMELIST/NAMTRANS/ LBARRIER_STATS, LBARRIER_STATS2, LDETAILED_STATS, &
                 & LUSERPNM, LKEEPRPNM, LUSEFLT, NQ, NLIN, &
                 & LNORMS, NITERS, ZMAXERR_CHECK
NAMELIST/NAMBIFFT/  NLATE, NFLD, NLON, NLAT, NPROMA
NAMELIST/NAMCT0/ NPRGPNS, NPRGPEW, NPRTRW, NPRTRV


REAL(KIND=JPRBT), ALLOCATABLE :: ZGTF(:,:)

!     ------------------------------------------------------------------

#include "setup_trans0.h"
#include "esetup_trans.h"
#include "einv_trans.h"
#include "edir_trans.h"
#include "edist_grid.h"
#include "etrans_inq.h"
#include "especnorm.h"
!#include "abor1.intfb.h"

! Initializations
NERR = 0
NULNAM = 4
NOUT = 6


! Message passing setup
CALL MPL_INIT()


! Initialize statistics
CALL GSTATS(0,0)
ZTINIT=TIMEF()

! Get number of procs and current proc
NPROC=MPL_NPROC()
MYPROC = MPL_MYRANK()
IF(MYPROC == 1) THEN
  WRITE(NERR,'(" NPROC=",I8)') NPROC
ENDIF

! Array with proc id's
ALLOCATE(NPRCIDS(NPROC))
DO JJ=1,NPROC
  NPRCIDS(JJ) = JJ
ENDDO

! Default statistics settings
LSTATS=.TRUE.
LDETAILED_STATS=.FALSE.
LSTATS_OMP=.FALSE.
LSTATS_COMMS=.FALSE.
LSTATS_MPL=.FALSE.
LBARRIER_STATS=.FALSE.
LBARRIER_STATS2=.FALSE.
LUSERPNM=.FALSE.
LKEEPRPNM=.FALSE.
LUSEFLT=.TRUE.
LNORMS=.FALSE.
NITERS=10

! Default BIFFT settings
!CDATAFILE='data_bifft_200x180.grb'
NLATE=0                                   ! matters for gridpoint distribution:
                                          ! if NLATE=0, all procs have +/- same number of gridpoints.
                                          ! if NLATE>0, procs treating the extension-zone have more gridpoints
                                          !             (useful because physics calculations are not necessary in E-zone)
NFLD=1
NLON=48
NLAT=40
NPROMA=16

! Default distribution properties
NPRGPNS=NPROC
NPRGPEW=1
NPRTRW=NPROC
NPRTRV=1
LSPLIT=.FALSE.

! Read settings from namelist
REWIND(NULNAM)
READ(NULNAM,NAMTRANS)
REWIND(NULNAM)
READ(NULNAM,NAMBIFFT)
REWIND(NULNAM)
READ(NULNAM,NAMCT0)

WRITE (*,*) 'JPRB = ',JPRB
WRITE (*,*) 'NLON = ',NLON,'; NLAT = ',NLAT


! Setup statistics
IF(MYPROC == 1) THEN
  WRITE(NOUT,'(" LDETAILED_STATS=",L1," LBARRIER_STATS=",L1," LBARRIER_STATS2=",L1)')&
   & LDETAILED_STATS,LBARRIER_STATS,LBARRIER_STATS2
ENDIF
IF(LDETAILED_STATS)THEN
  LSTATS_OMP=.TRUE.
  LSTATS_COMMS=.TRUE.
  LSTATS_MPL=.TRUE.
ENDIF
CALL GSTATS_SETUP(NPROC,MYPROC,NPRCIDS,&
 & LSTATS,LSTATSCPU,LSYNCSTATS,LDETAILED_STATS,LBARRIER_STATS,LBARRIER_STATS2,&
 & LLSTATS_OMP,LLSTATS_COMMS,LSTATS_MEM,NSTATS_MEM,LSTATS_ALLOC,&
 & LLTRACE_STATS,NTRACE_STATS,NPRNT_STATS,LXML_STATS,LSTATS_MPL)
CALL GSTATS_PSUT
!CALL GSTATS_LABEL_IFS

! Spectral truncation: linear -- maybe make this a namelist setting?
NMSMAX=(NLON-1)/2
NSMAX=(NLAT-1)/2

! Number of gridpoints w/o extension zone
NDGUX=NLAT-NLATE

! Setup distribution
ALLOCATE(NLOEN(NLAT))
NLOEN(:) = NLON
!write (*,*) __FILE__, __LINE__; call flush(6)
CALL SETUP_TRANS0(KOUT=NOUT,KERR=NERR,KPRINTLEV=0,KMAX_RESOL=1,&
 & KPRGPNS=NPRGPNS,KPRGPEW=NPRGPEW,KPRTRW=NPRTRW)

!write (0,*) __FILE__, __LINE__,'; cudaDeviceSynchronize returns ',cudaDeviceSynchronize()

! 1km resolution -- doesn't matter for FFT's
ZDX=1.e3
ZDY=1.e3
!write (*,*) __FILE__, __LINE__; call flush(6)
CALL ESETUP_TRANS(KMSMAX=NMSMAX,KSMAX=NSMAX,KDGL=NLAT,KDGUX=NDGUX, &
 & KLOEN=NLOEN,LDSPLIT=LSPLIT,PEXWN=ZDX,PEYWN=ZDY,LDUSEFFTW=.FALSE.)

! Find distribution properties: total number of waves, number of waves on this proc, total number of gridpoints, gridpoints on this proc
!write (*,*) __FILE__, __LINE__; call flush(6)
CALL ETRANS_INQ(KSPEC2=NSPEC2,KSPEC2G=NSPEC2G,KGPTOT=NGPTOT,KGPTOTG=NGPTOTG,KMYSETV=MYSETV)


! number of blocks
NBLK=(NGPTOT-1)/NPROMA+1

! Setup vertical distribution: IVSET contains on which proc each field is; NFLDSP is the number of fields on MYPROC
ALLOCATE(IVSET(NFLD))
IF (NPRTRV==1) THEN
  IVSET(:)=1
  NFLDSP=NFLD
ELSE
  NFLDSP=0
  JST=1
  DO JPRV=1,NPRTRV
    IF ( JPRV <= NFLD-NPRTRV*(NFLD/NPRTRV)) THEN
      JNUM=NFLD/NPRTRV+1
    ELSE
      JNUM=NFLD/NPRTRV
    ENDIF
    IVSET(JST:JST+JNUM-1)=JPRV
    IF (JPRV==MYSETV) NFLDSP=JNUM
    JST=JST+JNUM
  ENDDO

      
  !DO IFLD=1,NFLD
  !  IVSET(IFLD)=MODULO(IFLD-1,NPRTRV)+1   ! not optimal: clustered (1,1,...,2,2,...) maybe better than scattered (1,2,1,2,...)
  !  IF (IVSET(IFLD)==MYSETV) NFLDSP=NFLDSP+1
  !ENDDO

ENDIF
      
! Allocations for data
ALLOCATE(ZSPEC(NFLDSP,NSPEC2))   ! Spectral coefficients (all fields)
ALLOCATE(ZG(NGPTOTG,NFLD))      ! global gridpoint field
ALLOCATE(Z(NPROMA,NFLD,NBLK))     ! local gridpoint field array (all fields)

! Array necessary for distributing data from proc 1 to others
ALLOCATE(ITO(NFLD))

! first proc fills global array
IF (MYPROC==1) THEN
  KX=0
  KY=1
  
  ZG=0.
  DO JLAT=1,NLAT
    DO JLON=1,NLON
	  ZG((JLAT-1)*NLON+JLON,1) = COS( KX * 2*RPI*(REAL(JLON-1)/NLON) ) * COS ( KY * 2*RPI*(REAL(JLAT-1)/NLAT) )
	ENDDO
  ENDDO
ENDIF

! Synchronize processors
IF(NPROC > 1) THEN
	CALL MPL_BARRIER()
ENDIF

! Distribute gridpoint field
ITO(:)=1
!write (*,*) __FILE__, __LINE__; call flush(6)


CALL EDIST_GRID(KPROMA=NPROMA,PGPG=ZG,KFDISTG=NFLD,KFROM=ITO,PGP=Z)
!write (*,*) __FILE__, __LINE__; call flush(6)


write (*,*) 'before transform:'
write (*,*) ' Z = ',Z(1:6,1,1),'...'

! Deallocate arrays
DEALLOCATE(ZG)


! Allocate norm arrays
ALLOCATE(ZNORM(NFLD))

DO IITER=1,NITERS
! First spectral transform: we want to start the timeloop in spectral space: that's where the norms are calculated
!write (*,*) __FILE__, __LINE__; call flush(6)
CALL EDIR_TRANS(KPROMA=NPROMA,PSPSCALAR=ZSPEC,KVSETSC=IVSET,PGP=Z)
!write (*,*) __FILE__, __LINE__; call flush(6)

!write (*,*) 'ZSPEC = ',ZSPEC(1,1:6),'...'

! Calculate  and write initial norms
IF (LNORMS) THEN
	CALL ESPECNORM(PSPEC=ZSPEC,KVSET=IVSET,PNORM=ZNORM)
	IF(MYPROC == 1) THEN
		WRITE (NOUT,*)
		DO IFLD=1,NFLD
		  WRITE(NOUT,'(" ZNORM(",I4,")=",F20.15)') IFLD,ZNORM(IFLD)
		ENDDO
		WRITE (NOUT,*)
	ENDIF
ENDIF

! inverse transform
CALL EINV_TRANS(KPROMA=NPROMA,PSPSCALAR=ZSPEC,KVSETSC=IVSET,PGP=Z)

!write (*,*) 'Z = ',Z(1:5,1,1),'...'

ENDDO

write (*,*) 'spectral:'
write (*,*) ' ZSPEC = '
write (*,'(5X,4F12.4)') ZSPEC(1,:)

write (*,*) 'after transform:'
write (*,*) ' Z = ',Z(1:6,1,1),'...'


! Close down message passing
CALL MPL_BARRIER()


CALL GSTATS(0,1)
CALL GSTATS_PRINT(NOUT,ZAVEAVE,JPMAXSTAT)

CALL MPL_END()

END PROGRAM DRIVER_LAM
