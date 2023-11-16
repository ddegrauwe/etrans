program driver_lam
	use parkind1,only: jpim,jprb,jprd
	use mpl_module
	use yomgstats

	use cudafor

	implicit none

	integer(kind=jpim) :: nerr,nulnam,nproc,nlin,nsmax,nq,nout,nspec2,ngptot,ngptotg,ifld,&
		nproma,nblk,itag,nspec2g,myproc,jj,i,niters,nlon,nlat,nmsmax,ndgux,nlate,nprgpns,&
		nprgpew,nprtrv,nprtrw,mysetv,jst,jnum,jprv,ifile,nfld,ngrib,jf,nfldsp,igrib,jp,&
		jlat,kx,ky
	integer :: jstep
	integer(kind=jpim),allocatable :: nloen(:),ito(:),nprcids(:),ivset(:)
	real(kind=jprb) :: ztinit,ztloop,ztstep,timef
	real(kind=jprd) :: zaveave(0:jpmaxstat)
	real(kind=jprb) :: zmaxerr_check=0,zdx,zdy
	real(kind=jprb),parameter :: rpi=acos(-1._jprb)
	real(kind=jprb),allocatable :: zspec(:,:),znorm(:),znorm1(:),zg(:,:),z(:,:,:)
	logical :: lsplit,luserpnm,lkeeprpnm,luseflt,lltrace_stats,llstats_omp,llstats_comms,&
		lnorms=.false.
	character(len=256) :: cdatafile

	namelist/namtrans/lbarrier_stats,lbarrier_stats2,ldetailed_stats,luserpnm,lkeeprpnm,&
		luseflt,nq,nlin,lnorms,niters,zmaxerr_check
	namelist/nambifft/nlate,nfld,nlon,nlat,nproma
	namelist/namct0/nprgpns,nprgpew,nprtrw,nprtrv

#include "setup_trans0.h"
#include "esetup_trans.h"
#include "einv_trans.h"
#include "edir_trans.h"
#include "edist_grid.h"
#include "etrans_inq.h"
#include "especnorm.h"
#include "abor1.intfb.h"

	nerr = 0
	nulnam = 4
	nout = 6

	call mpl_init()

	call gstats(0,0)
	ztinit = timef()

	nproc = mpl_nproc()
	myproc = mpl_myrank()
	if (myproc /= 1) open(nout,file="/dev/null")

	write(nout,*) "NPROC:",nproc

	allocate(nprcids(nproc))
	do jj=1,nproc
		nprcids(jj) = jj
	end do

	lstats = .true.
	ldetailed_stats = .false.
	lstats_omp = .false.
	lstats_comms = .false.
	lstats_mpl = .false.
	lbarrier_stats = .false.
	lbarrier_stats2 = .false.
	luserpnm = .false.
	lkeeprpnm = .false.
	luseflt = .true.
	lnorms = .false.
	niters = 10

	nlate = 0
	nfld = 1
	nlon = 48
	nlat = 40
	nproma = 16

	nprgpns = nproc
	nprgpew = 1
	nprtrw = nproc
	nprtrv = 1
	lsplit = .false.

	write(nout,*) "Read settings from namelist"
	open(unit=nulnam,file="driver/driver.nam",form="formatted")
	read(nulnam,namtrans)
	rewind(nulnam)
	read(nulnam,nambifft)
	rewind(nulnam)
	read(nulnam,namct0)
	close(nulnam)

	write(nout,*) "Set up distribution: A/B sets"
   if (nproc == 1) then
		nprgpns = 1
		nprgpew = 1
	else
		call setab(nproc,nprgpns,nprgpew)
	   if (nprgpns*nprgpew /= nproc) call abor1("Error: no partition found for A/B sets")
	end if

	write(nout,*) "Distribution values (NS EW W V):",nprgpns,nprgpew,nprtrw,nprtrv

	ndgux = nlat-nlate

	write(nout,*) "NLAT/NLON/NDGUX:",nlat,nlon,ndgux

	write(nout,*) "Setup statistics"
	write(nout,'(" LDETAILED_STATS=",L1," LBARRIER_STATS=",L1," LBARRIER_STATS2=",&
		L1)') ldetailed_stats,lbarrier_stats,lbarrier_stats2
	if (ldetailed_stats) then
		lstats_omp = .true.
		lstats_comms = .true.
		lstats_mpl = .true.
	end if

	call gstats_setup(nproc,myproc,nprcids,lstats,lstatscpu,lsyncstats,ldetailed_stats,&
		lbarrier_stats,lbarrier_stats2,llstats_omp,llstats_comms,lstats_mem,nstats_mem,&
		lstats_alloc,lltrace_stats,ntrace_stats,nprnt_stats,lxml_stats,lstats_mpl)
	call gstats_psut

	nmsmax = (nlon-1)/2
	nsmax = (nlat-1)/2

	write(nout,*) "Setup distribution - part 0"
	allocate(nloen(nlat))
	nloen(:) = nlon

	call setup_trans0(kout=nout,kerr=nerr,kprintlev=0,kmax_resol=1,kprgpns=nprgpns,&
		kprgpew=nprgpew,kprtrw=nprtrw)

	write(nout,*) "Setup distribution - 1km resolution"
	zdx = 1.e3
	zdy = 1.e3
	call esetup_trans(kmsmax=nmsmax,ksmax=nsmax,kdgl=nlat,kdgux=ndgux,kloen=nloen,&
		ldsplit=lsplit,pexwn=zdx,peywn=zdy,ldusefftw=.false.)

	call etrans_inq(kspec2=nspec2,kspec2g=nspec2g,kgptot=ngptot,kgptotg=ngptotg,&
		kmysetv=mysetv)

	nblk = (ngptot-1)/nproma+1

	write(nout,*) "Runtime parameters:"
	write(nout,'("ngptot, ngptotg ",2(x,i0))') ngptot,ngptotg
	write(nout,'("nproma, nblk ",2(x,i0))') nproma,nblk
	write(nout,'("nspec2, nspec2g ",2(x,i0))') nspec2,nspec2g

	allocate(ivset(nfld))
	if (nprtrv == 1) then
		ivset(:) = 1
		nfldsp = nfld
	else
		nfldsp = 0
		jst = 1
		do jprv=1,nprtrv
			if (jprv <= nfld-nprtrv*(nfld/nprtrv)) then
				jnum = nfld/nprtrv+1
			else
				jnum = nfld/nprtrv
			end if
			ivset(jst:jst+jnum-1) = jprv
			if (jprv == mysetv) nfldsp = jnum
			jst = jst+jnum
		end do
	end if

	write(nout,*) "Initialize spectral/gridpoint data, 3D fields:",nfld
	allocate(zspec(nfldsp,nspec2))
	allocate(zg(ngptotg,nfld),z(nproma,nfld,nblk))

	allocate(ito(nfld))

	if (myproc == 1) then
		kx = 0
		ky = 1

		zg = 0
		do jlat=1,nlat
			do jp=1,nlon
				zg((jlat-1)*nlon+jp,1) = cos(2*rpi*kx*(jp-1)/nlon)*cos(2*rpi*ky*(jlat-1)/nlat)
			end do
		end do
	end if

	if (nproc > 1) call mpl_barrier()

	ito(:) = 1

	call edist_grid(kproma=nproma,pgpg=zg,kfdistg=nfld,kfrom=ito,pgp=z)

	write(nout,*) "Z(1:6):",z(1:6,1,1)

	deallocate(zg)

	allocate(znorm(nfld))

	do i=1,niters
		write(nout,"('. iteration',x,i0,'/',i0)") i,niters
		call edir_trans(kproma=nproma,pspscalar=zspec,kvsetsc=ivset,pgp=z)

		if (lnorms) then
			call especnorm(pspec=zspec,kvset=ivset,pnorm=znorm)

			do ifld=1,nfld
				write(nout,'("ZNORM(",I4,"):",F20.15)') ifld,znorm(ifld)
			end do
		end if

		write(nout,"('. inverse transforms')")
		call einv_trans(kproma=nproma,pspscalar=zspec,kvsetsc=ivset,pgp=z)
	end do

	write(nout,*) "ZSPEC:",zspec(1,1:5)
	write(nout,*) "Z:",minval(z),sum(z)/size(z),maxval(z)
	write(nout,*) "Z(1:6):",z(1:6,1,1)

	call mpl_barrier()

	call gstats(0,1)
	call gstats_print(nout,zaveave,jpmaxstat)

	call mpl_end()
contains
	subroutine setab(n,na,nb)
		integer(kind=jpim),intent(in) :: n
		integer(kind=jpim),intent(out) :: na,nb

		integer(kind=jpim) :: ja,ib,isqr

		isqr = sqrt(real(n))

		do ja=isqr,n
			ib = n/ja
			if (ja*ib /= n) cycle

			na = max(ja,ib)
			nb = min(ja,ib)
			exit
		end do
	end subroutine
end program
