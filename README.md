# etrans
LAM spectral transforms, related to ecmwf-ifs/ectrans

## Compilation/Installation/Running

### Environment settings

    # (load modules, I used NVHPC/21.2 ... )
    BASEDIR=${HOME}/etrans_test    # modify as you wish
    SOURCEDIR=${BASEDIR}/sources
    BUILDDIR=${BASEDIR}/build
    INSTALLDIR=${BASEDIR}/install
    mkdir -p ${BASEDIR} ${SOURCEDIR} ${BUILDDIR} ${INSTALLDIR}
    export CC=mpicc
    export CXX=mpicxx
    export FC=mpif90
    export PATH=${PATH}:${INSTALLDIR}/ecbuild/bin/
    TOOLCHAIN_FILE=${SOURCEDIR}/my_toolchain.cmake    # Note: a sample toolchain file is provided under etrans/cmake/toolchain/sample_toolchain.cmake
    
    cd ${BASEDIR}

### Prerequisites

* ecbuild, installed with

        cd ${SOURCEDIR}
        git clone git@github.com:ecmwf/ecbuild.git
        cd ecbuild
        git checkout master   # at time of writing, commit c39e3b0aaa1b9da7b5d6d9d419e8e25f37d70794
        rm -rf ${BUILDDIR}/ecbuild ${INSTALLDIR}/ecbuild
        mkdir -p ${BUILDDIR}/ecbuild
        cd ${BUILDDIR}/ecbuild
        ${SOURCEDIR}/ecbuild/bin/ecbuild --prefix=${INSTALLDIR}/ecbuild -DENABLE_INSTALL=ON ${SOURCEDIR}/ecbuild 
        make
        make install
    
* eccodes, installed with

        cd ${SOURCEDIR}
        git clone git@github.com:ecmwf/eccodes.git
        cd eccodes
        git checkout master    # at time of writing, commit 6efa9fc52899863d32311e5498f79abde8f303f3
        rm -rf ${BUILDDIR}/eccodes
        mkdir -p ${BUILDDIR}/eccodes
        cd ${BUILDDIR}/eccodes
        ecbuild --toolchain=${TOOLCHAIN_FILE} --prefix=${INSTALLDIR}/eccodes -DENABLE_MEMFS=ON -DENABLE_AEC=OFF ${SOURCEDIR}/eccodes
        make -j24
        rm -rf ${INSTALLDIR}/eccodes
        make install


* fiat, installed with

        cd ${SOURCEDIR}
        git clone git@github.com:ecmwf-ifs/fiat.git
        cd fiat
        git checkout main    # at time of writing, commit e214bef6aa3a67ebc89b675afc82ceb97e93c329
        rm -rf ${BUILDDIR}/fiat
        mkdir -p ${BUILDDIR}/fiat
        cd ${BUILDDIR}/fiat
        ecbuild --toolchain=${TOOLCHAIN_FILE} --prefix=${INSTALLDIR}/fiat ${SOURCEDIR}/fiat
        make -j8
        rm -rf ${INSTALLDIR}/fiat
        make install


* nvtx.mod; may be included in some NVHPC installations; otherwise compile it yourself from https://gist.github.com/jefflarkin/b64f63d79bdc978a2503 with

        nvfortran -c -fPIC nvtx.F90
        ar rcs libnvtx.a nvtx.o
        mkdir -p ${INSTALLDIR}/nvtx/include
        mv nvtx.mod ${INSTALLDIR}/nvtx/include
        mkdir -p ${INSTALLDIR}/nvtx/lib
        mv libnvtx.a ${INSTALLDIR}/nvtx/lib


### Compilation of ectrans

    cd ${SOURCEDIR}
    git clone git@github.com:ddegrauwe/ectrans.git
    cd ectrans
    git checkout gpu_daand_lam
    rm -rf ${BUILDDIR}/ectrans
    mkdir -p ${BUILDDIR}/ectrans
    cd ${BUILDDIR}/etrans
    ecbuild --toolchain=${TOOLCHAIN_FILE} --prefix=${INSTALLDIR}/ectrans -Dfiat_ROOT=${INSTALLDIR}/fiat ${SOURCEDIR}/ectrans
    make -j12
    rm -rf ${INSTALLDIR}/ectrans
    make install


### Compilation of etrans

    cd ${SOURCEDIR}
    git clone git@github.com:ddegrauwe/etrans.git
    cd etrans
    git checkout gpu_daand_lam
    rm -rf ${BUILDDIR}/etrans
    mkdir -p ${BUILDDIR}/etrans
    cd ${BUILDDIR}/etrans
    ecbuild --toolchain=${TOOLCHAIN_FILE} --prefix=${INSTALLDIR}/etrans  -Dectrans_ROOT=${INSTALLDIR}/ectrans -Dfiat_ROOT=${INSTALLDIR}/fiat -DNVTX_ROOT=${INSTALLDIR}/nvtx ${SOURCEDIR}/etrans
    make -j12
    rm -rf ${INSTALLDIR}/etrans
    make install

### Running test driver

    # path to shared libraries
    export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${INSTALLDIR}/fiat/lib64/:${INSTALLDIR}/ectrans/lib64/:${INSTALLDIR}/etrans/lib64/
    # create namelist file
    echo -e "&NAMTRANS\nNITERS=1\n/\n&NAMBIFFT\nNLON=6,NLAT=4,NFLD=1,NPROMA=24\n/\n&NAMCT0\n/\n" > fort.4
    # some environment settings
    export DR_HOOK=1
    export DR_HOOK_OPT=prof
    export OMP_NUM_THREADS=1
    export PGI_ACC_CUDA_HEAPSIZE=$((1024*1024*1024*8))
    export NPROC_PERNODE='1'
    # run with profiling
    precision=dp
    mpirun -np 1 nsys profile -t openacc,cuda -f true -o report_nsys ./dwarf_bifft_${precision} > log.out 2>log.err
    
    

