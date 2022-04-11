# etrans
LAM spectral transforms, related to ecmwf-ifs/ectrans

## Compilation/Installation/Running

### Prerequisites

* ecbuild
* eccodes: make sure to configure with -DENABLE_MEMFS=ON 
* fiat
* (nvtx.mod; may be included in some NVHPC installations; otherwise compile it yourself from https://gist.github.com/jefflarkin/b64f63d79bdc978a2503)

(all supposed to be installed under INSTALLDIR)

* sources of ectrans and etrans are supposed to be in SOURCEDIR

### Environment

    # (load modules)
    BASEDIR=${HOME}/etrans_test    # modify as you wish
    SOURCEDIR=${BASEDIR}/sources
    BUILDDIR=${BASEDIR}/build
    INSTALLDIR=${BASEDIR}/install
    mkdir -p ${BASEDIR} ${SOURCEDIR} ${BUILDDIR} ${INSTALLDIR}
    export CC=mpicc
    export CXX=mpicxx
    export FC=mpif90
    export PATH=${PATH}:${INSTALLDIR}/ecbuild/bin/

    export INSTALLDIR
    cd ${BASEDIR}

### Compilation of ectrans

Note: a sample toolchain file is provided under cmake/toolchain/

    cd ${SOURCEDIR}
    rm -rf ${BUILDDIR}/ectrans
    mkdir -p ${BUILDDIR}/ectrans
    cd ${BUILDDIR}/etrans
    ecbuild --toolchain=${TOOLCHAIN_FILE} --prefix=${INSTALLDIR}/ectrans -Dfiat_ROOT=${INSTALLDIR}/fiat ${SOURCEDIR}/ectrans
    make -j12
    rm -rf ${INSTALLDIR}/ectrans
    make install


### Compilation of etrans

    cd ${SOURCEDIR}
    rm -rf ${BUILDDIR}/etrans
    mkdir -p ${BUILDDIR}/etrans
    cd ${BUILDDIR}/etrans
    ecbuild --toolchain=${TOOLCHAIN_FILE} --prefix=${INSTALLDIR}/etrans  -Dectrans_ROOT=${INSTALLDIR}/ectrans -Dfiat_ROOT=${INSTALLDIR}/fiat ${SOURCEDIR}/etrans
    make -j12
    rm -rf ${INSTALLDIR}/etrans
    make install

