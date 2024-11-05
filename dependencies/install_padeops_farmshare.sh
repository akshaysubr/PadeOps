#!/bin/bash
CURRENT_DIR="${PWD}"

module purge
module load cmake intel
export COMPILER_ID=Intel
export FC=mpiifort
export CC=mpiicc
export CXX=mpiicpc
export I_MPI_CC=icc
export I_MPI_CXX=icpc
export I_MPI_F90=ifort
export I_MPI_F77=ifort



#####################################################
# Download PadeOps and configure the basic settings #
#####################################################

cd .. && export PADEOPS_PATH="${PWD}"

ls $PADEOPS_PATH/dependencies

cd $PADEOPS_PATH/dependencies \
    && tar -zvxf fftw-3.3.10.tar.gz \
    && tar -zvxf 2decomp_fft-1.5.847.tar.gz \
    && tar -zvxf Lib_VTK_IO.tar.gz \
    && tar -zvxf hdf5-1.8.18.tar.gz

export FFTW_PATH=$PADEOPS_PATH/dependencies/fftw-3.3.5
export DECOMP_PATH=$PADEOPS_PATH/dependencies/2decomp_fft
export VTK_IO_PATH=$PADEOPS_PATH/dependencies/Lib_VTK_IO
export HDF5_PATH=$PADEOPS_PATH/dependencies/hdf5-1.8.18


################
# INSTALL FFTW #
################
cd $FFTW_PATH
F77=mpiifort MPICC=mpiicc ./configure --prefix=$FFTW_PATH --enable-avx \
    && make && make install


###################
# INSTALL 2DECOMP #
###################
cd $DECOMP_PATH/src \
    && perl -pi -e "s/FFT=generic/FFT=fftw3_f03/ if $. == 25"        Makefile.inc.x86 \
    && perl -pi -e "s/.*// if $. == 32"                              Makefile.inc.x86 \
    && perl -pi -e 'print "FFTW_PATH=$ENV{'FFTW_PATH'}" if $. == 32' Makefile.inc.x86 \
    && perl -pi -e "s/F90=mpif90/F90=mpiifort/ if $. == 57"          Makefile.inc.x86 \
    && perl -pi -e "s/CC=mpicc/CC=mpiicc/ if $. == 71"               Makefile.inc.x86 \
    && perl -pi -e "s/-lfftw3f// if $. == 77"                        Makefile.inc.x86 \
    && ln -s Makefile.inc.x86 Makefile.inc \
    && cd .. \
    && make


###############
# INSTALL VTK #
###############
cd $VTK_IO_PATH \
    && mkdir -p build \
    && cd build \
    && rm -rf ./* \
    && CC=mpiicc CXX=mpiicpc FC=mpiifort cmake .. \
    && make


################
# INSTALL HDF5 #
################
cd $HDF5_PATH \
    && CC=mpiicc CXX=mpiicpc FC=mpiifort ./configure \
    --enable-parallel --enable-fortran --enable-build-mode=production \
    --prefix=$HDF5_PATH \
    && make \
    && make install


#############################################
# Add the configuration to "~/.bashrc" file #
#############################################
echo
echo
echo '# BEGIN OF AUTO-GENERATED SCRIPT BY PADEOPS_AUTO_INSTALL'
echo 'function setup_padeops() {'
echo '    module purge'
echo '    module load cmake intel'
echo '    export COMPILER_ID=Intel'
echo '    export FC=mpiifort'
echo '    export CC=mpiicc'
echo '    export CXX=mpiicpc'
echo '    export I_MPI_CC=icc'
echo '    export I_MPI_CXX=icpc'
echo '    export I_MPI_F90=ifort'
echo '    export I_MPI_F77=ifort'
echo "    export PADEOPS_PATH=$PADEOPS_PATH"
echo "    export FFTW_PATH=$FFTW_PATH"
echo "    export DECOMP_PATH=$DECOMP_PATH"
echo "    export VTK_IO_PATH=$VTK_IO_PATH/build"
echo "    export HDF5_PATH=$HDF5_PATH"
echo '}'
echo "# END OF AUTO-GENERATED SCRIPT BY PADEOPS_AUTO_INSTALL"

echo
echo
echo "[PADEOPS] INSTALLATION HAS COMPLETED!"
cd $CURRENT_DIR

