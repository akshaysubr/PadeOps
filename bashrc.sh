function setup_padeops(){
        module load gcc/7.3.1 spectrum-mpi cmake
        export COMPILER_ID=GNU
        export FC=mpif90
        export CC=mpicc
        export CXX=mpicxx
        export PADEOPS_ROOT_DIR=/g/g16/barbeau2/Codes/PadeOps
        export PADEOPS_TPL_DIR=$PADEOPS_ROOT_DIR/dependencies
        export PADEOPS_PATH=$PADEOPS_ROOT_DIR
        export PADEOPS_FFTW_PATH=$PADEOPS_TPL_DIR/fftw-3.3.5
        export DECOMP_PATH=$PADEOPS_TPL_DIR/2decomp_fft
        export VTK_IO_PATH=$PADEOPS_TPL_DIR/Lib_VTK_IO/build
        export HDF5_PATH=$PADEOPS_TBL_DIR/hdf5-1.8.18
}

