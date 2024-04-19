##
# Ubuntu Linux, gfortran Fortran compiler, FFTW 3.*, MPI
#
FCPP    = /usr/bin/cpp -P -traditional-cpp
CPPOPT  = -DUSEFFTW3 -DMPI -DDEBUG 

EXT     = .nvcc_openblas.mpi

F90s    = mpif90 # Serial compiler
F90     = $(F90s)

FFTW3_INC_DIR = /home/weiwei/Programs/FFTW3/fftw-3.3.7-nv/LIBS/include
FFTW3_LIB_DIR = /home/weiwei/Programs/FFTW3/fftw-3.3.7-nv/LIBS/lib
HDF5_INC_DIR  = /home/weiwei/Programs/hdf5/hdf5-1.8.16-nv/lib/include 
HDF5_LIB_DIR  = /home/weiwei/Programs/hdf5/hdf5-1.8.16-nv/lib/lib
EXTRA_DIR     = /home/weiwei/Programs/nanogw/LIB
LAPACK_DIR    = /home/weiwei/Programs/lapack-3.10.0-nv
SCALAPACK_DIR = /home/weiwei/Programs/scalapack-2.0.2-nv
BLAS_DIR      = /home/weiwei/Programs/BLAS-3.10.0-nv

OPTS    = -O2 -fast
OPTS2   = -I${HDF5_INC_DIR}

LIBFFT  = -I$(FFTW3_INC_DIR) -L$(FFTW3_LIB_DIR) -lfftw3

LIBLAPACK = -L${LAPACK_DIR} -L${SCALAPACK_DIR}/lib -llapack -lblas -lscalapack

LIBXC = -L${EXTRA_DIR}/libstring_f_nv -L${EXTRA_DIR}/libxc_nv/src -lxc -lstring_f 

LIBHDF5 = -L${HDF5_LIB_DIR} -lhdf5_fortran

AUX_SRC = aux_generic.f90

