##
# Ubuntu Linux, gfortran Fortran compiler, FFTW 3.*, MPI
#
FCPP    = /usr/bin/cpp -P -traditional-cpp
CPPOPT  = -DUSEFFTW3 -DMPI 

EXT     = .cpu.mpi

F90s    = mpif90  # Serial compiler
F90     = $(F90s)
DTK_PATH = /public/software/compiler/rocm/dtk-23.04
HIPCC   = $(DTK_PATH)/bin/hipcc
HIPFORT_PATH = /public/home/ghfund2_a28/src/hipfort-master

FFTW3_INC_DIR = /public/home/ghfund2_a28/src/fftw-3.3.4/include
FFTW3_LIB_DIR = /public/home/ghfund2_a28/src/fftw-3.3.4/lib
HDF5_INC_DIR  = /public/home/ghfund2_a28/src/test-compiles/gcc/hdf5-1.8.16/hdf5/include
HDF5_LIB_DIR  = /public/home/ghfund2_a28/src/test-compiles/gcc/hdf5-1.8.16/hdf5/lib
EXTRA_DIR     = /public/home/ghfund2_a28/src/test-compiles/gcc/nanogw-2022_4_7/LIB
#LAPACK_DIR    = /public/home/ghfund2_a28/src/lapack-3.2.2
LAPACK_DIR    = /public/software/mathlib/lapack/gnu/3.8.0/lib64
SCALAPACK_DIR = /public/home/ghfund2_a28/src/scalapack-2.0.2

OPTS       =  -ffree-line-length-500 -Ofast -fcheck=all
OPTS2      = -I${HDF5_INC_DIR} -I$(FFTW3_INC_DIR) 

LIBFFT     = -L$(FFTW3_LIB_DIR) -lfftw3

# LIBLAPACK  = $(SCALAPACK_DIR)/libscalapack.a $(LAPACK_DIR)/blas_LINUX.a $(LAPACK_DIR)/lapack_LINUX.a
LIBLAPACK  = -L$(LAPACK_DIR) -lblas -llapack $(SCALAPACK_DIR)/libscalapack.a 

LIBXC      = -L${EXTRA_DIR}/libstring_f/lib -L${EXTRA_DIR}/libxc/lib -lxc -lstring_f 

LIBHDF5    = -L${HDF5_LIB_DIR} -lhdf5hl_fortran -lhdf5_fortran

LIBHIP = -L$(HIPFORT_PATH)/lib -lhipfort-amdgcn \
  -L$(DTK_PATH)/lib -lamdhip64 -L$(DTK_PATH)/hipblas/lib -lhipblas

AUX_SRC = aux_generic.f90

