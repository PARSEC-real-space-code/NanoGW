##
# Ubuntu Linux, gfortran Fortran compiler, FFTW 3.*, MPI
#
FCPP    = /usr/bin/cpp -P -traditional-cpp
CPPOPT  = -DUSEFFTW3 -DMPI # -DDEBUG # -DHIPMAGMA -Dmagma_devptr_t="integer*8" 

EXT     = .mpi

F90s    = mpiifort # Serial compiler
F90     = $(F90s)

#FFTW3_INC_DIR = 
#FFTW3_LIB_DIR = 
HDF5_INC_DIR  = /public/software/mathlib/hdf5/1.8.20/intel/include 
HDF5_LIB_DIR  = /public/software/mathlib/hdf5/1.8.20/intel/lib
EXTRA_DIR     = /work/home/accflz0nn6/user/weiwei/NanoGW-dev/nanogw-2022_6_21/LIB
LAPACK_DIR    = /public/software/compiler/intel-compiler/2021.3.0/mkl

OPTS       = -O2 -mkl=sequential # -check all
OPTS2      = -I${HDF5_INC_DIR} 

LIBLAPACK  = -L${LAPACK_DIR}/lib/intel64 -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64 -lpthread -lm -ldl

LIBXC      = -L${EXTRA_DIR}/libstring_f/lib -L${EXTRA_DIR}/libxc/lib -lxc -lstring_f 

LIBHDF5    = -L${HDF5_LIB_DIR} -lhdf5_fortran

AUX_SRC = aux_generic.f90

