##
# stampede2.tacc.utexas.edu, Intel Fortran compiler, FFTW 3.*, MPI, HDF5
#
# Note: module load intel hdf5/1.8.16
#
FCPP    = /usr/bin/cpp -P -traditional-cpp
CPPOPT  = -DUSEFFTW3 -DMPI -DINTEL # -DDEBUG

EXT     = .mpi

F90s    = mpiifort # Serial compiler
F90     = $(F90s)

#FFTW3_INC_DIR = /public/home/users/app/lib/fftw/intel/double/include
#FFTW3_LIB_DIR = /public/home/users/app/lib/fftw/intel/double/lib
HDF5_INC_DIR  = /public/home/users/app/lib/hdf5.2017/include
HDF5_LIB_DIR  = /public/home/users/app/lib/hdf5.2017/lib
EXTRA_DIR     = /public/home/users/dlut010/source/nanogw/LIB
LAPACK_DIR    = /public/software/compiler/intel/intel-compiler-2017.5.239/mkl

OPTS    = -O2  -mkl=parallel # -check all # -g 
OPTS2   = -I${HDF5_INC_DIR}

#LIBFFT  = -I$(FFTW3_INC_DIR) -L$(FFTW3_LIB_DIR) -lfftw3_mpi -lfftw3

#LIBLAPACK = -L${LAPACK_DIR}  \
#       ${LAPACK_DIR}/lib/intel64/libmkl_scalapack_ilp64.a -Wl,--start-group \
#       ${LAPACK_DIR}/lib/intel64/libmkl_intel_ilp64.a \
#       ${LAPACK_DIR}/lib/intel64/libmkl_sequential.a \
#       ${LAPACK_DIR}/lib/intel64/libmkl_core.a \
#       ${LAPACK_DIR}/lib/intel64/libmkl_blacs_intelmpi_ilp64.a -Wl,--end-group -liomp5 -lpthread -lm -ldl
LIBLAPACK =   -L${LAPACK_DIR}/lib/intel64 -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64 -liomp5 -lpthread -lm -ldl  

LIBXC = -L${EXTRA_DIR}/libstring_f -L${EXTRA_DIR}/libxc/src -lxc -lstring_f 

LIBHDF5 = -L${HDF5_LIB_DIR} -lhdf5_fortran 

AUX_SRC = aux_generic.f90

##

