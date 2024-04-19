##
# Ubuntu Linux, gfortran Fortran compiler, FFTW 3.*, MPI
#
FCPP    = /usr/bin/cpp -P -traditional-cpp
CPPOPT  = -DUSEFFTW3 -DMPI 

EXT     = .mpi

F90s    = mpiifx # Serial compiler
F90     = $(F90s)

HDF5_INC_DIR  = /home/weiwei/Programs/hdf5-1.12.3/hdf5/include 
HDF5_LIB_DIR  = /home/weiwei/Programs/hdf5-1.12.3/hdf5/lib
EXTRA_DIR     = /home/weiwei/Programs/nanogw/LIB
LAPACK_DIR    = /opt/intel/oneapi/mkl/latest/

OPTS    = -qmkl=cluster # -check bounds -check shape #
OPTS2   = -I${HDF5_INC_DIR}


LIBLAPACK = 


LIBXC = -L${EXTRA_DIR}/libstring_f -L${EXTRA_DIR}/libxc/src -lxc -lstring_f 

LIBHDF5 = -L${HDF5_LIB_DIR} -lhdf5_fortran

AUX_SRC = aux_generic.f90

