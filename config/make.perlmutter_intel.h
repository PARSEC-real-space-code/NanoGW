##
# Ubuntu Linux, gfortran Fortran compiler, FFTW 3.*, MPI
#
FCPP    = /usr/bin/cpp -P -traditional-cpp
CPPOPT  = -DUSEFFTW3 -DMPI # -DDEBUG 

EXT     = .intel.mpi

F90s    = ftn
F90     = $(F90s)

#FFTW3_INC_DIR = 
#FFTW3_LIB_DIR = 
HDF5_INC_DIR  = /opt/cray/pe/hdf5/1.12.2.3/intel/19.0/include 
HDF5_LIB_DIR  = /opt/cray/pe/hdf5/1.12.2.3/intel/19.0/lib
EXTRA_DIR     = /global/u1/w/weiwei/src/new-nanogw/nanogw-2023_8_23/LIB_intel
LAPACK_DIR    = 

OPTS       = -O2 # -check all
OPTS2      = -I${HDF5_INC_DIR} 

LIBLAPACK  = -qmkl="cluster" 

LIBXC      = -L${EXTRA_DIR}/libstring_f/lib -L${EXTRA_DIR}/libxc/lib -lxc -lstring_f 

LIBHDF5    = -L${HDF5_LIB_DIR} -lhdf5_fortran_intel -lhdf5_fortran

AUX_SRC = aux_generic.f90

