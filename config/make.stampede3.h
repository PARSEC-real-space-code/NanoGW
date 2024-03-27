##
# stampede3.tacc.utexas.edu, Intel Fortran compiler ifx, FFTW 3.*, MPI, HDF5
#
# Note: module load intel/24.0 hdf5/1.14.3
#
FCPP    = /usr/bin/cpp -P -traditional-cpp
CPPOPT  = -DUSEFFTW3 -DMPI -DINTEL # -DDEBUG

EXT     = .mpi

F90s    = mpif90 # Serial compiler
F90     = $(F90s)

OPTS    = -O3 -axCORE-AVX512 -qmkl=cluster # -check bounds
OPTS2   = -I${TACC_HDF5_INC} 

FFTW_DIR = 

LIBLAPACK = 

EXTRA_DIR=../LIB
LIBXC = -L${EXTRA_DIR}/libstring_f -L${EXTRA_DIR}/libxc/src -lxc -lstring_f # -lxcf90

LIBHDF5 = -Wl,-rpath,${TACC_HDF5_LIB} -L${TACC_HDF5_LIB} -lhdf5_fortran -lz

AUX_SRC = aux_generic.f90

##

