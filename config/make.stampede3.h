##
# stampede3.tacc.utexas.edu, Intel Fortran compiler ifx, FFTW 3.*, MPI, HDF5
#
# Note: module load phdf5/1.14.3
# Currently Loaded Modules:
#  1) intel/24.0   3) autotools/1.4   5) xalt/3.0.1   7) phdf5/1.14.3
#  2) impi/21.11   4) cmake/3.28.1    6) TACC
#
FCPP    = /usr/bin/cpp -P -traditional-cpp
CPPOPT  = -DUSEFFTW3 -DMPI -DINTEL # -DDEBUG

EXT     = .mpi

F90s    = mpif90 # Serial compiler
F90     = $(F90s)

OPTS    = -O3 -xCORE-AVX512 -qmkl=cluster # -check bounds
OPTS2   = -I${TACC_HDF5_INC} 

FFTW_DIR = 

LIBLAPACK = 

EXTRA_DIR=../LIB
LIBXC = -L${EXTRA_DIR}/libstring_f -L${EXTRA_DIR}/libxc/src -lxc -lstring_f # -lxcf90

LIBHDF5 = -L${TACC_HDF5_LIB} -lhdf5_fortran

AUX_SRC = aux_generic.f90

##

