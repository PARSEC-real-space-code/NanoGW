##
# Ubuntu Linux, intel Fortran compiler ifx, MPI, MKL FFT with FFTW3 interface
#
# Variables HDF5_LIB and HDF5_INC should be set before compiling:
# 
# export HDF5_LIB=...
# export HDF5_INC=...
# 
# Note: Computing time may be unstable on WSL due to core distribution 
#       and competition with other applications running outside WSL.
#

FCPP    = /usr/bin/cpp -P -traditional-cpp
CPPOPT  = -DMPI -DUSEFFTW3 -DDEBUG

EXT     = .mpi

F90s    = mpiifx # Serial compiler
F90     = $(F90s)

EXTRA_DIR     = ../LIB

OPTS    = -O3 -qmkl=cluster -xHost 
OPTS2   = -I${HDF5_INC}

LIBLAPACK = 

LIBXC = -L${EXTRA_DIR}/libstring_f -L${EXTRA_DIR}/libxc/src -lxc -lstring_f 

LIBHDF5 = -L${HDF5_LIB} -lhdf5_fortran

AUX_SRC = aux_generic.f90

