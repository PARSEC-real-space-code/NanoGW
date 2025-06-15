##
#
# Ubuntu Linux, Nvidia compiler, OpenMPI, FFTW3, CUDA
#
# Following variables should be set before compiling:
#  export HDF5_LIB=...
#  export HDF5_INC=...
#  export FFTW3_LIB=...
#  export FFTW3_INC=...
# 
# Note: Computing time may be unstable on WSL due to core distribution 
#       and competition with other applications running outside WSL.
#

FCPP    = /usr/bin/cpp -P -traditional-cpp
CPPOPT  = -DMPI -DUSEFFTW3 -D_CUDA -DDEBUG

EXT     = .nvidia.mpi

F90s    = mpif90  # Serial compiler
F90     = $(F90s)

EXTRA_DIR = ../LIB

OPTS    = -fast -cuda -cudalib=cublas -gpu=ccnative
OPTS2   = -I${HDF5_INC} -I${FFTW3_INC}

LIBLAPACK = -Mscalapack
LIBFFT = -L${FFTW3_LIB} -lfftw3
LIBHDF5 = -L${HDF5_LIB} -lhdf5_fortran
LIBXC = -L${EXTRA_DIR}/libstring_f -L${EXTRA_DIR}/libxc/src -lxc -lstring_f 

AUX_SRC = aux_generic.f90

