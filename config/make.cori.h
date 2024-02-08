##
# cori.nersc.gov, Intel Fortran compiler, FFTW 3.*, MPI
# module load cray-hdf5 cray-fftw
#
FCPP    = /usr/bin/cpp -P -traditional-cpp
CPPOPT  = -DUSEFFTW3 -DMPI -DINTEL #-DDEBUG

EXT     = .mpi.cori

F90s    = ftn #-mkl # Serial compiler
F90     = $(F90s)

#OPTS   = -g
OPTS    = -no-ipo

FFTW_DIR =
LIBFFT  = -I$(FFTW_INC) -L$(FFTW_DIR) -lm

EXTRA_DIR=../LIB
LIBXC = -L${EXTRA_DIR}/libstring_f -L${EXTRA_DIR}/libxc/src -lxc -lstring_f 
# LIBXC = -L/global/homes/w/weiwei/lib/libstring_f -L/global/homes/w/weiwei/lib/libxc/src  -lxc -lstring_f   # -lxcf90
# LIBXC = -L/global/homes/w/weiwei/lib/libxc2.2.1/mylib/lib -lxcf90 -lxc

AUX_SRC = aux_generic.f90

##

