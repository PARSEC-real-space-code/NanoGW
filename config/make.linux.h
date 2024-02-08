##
# Ubuntu Linux, gfortran Fortran compiler, FFTW 3.*, MPI
#
FCPP    = /usr/bin/cpp -P -traditional-cpp
CPPOPT  = -DUSEFFTW3 -DMPI -DDEBUG -DGNU # DINTEL

EXT     = .mpi

F90s    = mpifort # Serial compiler
F90     = $(F90s)

FFTW3_INC_DIR = /usr/include/
FFTW3_LIB_DIR = /usr/lib/x86_64-linux-gnu/ 
HDF5_INC_DIR  = /usr/include/hdf5/serial
HDF5_LIB_DIR  = /usr/lib/x86_64-linux-gnu/
EXTRA_DIR     = /home/weiwei/Documents/NanoGW/NanoGW_optimized2-master/LIB
LAPACK_DIR    = /usr/lib/x86_64-linux-gnu/

OPTS    = -ffree-line-length-500 -O3 # -fcheck=all #-g 
OPTS2   = -I${HDF5_INC_DIR}

LIBFFT  = -I$(FFTW3_INC_DIR) -L$(FFTW3_LIB_DIR) -lfftw3f_mpi -lfftw3

LIBLAPACK = -L${LAPACK_DIR} -llapack -lblas -lscalapack-openmpi

LIBXC = -L${EXTRA_DIR}/libstring_f -L${EXTRA_DIR}/libxc/src -lxc -lstring_f 

LIBHDF5 = -L${HDF5_LIB_DIR} -lhdf5_serial_fortran 

AUX_SRC = aux_generic.f90

