##
# stampede2.tacc.utexas.edu, Intel Fortran compiler, FFTW 3.*, MPI, HDF5
#
# Note: module load intel hdf5/1.8.16
#
FCPP    = /usr/bin/cpp -P -traditional-cpp
CPPOPT  = -DUSEFFTW3 -DMPI -DINTEL # -DDEBUG

EXT     = .knl_skx.mpi

F90s    = mpif90 # Serial compiler
F90     = $(F90s)

#OPTS   = -g
#OPTS   = -xCORE-AVX2 -axCORE-AVX512,MIC-AVX512 -O2 -mkl=cluster # -check bounds
OPTS    = -O3 -mkl=cluster # -check bounds
OPTS2   = -I${TACC_HDF5_INC} 

FFTW_DIR = 

LIBLAPACK = 

EXTRA_DIR=../LIB
LIBXC = -L${EXTRA_DIR}/libstring_f -L${EXTRA_DIR}/libxc/src -lxc -lstring_f # -lxcf90

LIBHDF5 = -Wl,-rpath,${TACC_HDF5_LIB} -L${TACC_HDF5_LIB} -lhdf5_fortran -lz

AUX_SRC = aux_generic.f90

##

