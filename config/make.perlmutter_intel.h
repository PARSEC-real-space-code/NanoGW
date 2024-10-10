##
#
# NERSC Perlmutter, Intel compiler, FFTW 3, MPI
#
# module load cpu intel/2023.2.0 cray-hdf5-parallel/1.12.2.9
#
##

FCPP    = /usr/bin/cpp -P -traditional-cpp
CPPOPT  = -DUSEFFTW3 -DMPI # -DDEBUG

EXT     = .intel.mpi

F90s    = ftn
F90     = $(F90s)

HDF5_DIR   = /opt/cray/pe/hdf5/1.12.2.9/intel/2023.2
HDF5_INC   = ${HDF5_DIR}/include
HDF5_LIB   = ${HDF5_DIR}/lib
EXTRA_DIR  = ../LIB

OPTS       = -O3 -qmkl=cluster # -check all
OPTS2      = -I${HDF5_INC}

LIBXC      = -L${EXTRA_DIR}/libstring_f/lib -L${EXTRA_DIR}/libxc/lib -lxc -lstring_f

LIBHDF5    = -L${HDF5_LIB} -lhdf5_fortran_intel -lhdf5_fortran

AUX_SRC = aux_generic.f90

