##
#
# NERSC Perlmutter, Nvidia compiler, FFTW 3, MPI, CUDA
#
# module load PrgEnv-nvidia/8.5.0 cray-hdf5-parallel/1.14.3.1 cray-fftw/3.3.10.8
#
##

FCPP    = /usr/bin/cpp -P -traditional-cpp
CPPOPT  = -DUSEFFTW3 -DMPI -D_CUDA  # -DDEBUG

EXT     = .nvidia.mpi

F90s    = ftn
F90     = $(F90s)

HDF5_INC   = ${HDF5_DIR}/include
HDF5_LIB   = ${HDF5_DIR}/lib
EXTRA_DIR  = ../LIB

OPTS       = -fast -cuda -cudalib=cublas -gpu=cc80 # -Mbounds -Mchkptr -g -traceback
OPTS2      = -I${HDF5_INC}

LIBXC      = -L${EXTRA_DIR}/libstring_f/lib -L${EXTRA_DIR}/libxc/lib -lxc -lstring_f

LIBHDF5    = -L${HDF5_LIB} -lhdf5_fortran_nvidia -lhdf5_fortran

AUX_SRC    = aux_generic.f90

