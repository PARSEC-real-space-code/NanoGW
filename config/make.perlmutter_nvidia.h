##
#
# NERSC Perlmutter, Nvidia compiler, FFTW 3, MPI, CUDA
#
# module load nvhpc/24.5 cray-hdf5-parallel/1.12.2.9 cray-fftw/3.3.10.6
#
##

FCPP    = /usr/bin/cpp -P -traditional-cpp
CPPOPT  = -DUSEFFTW3 -DMPI -D_CUDA # -DDEBUG

EXT     = .nvidia.mpi

F90s    = ftn
F90     = $(F90s)

HDF5_DIR   = /opt/cray/pe/hdf5/1.12.2.9/nvidia/23.3
HDF5_INC   = ${HDF5_DIR}/include
HDF5_LIB   = ${HDF5_DIR}/lib
EXTRA_DIR  = ../LIB

OPTS       = -O3 -cuda -cudalib=cublas # -Mbounds -Mchkptr -g -traceback
OPTS2      = -I${HDF5_INC}

LIBXC      = -L${EXTRA_DIR}/libstring_f/lib -L${EXTRA_DIR}/libxc/lib -lxc -lstring_f

LIBHDF5    = -L${HDF5_LIB} -lhdf5_fortran_nvidia -lhdf5_fortran

AUX_SRC = aux_generic.f90

