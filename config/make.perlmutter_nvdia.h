##
# Ubuntu Linux, gfortran Fortran compiler, FFTW 3.*, MPI
#
FCPP    = /usr/bin/cpp -P -traditional-cpp
CPPOPT  = -DUSEFFTW3 -DMPI -D_CUDA # -DDEBUG 

EXT     = .nvdia.mpi

F90s    = ftn
F90     = $(F90s)

#FFTW3_INC_DIR = 
#FFTW3_LIB_DIR = 
HDF5_INC_DIR  = /opt/cray/pe/hdf5/1.12.2.3/nvidia/20.7/include
HDF5_LIB_DIR  = /opt/cray/pe/hdf5/1.12.2.3/nvidia/20.7/lib
EXTRA_DIR     = /global/u1/w/weiwei/src/new-nanogw/nanogw-2023_8_23/LIB_nv
# MKLROOT       = /opt/intel/oneapi/mkl/2023.1.0

OPTS       = -O3 -cuda -cudalib=cublas # -Mbounds -Mchkptr -g -traceback 
OPTS2      = -I${HDF5_INC_DIR}  

# LIBLAPACK  =  -L${MKLROOT}/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lmkl_blacs_intelmpi_lp64 -lpthread -lm -ldl

LIBXC      = -L${EXTRA_DIR}/libstring_f/lib -L${EXTRA_DIR}/libxc/lib -lxc -lstring_f 

LIBHDF5    = -L${HDF5_LIB_DIR} -lhdf5_fortran_nvidia -lhdf5_fortran

AUX_SRC = aux_generic.f90

