##
# Ubuntu Linux, gfortran Fortran compiler, FFTW 3.*, MPI
#
FCPP    = /usr/bin/cpp -P -traditional-cpp
CPPOPT  = -DUSEFFTW3 -DMPI -DDEBUG 

EXT     = .nvcc_mkl.mpi

F90s    = nvfortran # Serial compiler
F90     = nvfortran -I/opt/intel/oneapi/mpi/latest/include

#FFTW3_INC_DIR = /home/weiwei/Programs/FFTW3/fftw-3.3.7-nv/LIBS/include
#FFTW3_LIB_DIR = /home/weiwei/Programs/FFTW3/fftw-3.3.7-nv/LIBS/lib
HDF5_INC_DIR  = /home/weiwei/Programs/hdf5/hdf5-1.8.16-nv/lib/include 
HDF5_LIB_DIR  = /home/weiwei/Programs/hdf5/hdf5-1.8.16-nv/lib/lib
EXTRA_DIR     = /home/weiwei/Programs/nanogw/LIB
#LAPACK_DIR    = /opt/nvidia/hpc_sdk/Linux_x86_64/22.2/compilers
#SCALAPACK_DIR = /opt/nvidia/hpc_sdk/Linux_x86_64/22.2/comm_libs/openmpi/openmpi-3.1.5
MKLROOT=/opt/intel/oneapi/mkl/2022.0.2


OPTS    = -O2 -fast
OPTS2   = -I${HDF5_INC_DIR}  -I${MKLROOT}/include

LIBFFT  = 

# LIBLAPACK = ${MKLROOT}/lib/intel64/libmkl_scalapack_lp64.a -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_blacs_openmpi_lp64.a -Wl,--end-group -lpthread -lm -ldl
# LIBLAPACK =  ${MKLROOT}/lib/intel64/libmkl_scalapack_lp64.a -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_blacs_openmpi_lp64.a -Wl,--end-group -lpthread -lm -ldl
LIBLAPACK = -L/opt/intel/oneapi/mpi/latest/lib/release -L/opt/intel/oneapi/mpi/latest/lib -L${MKLROOT}/lib/intel64 -lmpi -lmpifort -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lmkl_blacs_intelmpi_lp64 -lpthread -lm -ldl

LIBXC = -L${EXTRA_DIR}/libstring_f_nv -L${EXTRA_DIR}/libxc_nv/src -lxc -lstring_f 

LIBHDF5 = -L${HDF5_LIB_DIR} -lhdf5_fortran

AUX_SRC = aux_generic.f90

