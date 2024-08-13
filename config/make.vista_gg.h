##
# TACC Vista gg nodes, Nvidia Fortran compiler, FFTW 3.*, MPI
#
# The Vista HDF5 module doesn`t have the Fortran version installed (???),
# so I installed the HDF5 with Fortran locally. You might want to modify the DEPS variable.
#
# module load fftw3 nvpl
#
# Currently Loaded Modules:
#  1) nvidia/24.7   2) ucx/1.17.0   3) openmpi/5.0.5   4) cmake/3.29.5   5) TACC   6) nvpl/24.7   7) fftw3/3.3.10
#
##

FCPP    = /usr/bin/cpp -P -traditional-cpp
CPPOPT  = -DUSEFFTW3 -DMPI # -D_CUDA -DDEBUG 

EXT     = .mpi

F90s    = mpif90
F90     = $(F90s)

EXTRA_DIR=../LIB

DEPS = /work/09183/zt3324/vista/gg/dependencies

OPTS       = -Ofast -mcpu=neoverse-v2 -g -traceback 
OPTS2      = -I${TACC_NVPL_INC} -I${TACC_FFTW3_INC} -I${DEPS}/include

LIBLAPACK  = -mp -L${TACC_NVPL_LIB} -lnvpl_blas_lp64_seq -lnvpl_lapack_lp64_seq -lnvpl_blacs_lp64_openmpi5 -lnvpl_scalapack_lp64
LIBFFT     = -L${TACC_FFTW3_LIB} -lfftw3
LIBHDF5    = -L${DEPS}/lib -lhdf5_fortran
LIBXC      = -L${EXTRA_DIR}/libstring_f/lib -L${EXTRA_DIR}/libxc/lib -lxc -lstring_f 

AUX_SRC = aux_generic.f90

