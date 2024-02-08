## Instruction for compiling the code on a personal computer running linux OS

Steps to compile this package on a PC with "Debian" based linux system.

1. Install libraries and GNU compilers

We need FFTW, lapack, blas, scalapack, and HDF5 libraries.
We also need serial fortran compiler (gfortran) and mpi compiler.
If you already have these packages, then skip the step 1 to step 2.

Run the following commands to install these softwares:

```bash
sudo apt-get install libblas-dev
sudo apt-get install liblapack-dev
sudo apt-get install libscalapack-mpi-dev
sudo apt-get install libfftw3-dev
sudo apt-get install libfftw3-mpi-dev
sudo apt-get install libhdf5-dev
sudo apt-get install gfortran libmpich-dev
```

If some libs or compilers are not found when installing, maybe try "sudo apt-get update" to update the database.

2. edit config files under the "config/" directory
For example, in config/make.linux.h, you need to define the following variables:

```bash
FFTW3_INC_DIR 
FFTW3_LIB_DIR 
HDF5_INC_DIR  
HDF5_LIB_DIR  
EXTRA_DIR     
LAPACK_DIR    
``` 
These are the locations where libraries are installed.
To find the location of the installed libraries, try:

```bash
dpkg -L [library_name]
```

3. compile libraries under `LIB/`

These libraries in LIB/ are required for calculating the exchange-correlation potentials.
Go to `LIB/libstring\_f/` and `LIB/libxc/`
Type the following commands the compile these two libraries.
```
./configure
make
```

4. Type the following commands to compile tdlda and sigma
```
make tdlda
make sigma
```
See Makefile for more descriptions.

