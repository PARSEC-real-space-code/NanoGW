## Information about the current version:

W. Gao implemented a symmetry-adapted interpolative separable density fitting method (ISDF) to 
accelerate the calculations of kernel matrix elements, which has the form of `<ij|K|nm>`.
The time spent in calculating kernel matrix elements can be reduced by two ~ three orders of
magnitudes. 
And the total computational cost for GW calculations can be reduced by a few times.
The calculation of interpolation vectors using ISDF method cost a lot of extra memory space.
Here the HDF5 package is used to temporarily store the interpolation vectors.
This version has been tested for some small molecules and silicon nanoclusters.

Things that is new in this version:
1. The ISDF method is only implemented for isolated systems for now.
The ISDF method seems to provide good accuracy and improvement of efficiency.
2. The parallelization of representation group seems to work correctly.

Issues to be improved:
1. The diagonalization (of Casida equation) step needs to redistribute data. 
This costs a lot of communication time, especially for large systems. We need to rethink the data layout.
2. Is it possible to avoid diagonalizing the whole Casida eq? Can we use approximate solutions for Casida eq?
3. Need more tests for periodic systems.
4. This version uses an old `libxc` library. Consider to update later.


