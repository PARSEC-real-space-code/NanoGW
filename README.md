## Brief intro of NanoGW

2020/05/20 Weiwei Gao

### Who are the developers
The first version of NanoGW was developed by Murilo Tiago and James R. Chelikowsky between 2004 and 2009. Later it is developed and maintained by Linda Hung and Weiwei Gao. Weiwei Gao implemented a symmetry-adapted interpolative separable density fitting method to drastically speed up the evaluation of kernel matrix elements in the calculations performed with NanoGW. 
NanoGW was originally named as RGWBS, which is more difficult to remember (so we changed the name to “NanoGW”). Since the main advantage of this package is for confined, nanoscale systems, so we call it NanoGW.

### What can NanoGW do:
This package can perform the following calculations:
1.	Linear-response time-dependent density functional theory (by solving Casida equation)
2.	Full-frequency GW calculation with or without LDA vertex function (does not support spin-orbit coupling)
3.	Construct and solve Bethe-salpeter equation
This package has been tested thoroughly and optimized for molecules and nanoclusters. It runs particularly efficient for small-size (less than 30 atoms) molecules or clusters. This package can also deal with crystalline systems. However, its functionality for dealing crystals has not been thoroughly tested yet. 
We have tested it one a few Linux/Unix machines, including NERSC cori, TACC stampede2, Workstation of Oden institute at UT Austin and personal computers that run Ubuntu or Mac OS. Now it can handles molecules (12 atoms or 50 electrons) with a few minutes and reasonable accuracy on personal computers with a 2-core/4-thread Intel processor. 

### What does NanoGW need as input:
This package requires wave functions and Kohn-Sham energies calculated with PARSEC as input. It also support plane-wave based DFT package PARATEC. NanoGW will convert the plane-wave based wave functions to real-space based wave functions. However, its compatibility with PARATEC has not been thoroughly tested yet. 

### How to cite this software:
If you find NanoGW is useful and use it for your publication, we appreciate if you can cite the following paper:
```
@article{PhysRevB.73.205334,
  title = {Optical excitations in organic molecules, clusters, and defects studied by first-principles Green's function methods},
  author = {Tiago, Murilo L. and Chelikowsky, James R.},
  journal = {Phys. Rev. B},
  volume = {73},
  issue = {20},
  pages = {205334},
  numpages = {19},
  year = {2006},
  month = {May},
  publisher = {American Physical Society},
  doi = {10.1103/PhysRevB.73.205334},
  url = {https://link.aps.org/doi/10.1103/PhysRevB.73.205334}
}
```
If you uses interpolative separable density fitting (ISDF) to speed up your calculation, please also cite the following reference:
```
@article{doi:10.1021/acs.jctc.9b01025,
author = {Gao, Weiwei and Chelikowsky, James R.},
title = {Accelerating Time-Dependent Density Functional Theory and GW Calculations for Molecules and Nanoclusters with Symmetry Adapted Interpolative Separable Density Fitting},
journal = {Journal of Chemical Theory and Computation},
volume = {16},
number = {4},
pages = {2216-2223},
year = {2020},
doi = {10.1021/acs.jctc.9b01025},
note ={PMID: 32074452},
URL = {https://doi.org/10.1021/acs.jctc.9b01025}
}
```

An incomplete list of papers that use NanoGW:
1.	Tiago, M. L. & Chelikowsky, J. R. First-principles GW–BSE excitations in organic molecules. Solid State Communications 136, 333-337, doi:https://doi.org/10.1016/j.ssc.2005.08.012 (2005).
2.	del Puerto, M. L., Tiago, M. L. & Chelikowsky, J. R. Excitonic Effects and Optical Properties of Passivated CdSe Clusters. Physical Review Letters 97, 096401, doi:10.1103/PhysRevLett.97.096401 (2006).
3.	Tiago, M. L. & Chelikowsky, J. R. Optical excitations in organic molecules, clusters, and defects studied by first-principles Green's function methods. Physical Review B 73, 205334, doi:10.1103/PhysRevB.73.205334 (2006).
4.	Lopez del Puerto, M., Tiago, M. L. & Chelikowsky, J. R. Ab initio methods for the optical properties of CdSe clusters. Physical Review B 77, 045404, doi:10.1103/PhysRevB.77.045404 (2008).
5.	Sai, N., Tiago, M. L., Chelikowsky, J. R. & Reboredo, F. A. Optical spectra and exchange-correlation effects in molecular crystals. Physical Review B 77, 161306, doi:10.1103/PhysRevB.77.161306 (2008).
6.	Tiago, M. L., Kent, P. R. C., Hood, R. Q. & Reboredo, F. A. Neutral and charged excitations in carbon fullerenes from first-principles many-body theories. The Journal of Chemical Physics 129, 084311, doi:10.1063/1.2973627 (2008).
7.	Chelikowsky, J. R. et al. Pseudopotentials on Grids: Application to the Electronic, Optical, and Vibrational Properties of Silicon Nanocrystals. Journal of Computational and Theoretical Nanoscience 6, 1247-1261, doi:10.1166/jctn.2009.1173 (2009).
8.	Frey, K., Idrobo, J. C., Tiago, M. L., Reboredo, F. & Öğüt, S. Quasiparticle gaps and exciton Coulomb energies in Si nanoshells: First-principles calculations. Physical Review B 80, 153411, doi:10.1103/PhysRevB.80.153411 (2009).
9.	Tiago, M. L., Idrobo, J. C., Öğüt, S., Jellinek, J. & Chelikowsky, J. R. Electronic and optical excitations in Agn clusters (n=1-8): Comparison of density-functional and many-body theories. Physical Review B 79, 155419, doi:10.1103/PhysRevB.79.155419 (2009).
10.	Tiago, M. L. & Reboredo, F. A. Controlling the gap of fullerene microcrystals by applying pressure: Role of many-body effects. Physical Review B 79, 195410, doi:10.1103/PhysRevB.79.195410 (2009).
11.	Hung, L., Baishya, K. & Öğüt, S. First-principles real-space study of electronic and optical excitations in rutile TiO2 nanocrystals. Physical Review B 90, 165424, doi:10.1103/PhysRevB.90.165424 (2014).
12.	Hung, L. et al. Excitation spectra of aromatic molecules within a real-space GW-BSE formalism: Role of self-consistency and vertex corrections. Physical Review B 94, 085125, doi:10.1103/PhysRevB.94.085125 (2016).
13.	Hung, L., Bruneval, F., Baishya, K. & Öğüt, S. Benchmarking the GW Approximation and Bethe–Salpeter Equation for Groups IB and IIB Atoms and Monoxides. Journal of Chemical Theory and Computation 13, 2135-2146, doi:10.1021/acs.jctc.7b00123 (2017).
14.	Gao, W., Hung, L., Ogut, S. & Chelikowsky, J. R. The stability, electronic structure, and optical absorption of boron-nitride diamondoids predicted with first-principles calculations. Physical Chemistry Chemical Physics 20, 19188-19194, doi:10.1039/C8CP02377H (2018).
15.	Gao, W. & Chelikowsky, J. R. Real-Space Based Benchmark of G0W0 Calculations on GW100: Effects of Semicore Orbitals and Orbital Reordering. Journal of Chemical Theory and Computation 15, 5299-5307, doi:10.1021/acs.jctc.9b00520 (2019).
16.	Gao, W. & Chelikowsky, J. R. Accelerating Time-Dependent Density Functional Theory and GW Calculations for Molecules and Nanoclusters with Symmetry Adapted Interpolative Separable Density Fitting. Journal of Chemical Theory and Computation 16, 2216-2223, doi:10.1021/acs.jctc.9b01025 (2020).
