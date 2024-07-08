## Brief intro of NanoGW

7/8/2024

### Updates
We implemented a Lanczos-based method to accelerate full-frequency GW calculations (named as LanczosGW). For now, we have demonstrate the great efficiency of LanczosGW method for finite systems. More details is discussed in [Physical Review Letters **132**, 126402](https://doi.org/10.1103/PhysRevLett.132.126402).

### developers
The first version of NanoGW was developed by Murilo Tiago and James R. Chelikowsky between 2004 and 2009. Later, it is developed and maintained by Linda Hung and Weiwei Gao. Weiwei Gao implemented a symmetry-adapted interpolative separable density fitting method to drastically speed up the evaluation of kernel matrix elements in the calculations performed with NanoGW.  
NanoGW was originally named as RGWBS, which is more difficult to remember (so we changed the name to “NanoGW”). Since the main advantage of this package is for confined, nanoscale systems, so we call it NanoGW.

### What can NanoGW do:
This package can perform the following calculations:
1.	Linear-response time-dependent density functional theory (by solving Casida equation)
2.	Full-frequency GW calculation with or without LDA vertex function (does not support spin-orbit coupling)
3.	Construct and solve Bethe-Salpeter equation
This package has been tested thoroughly and optimized for molecules and nanoclusters. It runs particularly efficient for small-size (less than 30 atoms) molecules or clusters. This package can also deal with crystalline systems. However, its functionality for dealing crystals has not been thoroughly tested yet. 
We have tested it one a few Linux/Unix machines, including NERSC perlmutter, TACC stampede3, Workstation of Oden institute at UT Austin and personal computers that run Ubuntu or Mac OS. Now it can handles molecules (12 atoms or 50 electrons) with a few minutes and reasonable accuracy on personal computers with a 2-core/4-thread Intel processor. 

### What does NanoGW need as input:
This package requires wave functions and Kohn-Sham energies calculated with PARSEC as input. It also support plane-wave based DFT package PARATEC. NanoGW will convert the plane-wave based wave functions to real-space based wave functions. However, its compatibility with PARATEC has not been thoroughly tested yet. 

### How to cite this software:
If you find NanoGW is useful and use it for your publication, we appreciate if you can cite the following paper:
```
@article{tiago2006optical,
  title={Optical excitations in organic molecules, clusters, and defects studied by first-principles Green’s function methods},
  author={Tiago, Murilo L and Chelikowsky, James R},
  journal={Physical Review B—Condensed Matter and Materials Physics},
  volume={73},
  number={20},
  pages={205334},
  year={2006},
  publisher={APS}
}
```
If you use interpolative separable density fitting (ISDF) to speed up your calculation, please also cite the following reference:
```
@article{gao2020accelerating,
  title={Accelerating time-dependent density functional theory and GW calculations for molecules and nanoclusters with symmetry adapted interpolative separable density fitting},
  author={Gao, Weiwei and Chelikowsky, James R},
  journal={Journal of Chemical Theory and Computation},
  volume={16},
  number={4},
  pages={2216--2223},
  year={2020},
  publisher={ACS Publications}
}
```
If you use LanczosGW to speed up your calculation, please also cite the following reference:
```
@article{gao2024efficient,
  title={Efficient Full-Frequency GW Calculations Using a Lanczos Method},
  author={Gao, Weiwei and Tang, Zhao and Zhao, Jijun and Chelikowsky, James R},
  journal={Physical Review Letters},
  volume={132},
  number={12},
  pages={126402},
  year={2024},
  publisher={APS}
}
```

### An incomplete list of papers that use NanoGW:
1.  Murilo L. Tiago and James R. Chelikowsky, *First-principles GW–BSE excitations in organic molecules*, [Solid State Communications **136**, 333](https://doi.org/10.1016/j.ssc.2005.08.012) (2005).
2.  Marie Lopez del Puerto, Murilo L. Tiago, and James R. Chelikowsky, *Excitonic effects and optical properties of passivated CdSe clusters*, [Physical Review Letters **97**, 096401](https://doi.org/10.1103/PhysRevLett.97.096401) (2006).
3.  Murilo L. Tiago and James R. Chelikowsky, *Optical excitations in organic molecules, clusters, and defects studied by first-principles Green's function methods*, [Physical Review B **73**, 205334](https://doi.org/10.1103/PhysRevB.73.205334) (2006).
4.  Marie Lopez del Puerto, Murilo L. Tiago, and James R. Chelikowsky, *Ab initio methods for the optical properties of CdSe clusters*, [Physical Review B **77**, 045404](https://doi.org/10.1103/PhysRevB.77.045404) (2008).
5.  Na Sai, Murilo L. Tiago, James R. Chelikowsky, and Fernando A. Reboredo, *Optical spectra and exchange-correlation effects in molecular crystals*, [Physical Review B **77**, 161306](https://doi.org/10.1103/PhysRevB.77.161306) (2008).
6.  Murilo L. Tiago, Paul R. C. Kent, Randolph Q. Hood, Fernando A. Reboredo, *Neutral and charged excitations in carbon fullerenes from first-principles many-body theories*, [The Journal of Chemical Physics **129**, 084311](https://doi.org/10.1063/1.2973627) (2008).
7.  James R. Chelikowsky, Yousef Saad, Tzu-Liang Chan, Murilo L. Tiago, Alexey T. Zayak, Yunkai Zhou, *Pseudopotentials on grids: application to the electronic, optical, and vibrational properties of silicon nanocrystals*, [Journal of Computational and Theoretical Nanoscience **6**, 1247](https://doi.org/10.1166/jctn.2009.1173) (2009).
8.  Kimberly Frey, Juan C. Idrobo, Murilo L. Tiago, Fernando Reboredo, and Serdar Öğüt, *Quasiparticle gaps and exciton Coulomb energies in Si nanoshells: first-principles calculations*, [Physical Review B **80**, 153411](https://doi.org/10.1103/PhysRevB.80.153411) (2009).
9.  Murilo L. Tiago, Juan C. Idrobo, Serdar Öğüt, Julius Jellinek, and James R. Chelikowsky, *Electronic and optical excitations in Ag<sub>n</sub> clusters (n=1-8): comparison of density-functional and many-body theories*, [Physical Review B **79**, 155419](https://doi.org/10.1103/PhysRevB.79.155419) (2009).
10. Murilo L. Tiago and Fernando A. Reboredo, *Controlling the gap of fullerene microcrystals by applying pressure: role of many-body effects*, [Physical Review B **79**, 195410](https://doi.org/10.1103/PhysRevB.79.195410) (2009).
11. Linda Hung, Kopinjol Baishya, and Serdar Öğüt, *First-principles real-space study of electronic and optical excitations in rutile TiO~2~ nanocrystals*, [Physical Review B **90**, 165424](https://doi.org/10.1103/PhysRevB.90.165424) (2014).
12. Linda Hung, Felipe H. da Jornada, Jaime Souto-Casares, James R. Chelikowsky, Steven G. Louie, and Serdar Öğüt, *Excitation spectra of aromatic molecules within a real-space GW-BSE formalism: role of self-consistency and vertex corrections*, [Physical Review B **94**, 085125](https://doi.org/10.1103/PhysRevB.94.085125) (2016).
13. Linda Hung, Fabien Bruneval, Kopinjol Baishya, and Serdar Öğüt *Benchmarking the GW approximation and Bethe–Salpeter equation for groups IB and IIB atoms and monoxides*, [Journal of Chemical Theory and Computation **13**, 2135](https://doi.org/10.1021/acs.jctc.7b00123) (2017).
14. Weiwei Gao, Linda Hung, Serdar Öğüt, and James R. Chelikowsky, *The stability, electronic structure, and optical absorption of boron-nitride diamondoids predicted with first-principles calculations*, [Physical Chemistry Chemical Physics **20**, 19188](https://doi.org/10.1039/C8CP02377H) (2018).
15. Weiwei Gao and James R. Chelikowsky, *Real-space based benchmark of G<sub>0</sub>W<sub>0</sub> calculations on GW100: effects of semicore orbitals and orbital reordering*, [Journal of Chemical Theory and Computation **15**, 5299](https://doi.org/10.1021/acs.jctc.9b00520) (2019).
16. Weiwei Gao and James R. Chelikowsky, *Accelerating time-dependent density functional theory and GW calculations for molecules and nanoclusters with symmetry adapted interpolative separable density fitting*, [Journal of Chemical Theory and Computation **16**, 2216](https://doi.org/10.1021/acs.jctc.9b01025) (2020).
17. Weiwei Gao, Zhao Tang, Jijun Zhao, and James R. Chelikowsky, *Efficient full-frequency GW calculations using a Lanczos method*, [Physical Review Letters **132**, 126402](https://doi.org/10.1103/PhysRevLett.132.126402) (2024).
