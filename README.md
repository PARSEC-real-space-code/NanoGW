# Brief introduction to NanoGW

## Updates
We implemented a Lanczos-based method to accelerate full-frequency GW calculations. Currently, we have demonstrated the high efficiency of this method for finite systems. More details are discussed in [Physical Review Letters **132**, 126402](https://doi.org/10.1103/PhysRevLett.132.126402).

## Developers
The first version of NanoGW was developed by Murilo L. Tiago and James R. Chelikowsky between 2004 and 2009. Later, it was further developed and maintained by Linda Hung and Weiwei Gao. Weiwei Gao implemented a symmetry-adapted interpolative separable density fitting (ISDF) method to drastically speed up the evaluation of kernel matrix elements in the NanoGW calculations.  

NanoGW was originally named "RGWBS", which was less memorable, prompting us to change its name. Given its primary advantage for confined, nanoscale systems, we adopted the name "NanoGW".

## Features of NanoGW
This package can perform the following calculations:
1.  Linear-response time-dependent density functional theory (by solving the Casida equation)
2.  Full-frequency GW calculation with or without LDA vertex function (does not support spin-orbit coupling)
3.  Constructing and solving the Bethe-Salpeter equation
	
This package has been thoroughly tested and optimized for molecules and nanoclusters. It runs particularly efficiently for small-size (less than 30 atoms) molecules or clusters. This package can also handle crystalline systems. However, its functionality for dealing with crystals has not been thoroughly tested yet.

We have tested it on a few Linux/Unix machines, including NERSC Perlmutter, TACC Stampede3, the workstation of Oden Institute at UT Austin, and personal computers that run Ubuntu or Mac OS. It can now handle molecules (12 atoms or 50 electrons) within a few minutes and with reasonable accuracy on personal computers with a 2-core/4-thread Intel processor.  

## How to compile
Before proceeding with these steps, you may need to load the required modules as specified in the configuration files. Compile the extra libraries before building the main NanoGW code:  
```bash
cd LIB/libstring_f
./configure --prefix=$(pwd)
make
make install

cd ../libxc
./configure --prefix=$(pwd)
make
make install
cd ../..
```

Then you may compile the NanoGW code for the sigma calculation:  
```bash
make sigma MACH=ubuntu_intel
```

Replace `ubuntu_intel` with the appropriate machine-dependent configuration name from the [config](config) directory. You may need to adjust some parameters to fit your programming environment.  

## Input requirements
This package requires wave functions and Kohn-Sham energies calculated with [PARSEC](https://github.com/PARSEC-real-space-code/PARSEC) as input. It also supports the plane-wave-based DFT package PARATEC. NanoGW converts plane-wave-based wave functions to real-space-based wave functions. However, its compatibility with PARATEC has not been thoroughly tested yet.

## How to cite NanoGW
If you find NanoGW useful and use it for your publication, we would appreciate it if you could cite the following paper:  
Murilo L. Tiago and James R. Chelikowsky, *Optical excitations in organic molecules, clusters, and defects studied by first-principles Green's function methods*, [Physical Review B **73**, 205334](https://doi.org/10.1103/PhysRevB.73.205334) (2006).
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
Weiwei Gao and James R. Chelikowsky, *Accelerating time-dependent density functional theory and GW calculations for molecules and nanoclusters with symmetry adapted interpolative separable density fitting*, [Journal of Chemical Theory and Computation **16**, 2216](https://doi.org/10.1021/acs.jctc.9b01025) (2020).
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
If you use the Lanczos method to speed up your calculation, please also cite the following reference:  
Weiwei Gao, Zhao Tang, Jijun Zhao, and James R. Chelikowsky, *Efficient full-frequency GW calculations using a Lanczos method*, [Physical Review Letters **132**, 126402](https://doi.org/10.1103/PhysRevLett.132.126402) (2024).
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

## An incomplete list of papers using NanoGW
1.  Murilo L. Tiago and James R. Chelikowsky, *First-principles GW–BSE excitations in organic molecules*, [Solid State Communications **136**, 333](https://doi.org/10.1016/j.ssc.2005.08.012) (2005).
2.  Marie Lopez del Puerto, Murilo L. Tiago, and James R. Chelikowsky, *Excitonic effects and optical properties of passivated CdSe clusters*, [Physical Review Letters **97**, 096401](https://doi.org/10.1103/PhysRevLett.97.096401) (2006).
3.  Murilo L. Tiago and James R. Chelikowsky, *Optical excitations in organic molecules, clusters, and defects studied by first-principles Green's function methods*, [Physical Review B **73**, 205334](https://doi.org/10.1103/PhysRevB.73.205334) (2006).
4.  Marie Lopez del Puerto, Murilo L. Tiago, and James R. Chelikowsky, *Ab initio methods for the optical properties of CdSe clusters*, [Physical Review B **77**, 045404](https://doi.org/10.1103/PhysRevB.77.045404) (2008).
5.  Na Sai, Murilo L. Tiago, James R. Chelikowsky, and Fernando A. Reboredo, *Optical spectra and exchange-correlation effects in molecular crystals*, [Physical Review B **77**, 161306](https://doi.org/10.1103/PhysRevB.77.161306) (2008).
6.  Murilo L. Tiago, Paul R. C. Kent, Randolph Q. Hood, and Fernando A. Reboredo, *Neutral and charged excitations in carbon fullerenes from first-principles many-body theories*, [The Journal of Chemical Physics **129**, 084311](https://doi.org/10.1063/1.2973627) (2008).
7.  James R. Chelikowsky, Yousef Saad, Tzu-Liang Chan, Murilo L. Tiago, Alexey T. Zayak, and Yunkai Zhou, *Pseudopotentials on grids: application to the electronic, optical, and vibrational properties of silicon nanocrystals*, [Journal of Computational and Theoretical Nanoscience **6**, 1247](https://doi.org/10.1166/jctn.2009.1173) (2009).
8.  Kimberly Frey, Juan C. Idrobo, Murilo L. Tiago, Fernando Reboredo, and Serdar Ogut, *Quasiparticle gaps and exciton Coulomb energies in Si nanoshells: first-principles calculations*, [Physical Review B **80**, 153411](https://doi.org/10.1103/PhysRevB.80.153411) (2009).
9.  Murilo L. Tiago, Juan C. Idrobo, Serdar Ogut, Julius Jellinek, and James R. Chelikowsky, *Electronic and optical excitations in Ag<sub>n</sub> clusters (n=1-8): comparison of density-functional and many-body theories*, [Physical Review B **79**, 155419](https://doi.org/10.1103/PhysRevB.79.155419) (2009).
10. Murilo L. Tiago and Fernando A. Reboredo, *Controlling the gap of fullerene microcrystals by applying pressure: role of many-body effects*, [Physical Review B **79**, 195410](https://doi.org/10.1103/PhysRevB.79.195410) (2009).
11. Linda Hung, Kopinjol Baishya, and Serdar Ogut, *First-principles real-space study of electronic and optical excitations in rutile TiO<sub>2</sub> nanocrystals*, [Physical Review B **90**, 165424](https://doi.org/10.1103/PhysRevB.90.165424) (2014).
12. Linda Hung, Felipe H. da Jornada, Jaime Souto-Casares, James R. Chelikowsky, Steven G. Louie, and Serdar Ogut, *Excitation spectra of aromatic molecules within a real-space GW-BSE formalism: role of self-consistency and vertex corrections*, [Physical Review B **94**, 085125](https://doi.org/10.1103/PhysRevB.94.085125) (2016).
13. Linda Hung, Fabien Bruneval, Kopinjol Baishya, and Serdar Ogut, *Benchmarking the GW approximation and Bethe–Salpeter equation for groups IB and IIB atoms and monoxides*, [Journal of Chemical Theory and Computation **13**, 2135](https://doi.org/10.1021/acs.jctc.7b00123) (2017).
14. Weiwei Gao, Linda Hung, Serdar Ogut, and James R. Chelikowsky, *The stability, electronic structure, and optical absorption of boron-nitride diamondoids predicted with first-principles calculations*, [Physical Chemistry Chemical Physics **20**, 19188](https://doi.org/10.1039/C8CP02377H) (2018).
15. Weiwei Gao and James R. Chelikowsky, *Real-space based benchmark of G<sub>0</sub>W<sub>0</sub> calculations on GW100: effects of semicore orbitals and orbital reordering*, [Journal of Chemical Theory and Computation **15**, 5299](https://doi.org/10.1021/acs.jctc.9b00520) (2019).
16. Weiwei Gao and James R. Chelikowsky, *Accelerating time-dependent density functional theory and GW calculations for molecules and nanoclusters with symmetry adapted interpolative separable density fitting*, [Journal of Chemical Theory and Computation **16**, 2216](https://doi.org/10.1021/acs.jctc.9b01025) (2020).
17. Weiwei Gao, Zhao Tang, Jijun Zhao, and James R. Chelikowsky, *Efficient full-frequency GW calculations using a Lanczos method*, [Physical Review Letters **132**, 126402](https://doi.org/10.1103/PhysRevLett.132.126402) (2024).
