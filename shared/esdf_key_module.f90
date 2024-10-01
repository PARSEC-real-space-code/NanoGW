!===============================================================
!
! Module to hold keyword list. this must be updated as
! new keywords are brought into existence.
!
! The 'label' is the label as used in calling the esdf routines
! 'typ' defines the type, with the following syntax. it is 3
! characters long.
! the first indicates:
!      i - integer
!      s - single
!      d - double
!      p - physical
!      t - string (text)
!      e - defined (exists)
!      l - boolean (logical)
!      b - block
! the second is always a colon (:)
! the third indicates the "level" of the keyword
!      b - basic
!      i - intermediate
!      e - expert
!      d - dummy
!
! 'Dscrpt' is a description of the variable. it should contain a
! (short) title enclosed between *! ... !*, and then a more detailed
! description of the variable.
!
! Copyright (C) 2009 Murilo L. Tiago, http://users.ices.utexas.edu/~mtiago
! This file is part of RGWBS. It is distributed under the GPL v1.
!
!---------------------------------------------------------------
module esdf_key

  implicit none
  ! maximum number of kewords
  integer, parameter :: NUMKW = 600

  character(len=80)   :: kw_label(NUMKW)
  character(len=3)    :: kw_typ(NUMKW)
  character(len=3000) :: kw_dscrpt(NUMKW)

  integer :: kw_index(NUMKW)

  ! now define the keywords

  data kw_label(1)/'tdlda_cutoff'/
  data kw_typ(1)/'P:E'/
  data kw_dscrpt(1)/'*! Energy cutoff in TDLDA !*'/

  data kw_label(2)/'buffer_size'/
  data kw_typ(2)/'P:E'/
  data kw_dscrpt(2)/'*! Buffer size in the output of kernel matrix elements !*'/

  data kw_label(3)/'cache_size'/
  data kw_typ(3)/'I:E'/
  data kw_dscrpt(3)/'*! Cache size in evaluation of kernel matrix elements !*'/

  data kw_label(4)/'qp_mixing_param'/
  data kw_typ(4)/'D:E'/
  data kw_dscrpt(4)/'*! Amount of mixing between old and new QP self-energy!*'/

  data kw_label(5)/'correlation_type'/
  data kw_typ(5)/'T:D'/
  data kw_dscrpt(5)/'*! Correlation type !*'/

  data kw_label(6)/'dft_program'/
  data kw_typ(6)/'T:E'/
  data kw_dscrpt(6)/'*! name of DFT program !*'/

  data kw_label(7)/'no_lda_kernel'/
  data kw_typ(7)/'E:E'/
  data kw_dscrpt(7)/'*! TDLDA kernel !*'/

  data kw_label(8)/'tdlda_triplet_kernel'/
  data kw_typ(8)/'E:E'/
  data kw_dscrpt(8)/'*! TDLDA triplet kernel !*'/

  data kw_label(9)/'no_exchange'/
  data kw_typ(9)/'E:E'/
  data kw_dscrpt(9)/'*! Exchange kernel !*'/

  data kw_label(10)/'tamm_dancoff'/
  data kw_typ(10)/'E:E'/
  data kw_dscrpt(10)/'*! Use Tamm-Dancoff approximation !*'/

  data kw_label(11)/'truncate_coulomb'/
  data kw_typ(11)/'E:E'/
  data kw_dscrpt(11)/'*! Truncation flag for Coulomb interaction !*'/

  data kw_label(12)/'tdlda_valence'/
  data kw_typ(12)/'B:E'/
  data kw_dscrpt(12)/'*! Valence(occupied) levels !*'/

  data kw_label(13)/'tdlda_valence_up'/
  data kw_typ(13)/'B:E'/
  data kw_dscrpt(13)/'*! Valence(occupied) levels, spin up !*'/

  data kw_label(14)/'tdlda_valence_down'/
  data kw_typ(14)/'B:E'/
  data kw_dscrpt(14)/'*! Valence(occupied) levels, spin down !*'/

  data kw_label(15)/'tdlda_conduction'/
  data kw_typ(15)/'B:E'/
  data kw_dscrpt(15)/'*! Conduction(unoccupied) levels !*'/

  data kw_label(16)/'tdlda_conduction_up'/
  data kw_typ(16)/'B:E'/
  data kw_dscrpt(16)/'*! Conduction(unoccupied) levels, spin up !*'/

  data kw_label(17)/'tdlda_conduction_down'/
  data kw_typ(17)/'B:E'/
  data kw_dscrpt(17)/'*! Conduction(unoccupied) levels, spin down !*'/

  data kw_label(18)/'bse_valence'/
  data kw_typ(18)/'B:E'/
  data kw_dscrpt(18)/'*! Valence(occupied) levels !*'/

  data kw_label(19)/'bse_valence_up'/
  data kw_typ(19)/'B:E'/
  data kw_dscrpt(19)/'*! Valence(occupied) levels, spin up !*'/

  data kw_label(20)/'bse_valence_down'/
  data kw_typ(20)/'B:E'/
  data kw_dscrpt(20)/'*! Valence(occupied) levels, spin down !*'/

  data kw_label(21)/'bse_conduction'/
  data kw_typ(21)/'B:E'/
  data kw_dscrpt(21)/'*! Conduction(unoccupied) levels !*'/

  data kw_label(22)/'bse_conduction_up'/
  data kw_typ(22)/'B:E'/
  data kw_dscrpt(22)/'*! Conduction(unoccupied) levels, spin up !*'/

  data kw_label(23)/'bse_conduction_down'/
  data kw_typ(23)/'B:E'/
  data kw_dscrpt(23)/'*! Conduction(unoccupied) levels, spin down !*'/

  data kw_label(24)/'energy_data_points'/
  data kw_typ(24)/'I:E'/
  data kw_dscrpt(24)/'*! Number of energy values where self-energy is calculated !*'/

  data kw_label(25)/'scratch_disk_size'/
  data kw_typ(25)/'P:E'/
  data kw_dscrpt(25)/'*! Amount of disk space used for scratch !*'/

  data kw_label(26)/'max_number_states'/
  data kw_typ(26)/'I:E'/
  data kw_dscrpt(26)/'*! Maximum number of states in Green function  !*'/

  data kw_label(27)/'energy_reference'/
  data kw_typ(27)/'P:E'/
  data kw_dscrpt(27)/'*! Energy reference in BSE !*'/

  data kw_label(28)/'renormalize_sumrule'/
  data kw_typ(28)/'E:E'/
  data kw_dscrpt(28)/'*! Renormalization of sum rule !*'/

  data kw_label(29)/'energy_range'/
  data kw_typ(29)/'P:E'/
  data kw_dscrpt(29)/'*! Energy range where self-energy is calculated !*'/

  data kw_label(30)/'diag'/
  data kw_typ(30)/'B:E'/
  data kw_dscrpt(30)/'*! Diagonal matrix elements in self-energy !*'/

  data kw_label(31)/'offdiag'/
  data kw_typ(31)/'B:E'/
  data kw_dscrpt(31)/'*! Off-diagonal matrix elements in self-energy !*'/

  data kw_label(32)/'bse_cutoff'/
  data kw_typ(32)/'P:E'/
  data kw_dscrpt(32)/'*! Energy cutoff in BSE !*'/

  data kw_label(33)/'rpa_spectrum_only'/
  data kw_typ(33)/'E:E'/
  data kw_dscrpt(33)/'*! Write out file eigenvalues_rpa only !*'/

  data kw_label(34)/'no_eigenvectors'/
  data kw_typ(34)/'E:E'/
  data kw_dscrpt(34)/'*! Do not write BSE eigenvectors !*'/

  data kw_label(35)/'bse_triplet_kernel'/
  data kw_typ(35)/'E:E'/
  data kw_dscrpt(35)/'*! BSE triplet kernel !*'/

  data kw_label(36)/'use_mixing'/
  data kw_typ(36)/'E:E'/
  data kw_dscrpt(36)/'*! Use mixing in BSE !*'/

  data kw_label(37)/'write_qp_wavefunctions'/
  data kw_typ(37)/'E:E'/
  data kw_dscrpt(37)/'*! Write QP wavefunctions in file parsec_qp.dat !*'/

  data kw_label(38)/'read_vxc_matrix_elements'/
  data kw_typ(38)/'E:E'/
  data kw_dscrpt(38)/'*! Read VXC matrix elements from sigma_mtxel.dat !*'/

  data kw_label(39)/'distribute_representations'/
  data kw_typ(39)/'I:E'/
  data kw_dscrpt(39)/'*! Distribute representations among MPI groups !*'/

  data kw_label(40)/'dynamic_energy_resolution'/
  data kw_typ(40)/'P:E'/
  data kw_dscrpt(40)/'*! Resolution in energy denominators, dynamical self-energy !*'/

  data kw_label(41)/'cohsex_approximation'/
  data kw_typ(41)/'E:E'/
  data kw_dscrpt(41)/'*! Calculate Sigma under the COHSEX (static) approximation !*'/

  data kw_label(42)/'read_orbital_occupancies'/
  data kw_typ(42)/'E:E'/
  data kw_dscrpt(42)/'*! Override level occupancies in parsec.dat with data in file occup.in !*'/

  data kw_label(43)/'scissors'/
  data kw_typ(43)/'B:E'/
  data kw_dscrpt(43)/'*! Scissors operator !*'/

  data kw_label(44)/'distribute_wavefunctions'/
  data kw_typ(44)/'I:E'/
  data kw_dscrpt(44)/'*! Distribute wavefunctions among MPI groups !*'/

  data kw_label(45)/'bse_energy_resolution'/
  data kw_typ(45)/'P:E'/
  data kw_dscrpt(45)/'*! Resolution in energy denominators, BSE kernel !*'/

  data kw_label(46)/'exchange_correlation'/
  data kw_typ(46)/'T:E'/
  data kw_dscrpt(46)/'*! Definition of the exchange-correlation model!*'/

  data kw_label(47)/'max_number_states_cohsex'/
  data kw_typ(47)/'I:E'/
  data kw_dscrpt(47)/'*! Maximum number of states used in the COHSEX !*'/

  data kw_label(48)/'self_energy_cutoff'/
  data kw_typ(48)/'P:E'/
  data kw_dscrpt(48)/'*! Cut-off in self-energy matrix elements !*'/

  data kw_label(49)/'number_iterations'/
  data kw_typ(49)/'I:E'/
  data kw_dscrpt(49)/'*! Number of self-consistent iterations !*'/

  data kw_label(50)/'convergence_potential'/
  data kw_typ(50)/'P:E'/
  data kw_dscrpt(50)/'*! Maximum component of perturbation potential !*'/

  data kw_label(51)/'read_checkpoint'/
  data kw_typ(51)/'E:E'/
  data kw_dscrpt(51)/'*! Read checkpoint files if they exist !*'/

  data kw_label(52)/'atomic_numbers'/
  data kw_typ(52)/'B:E'/
  data kw_dscrpt(52)/'*! Block of atomic numbers !*'/

  data kw_label(53)/'number_localized_orbitals'/
  data kw_typ(53)/'I:E'/
  data kw_dscrpt(53)/'*! Number of localized orbitals for BLIP transformation !*'/

  data kw_label(54)/'no_hqp_symmetrize'/
  data kw_typ(54)/'E:E'/
  data kw_dscrpt(54)/'*! do not symmetrize QP Hamiltonian !*'/

  data kw_label(55)/'point_group_tables'/
  data kw_typ(55)/'B:E'/
  data kw_dscrpt(55)/'*! Names of files with point group tables !*'/

  data kw_label(56)/'qpoints'/
  data kw_typ(56)/'B:E'/
  data kw_dscrpt(56)/'*! List of q-points in polarizability !*'/

  data kw_label(57)/'sigma_kpoints'/
  data kw_typ(57)/'B:E'/
  data kw_dscrpt(57)/'*! List of k-points in self-energy !*'/

  data kw_label(58)/'qpoints_bse'/
  data kw_typ(58)/'B:E'/
  data kw_dscrpt(58)/'*! List of q-points in BSE polarizability !*'/

  data kw_label(59)/'sigma_energy'/
  data kw_typ(59)/'T:E'/
  data kw_dscrpt(59)/'*! energy used to calculate dynamic sigma !*'/

  data kw_label(60)/'bwfn_filename'/
  data kw_typ(60)/'T:D'/
  data kw_dscrpt(60)/'*! Bwfn file name !*'/

  data kw_label(61)/'fft_expansion_factor'/
  data kw_typ(61)/'D:E'/
  data kw_dscrpt(61)/'*! Factor used to expand the FFT grid !*'/

  data kw_label(62)/'maximum_orbital'/
  data kw_typ(62)/'I:E'/
  data kw_dscrpt(62)/'*! Maximum orbital to be writen in BLIP file !*'/

  data kw_label(63)/'static_type'/
  data kw_typ(63)/'I:E'/
  data kw_dscrpt(63)/'*! Type of static correction for correlation self energy !*'/

  data kw_label(65)/'only_diagonal'/
  data kw_typ(65)/'E:E'/
  data kw_dscrpt(65)/'*! Calculate only diagonal of SIGMA matrix (no off-diagonal) !*'/

  data kw_label(66)/'qp_min'/
  data kw_typ(66)/'I:E'/
  data kw_dscrpt(66)/'*! Lowest index quasiparticle in SIGMA matrix !*'/

  data kw_label(67)/'qp_max'/
  data kw_typ(67)/'I:E'/
  data kw_dscrpt(67)/'*! Highest index quasiparticle in SIGMA matrix !*'/

  data kw_label(68)/'qp_estimate'/
  data kw_typ(68)/'B:E'/
  data kw_dscrpt(68)/'*! Estimated qp energy !*'/

  data kw_label(69)/'read_im_sigma'/
  data kw_typ(69)/'E:E'/
  data kw_dscrpt(69)/'*! Read imaginary part of sigma from sigma_mtxel.dat !*'/

  data kw_label(201)/'periodic_system'/
  data kw_typ(201)/'L:B'/
  data kw_dscrpt(201)/'*!System type (confined or periodic)!*'/

  data kw_label(202)/'expansion_order'/
  data kw_typ(202)/'I:B'/
  data kw_dscrpt(202)/'*!finite difference expansion order!*'/

  data kw_label(203)/'atom_types_num'/
  data kw_typ(203)/'I:E'/
  data kw_dscrpt(203)/'*! Number of atom types !*'/

  data kw_label(204)/'atom_type'/
  data kw_typ(204)/'T:E'/
  data kw_dscrpt(204)/'*! Type of the atom !*'/

  data kw_label(208)/'core_cutoff_radius'/
  data kw_typ(208)/'P:E'/
  data kw_dscrpt(208)/'*!pseudopotential core cutoff radius!*'/

  data kw_label(209)/'potential_num'/
  data kw_typ(209)/'I:E'/
  data kw_dscrpt(209)/'*! Number of potentials !*'/

  data kw_label(210)/'local_component'/
  data kw_typ(210)/'T:B'/
  data kw_dscrpt(210)/'*! The local component (s,p or d) !*'/

  data kw_label(211)/'lattice_vector_scale'/
  data kw_typ(211)/'P:D'/
  data kw_dscrpt(211)/'*! Unit for lattice vectors !*'/

  data kw_label(212)/'atom_coord'/
  data kw_typ(212)/'B:E'/
  data kw_dscrpt(212)/'*! coordinates of atom!*'/

  data kw_label(213)/'coordinate_unit'/
  data kw_typ(213)/'T:D'/
  data kw_dscrpt(213)/'*! Unit for atomic coordinates !*'/

  data kw_label(214)/'cell_shape'/
  data kw_typ(214)/'B:E'/
  data kw_dscrpt(214)/'*! length of sides of periodic box !*'/

  data kw_label(215)/'pseudopotential_format'/
  data kw_typ(215)/'T:I'/
  data kw_dscrpt(215)/'*! Format of pseudopotential files !*'/

  data kw_label(216)/'boundary_conditions'/
  data kw_typ(216)/'T:B'/
  data kw_dscrpt(216)/'*!Type of boundary conditions: cluster, tube, slab, crystal!*'/

  data kw_label(217)/'coordinate_scale'/
  data kw_typ(217)/'P:E'/
  data kw_dscrpt(217)/'*! Scale in atomic coordinates !*'/

  data kw_label(401)/'latticevecs'/
  data kw_typ(401)/'B:E'/
  data kw_dscrpt(401)/'*! length of sides of periodic box !*'/

  data kw_label(402)/'coordinates'/
  data kw_typ(402)/'B:E'/
  data kw_dscrpt(402)/'*! coordinates of atom !*'/

  data kw_label(403)/'coordinates_absolute'/
  data kw_typ(403)/'E:E'/
  data kw_dscrpt(403)/'*! specify atom coordinates in Cartesian coord. !*'/

  data kw_label(404)/'pp_format'/
  data kw_typ(404)/'I:I'/
  data kw_dscrpt(404)/'*! Format of pseudopotential files !*'/

  data kw_label(405)/'pseudopotential'/
  data kw_typ(405)/'B:E'/
  data kw_dscrpt(405)/'*! specifications of pseudopotentials !*'/

  ! Parameter for Interpolative Separable Density Fitting Method (ISDF)
  data kw_label(406)/'doisdf'/
  data kw_typ(406)/'E:E'/
  data kw_dscrpt(406)/'*! Do interpolative separable density fitting !*'/

  ! Number of the interpolation points for ISDF method
  data kw_label(407)/'num_isdf_points'/
  data kw_typ(407)/'I:E'/
  data kw_dscrpt(407)/'*! The number of interpolation points for isdf !*'/

  ! Method for finding the interpolation points for ISDF method
  data kw_label(408)/'intp_type'/
  data kw_typ(408)/'I:E'/
  data kw_dscrpt(408)/'*! &
&    1. cvt method using charge density as weight; &
&    2. cvt method using the modular square of all wavefunctions as weight  *!'/

  ! Method for finding the interpolation points for ISDF method
  data kw_label(409)/'isdf_type'/
  data kw_typ(409)/'I:E'/
  data kw_dscrpt(409)/'*! different implementations of isdf. &
&    1. The old paralleization. 2. New implementation with P and Q  *!'/

  data kw_label(410)/'c_intp_filter'/
  data kw_typ(410)/'D:E'/
  data kw_dscrpt(410)/'*!  *!'/

  ! Parameter for using less memory
  data kw_label(411)/'lessmemory'/
  data kw_typ(411)/'I:E'/
  data kw_dscrpt(411)/'*! Use lessmemory. Call isdf_parallel_sym_lessmem !*'/

  ! Parameter for faster matrix elements selection
  data kw_label(412)/'fastselect'/
  data kw_typ(412)/'E:E'/
  data kw_dscrpt(412)/'*! Use fastselect. Call k_integrate_isdf_fastselect !*'/

  ! Parameter for eigensolver
  data kw_label(413)/'eigsolver'/
  data kw_typ(413)/'T:E'/
  data kw_dscrpt(413)/'*! call which eigsolver to diagonalize Casida equation !*'/

  ! Parameter for linear algebra subroutines (such as dgemm)
  data kw_label(414)/'linear_algebra'/
  data kw_typ(414)/'T:E'/
  data kw_dscrpt(414)/'*! use gpu or cpu to perform matrix multiplications? !*'/

  ! Parameter for lanczos GW method
  data kw_label(415)/'lanczos_gw'/
  data kw_typ(415)/'E:E'/
  data kw_dscrpt(415)/'*! Use lanczos method to speed up GW calculation !*'/

  ! Parameter for using less memory
  data kw_label(416)/'nbl'/
  data kw_typ(416)/'I:E'/
  data kw_dscrpt(416)/'*! The scalapack blocksize for linear solver pdgesv !*'/

  ! Parameter for the dgree of polynomial
  data kw_label(417)/'lanczos_npoly'/
  data kw_typ(417)/'I:E'/
  data kw_dscrpt(417)/'*! The degree of fitting polynomial for lanczos method !*'/

  ! Parameter for the number of iterations for lanczos method
  data kw_label(418)/'lanczos_niter'/
  data kw_typ(418)/'I:E'/
  data kw_dscrpt(418)/'*! The number of iterations for lanczos method !*'/

  data kw_label(419)/'jnblock'/
  data kw_typ(419)/'I:E'/
  data kw_dscrpt(419)/'*! The number of bands treated in each iteration for the summation over bands !*'/

end module esdf_key
!===============================================================
