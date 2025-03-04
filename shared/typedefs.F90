!===================================================================
!
! Definition of derived types.
!
! Copyright (C) 2009 Murilo L. Tiago, http://users.ices.utexas.edu/~mtiago
! This file is part of RGWBS. It is distributed under the GPL v1.
!
!-------------------------------------------------------------------
module typedefs

  use myconstants
  !
  type options
    !
    integer :: eigsolver
    integer :: linear_algebra
    integer :: pdgesv_nbl2d
    integer :: lanczos_npoly
    integer :: lanczos_niter
    integer :: jnblock
    !
  end type options

  type ISDF
    !
    logical :: fastselect
    integer :: lessmemory
    integer :: n_slice
    integer :: maxicc
    integer :: maxivv ! this is not the highest occupied state, for GW
    ! calculation, maxivv would be the largest index of
    ! state for which GW energy is calculated
    integer :: maxncv
    integer :: maxnv
    integer :: maxnc
    ! number of interpolation points in irreducible domain
    integer :: n_intp_r
    ! number of interpolation points in irreducible domain stored locally in
    ! the current process
    integer :: n_intp_r_loc
    ! the index of interpolation points in irreducible domain
    integer, pointer :: intp_r(:)
    ! the number orbital pairs for each (kpt, spin and representation)
    integer, pointer :: ncv(:, :, :)
    integer, pointer :: invpairmap(:, :, :, :, :)
    integer, pointer :: pairmap(:, :, :, :, :)
    integer, pointer :: nv(:, :, :)
    integer, pointer :: ivlist(:, :, :, :)
    integer, pointer :: nc(:, :, :)
    integer, pointer :: iclist(:, :, :, :)
    real(dp), pointer :: Cmtrx(:, :, :, :, :)
    real(dp), pointer :: Mmtrx(:, :, :, :, :, :, :)
    real(dp), pointer :: Psi_intp(:, :, :, :)
    real(dp), pointer :: Mmtrx_loc(:, :, :, :, :, :, :)
    real(dp), pointer :: Psi_intp_loc(:, :, :, :)
    real(dp), pointer :: PsiC_intp_bl(:, :, :, :), PsiV_intp_bl(:, :, :, :)
    integer :: xvgrp, xcgrp, myvgrp, mycgrp, mynv, mync, &
               ldv, ldc, myvstart, myvend, mycstart, mycend

    ! We assume nkp = nspin = 1 for now, kp and spin indices
    ! are needed later
    !
    ! Index arrays for calculate_sigma_lanczos_lowcomsym
    ! Arrays map to the real index of wave functions
    integer, pointer :: lrep(:), lrep_bound(:, :)
    ! Index arrays for PsiV_intp_bl and PsiC_intp_bl
    ! Arrays map to the index of local arrays PsiV_intp_bl and
    ! PsiC_intp_bl
    integer, pointer :: lvrep(:), lvrep_bound(:, :), &
                        lcrep(:), lcrep_bound(:, :)
    integer, pointer :: myncv_sym(:), mync_sym(:), mynv_sym(:)
    integer :: maxmync_sym, maxmynv_sym, maxmyncv_sym
  end type ISDF
  !
  ! Parameters for symmetry groups
  !
  type symmetries
    ! number of symmetry groups, excluding the Abelian group
    integer :: ngr
    ! label of each symmetry group, excluding Abelian
    character(len=80), pointer :: grfilename(:)
    ! number of representations in each symmetry groups, excluding Abelian
    integer, pointer :: nrep(:)
    ! product table of representations in each symmetry group, excluding
    ! Abelian. If r_prod(i,j,l,g) = n then the product of representations
    ! 'i' and 'j' of group 'g' contains 'n' replicas of representation 'l'.
    integer, pointer :: g_prod(:, :, :, :)
    ! multiplication table of Abelian group
    integer, pointer :: prod(:, :)
    ! number of representations in Abelian group = order of group
    integer :: ntrans
    ! rotation matrices for all operations in Abelian group, cartesian units
    real(dp), pointer :: trans(:, :, :)
    ! character table: chi(irep,iop) is character of representation irep
    ! under operation iop in Abelian group
    integer, pointer :: chi(:, :)
  end type symmetries
  !
  ! Real space and reciprocal space variables
  !
  type gspace
    ! nr : number of points in irreducible wedge
    integer :: nr
    ! number of periodic directions
    integer :: per
    ! unit lattice vectors (identity matrix for non-periodic system)
    real(dp), dimension(3, 3) :: avec
    ! normalized unit lattice vectors (identity matrix for non-periodic system).
    real(dp), dimension(3, 3) :: avec_norm
    ! reciprocal unit lattice vectors
    real(dp), dimension(3, 3) :: bvec
    ! metric matrix in reciprocal space,
    ! bdot(i,j) = dot_product(bvec(:,i),bvec(:,j))
    real(dp), dimension(3, 3) :: bdot
    ! boundary radius
    real(dp) :: rmax
    ! grid spacing along each direction, in units of Bohr radius
    real(dp) :: step(3)
    ! grid shift, the cartesian coordinates of grid point (i,j,k) are
    ! defined as x = ( i + shift )*step, and similarly for y,z
    real(dp) :: shift(3)
    ! hcub = volume element in the grid, units (bohr)^3
    real(dp) :: hcub
    ! 1-D (wire) systems: length of periodic cell
    ! 2-D (slab) systems: area of periodic cell
    ! 3-D (bulk) systems: volume of periodic cell
    real(dp) :: celvol
    ! coordinates of real space points
    integer, pointer :: r(:, :)
    ! irreducible wedge parameters
    integer, dimension(:), pointer :: rindex
    integer, dimension(:), pointer :: rtrans
    ! long-wavelength part of Coulomb potential. Since this potential is
    ! 1/r, its Fourier transform diverges at q = 0, with a different
    ! divergence behavior for each dimensionality (wire, slab, bulk)
    real(dp) :: long

    ! symmetry operations in the Abelian group
    type(symmetries) :: syms

  end type gspace
  !
  ! Electron wavefunctions
  !
  type wavefunction
    ! number of LDA states, read from wfn.dat
    integer :: nstate
    ! number of LDA states for which kernel matrix elements are computed
    integer :: nmem
    ! occupancy of DFT states
    real(dp), pointer :: occ0(:)
    ! occupancy of DFT/QP states (read from occup.in or, if absent, from wfn.dat)
    real(dp), pointer :: occ1(:)
    ! E_lda energy eigenvalues, in Ry
    real(dp), pointer :: e0(:)
    ! E_lda energy eigenvalues + scissors operator (if defined), in Ry
    real(dp), pointer :: e1(:)
    ! irreducible representation of electronic states, Abelian group
    integer, pointer :: irep(:)
    ! irreducible representation of electronic states, generic groups
    integer, pointer :: jrep(:, :)
    ! projections of electron wavefunctions in symmetry representations
    real(dp), pointer :: proj(:, :)
    ! map: states used in kernel -> states read from wfn.dat
    ! imap: inverse map
    integer, pointer :: map(:), imap(:)
    ! number of dipole matrix elements
    integer :: ndip
    ! map for dipole matrix elements
    integer, pointer :: mapd(:, :)
    ! dipole matrix elements
    real(dp), pointer :: ddipole(:, :)
    complex(dpc), pointer :: zdipole(:, :)
    ! electron wavefunctions
    real(dp), pointer :: dwf(:, :)
    complex(dpc), pointer :: zwf(:, :)
    ! inverse maps
    integer, pointer :: cmapi(:), vmapi(:)
    ! quasi-particle Hamiltonian
    real(dp), pointer :: hqpvv(:, :), hqpcc(:, :)
    ! quasi-particle energy of occupied and unoccupied
    ! states (after diagonalizing QP Hamiltonian)
    real(dp), pointer :: eqpv(:), eqpc(:)
  end type wavefunction
  !
  ! k-point structure (k-point = generic crystal momentum)
  !
  type kptinfo
    ! flag for complex wave-functions: true if they are complex;
    ! false otherwise
    logical :: lcplx
    ! number of k-points in the full Brillouin zone
    integer :: nk
    ! coordinates of k-points in units of reciprocal lattice vectors
    real(dp), dimension(:, :), pointer :: fk
    ! weights of k-points
    real(dp), dimension(:), pointer :: weight
    ! electron density for each spin channel (units of bohr^-3)
    real(dp), pointer :: rho(:, :)
    ! electron wave-functions
    type(wavefunction), dimension(:, :), pointer :: wfn
  end type kptinfo
  !
  ! q-points structure (q-point = generic photon momentum)
  !
  type qptinfo
    ! number of k-points in the full Brillouin zone
    integer :: nk
    ! coordinates of k-points in units of reciprocal lattice vectors
    real(dp), dimension(:, :), pointer :: fk
    ! true if this k-point has zero length, false otherwise
    logical, dimension(:), pointer :: zerok
    ! match of k-points
    integer, dimension(:, :), pointer :: match
  end type qptinfo
  !
  ! Kernels used in sigma, tdlda, and bsesolv calculations
  ! matrix elements are defined in double basis as:
  ! kernel(ii,jj) == int_dx dx' [ phi_i1(x)  phi_i2(x)  ] *
  !                       [ phi_i3(x') phi_i4(x') ] * V_int(x,x')
  !       i1 = map1( row(ii,1) )
  !       i2 = map2( row(ii,2) )
  !       i3 = map3( col(jj,1) )
  !       i4 = map4( col(jj,2) )
  !
  type kernelinfo
    logical :: isdf
    ! number of rows and columns of full kernel matrix
    integer :: nrow, ncol
    ! number of rows and columns for spin up kernel
    integer :: nrow_up, ncol_up
    ! indices of pairs for rows and columns of kernel matrix
    integer, dimension(:, :), pointer :: row, col
    ! number of columns or rows of distributed kernel
    integer :: nn
    ! length of buffer space
    integer :: nbuff
    ! length of integration boxes
    integer :: lcache
    ! distributed interaction kernel (local array)
    ! K_ij^mn = < i j | K | m n> = m(u,v) where
    !   row(u,1)=i   row(u,2)=j  row(u,3)=k_i  row(u,4)=k_j, 1 <= u <= ncol
    !   col(v,1)=m   col(v,2)=n  col(v,3)=k_m  col(v,4)=k_n, 1 <= v <= nrow
    real(dp), pointer :: dm(:, :)
    complex(dpc), pointer :: zm(:, :)
  end type kernelinfo
  !
  ! Scratch structure used in kernel and k_integrate.
  !
  type kernel_scratch
    integer :: nblock, nblock_max
    integer, pointer :: quadr(:, :)
    real(dp), pointer :: dm_f(:)
    complex(dpc), pointer :: zm_f(:)
  end type kernel_scratch
  !
  ! Polarizability (TDLDA or BSE)
  !
  type polinfo
    ! number of occupied states in
    integer :: nval
    ! number of unoccupied states in
    integer :: ncond
    ! absolute (i.e. input) index of occupied (valence) states
    integer, pointer :: vmap(:)
    ! absolute (i.e. input) index of unoccupied (conduction) states
    integer, pointer :: cmap(:)
    ! total number of pair transitions
    integer :: ntr
    ! number of pair transitions with spin up (equal to ntr if nspin = 1)
    integer :: n_up
    ! number of eigenstates stored on each processor
    integer :: nn
    ! pair transition
    integer, pointer :: tr(:, :)
    ! dipole matrix elements
    real(dp), pointer :: ddipole(:, :)
    complex(dpc), pointer :: zdipole(:, :)
    ! eigenvalues
    real(dp), pointer :: eig(:)
    ! oscillator strength matrix elements
    real(dp), pointer :: ostr(:, :)
    ! true if v is allocated
    logical :: lv
    ! true if tv is allocated
    logical :: ltv
    ! eigenvectors
    ! v(i,j) = amplitude of transition between occupied LDA orbital tr(i,1)
    !          and unnocuppied LDA orbital tr(i,2) in j-th eigenstate.
    real(dp), pointer :: dv(:, :)
    complex(dpc), pointer :: zv(:, :)
    ! conjugate transpose of eigenvectors: tv(j,i) = conjg( v(i,j) )
    real(dp), pointer :: dtv(:, :)
    complex(dpc), pointer :: ztv(:, :)
  end type polinfo
  !
  ! Self-energy, sigma calculation
  !
  type siginfo
    ! type of exchange-correlation model (GW by default)
    integer :: xc
    ! maximum number of DFT orbitals included in Green's function (correlation)
    integer :: nmax_c
    ! energy range in derivative of self-energy
    real(dp) :: deltae
    ! number of data points in derivative of self-energy
    integer :: nen

    ! index of this k-point in the list kpt%fk(:,1:nk)
    integer :: indxk
    ! number of sigma diagonal matrix elements to be computed within COHSEX
    integer :: ndiag_s
    ! number of sigma diagonal matrix elements to be computed
    integer :: ndiag
    ! number of off-diagonal matrix elements to be computed within COHSEX
    integer :: noffd_s
    ! number of off-diagonal matrix elements to be computed
    integer :: noffd
    ! total number of states for which self-energy is computed
    integer :: nmap
    ! absolute (i.e. input) index of states for which self-energy is computed
    integer, pointer :: map(:)
    ! position of states for diagonal elements in map
    integer, pointer :: diag(:)
    ! position of states for off-diagonal matrix elements of sigma in map
    integer, pointer :: off1(:)
    integer, pointer :: off2(:)
    ! V_xc matrix element (LDA)
    complex(dpc), pointer :: vxcdiag(:)
    complex(dpc), pointer :: vxcoffd(:)
    ! self-energy components
    complex(dpc), pointer :: xdiag(:)
    complex(dpc), pointer :: xoffd(:)
    complex(dpc), pointer :: scdiag(:, :)
    complex(dpc), pointer :: scoffd(:, :)
    complex(dpc), pointer :: sexdiag(:, :)
    complex(dpc), pointer :: scsdiag(:)
    complex(dpc), pointer :: scsoffd(:)
    complex(dpc), pointer :: sgdiag(:, :)
    complex(dpc), pointer :: sgoffd(:, :)
    complex(dpc), pointer :: sgsdiag(:)
    complex(dpc), pointer :: sgsoffd(:)

    real(dp), pointer :: sigmai(:)
    ! self-energy operator in the basis of QP orbitals.
    ! This operator is defined as the difference between (Sigma + V_hartree)
    ! and (Vxc + V_hartree). If the input wavefunctions were
    ! obtained from a previous GW iteration, then the self-energy operator is
    ! the difference between (Sigma + V_hartree) and (Sigma + V_hartree) from
    ! the previous iteration. Potential V_hartree is included because the
    ! electron density could change after a orthogonalization of the QP
    ! Hamiltonian.
    complex(dpc), pointer :: sig_mat(:, :, :)

    ! whether or not using lanczos algorithm for GW calculation
    logical :: lanczos

  end type siginfo
  !
  !  Quasi-particle orbitals
  !
  type qpinfo
    ! number of QP orbitals
    integer :: neig
    ! order of this orbital in the set of QP
    integer, pointer :: jrep(:)
    ! QP hamiltonian in the original basis
    complex(dpc), pointer :: hqp(:, :)
    ! Self-energy in the original basis, real part
    complex(dpc), pointer :: sigmaqp(:, :)
    ! Self-energy in the original basis, imaginary part, diagonal part only
    real(dp), pointer :: sigmai(:)
    ! QP eigenvalues
    real(dp), pointer :: eqp(:)
    ! QP eigenvectors in the original basis, this defines the unitary
    ! rotation of basis.
    complex(dpc), pointer :: vqp(:, :)
  end type qpinfo

end module typedefs
!===================================================================
