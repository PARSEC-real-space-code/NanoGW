!===================================================================
!
! Definitions of numerical constants, physical constants and various
! internal constants.
!
! Copyright (C) 2009 Murilo L. Tiago, http://users.ices.utexas.edu/~mtiago
! This file is part of RGWBS. It is distributed under the GPL v1.
!
!-------------------------------------------------------------------
module myconstants

  implicit none
  !
  ! variable types
  !
  integer, parameter :: dp = kind(1.0d0)
  integer, parameter :: dpc = kind((1.0d0,1.0d0))
  !
  ! numerical constants
  !
  real(dp), parameter :: zero = 0.0_dp
  real(dp), parameter :: one = 1.0_dp
  complex(dpc), parameter :: zzero = (0.0_dp,0.0_dp)
  complex(dpc), parameter :: zone = (1.0_dp,0.0_dp)
  complex(dpc), parameter :: zi = (0.0_dp,1.0_dp)
  real(dp), parameter :: two = 2.0_dp
  real(dp), parameter :: three = 3.d0
  real(dp), parameter :: four = 4.0_dp
  real(dp), parameter :: five = 5.d0
  real(dp), parameter :: six = 6.0_dp
  real(dp), parameter :: eight = 8.d0
  real(dp), parameter :: nine = 9.0_dp
  real(dp), parameter :: half = 0.5_dp
  real(dp), parameter :: third = 1.0_dp/3.0_dp
  real(dp), parameter :: mone = -1.0_dp
  real(dp), parameter :: pi = 3.1415926535897932384626433832795_dp
  !
  ! physical constants
  !
  ! Rydberg constant
  real(dp), parameter :: ryd = 13.60569193_dp
  ! Convert length units from atomic units (bohrs) to angstroms
  real(dp), parameter :: angs = 0.52917720859_dp
  !
  ! numerical parameters
  !
  ! minimum deviation of occupancy factors from 0.0 or 1.0
  !      (this helps searching for occupied/empty LDA orbitals)
  ! if occ(i) < tol, i is assumed completely empty
  ! if occ(i) > 1.0 - tol, i is assued completely occupied
  real(dp), parameter :: tol_occ = 0.01_dp
  ! maximum size of arrays
  integer, parameter :: maxdata = 400000
  ! maximum length of buffers in MPI_ALLREDUCE( dpsum function)
  integer, parameter :: MAXSIZE_MPI = 1000000
  ! constants for DFT code
  integer, parameter :: PARSEC = 1, PARATEC = 2
  ! constants for the energy dependence in self-energy (sigma only)
  integer, parameter :: SIG_LEFT = 1, SIG_RIGHT = 2, SIG_AV = 3
  ! exchange-correlation models, used by "sigma" only.
  integer, parameter :: &
       XC_GW = 1, &               ! GW exchange-correlation
       XC_HF = 2, &               ! Hartree-Fock exchange, no correlation
       XC_B3LYP = 3, &            ! B3LYP (model) functional
       ! see: A.D. Becke, J. Chem. Phys. 98, 5648 (1993)
       XC_LDA_CA = 4, &           ! LDA Ceperley-Alder Perdew-Zunger functional
       XC_GGA_PBE = 5, &          ! GGA Perdew-Burke-Ernzerhof functional
       XC_GGA_BLYP = 6            ! GGA BLYP functional

end module myconstants
!===================================================================
