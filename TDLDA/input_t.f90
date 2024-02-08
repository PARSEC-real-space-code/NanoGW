!===================================================================
!
! Read input parameters from file rgwbs.in
! Only search for parameters specific to TDLDA.
!
! Copyright (C) 2009 Murilo L. Tiago, http://users.ices.utexas.edu/~mtiago
! This file is part of RGWBS. It is distributed under the GPL v1.
!
!-------------------------------------------------------------------
subroutine input_t(tamm_d,rpaonly,trip_flag,noxchange,trunc_c)
  !
  ! Read input parameters from file rgwbs.in
  !
  use myconstants
  use esdf
  implicit none

  ! arguments
  ! true if Tamm-Dancof approximation is used
  logical, intent(inout) :: tamm_d
  ! true if only RPA excitations are calculated (useful to get the sum rule quickly)
  logical, intent(out) :: rpaonly
  ! true if triplet excitations are calculated, false otherwise
  logical, intent(out) :: trip_flag
  ! true if exchange part of kernel is ignored (debug purposes only)
  logical, intent(out) :: noxchange
  ! true if we want to truncate the Coulomb interaction (important if we calculate the macroscopic dielectric function of a periodic system)
  logical, intent(out) :: trunc_c
  !-----------------------------------------------------------------------
  ! Start reading input parameters until the end of file is reached.
  !
  call esdf_init('rgwbs.in')

  rpaonly =  (esdf_defined('rpa_spectrum_only'))

  trip_flag =  (esdf_defined('tdlda_triplet_kernel'))

  noxchange =  (esdf_defined('no_exchange'))

  trunc_c =  (esdf_defined('truncate_coulomb'))

  call esdf_close

  if ( trip_flag ) noxchange = .true.
  if ( trunc_c ) tamm_d = .true.

  return

end subroutine input_t
!===================================================================
