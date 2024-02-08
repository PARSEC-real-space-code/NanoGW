!===================================================================
!
! This is a simple wrapper to the subroutines that calculate laplacian
! with finite differences or with FFT.
!
! Copyright (C) 2009 Murilo L. Tiago, http://users.ices.utexas.edu/~mtiago
! This file is part of RGWBS. It is distributed under the GPL v1.
!
!-------------------------------------------------------------------
subroutine get_lap_sym(gvec,given_func,laplf) 

  use myconstants
  use typedefs
  use fft_module
  use fd_module
  implicit none

  ! arguments
  ! grid
  type (gspace), intent(inout) :: gvec
  ! input function
  real(dp), dimension (gvec%nr), intent(in) :: given_func
  ! output function, lapl = -Laplacian(given_func)
  real(dp), dimension (gvec%nr), intent(out) :: laplf

  !-------------------------------------------------------------------

  if (fd%norder < 0) then
     call dget_lap_FFT(gvec,given_func,laplf,1)
  else
     call dget_lap_fd(gvec%syms,given_func,laplf,1)
  endif

end subroutine get_lap_sym
!===============================================================
