!===================================================================
!
! Calculates the gradient of a function specified on the
! real-space grid using finite differences. Input function is assumed
! to be completely symmetric (trivial representation).
!
! WARNING: Array gradf returns with the gradient on the irreducible
! wedge first, and then the rest of the grid.
!
! Copyright (C) 2009 Murilo L. Tiago, http://users.ices.utexas.edu/~mtiago
! This file is part of RGWBS. It is distributed under the GPL v1.
!
!-------------------------------------------------------------------
subroutine get_grad_sym(gvec,given_func,gradf) 

  use myconstants
  use typedefs
  use fft_module
  use fd_module
  implicit none

  ! arguments
  ! grid
  type (gspace), intent(inout) :: gvec
  ! the function
  real(dp), dimension (gvec%nr), intent(in) :: given_func
  ! gradient
  real(dp), dimension (3,gvec%nr*gvec%syms%ntrans), intent(out) :: gradf

!-------------------------------------------------------------------

  if (fd%norder < 0) then
     call dget_grad_FFT(gvec,given_func,gradf,1)
  else
     call dget_grad_fd(gvec%syms,given_func,gradf,1)
  endif

end subroutine get_grad_sym
!===============================================================
