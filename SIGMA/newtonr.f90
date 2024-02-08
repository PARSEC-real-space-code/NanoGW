!===================================================================
!
! Newton-Raphson node searcher adapted from Numerical Recipes.
!
! Search for root Re{f(x)} = 0 using Newton-Raphson algorithm.
! Input is function f(x) tabulated at set of regularly spaced points:
!   f(x) = ftable[1] for x = xtable[1]
!   f(x) = ftable[2] for x = xtable[2]
!     . . .
!   f(x) = ftable[ntable] for x = xtable[ntable]
!
! Output is the estimate of root, Re{f(x)} = 0 in the interval
!             x_zero - xacc < x < x_zero + xacc
! At input, x_zero is a first estimate of the root, and xacc is the
! desired accuracy.
! 
!-------------------------------------------------------------------
subroutine newtonr(ntable,xtable,ftable,x_zero,xacc,fx,ierr)

  use myconstants
  implicit none

  ! arguments
  ! number of tabulated data values
  integer, intent(in) :: ntable
  ! tabulated values of the argument
  real(dp), intent(in) :: xtable(ntable)
  ! tabulated values of the function, ftable(i) = F[xtable(i)]
  complex(dpc), intent(in) :: ftable(ntable)
  ! input: initial value of zero of F
  ! output: final value of zero of F, F[x_zero] = 0
  real(dp), intent(inout) :: x_zero
  ! accuracy in determination of root
  real(dp), intent(in) :: xacc
  ! value of function at the calculated root
  complex(dpc), intent(out) :: fx
  ! error flag
  integer, intent(out) :: ierr

  ! local variables
  integer :: jj
  real(dp) :: df, dx
  ! maximum number of Newton-Raphson iterations
  integer, parameter :: MAX_IT = 20

  ierr = 0
  do jj = 1, MAX_IT
     call spline_evaluate(ntable,xtable,ftable,x_zero,fx,df)
     dx = real(fx,dp) / df
     x_zero = x_zero - dx
     if( abs(dx) < xacc ) then
        call spline_evaluate(ntable,xtable,ftable,x_zero,fx,df)
        return
     endif
  enddo
  ierr = jj

end subroutine newtonr
!===================================================================
