!===================================================================
!
! Wraper for the cubic spline interpolator.
! Input is ntable, xtable, ftable with the same meaning as in newtonr
! Additional input:
!      xx : value at which f(x) will be calculated
! Output is fx = f(x) and its derivative df = Re{df/dx} at x.
!
! Copyright (C) 2009 Murilo L. Tiago, http://users.ices.utexas.edu/~mtiago
! This file is part of RGWBS. It is distributed under the GPL v1.
!
!-------------------------------------------------------------------
subroutine spline_evaluate(ntable,xtable,ftable,xx,fx,df)

  use myconstants
  implicit none

  ! arguments
  ! number of tabulated data values
  integer, intent(in) :: ntable
  ! tabulated values of the argument
  real(dp), intent(in) :: xtable(ntable)
  ! tabulated values of the function, ftable(i) = F[xtable(i)]
  complex(dpc), intent(in) :: ftable(ntable)
  ! value of argument where function is calculated
  real(dp), intent(in) :: xx
  ! interpolated values of function and its derivative at xx
  complex(dpc), intent(out) :: fx
  real(dp), intent(out) :: df

  ! local variables
  real(dp) :: xm, xp
  complex(dpc) :: yp1, ypn, ym, yp, y2(ntable-2)
  ! delta_x increment used to calculate the numerical derivative
  real(dp), parameter :: DELTAX = 1.d-2

  yp1 = (ftable(2) - ftable(1))/(xtable(2) - xtable(1))
  ypn = (ftable(ntable) - ftable(ntable-1))/(xtable(ntable) - xtable(ntable-1))

  call spline(xtable(2),ftable(2),ntable-2,yp1,ypn,y2)
  call splint(xtable(2),ftable(2),y2,ntable-2,xx,fx)

  xm = xx - DELTAX
  call splint(xtable(2),ftable(2),y2,ntable-2,xm,ym)
  xp = xx + DELTAX
  call splint(xtable(2),ftable(2),y2,ntable-2,xp,yp)

  df = ( real(yp - ym,dp) )/ DELTAX * half

  return
end subroutine spline_evaluate
!===================================================================
