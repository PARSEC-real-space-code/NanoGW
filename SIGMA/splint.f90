!===================================================================
!
!   Cubic spline adapted from Numerical Recipes.
!
!   Given a value x_0, given the bounds xx(1) and xx(2) with
!   xx(1) < x_0 < xx(2), given the corresponding values of a function
!   y(xhi) = yy(1), y(xlo) = yy(2), and given the array y2a(1:nn),
!   which is the output from spline, this routine returns a
!   cubic-spline interpolated value y_0 = y(x_0).
!
!---------------------------------------------------------------
subroutine splint(xx,yy,y2a,nn,x_0,y_0)

  use myconstants
  implicit none

  ! arguments
  integer, intent(in) :: nn
  real(dp), intent(in) :: x_0, xx(nn)
  complex(dpc), intent(in) :: yy(nn), y2a(nn)
  complex(dpc), intent(out) :: y_0

  ! local variables
  integer :: kk, khi, klo
  real(dp) :: aa, bb, hh

  !---------------------------------------------------------------

  klo = 1
  khi = nn

  do
     if (khi-klo > 1) then 
        kk = (khi + klo)/2
        if(xx(kk) > x_0)then 
           khi = kk 
        else 
           klo = kk 
        endif
     else
        exit
     endif
  enddo

  hh = xx(khi) - xx(klo) 
  if (hh == zero) call die('bad xx input in splint')
  aa = (xx(khi) - x_0)/hh
  bb = (x_0 - xx(klo))/hh 
  y_0 = aa*yy(klo)+bb*yy(khi) + ((aa**3-aa)*y2a(klo) + &
       (bb**3-bb)*y2a(khi))*(hh**2)/six

end subroutine splint
!===============================================================
