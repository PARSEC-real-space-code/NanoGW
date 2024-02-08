!===================================================================
!
!  Cubic spline adapted from Numerical Recipes.
!
!  Given arrays xx(1:nn) and yy(1:nn) containing a tabulated
!  function, i.e., yi = f(xi), with x1 < x2 < .. . < xN, and
!  given values yp1 and ypn for the first derivative of the
!  interpolating function at points 1 and n, respectively, this
!  routine returns an array y2(1:n) of length n which contains
!  the second derivatives of the interpolating  function at the
!  tabulated points xi. If yp1 and/or ypn are equal to 1 x 10^30
!  or larger, the routine is signaled to set the corresponding
!  boundary condition for a natural spline, with zero second
!  derivative on that boundary.
!
!  Notice that yy, yp1 and ypn have two components, to account for real
!  and imaginary parts of the complex variables in parent routine.
!
!---------------------------------------------------------------
subroutine spline(xx,yy,nn,yp1,ypn,y2)

  use myconstants
  implicit none

  ! arguments
  integer, intent(in) :: nn
  real(dp), intent(in) :: xx(nn), yy(2,nn), yp1(2), ypn(2)
  real(dp), intent(out) :: y2(2,nn)

  ! local variables
  integer :: ii, kk, ll 
  real(dp) :: p, qn, sig, un, uu(2,nn) 
  real, parameter :: HUGE = 1.d10

  do ll = 1, 2
     !---------------------------------------------------------------
     ! The lower boundary condition is set either to be "natural"
     ! or else to have a specified first derivative. 
     if (abs(yp1(ll)) > HUGE) then
        y2(ll,1)= zero
        uu(ll,1) = zero
     else
        y2(ll,1) = -half
        uu(ll,1) = (three/(xx(2)-xx(1)))* &
             ((yy(ll,2)-yy(ll,1))/(xx(2)-xx(1))-yp1(ll))
     endif

     ! This is the decomposition loop of the tridiagonal algorithm.
     ! y2 and u are used for temporary storage of the decomposed factors. 
     do ii = 2, nn - 1
        sig = (xx(ii)-xx(ii-1))/(xx(ii+1)-xx(ii-1))
        p = sig*y2(ll,ii-1) + two
        y2(ll,ii)= (sig-one)/p
        uu(ll,ii) = (six*((yy(ll,ii+1)-yy(ll,ii))/(xx(ii+1)-xx(ii))- &
             (yy(ll,ii)-yy(ll,ii-1))/(xx(ii)-xx(ii-1)))/(xx(ii+1)-xx(ii-1))- &
             sig*uu(ll,ii-1))/p
     enddo

     ! The upper boundary condition is set either to be "natural"
     ! or else to have a specified first derivative. 
     if (ypn(ll) >  HUGE) then
        qn = zero
        un = zero
     else
        qn = half
        un = (three/(xx(nn)-xx(nn-1)))* &
             (ypn(ll)-(yy(ll,nn)-yy(ll,nn-1))/(xx(nn)-xx(nn-1)))
     endif

     ! This is the backsubstitution loop of the tridiagonal algorithm. 
     y2(ll,nn) = (un-qn*uu(ll,nn-1))/(qn*y2(ll,nn-1)+one) 
     do kk = nn - 1, 1, -1
        y2(ll,kk) = y2(ll,kk)*y2(ll,kk+1)+uu(ll,kk)
     enddo
  enddo

end subroutine spline
!===============================================================
