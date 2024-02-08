!===================================================================
!
! Takes a vector gw(1:3) in irreducible wedge and operates a rotation
! trans(3,3). Returns the rotated vector gr(1:3), properly shifted.
!
! Copyright (C) 2009 Murilo L. Tiago, http://users.ices.utexas.edu/~mtiago
! This file is part of RGWBS. It is distributed under the GPL v1.
!
!-------------------------------------------------------------------
subroutine unfold(gw,trans,shift,gr)

  use myconstants
  implicit none

  ! arguments
  integer, dimension(3), intent(in) :: gw
  real(dp), intent(in) :: trans(3,3), shift(3)
  integer, dimension(3), intent(out) :: gr

  ! local variables
  integer :: ii
  real(dp) :: xx

  do ii = 1, 3
     xx = trans(ii,ii) * gw(ii) + trans(ii,ii) * shift(ii) - shift(ii)
     gr(ii) = nint(xx)
  enddo

end subroutine unfold
!===============================================================
