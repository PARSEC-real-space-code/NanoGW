!===================================================================
!
! Inverts a 3x3 matrix m and returns the result in m. It also calculates
! determinant of the input matrix. Matrix inversion is aborted if
! det < del.
!
! INPUT:
!    m : input matrix
!
! OUTPUT:
!    m : inverse of input matrix
!    det : determinant of input matrix
!    ierr : error flag (0 for successful exit; 1 if m has zero determinant)
!
! Copyright (C) 2009 Murilo L. Tiago, http://users.ices.utexas.edu/~mtiago
! This file is part of RGWBS. It is distributed under the GPL v1.
!
!-------------------------------------------------------------------
subroutine mtrxin(m,det,ierr)

  use myconstants
  implicit none

  ! arguments
  real(dp), intent(inout) :: m(3,3)
  real(dp), intent(out) :: det
  integer, intent(out) :: ierr

  ! local variables
  real(dp) :: a(3,3)
  real(dp), parameter :: DEL = 1.d-5

  !-------------------------------------------------------------------
  ierr = 0

  ! Compute matrix of cofactors.
  a(1,1) = m(2,2)*m(3,3) - m(2,3)*m(3,2)
  a(2,1) = m(2,3)*m(3,1) - m(2,1)*m(3,3)
  a(3,1) = m(2,1)*m(3,2) - m(2,2)*m(3,1)
  a(1,2) = m(1,3)*m(3,2) - m(1,2)*m(3,3)
  a(2,2) = m(1,1)*m(3,3) - m(1,3)*m(3,1)
  a(3,2) = m(1,2)*m(3,1) - m(1,1)*m(3,2)
  a(1,3) = m(1,2)*m(2,3) - m(1,3)*m(2,2)
  a(2,3) = m(1,3)*m(2,1) - m(1,1)*m(2,3)
  a(3,3) = m(1,1)*m(2,2) - m(1,2)*m(2,1)

  ! Compute determinant.
  det = m(1,1) * a(1,1) + m(1,2) * a(2,1) + m(1,3) * a(3,1)
  if (abs(det) < DEL) then
     write(6,*) 'ERROR in mtrxin: zero determinant'
     ierr = 1
     return
  endif

  ! Compute inverse matrix.
  m = a/det

end subroutine mtrxin
!===================================================================
