!===================================================================
!
! For each vector k3 in list kpt_3 and each vector k2 in list
! kpt_2, find the vector k1 = k2 - k3 that is present in list kpt_1.
! Output is array match(i,j) = l, such that k1 = k2 - k3, where:
! k1 = l-th vector in kpt_1; k2 = i-th vector in kpt_2; k3 = j-th vector
! in kpt_3. If no match is found for a pair (i,j), return match(i,j) = 0.
!
! INPUT:
!    kpt_1, kpt_2, kpt_3: k-point types
!    bdot : reciprocal space metric,
!           bdot(i,j) = inner product between vectors b_i and b_j
!
! OUTPUT:
!    qmatch : matching array (see above)
!
! Copyright (C) 2009 Murilo L. Tiago, http://users.ices.utexas.edu/~mtiago
! This file is part of RGWBS. It is distributed under the GPL v1.
!
!-------------------------------------------------------------------
subroutine qpt_matching(nkpt1,nkpt2,nkpt3,kpt_1,kpt_2,kpt_3,bdot,qmatch)

  use typedefs
  implicit none

  ! arguments
  ! length of vector lists
  integer, intent(in) :: nkpt1, nkpt2, nkpt3
  ! coordinates of vectors, assumed periodic with period one
  real(dp), intent(in) :: kpt_1(3,nkpt1), kpt_2(3,nkpt2), kpt_3(3,nkpt3)
  ! metric
  real(dp), intent(in) :: bdot(3,3)
  integer, intent(out) :: qmatch(nkpt2,nkpt3)

  ! local variables
  integer :: ik, iq, ikp
  real(dp) :: kmq(3), diff(3), length
  real(dp), parameter :: tol_q = 1.d-8

!-------------------------------------------------------------------

  qmatch = 0
  do iq = 1, nkpt3
     do ik = 1, nkpt2
        kmq = kpt_2(:,ik) - kpt_3(:,iq)
        do ikp = 1, nkpt1
           diff = kmq - kpt_1(:,ikp) + 4.d0
           diff = diff - anint(diff)
           length = sqrt(dot_product(diff,matmul(bdot,diff)))
           if (length < tol_q) then
              qmatch(ik,iq) = ikp
              exit
           endif
        enddo
     enddo
  enddo

  return

end subroutine qpt_matching
!===================================================================
