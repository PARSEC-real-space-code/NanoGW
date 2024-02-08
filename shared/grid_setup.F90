!===============================================================
!
!  Define the full real-space grid from the irreducible wedge.
!
!   If a point r_i in the full grid is equivalent to point r_j in
!   the irreducible wedge by applying symmetry operation k, then
!        gvec%rindex(i) = j
!        gvec%rtrans(i) = k
!   The grid is built so that the first Nr points (Nr: size of the
!   irreducible wedge) are in the irreducible wedge. The next set of
!   Nr points are obtained by applying the first symmetry operation
!   on the wedge, and so on.
!
! INPUT:
!    gvec : definition of the irreducible wedge.
!
! OUTPUT:
!    gvec%avec_norm : normalized lattice vectors
!    gvec%rindex : index of the grid point in irreducible wedge that
!       is equivalent to each point in the full grid.
!    gvec%rtrans : symmetry operation needed to bring each point
!       in the full grid to the irreducible wedge.
!    kr : indices of grid points in the full grid.
!
! Copyright (C) 2009 Murilo L. Tiago, http://users.ices.utexas.edu/~mtiago
! This file is part of RGWBS. It is distributed under the GPL v1.
!
!---------------------------------------------------------------
subroutine grid_setup(gvec,kr)

  use typedefs

  implicit none

  ! arguments
  type (gspace), intent(inout):: gvec
  ! Coordinates of grid points in real space: rr_i = (shift + kr_i)*h
  integer, intent(out) :: kr(3,gvec%nr*gvec%syms%ntrans)

  ! local variables
  integer :: ii, jj, kk, info, ndimp1
  real(dp) :: rvec(3), tvec(3)

  !---------------------------------------------------------------
  !  Calculate reciprocal space metric.
  !
  gvec%avec_norm = gvec%avec
  do ii = 1, 3
     rvec(1) = sqrt(sum(gvec%avec_norm(:,ii)**2))
     gvec%avec_norm(:,ii) = gvec%avec_norm(:,ii)/rvec(1)
  enddo

  ndimp1 = gvec%nr * gvec%syms%ntrans + 1

  allocate(gvec%rindex(ndimp1),stat=info)
  call alccheck('rindex','grid_setup',ndimp1,info)
  allocate(gvec%rtrans(ndimp1),stat=info)
  call alccheck('rtrans','grid_setup',ndimp1,info)

  ! Define full grid. The first part of the grid is the irreducible wedge,
  ! followed by images of the wedge.
  do kk = 1, gvec%syms%ntrans
     do jj = 1, gvec%nr
        ii = jj + (kk-1)*gvec%nr
        rvec = gvec%r(:,jj) + gvec%shift
        tvec = matmul(rvec,gvec%syms%trans(:,:,kk)) - gvec%shift
        kr(:,ii) = nint(tvec(:))
        gvec%rindex(ii) = jj
        gvec%rtrans(ii) = kk
     enddo
  enddo

  gvec%rindex(ndimp1) = gvec%nr + 1
  gvec%rtrans(ndimp1) = 1

end subroutine grid_setup
!===============================================================
