
! Compute the spectra function of Hermitian matrix using
! lanczos algorithm
! what we calculate is
!
!  <vec| 1/(z-H) |vec>
!
subroutine lanczos_spectra_isdf(v0, ndim, zz, nz, niter, H, spectra)

  use typedefs

  integer, intent(in) :: ndim, niter, nz
  complex(8), intent(in) :: zz(nz)
  real(8), intent(in) :: v0(ndim)
  type(dmat), intent(in) :: H
  complex(8), intent(out) :: spectra(nz)
  complex(8) :: enum, deno
  real(8) :: v0_norm, vec(ndim), wvec(ndim), tmpvec(ndim)
  real(8), allocatable :: a_elmt(:), b_elmt(:)
  real(8), parameter :: eps = 1.0d-8
  integer :: i, j, k, iz
  real(8), external :: ddot

  ! normalize the vector
  v0_norm = sqrt(ddot(ndim, v0, 1, v0, 1)) ! compute the norm of the input vector
  !print *, v0_norm
  vec = v0/v0_norm
  !print *, vec

  allocate (a_elmt(niter))
  allocate (b_elmt(niter))

  b_elmt(1) = 1.d0

  ! compute wvec = H*vec
  ! wvec = polynomial_matvec(vec, ndim)
  call dgemv('n', ndim, ndim, 1.d0, H%h, ndim, vec, 1, 0.d0, wvec, 1)
  ! dot(vec, wvec)
  a_elmt(1) = ddot(ndim, wvec, 1, vec, 1)
  ! wvec = wvec - a_elmt(1) * vec
  call daxpy(ndim, -a_elmt(1), vec, 1, wvec, 1)
  b_elmt(2) = sqrt(ddot(ndim, wvec, 1, wvec, 1))
  do k = 2, niter - 1
    if (abs(b_elmt(k)) < eps) exit
    wvec = wvec/b_elmt(k)
    vec = -vec*b_elmt(k)
    ! swap wvec and vec
    tmpvec = vec
    vec = wvec
    wvec = tmpvec
    ! wvec = wvec + H*vec
    ! wvec = wvec + polynomial_matvec(vec, ndim)
    call dgemv('n', ndim, ndim, 1.d0, H%h, ndim, vec, 1, 1.d0, wvec, 1)
    a_elmt(k) = ddot(ndim, wvec, 1, vec, 1)
    call daxpy(ndim, -a_elmt(k), vec, 1, wvec, 1)
    b_elmt(k + 1) = sqrt(ddot(ndim, wvec, 1, wvec, 1))
  end do
  spectra = 0.d0
  do iz = 1, nz
    do j = k, 1, -1
      enum = b_elmt(k)*b_elmt(k)
      if (j == 1) enum = (1.d0, 0.d0)
      deno = zz(iz) - a_elmt(j) - spectra(iz)
      spectra(iz) = enum/deno
    end do ! j
  end do ! iz
  !print *, "a_elmt ", a_elmt
  !print *, "a_elmt ", b_elmt
  deallocate (a_elmt)
  deallocate (b_elmt)

end subroutine

