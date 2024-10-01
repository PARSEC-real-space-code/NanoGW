! Compute the spectra function of Hermitian matrix using
! lanczos algorithm
! what we calculate is
!
!  <vec| 1/(z-H) |vec>
!
subroutine lanczos_spectra_isdf(v0, v0_norm, ncv_loc, zz, nz, niter, isdf_in, &
                                wfn, tmp_Cmtrx, sqrtR, polynomial, npoly, spectra)

  use typedefs
  use mpi_module
  implicit none
  include 'mpif.h'

  integer, intent(in) :: ncv_loc, niter, nz, npoly
  complex(dp), intent(in) :: zz(nz)
  type(ISDF), intent(in) :: isdf_in
  real(dp), intent(in) :: v0(ncv_loc), polynomial(npoly + 1), v0_norm, &
                          tmp_Cmtrx(isdf_in%n_intp_r, ncv_loc), sqrtR(ncv_loc)
  complex(dp), intent(out):: spectra(nz)
  type(wavefunction), intent(in) :: wfn
  real(dp), allocatable :: a_elmt(:), b_elmt(:)
  complex(8) :: enum, deno
  real(dp), parameter :: eps = 1.0d-8
  integer :: i, j, k, iz, info
  real(dp), external :: ddot
  real(dp) :: vec(ncv_loc), wvec(ncv_loc), tmpvec(ncv_loc)

  ! normalize the vector
  vec = v0/v0_norm

  allocate (a_elmt(niter))
  allocate (b_elmt(niter))

  b_elmt(1) = 1.d0

  ! compute wvec = H*vec
  if (peinf%master) print *, "vec = ", vec(1:10)
  call polynomial_matvec_isdf(vec, wvec, ncv_loc, polynomial, npoly, &
                              isdf_in, wfn, tmp_Cmtrx, sqrtR)
  if (peinf%master) print *, "wvec = ", wvec(1:10)
  ! stop
  ! dot(vec, wvec)
  a_elmt(1) = ddot(ncv_loc, wvec, 1, vec, 1)
  !if (peinf%master) print *, "before allreduce a_elmt ", a_elmt(1)
  call mpi_allreduce(MPI_IN_PLACE, a_elmt(1), 1, &
                     MPI_DOUBLE, MPI_SUM, peinf%comm, info)
  !if (peinf%master) print *, "after allreduce a_elmt ", a_elmt(1)
  ! wvec = wvec - a_elmt(1) * vec
  call daxpy(ncv_loc, -a_elmt(1), vec, 1, wvec, 1)
  b_elmt(2) = ddot(ncv_loc, wvec, 1, wvec, 1)
  !if (peinf%master) print *, "before allreduce b_elmt ", b_elmt(2)
  call mpi_allreduce(MPI_IN_PLACE, b_elmt(2), 1, &
                     MPI_DOUBLE, MPI_SUM, peinf%comm, info)
  !if (peinf%master) print *, "before allreduce b_elmt ", b_elmt(2)
  b_elmt(2) = sqrt(b_elmt(2))
  do k = 2, niter
    if (abs(b_elmt(k)) < eps) exit
    wvec = wvec/b_elmt(k)
    vec = -vec*b_elmt(k)
    ! swap wvec and vec
    tmpvec = vec
    vec = wvec
    wvec = tmpvec
    ! wvec = wvec + H*vec
    call polynomial_matvec_isdf(vec, tmpvec, ncv_loc, polynomial, npoly, isdf_in, &
                                wfn, tmp_Cmtrx, sqrtR)
    wvec = wvec + tmpvec
    a_elmt(k) = ddot(ncv_loc, wvec, 1, vec, 1)
    call mpi_allreduce(MPI_IN_PLACE, a_elmt(k), 1, &
                       MPI_DOUBLE, MPI_SUM, peinf%comm, info)
    if (k + 1 <= niter) then
      call daxpy(ncv_loc, -a_elmt(k), vec, 1, wvec, 1)
      b_elmt(k + 1) = ddot(ncv_loc, wvec, 1, wvec, 1)
      call mpi_allreduce(MPI_IN_PLACE, b_elmt(k + 1), 1, &
                         MPI_DOUBLE, MPI_SUM, peinf%comm, info)
      b_elmt(k + 1) = sqrt(b_elmt(k + 1))
    end if
  end do
  !if (peinf%master) then
  !  print *, "j  a_elmt  b_elmt "
  !  do j = 1, niter
  !    print *, j, a_elmt(j), b_elmt(j)
  !  enddo
  !endif
  spectra = 0.d0
  do iz = 1, nz
    do j = niter, 1, -1
      enum = b_elmt(j)*b_elmt(j)
      if (j == 1) enum = (1.d0, 0.d0)
      deno = zz(iz) - a_elmt(j) - spectra(iz)
      spectra(iz) = enum/deno
      !if (nz .eq. 1) print *, "check ", j, enum, deno
    end do ! j
    !print *, "spectra ", iz, zz(iz), spectra(iz)
  end do ! iz
  deallocate (a_elmt)
  deallocate (b_elmt)

end subroutine

