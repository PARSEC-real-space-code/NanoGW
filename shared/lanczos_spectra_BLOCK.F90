! Compute the spectra function of Hermitian matrix using
! lanczos algorithm
! what we calculate is
!
!  <vec| 1/(z-H) |vec>
!
subroutine lanczos_spectra_isdf_BLOCK(v0, v0_norm, ncv_loc, zz, nz, niter, isdf_in, &
                                      wfn, cvpair, ncv, sqrtR, polynomial, npoly, spectra, &
                                      MC_vec, CMC_vec, tmp_vec, tmp_Hvec, tmp_cvpair, blksz)

  use typedefs
  use mpi_module
  implicit none
  include 'mpif.h'

  integer, intent(in) :: ncv_loc, niter, nz, npoly, blksz, ncv
  complex(dp), intent(in) :: zz(nz, blksz)
  type(ISDF), intent(in) :: isdf_in
  real(dp), intent(in) :: v0(ncv_loc, blksz), &
                          polynomial(npoly + 1), v0_norm(blksz), &
                          sqrtR(ncv_loc)
  integer, intent(in) :: cvpair(4, ncv)
  complex(dp), intent(out):: spectra(nz, blksz)
  type(wavefunction), intent(in) :: wfn
  real(dp), intent(inout) :: MC_vec(w_grp%myn_intp_r, blksz), &
                             CMC_vec(w_grp%ldncv, blksz), &
                             tmp_vec(isdf_in%n_intp_r, blksz), &
                             tmp_Hvec(w_grp%ldncv, blksz)
  integer, intent(inout) :: tmp_cvpair(2, w_grp%ldncv)
  real(dp)   :: a_elmt(niter, blksz), b_elmt(niter, blksz)
  complex(8) :: enum, deno
  real(dp), parameter :: eps = 1.0d-8
  integer :: ii, i, j, k, iz, info
  real(dp), external :: ddot
  real(dp) :: vec(ncv_loc, blksz), &
              wvec(ncv_loc, blksz), tmpvec(ncv_loc, blksz)

  !print *, "inside lanczos_spectra_BLOCK"
  ! normalize the vector
  do ii = 1, blksz
    vec(1:ncv_loc, ii) = v0(1:ncv_loc, ii)/v0_norm(ii)
  end do

! allocate(a_elmt(niter))
! allocate(b_elmt(niter))

  b_elmt(1, 1:blksz) = 1.d0

  ! compute wvec = H*vec
  !print *, "vec = ", vec(1:5,1), "...", &
  !    vec(ncv_loc-4:ncv_loc,1)

  call stopwatch(peinf%master, "Call matvec_isdf_BLOCK")
  call polynomial_matvec_isdf_BLOCK(vec, wvec, ncv_loc, polynomial, npoly, &
                                    isdf_in, wfn, cvpair, ncv, sqrtR, &
                                    MC_vec, CMC_vec, tmp_vec, tmp_Hvec, tmp_cvpair, blksz)
  !call stopwatch(peinf%master, "Done matvec_isdf")
  !if (peinf%master) print *, "wvec = ",wvec(1:10, 2), "...", &
  !    vec(ncv_loc-4:ncv_loc,2)
  ! stop
  ! dot(vec, wvec)
  !
  do ii = 1, blksz
    a_elmt(1, ii) = ddot(ncv_loc, wvec(1, ii), 1, vec(1, ii), 1)
    !if (peinf%master) print *, "before allreduce a_elmt ", a_elmt(1,2)
    call MPI_ALLREDUCE(MPI_IN_PLACE, a_elmt(1, ii), 1, &
                       MPI_DOUBLE, MPI_SUM, peinf%comm, info)
    !if (peinf%master) print *, "after allreduce a_elmt ", a_elmt(1,2)
    ! wvec = wvec - a_elmt(1) * vec
    call daxpy(ncv_loc, -a_elmt(1, ii), vec(1, ii), 1, wvec(1, ii), 1)
    b_elmt(2, ii) = ddot(ncv_loc, wvec(1, ii), 1, wvec(1, ii), 1)
    !if (peinf%master) print *, "before allreduce b_elmt ", b_elmt(2,2)
    call MPI_ALLREDUCE(MPI_IN_PLACE, b_elmt(2, ii), 1, &
                       MPI_DOUBLE, MPI_SUM, peinf%comm, info)
    !if (peinf%master) print *, "before allreduce b_elmt ", b_elmt(2,2)
    b_elmt(2, ii) = sqrt(b_elmt(2, ii))
  end do ! ii
  !
  do k = 2, niter
    !if (abs(b_elmt(k)) .lt. eps) exit
    do ii = 1, blksz
      wvec(1:ncv_loc, ii) = wvec(1:ncv_loc, ii)/b_elmt(k, ii)
      vec(1:ncv_loc, ii) = -vec(1:ncv_loc, ii)*b_elmt(k, ii)
    end do ! ii loop
    ! swap wvec and vec
    tmpvec = vec
    vec = wvec
    wvec = tmpvec
    ! wvec = wvec + H*vec
    !call stopwatch(peinf%master, "Call matvec_isdf")
    call polynomial_matvec_isdf_BLOCK(vec, tmpvec, ncv_loc, polynomial, npoly, isdf_in, &
                                      wfn, cvpair, ncv, sqrtR, &
                                      MC_vec, CMC_vec, tmp_vec, tmp_Hvec, tmp_cvpair, blksz)
    !call stopwatch(peinf%master, "Done matvec_isdf")
    wvec = wvec + tmpvec
    do ii = 1, blksz
      a_elmt(k, ii) = ddot(ncv_loc, wvec(1, ii), 1, vec(1, ii), 1)
      call MPI_ALLREDUCE(MPI_IN_PLACE, a_elmt(k, ii), 1, &
                         MPI_DOUBLE, MPI_SUM, peinf%comm, info)
      if (k + 1 <= niter) then
        call daxpy(ncv_loc, -a_elmt(k, ii), vec(1, ii), 1, wvec(1, ii), 1)
        b_elmt(k + 1, ii) = ddot(ncv_loc, wvec(1, ii), 1, wvec(1, ii), 1)
        call MPI_ALLREDUCE(MPI_IN_PLACE, b_elmt(k + 1, ii), 1, &
                           MPI_DOUBLE, MPI_SUM, peinf%comm, info)
        b_elmt(k + 1, ii) = sqrt(b_elmt(k + 1, ii))
      end if
    end do ! ii loop
  end do
  !if (peinf%master) then
  !  print *, "j  a_elmt  b_elmt "
  !  do ii = 1, blksz
  !  do j = 1, niter
  !    print *, ii, j, a_elmt(j, ii), b_elmt(j, ii)
  !  enddo
  !  enddo
  !endif
  spectra = 0.d0

  do ii = 1, blksz
  do iz = 1, nz
    do j = niter, 1, -1
      enum = b_elmt(j, ii)*b_elmt(j, ii)
      if (j == 1) enum = (1.d0, 0.d0)
      deno = zz(iz, ii) - a_elmt(j, ii) - spectra(iz, ii)
      spectra(iz, ii) = enum/deno
      !if (nz .eq. 1) print *, "check ", j, enum, deno
    end do ! j
    !print *, "spectra ", iz, zz(iz), spectra(iz)
  end do ! iz
  end do ! ii
! deallocate(a_elmt)
! deallocate(b_elmt)

end subroutine

