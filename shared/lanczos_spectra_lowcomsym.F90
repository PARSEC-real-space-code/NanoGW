! Compute the spectra function of Hermitian matrix using
! lanczos algorithm
! what we calculate is
!
!  <vec| 1/(z-H) |vec>
!
! ncv_loc = isdf_in%ncv_sym(nmrep)
subroutine lanczos_spectra_isdf_lowcomsym(v0, v0_norm, ncv_loc, zz, nz, niter, isdf_in, &
                                          ncv, sqrtR, polynomial, npoly, spectra, &
#if defined DCU
                                          !MC_vec, CMC_vec, tmp_vec,
                                          blksz, nmrep, gvec, d_PsiV, d_PsiC, d_Cmtrx, &
                                          d_pvec, d_cvec, d_lcrep, hipblasHandle &
#elif defined _CUDA
                                          !MC_vec, CMC_vec, tmp_vec,
                                          blksz, nmrep, gvec, d_PsiV, d_PsiC, d_Cmtrx, &
                                          d_pvec, d_cvec, d_lcrep, vstep, tmp_vec, streamid &
#else
                                          blksz, nmrep, gvec &
#endif
                                          )

  use typedefs
  use mpi_module
#ifdef DCU
  use hipfort_types
#endif
#ifdef _CUDA
! Need to include a couple of modules:
!   cublas: required to use generic BLAS interface
!   cudafor: required to use CUDA runtime API routines (e.g. cudaDeviceSynchronize)
!            not explicitly required if file has *.cuf suffix
  use cublas
  use cudafor
#endif
  implicit none
  include 'mpif.h'

  integer, intent(in) :: ncv_loc, niter, nz, npoly, blksz, ncv, nmrep
  complex(dp), intent(in) :: zz(nz, blksz)
  type(ISDF), intent(in) :: isdf_in
  real(dp), intent(in) :: v0(ncv_loc, blksz), &
                          polynomial(npoly + 1), v0_norm(blksz), &
                          sqrtR(ncv_loc)
  type(gspace), intent(in) :: gvec
  complex(dp), intent(out):: spectra(nz, blksz)
  !real(dp), target, intent(inout) :: MC_vec(w_grp%myn_intp_r, blksz), &
  !  CMC_vec(ncv_loc, blksz), tmp_vec(isdf_in%n_intp_r, blksz)
#if defined DCU
  type(c_ptr), intent(inout) :: d_PsiV, d_PsiC, d_Cmtrx, d_cvec, &
                                d_pvec, d_lcrep
  type(c_ptr), intent(in) :: hipblasHandle
#elif defined _CUDA
  real(dp), device, intent(in) :: d_PsiV(isdf_in%n_intp_r, *), &
                                  d_PsiC(isdf_in%n_intp_r, *)
  real(dp), device, intent(inout) :: d_Cmtrx(isdf_in%maxmync_sym*vstep, isdf_in%n_intp_r), &
                                     d_cvec(isdf_in%n_intp_r, blksz), &
                                     d_pvec(isdf_in%maxmync_sym*vstep, blksz)
  real(dp), intent(inout) :: tmp_vec(isdf_in%n_intp_r, blksz)
  integer, device, intent(in) :: d_lcrep(*)
  integer, intent(in) :: vstep
  integer(kind=cuda_stream_kind), intent(in) :: streamid
#endif
  real(dp)   :: a_elmt(niter, blksz), b_elmt(niter, blksz)
  complex(8) :: enum, deno
  real(dp), parameter :: eps = 1.0d-8
  integer :: ii, i, j, k, iz, info
#ifndef _CUDA
  real(dp), external :: ddot
#endif
  real(dp) :: vec(ncv_loc, blksz), &
              wvec(ncv_loc, blksz), tmpvec(ncv_loc, blksz)

  !print *, "inside lanczos_spectra_lowcom "
  ! normalize the vector
  do ii = 1, blksz
    vec(1:ncv_loc, ii) = v0(1:ncv_loc, ii)/v0_norm(ii)
  end do

! allocate(a_elmt(niter))
! allocate(b_elmt(niter))

  b_elmt(1, 1:blksz) = 1.d0

  ! compute wvec = H*vec
  !if (peinf%master) print *, "vec = ", vec(1:5,1), "...", &
  !    vec(ncv_loc-4:ncv_loc,1)
  call polynomial_matvec_isdf_lowcomsym(vec, wvec, ncv_loc, &
                                        polynomial, npoly, isdf_in, ncv, sqrtR, &
#if defined DCU
                                        blksz, nmrep, gvec, d_PsiV, d_PsiC, d_Cmtrx, &
                                        d_pvec, d_cvec, d_lcrep, hipblasHandle &
#elif defined _CUDA
                                        blksz, nmrep, gvec, d_PsiV, d_PsiC, d_Cmtrx, &
                                        d_pvec, d_cvec, d_lcrep, vstep, tmp_vec, streamid &
#else
                                        blksz, nmrep, gvec &
#endif
                                        )
  !if (peinf%master) print *, "wvec = ",wvec(1:10, 2), "...", &
  !    vec(ncv_loc-4:ncv_loc,2)
  ! stop
  ! dot(vec, wvec)
  !
  !call stopwatch(peinf%master,'mpi_allreduce 0 start')
  do ii = 1, blksz
    a_elmt(1, ii) = ddot(ncv_loc, wvec(1, ii), 1, vec(1, ii), 1)
    !if (peinf%master) print *, "before allreduce a_elmt ", a_elmt(1,2)
    call mpi_allreduce(MPI_IN_PLACE, a_elmt(1, ii), 1, &
                       MPI_DOUBLE, MPI_SUM, peinf%comm, info)
    !if (peinf%master) print *, "after allreduce a_elmt ", a_elmt(1,2)
    ! wvec = wvec - a_elmt(1) * vec
    call daxpy(ncv_loc, -a_elmt(1, ii), vec(1, ii), 1, wvec(1, ii), 1)
    b_elmt(2, ii) = ddot(ncv_loc, wvec(1, ii), 1, wvec(1, ii), 1)
    !if (peinf%master) print *, "before allreduce b_elmt ", b_elmt(2,2)
    call mpi_allreduce(MPI_IN_PLACE, b_elmt(2, ii), 1, &
                       MPI_DOUBLE, MPI_SUM, peinf%comm, info)
    !if (peinf%master) print *, "before allreduce b_elmt ", b_elmt(2,2)
    b_elmt(2, ii) = sqrt(b_elmt(2, ii))
  end do ! ii
  !call stopwatch(peinf%master,'mpi_allreduce 0 finish')
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
    call polynomial_matvec_isdf_lowcomsym(vec, tmpvec, ncv_loc, &
                                          polynomial, npoly, isdf_in, ncv, sqrtR, &
#if defined DCU
                                          blksz, nmrep, gvec, d_PsiV, d_PsiC, d_Cmtrx, &
                                          d_pvec, d_cvec, d_lcrep, hipblasHandle &
#elif defined _CUDA
                                          blksz, nmrep, gvec, d_PsiV, d_PsiC, d_Cmtrx, &
                                          d_pvec, d_cvec, d_lcrep, vstep, tmp_vec, streamid &
#else
                                          blksz, nmrep, gvec &
#endif
                                          )
    wvec = wvec + tmpvec
    !call stopwatch(peinf%master,'mpi_allreduce k start')
    do ii = 1, blksz
      a_elmt(k, ii) = ddot(ncv_loc, wvec(1, ii), 1, vec(1, ii), 1)
      call mpi_allreduce(MPI_IN_PLACE, a_elmt(k, ii), 1, &
                         MPI_DOUBLE, MPI_SUM, peinf%comm, info)
      if (k + 1 <= niter) then
        call daxpy(ncv_loc, -a_elmt(k, ii), vec(1, ii), 1, wvec(1, ii), 1)
        b_elmt(k + 1, ii) = ddot(ncv_loc, wvec(1, ii), 1, wvec(1, ii), 1)
        call mpi_allreduce(MPI_IN_PLACE, b_elmt(k + 1, ii), 1, &
                           MPI_DOUBLE, MPI_SUM, peinf%comm, info)
        b_elmt(k + 1, ii) = sqrt(b_elmt(k + 1, ii))
      end if
    end do ! ii loop
    !call stopwatch(peinf%master,'mpi_allreduce k finish')
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

