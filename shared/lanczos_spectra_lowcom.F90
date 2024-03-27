! Compute the spectra function of Hermitian matrix using 
! lanczos algorithm
! what we calculate is
!   
!  <vec| 1/(z-H) |vec> 
!
subroutine lanczos_spectra_isdf_lowcom ( v0, v0_norm, ncv_loc, zz, nz, niter, isdf_in, &
  ncv, sqrtR, polynomial, npoly, spectra, &
#ifdef DCU
  !MC_vec, CMC_vec, tmp_vec, 
  blksz, d_PsiV, d_PsiC, d_Cmtrx, d_pvec, d_cvec, &
  hipblasHandle)
#else
  !MC_vec, CMC_vec, tmp_vec, 
  blksz)
#endif
  
  use typedefs
  use mpi_module
#ifdef DCU
  use hipfort_types
#endif
  implicit none
  include 'mpif.h'

  integer, intent(in) :: ncv_loc, niter, nz, npoly, blksz, ncv
  complex(dp), intent(in) :: zz(nz, blksz)
  type(ISDF), intent(in) :: isdf_in
  real(dp), intent(in) :: v0(ncv_loc, blksz), &
   polynomial(npoly+1), v0_norm(blksz), &
   sqrtR(ncv_loc)
  complex(dp), intent(out):: spectra(nz, blksz)
  !real(dp), target, intent(inout) :: MC_vec(w_grp%myn_intp_r, blksz), &
  !  CMC_vec(ncv_loc, blksz), tmp_vec(isdf_in%n_intp_r, blksz)
#ifdef DCU
  type(c_ptr), intent(inout) :: d_PsiV, d_PsiC, d_Cmtrx, d_cvec, d_pvec
  type(c_ptr), intent(in) :: hipblasHandle
#endif
  real(dp)   :: a_elmt(niter, blksz), b_elmt(niter, blksz)
  complex(8) :: enum, deno
  real(dp), parameter :: eps = 1.0d-8
  integer :: ii, i, j, k, iz, info
  real(dp), external :: ddot
  real(dp) :: vec(ncv_loc, blksz), &
    wvec(ncv_loc, blksz), tmpvec(ncv_loc, blksz)

  !print *, "inside lanczos_spectra_lowcom "
  ! normalize the vector 
  do ii = 1, blksz
  vec(1:ncv_loc, ii) = v0(1:ncv_loc, ii)/v0_norm(ii)
  enddo

! allocate(a_elmt(niter))
! allocate(b_elmt(niter))

  b_elmt(1, 1:blksz) = 1.d0

  ! compute wvec = H*vec
  !if (peinf%master) print *, "vec = ", vec(1:5,1), "...", &
  !    vec(ncv_loc-4:ncv_loc,1)

  call stopwatch(peinf%master, "Call matvec_isdf_lowcom")
  call polynomial_matvec_isdf_lowcom (vec, wvec, ncv_loc, polynomial, npoly, &
     isdf_in, ncv, sqrtR, &
#ifdef DCU
     !MC_vec, CMC_vec, tmp_vec, &
     blksz, d_PsiV, d_PsiC, d_Cmtrx, d_pvec, d_cvec, hipblasHandle &
#else
     !MC_vec, CMC_vec, tmp_vec, &
     blksz &
#endif
   )
  !call stopwatch(peinf%master, "Done matvec_isdf")
  !if (peinf%master) print *, "wvec = ",wvec(1:10, 2), "...", &
  !    vec(ncv_loc-4:ncv_loc,2)
  ! stop
  ! dot(vec, wvec)
  !
  do ii = 1, blksz
  a_elmt(1, ii) = ddot(ncv_loc, wvec(1,ii), 1, vec(1,ii), 1)
  !if (peinf%master) print *, "before allreduce a_elmt ", a_elmt(1,2)
  call mpi_allreduce(MPI_IN_PLACE, a_elmt(1, ii), 1, &
    MPI_DOUBLE, MPI_SUM, peinf%comm, info)
  !if (peinf%master) print *, "after allreduce a_elmt ", a_elmt(1,2)
  ! wvec = wvec - a_elmt(1) * vec
  call daxpy(ncv_loc, -a_elmt(1, ii), vec(1, ii), 1, wvec(1, ii), 1)
  b_elmt(2, ii) = ddot(ncv_loc , wvec(1, ii), 1, wvec(1, ii), 1)
  !if (peinf%master) print *, "before allreduce b_elmt ", b_elmt(2,2)
  call mpi_allreduce(MPI_IN_PLACE, b_elmt(2, ii), 1, &
    MPI_DOUBLE, MPI_SUM, peinf%comm, info)
  !if (peinf%master) print *, "before allreduce b_elmt ", b_elmt(2,2)
  b_elmt(2, ii) = sqrt(b_elmt(2, ii))
  enddo ! ii
  !
  do k = 2, niter
    !if (abs(b_elmt(k)) .lt. eps) exit
    do ii = 1, blksz
    wvec(1:ncv_loc, ii) = wvec(1:ncv_loc, ii)/b_elmt(k, ii)
    vec(1:ncv_loc, ii)  = -vec(1:ncv_loc, ii)*b_elmt(k, ii)
    enddo ! ii loop
    ! swap wvec and vec
    tmpvec = vec 
    vec = wvec
    wvec = tmpvec
    ! wvec = wvec + H*vec
    !call stopwatch(peinf%master, "Call matvec_isdf") 
    call polynomial_matvec_isdf_lowcom (vec, tmpvec, ncv_loc, &
      polynomial, npoly, isdf_in, ncv, sqrtR, &
#ifdef DCU
      !MC_vec, CMC_vec, tmp_vec, &
      blksz, d_PsiV, d_PsiC, d_Cmtrx, d_pvec, d_cvec, hipblasHandle &
#else
      !MC_vec, CMC_vec, tmp_vec, &
      blksz &
#endif
    )
    !call stopwatch(peinf%master, "Done matvec_isdf")
    wvec = wvec + tmpvec
    do ii = 1, blksz
    a_elmt(k, ii) = ddot(ncv_loc, wvec(1, ii), 1, vec(1, ii), 1)
    call mpi_allreduce(MPI_IN_PLACE, a_elmt(k, ii), 1,&
      MPI_DOUBLE, MPI_SUM, peinf%comm, info)
    if (k+1 .le. niter) then
      call daxpy(ncv_loc, -a_elmt(k, ii), vec(1, ii), 1, wvec(1, ii), 1)
      b_elmt(k+1, ii) = ddot(ncv_loc, wvec(1, ii), 1, wvec(1, ii), 1)
      call mpi_allreduce(MPI_IN_PLACE, b_elmt(k+1, ii), 1, &
        MPI_DOUBLE, MPI_SUM, peinf%comm, info)
      b_elmt(k+1, ii) = sqrt(b_elmt(k+1, ii))
    endif
    enddo ! ii loop
  enddo
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
      if (j .eq. 1) enum = (1.d0, 0.d0)
      deno = zz(iz, ii)-a_elmt(j, ii)-spectra(iz, ii)
      spectra(iz, ii) = enum/deno
      !if (nz .eq. 1) print *, "check ", j, enum, deno
    enddo ! j
    !print *, "spectra ", iz, zz(iz), spectra(iz)
  enddo ! iz
  enddo ! ii
! deallocate(a_elmt)
! deallocate(b_elmt)

end subroutine 

