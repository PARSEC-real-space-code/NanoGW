! The casida equation has the structure
!  C Z^s = omega^2 Z^s, where
!  C = (A-B)^(1/2) @ (A+B) @ (A-B)^1/2
! Here "@" denotes matrix-matrix multiplication.
!
! Note R = (A-B) is a diagonal matrix, the matrix elements is
!   R_{cv, c'v'} = delta_cc' delta_vv' (E_c - E_v)
!
! And the matrix A+B = R+4*K, where K is a low-rank matrix constructed using ISDF method
!   K = Cmtrx^T @ Mmtrx @ Cmtrx
! Pay attention that Cmtrx should not be confused with C.
! Cmtrx is the coefficients found in ISDF method.
!   Cmtrx(n_intp_r, icv) = Psi_intp(n_intp_r, ic) x Psi_intp(n_intp_r, iv)
! Matrix dimensions:
!   Mmtrx: n_intp_r * n_intp_r
!   Cmtrx: n_intp_r * (nv*nc)
!   K:      (nv*nc) * (nv*nc)
!   C:      (nv*nc) * (nv*nc)
!   R:      (nv*nc) * (nv*nc)
! In this subroutine, what we compute is
!   Hvec = C @ vec
!        = R^0.5 (R + 4 * Cmtrx^T @ Mmtrx @ Cmtrx) R^0.5 @ vec
!        = R^1.5 @ ( R^0.5 @ vec ) + 4 * ( R^0.5 @ Cmtrx^T @ Mmtrx @ Cmtrx ) @ ( R^0.5 @ vec )
!
!-----------------------------------------------------------------------------------------
!
! With the TDA (tamm_d), we do calculation directly from A, instead of C
!   Hvec = A @ vec
!        = R @ vec + 2 * Cmtrx^T @ Mmtrx @ Cmtrx @ vec

subroutine matvec_isdf_lowcomsym(vec, Hvec, ncv_loc, isdf_in, sqrtR, tamm_d, &
                                 blksz, nmrep, gvec &
#if defined DCU
                                 , d_PsiV, d_PsiC, d_Cmtrx, d_pvec, d_cvec, d_lcrep, &
                                 hipblasHandle &
#elif defined _CUDA
                                 , d_PsiV, d_PsiC, d_Cmtrx, d_pvec, d_cvec, d_lcrep, &
                                 vstep, tmp_vec, streamid &
#endif
                                 )

  use typedefs
  use mpi_module

#ifdef DCU
  use sub_module
  use hipfort
  use hipfort_types
  use hipfort_check
  use hipfort_hipblas
#elif defined _CUDA
  use cublas  ! required to use generic BLAS interface
  use cudafor ! CUDA runtime API routines (e.g. cudaDeviceSynchronize)
#endif

  implicit none

#ifdef MPI
  include 'mpif.h'
#endif

  integer, intent(in) :: ncv_loc, blksz, nmrep
  type(gspace), intent(in) :: gvec
  real(dp), intent(in) :: vec(ncv_loc, blksz)
  type(ISDF), intent(in) :: isdf_in
  real(dp), intent(out) :: Hvec(ncv_loc, blksz)
  real(dp), intent(in) :: sqrtR(ncv_loc)
  logical, intent(in) :: tamm_d
  real(dp), target :: MC_vec(w_grp%myn_intp_r, blksz), &
                      CMC_vec(ncv_loc, blksz)
#ifdef _CUDA
  real(dp), intent(inout) :: tmp_vec(isdf_in%n_intp_r, blksz)
#else
  real(dp), target :: tmp_vec(isdf_in%n_intp_r, blksz)
#endif

#ifdef DCU
  type(c_ptr), intent(inout) :: d_PsiV, d_PsiC, d_Cmtrx, d_cvec, d_pvec, d_lcrep
  type(c_ptr), intent(in) :: hipblasHandle
  real(dp), target :: tmparray(isdf_in%maxmync_sym, blksz)
  real(dp) :: beta
  real(dp), parameter :: d_one = 1.d0, d_zero = 0.d0
  integer(kind(HIPBLAS_OP_N)), parameter :: transt = HIPBLAS_OP_T, transn = HIPBLAS_OP_N
  integer(c_size_t), parameter :: Bytes_double = 8
  integer(c_size_t) :: Nbytes
  real(dp) :: tmp_Cmtrx(isdf_in%n_intp_r, isdf_in%maxmync_sym)
#elif defined _CUDA
  real(dp), device, intent(in) :: d_PsiV(isdf_in%n_intp_r, *), d_PsiC(isdf_in%n_intp_r, *)
  real(dp), device, intent(inout) :: d_Cmtrx(isdf_in%maxmync_sym*vstep, isdf_in%n_intp_r), &
                                     d_cvec(isdf_in%n_intp_r, blksz), d_pvec(isdf_in%maxmync_sym*vstep, blksz)
  integer, device, intent(in) :: d_lcrep(*)
  integer, intent(in) :: vstep
  integer(kind=cuda_stream_kind), intent(in) :: streamid
  real(dp) :: tmparray(isdf_in%maxmync_sym*vstep, blksz)
  integer :: icstart, icend, istat, v_incr, v_incr_end, voffset, ir
  real(dp) :: beta
  integer, device :: d_lvrep(isdf_in%mynv), d_mync_sym
#else
  real(dp) :: tmp_Cmtrx(isdf_in%n_intp_r, isdf_in%maxmync_sym)
#endif

  real(dp) :: tsec(2)
  integer :: iv, ic, intp_end, intp_start, err, n_intp_r, ii, cvstart, jc, jv, jcrep, jvrep
  ! external functions

#ifndef _CUDA
  real(dp), external :: ddot
#endif

  call timacc(58, 1, tsec)

  if (tamm_d) then
    ! if TDA, apply vec <-- vec
    ! if (w_grp%master) print *, "1 vec", vec(1:10,1)
    call timacc(68, 1, tsec)
    Hvec(1:ncv_loc, 1:blksz) = vec(1:ncv_loc, 1:blksz)
    call timacc(68, 2, tsec)
    ! if (w_grp%master) print *, "2 Hvec", Hvec(1:10,1)
  else
    ! if not TDA, apply vec <-- R^1/2 vec
    ! if (w_grp%master) print *, "1 vec", vec(1:10,1)
    call timacc(68, 1, tsec)
    do ii = 1, blksz
      Hvec(1:ncv_loc, ii) = sqrtR(1:ncv_loc)*vec(1:ncv_loc, ii)
    end do
    call timacc(68, 2, tsec)
    ! if (w_grp%master) print *, "2 Hvec", Hvec(1:10,1)
  end if

  ! Todo: implement dependence on spin and kpt index
  ! ----------------------------------------
  !  subroutine dgemv  (   character   TRANS,
  !   integer   M,
  !   integer   N,
  !   double precision  ALPHA,
  !   double precision, dimension(lda,*)    A,
  !   integer   LDA,
  !   double precision, dimension(*)    X,
  !   integer   INCX,
  !   double precision  BETA,
  !   double precision, dimension(*)    Y,
  !   integer   INCY
  !  )
  ! ----------------------------------------
  ! compute C_vec = Cmtrx @ vec
  ! call dgemv('N', isdf_in%n_intp_r, ncv_loc, 1.0d0, &
  !   tmpCmtrx, isdf_in%n_intp_r, Hvec, 1, 0.0d0, tmp_vec, 1)
  n_intp_r = isdf_in%n_intp_r
  intp_start = w_grp%n_intp_start(w_grp%inode)
  intp_end = w_grp%n_intp_end(w_grp%inode)
  tmp_vec(1:n_intp_r, 1:blksz) = 0.d0
  MC_vec = 0.d0
  CMC_vec = 0.d0
  !if (w_grp%master) print *, "breakpoint 1"
  ! d_cvec ( n_intp_r, jnblock )
  ! d_pvec ( mync, jnblock )
  cvstart = 0
  do jvrep = 1, gvec%syms%ntrans
    do jcrep = 1, gvec%syms%ntrans
      if (gvec%syms%prod(jcrep, jvrep) == nmrep) exit
    end do ! jcrep

#ifdef _CUDA
    icstart = isdf_in%lcrep_bound(jcrep, 1)
    icend = isdf_in%lcrep_bound(jcrep, 2)
    d_lvrep = isdf_in%lvrep
    d_mync_sym = isdf_in%mync_sym(jcrep)
    do iv = isdf_in%lvrep_bound(jvrep, 1), isdf_in%lvrep_bound(jvrep, 2), vstep
      if (iv + vstep - 1 <= isdf_in%lvrep_bound(jvrep, 2)) then
        v_incr_end = vstep
      else
        v_incr_end = isdf_in%lvrep_bound(jvrep, 2) - iv + 1
      end if
#else
    do iv = isdf_in%lvrep_bound(jvrep, 1), isdf_in%lvrep_bound(jvrep, 2)
      jv = isdf_in%lvrep(iv)
#endif

#ifdef DCU
      call timacc(69, 1, tsec)
      Nbytes = isdf_in%maxmync_sym*blksz*Bytes_double
      tmparray = 0.d0
      tmparray(1:isdf_in%mync_sym(jcrep), 1:blksz) = Hvec(cvstart + 1:cvstart + isdf_in%mync_sym(jcrep), 1:blksz)
      call hipCheck(hipMemcpy(d_pvec, c_loc(tmparray), &
                              Nbytes, hipMemcpyHostToDevice))
      call hadamard_prod_mat_sym(d_PsiV, d_PsiC, d_Cmtrx, d_lcrep, isdf_in%n_intp_r, jv - 1, isdf_in%mynv, &
                                 isdf_in%mync, isdf_in%maxmync_sym, isdf_in%lcrep_bound(jcrep, 1), &
                                 isdf_in%lcrep_bound(jcrep, 2))
      call timacc(69, 2, tsec)
      call timacc(71, 1, tsec)
      if (cvstart == 0) then
        beta = d_zero
      else
        beta = d_one
      end if
      call hipblasCheck(hipblasDgemm(hipblasHandle, transn, transn, n_intp_r, blksz, isdf_in%mync_sym(jcrep), d_one, &
                                     d_Cmtrx, n_intp_r, d_pvec, isdf_in%maxmync_sym, beta, d_cvec, n_intp_r))
      call timacc(71, 2, tsec)
#elif defined _CUDA
      istat = cudaDeviceSynchronize()
      call timacc(69, 1, tsec)
!$cuf kernel do(3) <<< (1,*,*), (*,16,128), stream=streamid >>>
      do v_incr = 1, v_incr_end
        do ic = icstart, icend
          do ir = 1, n_intp_r
            jc = d_lcrep(ic)
            voffset = (v_incr - 1)*d_mync_sym
            jv = d_lvrep(iv + v_incr - 1)
            d_Cmtrx(voffset + ic - icstart + 1, ir) = &
              d_PsiC(ir, jc)*d_PsiV(ir, jv)
          end do ! ir
        end do ! ic
      end do ! v_incr
      !
      if (cvstart == 0) then
        beta = 0.d0
      else
        beta = 1.d0
      end if
      istat = cudaDeviceSynchronize()
      call timacc(69, 2, tsec)
      call timacc(71, 1, tsec)
      tmparray(1:isdf_in%mync_sym(jcrep)*v_incr_end, 1:blksz) = &
        Hvec(cvstart + 1:cvstart + isdf_in%mync_sym(jcrep)*v_incr_end, 1:blksz)
      d_pvec(1:isdf_in%maxmync_sym*v_incr_end, 1:blksz) = &
        tmparray(1:isdf_in%maxmync_sym*v_incr_end, 1:blksz)
      istat = cudaDeviceSynchronize()
      call cublasDGEMM('T', 'N', n_intp_r, blksz, isdf_in%mync_sym(jcrep)*v_incr_end, 1.d0, d_Cmtrx, &
                       isdf_in%maxmync_sym*vstep, d_pvec, isdf_in%maxmync_sym*vstep, beta, d_cvec, n_intp_r)
      cvstart = cvstart + isdf_in%mync_sym(jcrep)*v_incr_end
      istat = cudaDeviceSynchronize()
      call timacc(71, 2, tsec)
#else
      call timacc(69, 1, tsec)
      do ic = isdf_in%lcrep_bound(jcrep, 1), isdf_in%lcrep_bound(jcrep, 2)
        jc = isdf_in%lcrep(ic)
        tmp_Cmtrx(1:n_intp_r, ic - isdf_in%lcrep_bound(jcrep, 1) + 1) = isdf_in%PsiC_intp_bl(1:n_intp_r, jc, 1, 1)* &
                                                                        isdf_in%PsiV_intp_bl(1:n_intp_r, jv, 1, 1)
      end do ! ic
      call timacc(69, 2, tsec)
      call timacc(71, 1, tsec)
      call dgemm('N', 'N', n_intp_r, blksz, isdf_in%mync_sym(jcrep), 1.d0, tmp_Cmtrx, n_intp_r, &
                 Hvec(cvstart + 1, 1), ncv_loc, 1.d0, tmp_vec(1, 1), n_intp_r)
      call timacc(71, 2, tsec)
#endif

#ifndef _CUDA
      cvstart = cvstart + isdf_in%mync_sym(jcrep)
#endif

    end do ! iv
  end do ! jvrep
#ifdef DCU
  call hipCheck(hipDeviceSynchronize())
  Nbytes = n_intp_r*blksz*Bytes_double
  call hipCheck(hipMemcpy(c_loc(tmp_vec), d_cvec, Nbytes, hipMemcpyDeviceToHost))
#elif defined _CUDA
  tmp_vec(1:n_intp_r, 1:blksz) = d_cvec(1:n_intp_r, 1:blksz)
#endif

  call timacc(70, 1, tsec)
  !call stopwatch(peinf%master,'matvec: mpi_allreduce start')
  !if ( w_grp%npes .gt. 1) then
  call MPI_allreduce(MPI_IN_PLACE, tmp_vec, isdf_in%n_intp_r*blksz, MPI_DOUBLE, MPI_SUM, peinf%comm, err)
  !endif
  !call stopwatch(peinf%master,'matvec: mpi_allreduce finish')
  call timacc(70, 2, tsec)
  !if (w_grp%master) print *, "Inside matvec_isdf_UltraLowMem, tmp_vec", tmp_vec(1:5)
  !stop
  ! compute MC_vec = Mmtrx @ tmp_vec
  !  Mtrx is distributed in columns (n_intp_r, n_intp_loc)
  !  Mtrx shape:
  !   allocate(isdf_in%Mmtrx( n_intp_r, n_intp_r, nspin, nspin, kpt%nk, 2, gvec%syms%ntrans ))
  !  Calculate Mmtrx(1) = <zeta(i)| V_coul(r,r') | zeta(j)> and
  !            Mmtrx(2) = <zeta(i)| f^LDA(r) | zeta(j)>
  !  print *, "Cmtrx R^1/2 vec", tmp_vec(1:10), sum(tmp_vec)

  !call dgemv('N', w_grp%myn_intp_r, isdf_in%n_intp_r, 1.0d0, &
  !  isdf_in%Mmtrx_loc(1,1,1,1,1,1,1), &
  !  w_grp%myn_intp_r, tmp_vec, 1, 0.0d0, MC_vec, 1)
  call timacc(72, 1, tsec)
  !if (w_grp%master) print *, "breakpoint 3"
  call dgemm('N', 'N', w_grp%myn_intp_r, blksz, isdf_in%n_intp_r, 1.d0, isdf_in%Mmtrx_loc(1, 1, 1, 1, 1, 1, nmrep), &
             w_grp%myn_intp_r, tmp_vec, isdf_in%n_intp_r, 0.d0, MC_vec, w_grp%myn_intp_r)
  call timacc(72, 2, tsec)
  !print *, "Mmtrx(1,1:n_intp_r)", isdf_in%Mmtrx(1,1:isdf_in%n_intp_r,1,1,1,1,1)
  !print *, "tmp_vec(1:n_intp_r)", tmp_vec(1:isdf_in%n_intp_r)
  !print *, "Mmtrx(1,1:n_intp_r)*tmp_vec(1:n_intp_r)", &
  ! sum(isdf_in%Mmtrx(1,1:isdf_in%n_intp_r,1,1,1,1,1)*tmp_vec(1:isdf_in%n_intp_r))
  call timacc(70, 1, tsec)
  !call mpi_barrier(peinf%comm, err)
  !call stopwatch(peinf%master,'matvec: mpi_allreduce start')
  tmp_vec(1:n_intp_r, 1:blksz) = 0.d0
  tmp_vec(intp_start:intp_end, 1:blksz) = MC_vec(1:w_grp%myn_intp_r, 1:blksz)
  call MPI_allreduce(MPI_IN_PLACE, tmp_vec, isdf_in%n_intp_r*blksz, MPI_DOUBLE, MPI_SUM, peinf%comm, err)
  !call stopwatch(peinf%master,'matvec: mpi_allreduce finish')
  call timacc(70, 2, tsec)
  !if (w_grp%master) print *, "Mmtrx Cmtrx R^1/2 vec", MC_vec(1:10)
  ! apply CMC_vec = Cmtrx^T MC_vec
  ! d_cvec ( n_intp_r, jnblock )
  ! d_pvec ( mync, jnblock )

#ifdef DCU
  Nbytes = n_intp_r*blksz*Bytes_double
  call hipCheck(hipMemcpy(d_cvec, c_loc(tmp_vec), Nbytes, hipMemcpyHostToDevice))
#elif defined _CUDA
  d_cvec(1:n_intp_r, 1:blksz) = tmp_vec(1:n_intp_r, 1:blksz)
#endif

  cvstart = 0
  do jvrep = 1, gvec%syms%ntrans
    do jcrep = 1, gvec%syms%ntrans
      if (gvec%syms%prod(jcrep, jvrep) == nmrep) exit
    end do ! jcrep

#ifdef _CUDA
    icstart = isdf_in%lcrep_bound(jcrep, 1)
    icend = isdf_in%lcrep_bound(jcrep, 2)
    d_lvrep = isdf_in%lvrep
    d_mync_sym = isdf_in%mync_sym(jcrep)
    do iv = isdf_in%lvrep_bound(jvrep, 1), isdf_in%lvrep_bound(jvrep, 2), vstep
      if (iv + vstep - 1 <= isdf_in%lvrep_bound(jvrep, 2)) then
        v_incr_end = vstep
      else
        v_incr_end = isdf_in%lvrep_bound(jvrep, 2) - iv + 1
      end if
#else
    do iv = isdf_in%lvrep_bound(jvrep, 1), isdf_in%lvrep_bound(jvrep, 2)
      jv = isdf_in%lvrep(iv)
#endif

#ifdef DCU
      call timacc(69, 1, tsec)
      call hadamard_prod_mat_sym(d_PsiV, d_PsiC, d_Cmtrx, d_lcrep, isdf_in%n_intp_r, jv - 1, isdf_in%mynv, &
                                 isdf_in%mync, isdf_in%maxmync_sym, isdf_in%lcrep_bound(jcrep, 1), &
                                 isdf_in%lcrep_bound(jcrep, 2))
      call timacc(69, 2, tsec)
      call timacc(71, 1, tsec)
      call hipblasCheck(hipblasDgemm(hipblasHandle, transt, transn, isdf_in%mync_sym(jcrep), blksz, n_intp_r, d_one, &
                                     d_Cmtrx, n_intp_r, d_cvec, n_intp_r, d_zero, d_pvec, isdf_in%maxmync_sym))
      call hipCheck(hipDeviceSynchronize())
      Nbytes = isdf_in%maxmync_sym*blksz*Bytes_double
      call hipCheck(hipMemcpy(c_loc(tmparray), d_pvec, Nbytes, hipMemcpyDeviceToHost))
      CMC_vec(cvstart + 1:cvstart + isdf_in%mync_sym(jcrep), 1:blksz) = tmparray(1:isdf_in%mync_sym(jcrep), 1:blksz)
      call timacc(71, 2, tsec)
#elif defined _CUDA
      istat = cudaDeviceSynchronize()
      call timacc(69, 1, tsec)
!$cuf kernel do(3) <<< (1,*,*), (*,16,128), stream=streamid >>>
      do v_incr = 1, v_incr_end
        do ic = icstart, icend
          do ir = 1, n_intp_r
            jc = d_lcrep(ic)
            voffset = (v_incr - 1)*d_mync_sym
            jv = d_lvrep(iv + v_incr - 1)
            d_Cmtrx(voffset + ic - icstart + 1, ir) = d_PsiC(ir, jc)* d_PsiV(ir, jv)
          end do ! ir
        end do ! ic
      end do ! v_incr
      istat = cudaDeviceSynchronize()
      call timacc(69, 2, tsec)
      !
      call timacc(71, 1, tsec)
      istat = cudaDeviceSynchronize()
      call cublasDGEMM('N', 'N', v_incr_end*isdf_in%mync_sym(jcrep), blksz, n_intp_r, 1.d0, d_Cmtrx, &
                       vstep*isdf_in%maxmync_sym, d_Cvec, n_intp_r, 0.d0, d_pvec, vstep*isdf_in%maxmync_sym)
      istat = cudaDeviceSynchronize()
      CMC_vec(cvstart + 1:cvstart + isdf_in%mync_sym(jcrep)*v_incr_end, 1:blksz) = &
        d_pvec(1:isdf_in%mync_sym(jcrep)*v_incr_end, 1:blksz)
      call timacc(71, 2, tsec)
      cvstart = cvstart + isdf_in%mync_sym(jcrep)*v_incr_end
#else
      call timacc(69, 1, tsec)
      do ic = isdf_in%lcrep_bound(jcrep, 1), isdf_in%lcrep_bound(jcrep, 2)
        jc = isdf_in%lcrep(ic)
        tmp_Cmtrx(1:n_intp_r, ic - isdf_in%lcrep_bound(jcrep, 1) + 1) = &
          isdf_in%PsiC_intp_bl(1:n_intp_r, jc, 1, 1)* &
          isdf_in%PsiV_intp_bl(1:n_intp_r, jv, 1, 1)
      end do
      call timacc(69, 2, tsec)
      call timacc(71, 1, tsec)
      call dgemm('T', 'N', isdf_in%mync_sym(jcrep), blksz, n_intp_r, 1.d0, tmp_Cmtrx, n_intp_r, tmp_vec, n_intp_r, &
                 0.d0, CMC_vec(cvstart + 1, 1), ncv_loc)
      call timacc(71, 2, tsec)
#endif

#ifndef _CUDA
      cvstart = cvstart + isdf_in%mync_sym(jcrep)
#endif

    end do ! iv
  end do ! jvrep

  !if (w_grp%master) print *, "breakpoint 4"
  call timacc(73, 1, tsec)
  if (tamm_d) then
    ! Hvec = R @ vec + 2 * Cmtrx^T @ Mmtrx @ Cmtrx @ vec
    do ii = 1, blksz
      !   Hvec <- R * Hvec
      Hvec(1:ncv_loc, ii) = sqrtR(1:ncv_loc)**2*Hvec(1:ncv_loc, ii)
      !   Hvec <- Hvec + 2 * CMC_vec
      Hvec(1:ncv_loc, ii) = Hvec(1:ncv_loc, ii) + 2.d0*CMC_vec(1:ncv_loc, ii)
    end do ! blksz
  else
    ! Hvec = R^1.5 @ ( R^0.5 @ vec ) + 4 * R^0.5 @ ( Cmtrx^T @ Mmtrx @ Cmtrx @ R^0.5 @ vec )
    do ii = 1, blksz
      !   Hvec <- R^3/2 * Hvec
      Hvec(1:ncv_loc, ii) = sqrtR(1:ncv_loc)**3*Hvec(1:ncv_loc, ii)
      !   Hvec <- Hvec + 4 * R^1/2 * CMC_vec
      Hvec(1:ncv_loc, ii) = Hvec(1:ncv_loc, ii) + 4.d0*sqrtR(1:ncv_loc)*CMC_vec(1:ncv_loc, ii)
    end do ! blksz
  end if
  call timacc(73, 2, tsec)
  call timacc(58, 2, tsec)

end subroutine matvec_isdf_lowcomsym


subroutine polynomial_matvec_isdf_lowcomsym(vec, Hvec, ncv_loc, polynomial, npoly, isdf_in, sqrtR, &
                                            blksz, nmrep, gvec &
#ifdef DCU
                                            , d_PsiV, d_PsiC, d_Cmtrx, d_pvec, d_cvec, d_lcrep, &
                                            hipblasHandle &
#elif defined _CUDA
                                            , d_PsiV, d_PsiC, d_Cmtrx, d_pvec, d_cvec, d_lcrep, &
                                            vstep, tmp_vec, streamid &
#endif
                                            )

  use typedefs
  use mpi_module

#ifdef DCU
  use hipfort_types
#elif defined _CUDA
  ! Need to include a couple of modules:
  !   cublas: required to use generic BLAS interface
  !   cudafor: required to use CUDA runtime API routines (e.g. cudaDeviceSynchronize)
  !            not explicitly required if file has *.cuf suffix
  use cublas
  use cudafor
#endif

  implicit none
#ifdef MPI
  include 'mpif.h'
#endif

  integer, intent(in) :: ncv_loc, npoly, blksz, nmrep
  type(ISDF), intent(in) :: isdf_in
  real(dp), intent(in) :: polynomial(npoly + 1), sqrtR(ncv_loc), vec(ncv_loc, blksz)
  type(gspace), intent(in) :: gvec
  ! polynomial is an array that stores the coefficients from higher order to lower order
  ! polynomial = [a_n, a_{n-1}, a_{n-2}, ..., a_0]
  ! polynomial(1) = a_n
  ! polynomial(2) = a_{n-1}
  ! ...
  ! polynomial(npoly+1) = a_0
  real(dp), intent(out) :: Hvec(ncv_loc, blksz)
  !real(dp), target, intent(inout) :: MC_vec(w_grp%myn_intp_r, blksz), &
  !  CMC_vec(ncv_loc, blksz), &
  !  tmp_vec(isdf_in%n_intp_r, blksz)

#ifdef DCU
  type(c_ptr), intent(inout) :: d_PsiV, d_PsiC, d_Cmtrx, d_cvec, d_pvec, d_lcrep
  type(c_ptr), intent(in) :: hipblasHandle
#elif defined _CUDA
  real(dp), device, intent(in) :: d_PsiV(isdf_in%n_intp_r, *), d_PsiC(isdf_in%n_intp_r, *)
  real(dp), device, intent(inout) :: d_Cmtrx(isdf_in%maxmync_sym*vstep, isdf_in%n_intp_r), &
                                     d_cvec(isdf_in%n_intp_r, blksz), d_pvec(isdf_in%maxmync_sym*vstep, blksz)
  real(dp), intent(inout) :: tmp_vec(isdf_in%n_intp_r, blksz)
  integer, device, intent(in) :: d_lcrep(*)
  integer, intent(in) :: vstep
  integer(kind=cuda_stream_kind), intent(in) :: streamid
#endif

  real(dp) :: Hvec_old(ncv_loc, blksz), tsec(2)
  integer :: k

  ! isdf_in%n_intp_r == w_grp%n_intp_r
  ! allocate(MC_vec(w_grp%myn_intp_r))
  ! allocate(CMC_vec(w_grp%myncv))
  ! allocate(tmp_vec(isdf_in%n_intp_r))
  ! allocate(tmp_Hvec(w_grp%ldncv))
  !   (a_n H^n + a_{n-1} H^{n-1} + ... + a_1 H + a_0 ) vec
  ! = (H*...(H*(H*(a_n H*vec + a_{n-1} vec) + a_{n-2} vec) + a_{n-3} vec)...a_1 vec) + a_0 vec
  ! Recursively:
  !   Hvec <- a_k H*vec + a_{k-1}

  Hvec_old = vec
  do k = 1, npoly
    call matvec_isdf_lowcomsym(Hvec_old, Hvec, ncv_loc, isdf_in, sqrtR, .false., &
                               blksz, nmrep, gvec &
#ifdef DCU
                               , d_PsiV, d_PsiC, d_Cmtrx, d_pvec, d_cvec, d_lcrep, &
                               hipblasHandle &
#elif defined _CUDA
                               , d_PsiV, d_PsiC, d_Cmtrx, d_pvec, d_cvec, d_lcrep, &
                               vstep, tmp_vec, streamid &
#endif
                               )
    !if (peinf%master) print *, "k =", k, polynomial(k)
    call timacc(57, 1, tsec)
    if (k == 1) then
      Hvec = polynomial(k)*Hvec + polynomial(k + 1)*vec
      !if (peinf%master) print *, "In polynomial_matvec_isdf: Hvec", Hvec(1:10)
    else
      Hvec = Hvec + polynomial(k + 1)*vec
      !if (peinf%master) print *, "In polynomial_matvec_isdf: Hvec", Hvec(1:10)
      !stop
    end if
    Hvec_old = Hvec
    call timacc(57, 2, tsec)
  end do ! k

end subroutine polynomial_matvec_isdf_lowcomsym
