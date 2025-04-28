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
!        = R^1/2 (R+4*Cmtrx^T @ Mmtrx @ Cmtrx) R^1/2 @ vec
!        = R^2 @ vec + 4*R^1/2 @ Cmtrx^T @ Mmtrx @ Cmtrx @ R^1/2 @ vec
!
subroutine matvec_isdf_UltraLowMem(vec, Hvec, ncv_loc, isdf_in, wfn, cvpair, sqrtR, &
                                   MC_vec, CMC_vec, tmp_vec, tmp_Hvec, tmp_cvpair)

  use typedefs
  use mpi_module
  implicit none
  include 'mpif.h'

  integer, intent(in) :: ncv_loc
  real(dp), intent(in) :: vec(ncv_loc)
  type(ISDF), intent(in) :: isdf_in
  real(dp), intent(out) :: Hvec(ncv_loc)
  real(dp), intent(in) :: sqrtR(ncv_loc)
  integer, intent(in) :: cvpair(2, ncv_loc)
  type(wavefunction), intent(in) :: wfn
  real(dp), intent(inout) :: MC_vec(w_grp%myn_intp_r), &
                             CMC_vec(w_grp%myncv), &
                             tmp_vec(isdf_in%n_intp_r), tmp_Hvec(w_grp%ldncv)
  integer, intent(inout) :: tmp_cvpair(2, w_grp%ldncv)
  integer :: iv, ic, icv, incr, intp_end, intp_start, &
             res, err, n_intp_loc, ncvdim, ipe, ncvpair
  ! external functions
  real(dp), external :: ddot

! apply vec <-- R^1/2 vec
  !if (w_grp%master) print *, "1 vec", vec(1:10)
  Hvec = sqrtR*vec
  !if (w_grp%master) print *, "2 Hvec", Hvec(1:10)

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
  intp_start = w_grp%n_intp_start(w_grp%inode)
  intp_end = w_grp%n_intp_end(w_grp%inode)
  tmp_vec = 0.d0
  MC_vec = 0.d0
  CMC_vec = 0.d0
  do ipe = 0, w_grp%npes - 1
    ncvdim = w_grp%ncv_end(ipe) - w_grp%ncv_start(ipe) + 1
    if (w_grp%inode == ipe) then
      tmp_Hvec(1:ncvdim) = Hvec(1:ncvdim)
      tmp_cvpair(1:2, 1:ncvdim) = cvpair(1:2, 1:ncvdim)
    end if
    if (w_grp%npes > 1) then
      call MPI_BCAST(tmp_Hvec, w_grp%ldncv, MPI_DOUBLE, ipe, &
                     w_grp%comm, err)
      call MPI_BCAST(tmp_cvpair, 2*w_grp%ldncv, MPI_INTEGER, ipe, &
                     w_grp%comm, err)
    end if
    do icv = 1, ncvdim
      iv = tmp_cvpair(1, icv)
      ic = tmp_cvpair(2, icv)
      tmp_vec(intp_start:intp_end) = &
        tmp_vec(intp_start:intp_end) + &
        isdf_in%dPsi_intp_loc(1:w_grp%myn_intp_r, ic, 1, 1)* &
        isdf_in%dPsi_intp_loc(1:w_grp%myn_intp_r, iv, 1, 1)* &
        tmp_Hvec(icv)
      ! if (w_grp%master) print *, iv, ic, tmp_vec(1:5)
    end do ! icv
  end do ! ipe
  if (w_grp%npes > 1) then
    call MPI_ALLREDUCE(MPI_IN_PLACE, tmp_vec, isdf_in%n_intp_r, MPI_DOUBLE, MPI_SUM, w_grp%comm, err)
  end if
  !if (w_grp%master) print *, "Inside matvec_isdf_UltraLowMem, tmp_vec", tmp_vec(1:5)
  !stop
! compute MC_vec = Mmtrx @ tmp_vec
!  Mtrx is distributed in columns (n_intp_r, n_intp_loc)
!  Mtrx shape:
!   allocate(isdf_in%Mmtrx( n_intp_r, n_intp_r, nspin, nspin, kpt%nk, 2, gvec%syms%ntrans ))
!  Calculate Mmtrx(1) = <zeta(i)| V_coul(r,r') | zeta(j)> and
!            Mmtrx(2) = <zeta(i)| f^LDA(r) | zeta(j)>
!  print *, "Cmtrx R^1/2 vec", tmp_vec(1:10), sum(tmp_vec)

  call dgemv('N', w_grp%myn_intp_r, isdf_in%n_intp_r, 1.0d0, &
             isdf_in%dMmtrx_loc(1, 1, 1, 1, 1, 1, 1), &
             w_grp%myn_intp_r, tmp_vec, 1, 0.0d0, MC_vec, 1)
  !print *, "Mmtrx(1,1:n_intp_r)", isdf_in%Mmtrx(1,1:isdf_in%n_intp_r,1,1,1,1,1)
  !print *, "tmp_vec(1:n_intp_r)", tmp_vec(1:isdf_in%n_intp_r)
  !print *, "Mmtrx(1,1:n_intp_r)*tmp_vec(1:n_intp_r)", &
  ! sum(isdf_in%Mmtrx(1,1:isdf_in%n_intp_r,1,1,1,1,1)*tmp_vec(1:isdf_in%n_intp_r))
  !tmp_vec = 0.d0
  !tmp_vec(intp_start:intp_end) = MC_vec(1:n_intp_loc)
  !call MPI_ALLREDUCE( MPI_IN_PLACE, tmp_vec, isdf_in%n_intp_r, MPI_DOUBLE, MPI_SUM, &
  !  peinf%comm, err )
  !if (w_grp%master) print *, "Mmtrx Cmtrx R^1/2 vec", MC_vec(1:10)

! apply CMC_vec = Cmtrx^T MC_vec
  do ipe = 0, w_grp%npes - 1
    !call dgemv('T', isdf_in%n_intp_r_loc, ncv_loc, 1.0d0, tmpCmtrx, isdf_in%n_intp_r, tmp_vec, 1, 0.0d0, CMC_vec, 1)
    !print *, "Cmtrx^T Mmtrx Cmtrx R^1/2 vec", CMC_vec(1:5)
    ncvdim = w_grp%ncv_end(ipe) - w_grp%ncv_start(ipe) + 1
    if (w_grp%inode == ipe) then
      tmp_cvpair(1:2, 1:ncvdim) = cvpair(1:2, 1:ncvdim)
    end if
    if (w_grp%npes > 1) then
      call MPI_BCAST(tmp_cvpair, 2*ncvdim, MPI_INTEGER, &
                     ipe, w_grp%comm, err)
    end if
    tmp_Hvec = 0.d0
    do icv = 1, ncvdim
      iv = tmp_cvpair(1, icv)
      ic = tmp_cvpair(2, icv)
      tmp_vec(1:w_grp%myn_intp_r) = &
        isdf_in%dPsi_intp_loc(1:w_grp%myn_intp_r, ic, 1, 1)* &
        isdf_in%dPsi_intp_loc(1:w_grp%myn_intp_r, iv, 1, 1)
      tmp_Hvec(icv) = &
        ddot(w_grp%myn_intp_r, tmp_vec(1), 1, MC_vec(1), 1)
    end do ! icv
    call MPI_REDUCE(tmp_Hvec(1), CMC_vec(1), ncvdim, MPI_DOUBLE, MPI_SUM, &
                    ipe, w_grp%comm, err)
    !if (w_grp%master) print *, w_grp%inode, "tmp_Hvec", tmp_Hvec(1:5)
    !if (w_grp%master) print *, w_grp%inode, "CMC_vec",  CMC_vec(1:5)
  end do ! ipe
  !if (w_grp%master) print *, "Inside apply_matvec_UltraLowMem, CMC_vec ", CMC_vec(1:5)
  !stop

! Hvec <- R^3/2 * Hvec
  Hvec = sqrtR**3*Hvec
  !print *, "4 Hvec", Hvec(1:5)

! Hvec <- Hvec + R^1/2 * CMC_vec
  Hvec = Hvec + 4.d0*sqrtR*CMC_vec
  !print *, "5 Hvec", Hvec(1:5)

end subroutine matvec_isdf_UltraLowMem

subroutine polynomial_matvec_isdf_UltraLowMem(vec, Hvec, ncv_loc, polynomial, npoly, isdf_in, wfn, cvpair, sqrtR, &
                                              MC_vec, CMC_vec, tmp_vec, tmp_Hvec, tmp_cvpair)

  use typedefs
  use mpi_module
  implicit none

  integer, intent(in) :: ncv_loc, npoly
  type(ISDF), intent(in) :: isdf_in
  real(dp), intent(in) :: polynomial(npoly + 1), &
                          sqrtR(ncv_loc)
  real(dp), intent(in) :: vec(ncv_loc)
  type(wavefunction), intent(in) :: wfn
  integer, intent(in) :: cvpair(2, ncv_loc)
  ! polynomial is an array that stores the coefficients from higher order to lower order
  ! polynomial = [a_n, a_{n-1}, a_{n-2}, ..., a_0]
  ! polynomial(1) = a_n
  ! polynomial(2) = a_{n-1}
  ! ...
  ! polynomial(npoly+1) = a_0
  real(dp), intent(out) :: Hvec(ncv_loc)
  real(dp), intent(inout) :: MC_vec(w_grp%myn_intp_r), CMC_vec(w_grp%myncv), &
                             tmp_vec(isdf_in%n_intp_r), tmp_Hvec(w_grp%ldncv)
  integer, intent(inout) :: tmp_cvpair(2, w_grp%ldncv)
  real(dp) :: Hvec_old(ncv_loc)
  integer  :: k

! isdf_in%n_intp_r == w_grp%n_intp_r
! allocate(MC_vec(w_grp%myn_intp_r))
! allocate(CMC_vec(w_grp%myncv))
! allocate(tmp_vec(isdf_in%n_intp_r))
! allocate(tmp_Hvec(w_grp%ldncv))
! allocate(tmp_cvpair(2, w_grp%ldncv))
  !   (a_n H^n + a_{n-1} H^{n-1} + ... + a_1 H + a_0 ) vec
  ! = (H*...(H*(H*(a_n H*vec + a_{n-1} vec) + a_{n-2} vec) + a_{n-3} vec)...a_1 vec) + a_0 vec
  ! Recursively:
  !   Hvec <- a_k H*vec + a_{k-1}
  do k = 1, npoly
    !if (peinf%master) print *, "k =", k, polynomial(k)
    if (k == 1) then
      call matvec_isdf_UltraLowMem(vec, Hvec, ncv_loc, isdf_in, wfn, cvpair, &
                                   sqrtR, MC_vec, CMC_vec, tmp_vec, tmp_Hvec, tmp_cvpair)
      Hvec = polynomial(k)*Hvec + polynomial(k + 1)*vec
      !if (peinf%master) print *, "In polynomial_matvec_isdf: Hvec", Hvec(1:10)
    else
      call matvec_isdf_UltraLowMem(Hvec_old, Hvec, ncv_loc, isdf_in, wfn, cvpair, &
                                   sqrtR, MC_vec, CMC_vec, tmp_vec, tmp_Hvec, tmp_cvpair)
      Hvec = Hvec + polynomial(k + 1)*vec
      !if (peinf%master) print *, "In polynomial_matvec_isdf: Hvec", Hvec(1:10)
      !stop
    end if
    Hvec_old = Hvec
  end do ! k

end subroutine polynomial_matvec_isdf_UltraLowMem

