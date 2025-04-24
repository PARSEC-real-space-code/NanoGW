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
subroutine matvec_isdf_BLOCK(vec, Hvec, ncv_loc, isdf_in, wfn, cvpair, ncv, sqrtR, &
                             MC_vec, CMC_vec, tmp_vec, tmp_Hvec, tmp_cvpair, blksz)

  use typedefs
  use mpi_module
  implicit none
  include 'mpif.h'

  integer, intent(in) :: ncv_loc, blksz, ncv
  real(dp), intent(in) :: vec(ncv_loc, blksz)
  type(ISDF), intent(in) :: isdf_in
  real(dp), intent(out) :: Hvec(ncv_loc, blksz)
  real(dp), intent(in) :: sqrtR(ncv_loc)
  integer, intent(in) :: cvpair(4, ncv)
  type(wavefunction), intent(in) :: wfn
  real(dp), intent(inout) :: MC_vec(w_grp%myn_intp_r, blksz), &
                             CMC_vec(w_grp%ldncv, blksz), &
                             tmp_vec(isdf_in%n_intp_r, blksz), tmp_Hvec(w_grp%ldncv, blksz)
  integer, intent(inout) :: tmp_cvpair(2, w_grp%ldncv)
  real(dp) :: tmp_Cmtrx(w_grp%myn_intp_r, w_grp%ldncv), tsec(2)
  integer :: iv, ic, icv, incr, intp_end, intp_start, &
             res, err, n_intp_loc, ncvdim, ipe, ncvpair, ii
  ! external functions
  real(dp), external :: ddot

! apply vec <-- R^1/2 vec
  !if (w_grp%master) print *, "1 vec", vec(1:10,1)
  call timacc(68, 1, tsec)
  do ii = 1, blksz
    Hvec(1:ncv_loc, ii) = sqrtR(1:ncv_loc)*vec(1:ncv_loc, ii)
  end do
  call timacc(68, 2, tsec)
  !if (w_grp%master) print *, "2 Hvec", Hvec(1:10,1)

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
  !if (w_grp%master) print *, "breakpoint 1"
  do ipe = 0, w_grp%npes - 1
    ncvdim = w_grp%ncv_end(ipe) - w_grp%ncv_start(ipe) + 1
    call timacc(69, 1, tsec)
    if (w_grp%inode == ipe) then
      tmp_Hvec(1:ncvdim, 1:blksz) = Hvec(1:ncvdim, 1:blksz)
      !tmp_cvpair(1:2, 1:ncvdim) = cvpair(1:2, &
      !  w_grp%ncv_start(ipe):w_grp%ncv_end(ipe))
    end if
    if (w_grp%npes > 1) then
      call MPI_BCAST(tmp_Hvec, w_grp%ldncv*blksz, MPI_DOUBLE, ipe, &
                     w_grp%comm, err)
      !call MPI_BCAST( tmp_cvpair, 2*w_grp%ldncv, MPI_INTEGER, ipe, &
      !  w_grp%comm, err )
    end if
    call timacc(69, 2, tsec)
    call timacc(70, 1, tsec)
    do icv = 1, ncvdim
      !iv = tmp_cvpair(1,icv)
      !ic = tmp_cvpair(2,icv)
      iv = cvpair(1, w_grp%ncv_start(ipe) + icv - 1)
      ic = cvpair(2, w_grp%ncv_start(ipe) + icv - 1)
      tmp_Cmtrx(1:w_grp%myn_intp_r, icv) = &
        isdf_in%dPsi_intp_loc(1:w_grp%myn_intp_r, ic, 1, 1)* &
        isdf_in%dPsi_intp_loc(1:w_grp%myn_intp_r, iv, 1, 1)
    end do ! icv
    call timacc(70, 2, tsec)
    call timacc(71, 1, tsec)
    call dgemm('n', 'n', w_grp%myn_intp_r, blksz, ncvdim, 1.d0, &
               tmp_Cmtrx(1, 1), w_grp%myn_intp_r, tmp_Hvec(1, 1), w_grp%ldncv, &
               1.d0, tmp_vec(intp_start, 1), isdf_in%n_intp_r)
    call timacc(71, 2, tsec)
    ! if (w_grp%master) print *, iv, ic, tmp_vec(1:5)
  end do ! ipe
  !if (w_grp%master) print *, "breakpoint 2"
  call timacc(72, 1, tsec)
  if (w_grp%npes > 1) then
    call MPI_allreduce(MPI_IN_PLACE, tmp_vec, isdf_in%n_intp_r*blksz, &
                       MPI_DOUBLE, MPI_SUM, &
                       w_grp%comm, err)
  end if
  call timacc(72, 2, tsec)
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
  call timacc(73, 1, tsec)
  !if (w_grp%master) print *, "breakpoint 3"
  call dgemm('N', 'N', w_grp%myn_intp_r, blksz, isdf_in%n_intp_r, 1.d0, &
             isdf_in%dMmtrx_loc(1, 1, 1, 1, 1, 1, 1), &
             w_grp%myn_intp_r, tmp_vec, isdf_in%n_intp_r, 0.d0, &
             MC_vec, w_grp%myn_intp_r)
  call timacc(73, 2, tsec)
  !print *, "Mmtrx(1,1:n_intp_r)", isdf_in%Mmtrx(1,1:isdf_in%n_intp_r,1,1,1,1,1)
  !print *, "tmp_vec(1:n_intp_r)", tmp_vec(1:isdf_in%n_intp_r)
  !print *, "Mmtrx(1,1:n_intp_r)*tmp_vec(1:n_intp_r)", &
  ! sum(isdf_in%Mmtrx(1,1:isdf_in%n_intp_r,1,1,1,1,1)*tmp_vec(1:isdf_in%n_intp_r))
  !tmp_vec = 0.d0
  !tmp_vec(intp_start:intp_end) = MC_vec(1:n_intp_loc)
  !call MPI_allreduce( MPI_IN_PLACE, tmp_vec, isdf_in%n_intp_r, MPI_DOUBLE, MPI_SUM, &
  !  peinf%comm, err )
  !if (w_grp%master) print *, "Mmtrx Cmtrx R^1/2 vec", MC_vec(1:10)
! apply CMC_vec = Cmtrx^T MC_vec
  do ipe = 0, w_grp%npes - 1
    !call dgemv('T', isdf_in%n_intp_r_loc, ncv_loc, 1.0d0, tmpCmtrx, isdf_in%n_intp_r, tmp_vec, 1, 0.0d0, CMC_vec, 1)
    !print *, "Cmtrx^T Mmtrx Cmtrx R^1/2 vec", CMC_vec(1:5)
    call timacc(74, 1, tsec)
    ncvdim = w_grp%ncv_end(ipe) - w_grp%ncv_start(ipe) + 1
    !if (w_grp%inode .eq. ipe) then
    !  tmp_cvpair(1:2, 1:ncvdim) = &
    !    cvpair(1:2, w_grp%ncv_start(ipe):w_grp%ncv_end(ipe))
    !endif
    !if (w_grp%npes .gt. 1) then
    !call MPI_BCAST(tmp_cvpair, 2*ncvdim, MPI_INTEGER, &
    !     ipe, w_grp%comm, err)
    !endif
    tmp_Hvec = 0.d0
    call timacc(74, 2, tsec)
    call timacc(75, 1, tsec)
    do icv = 1, ncvdim
      !iv = tmp_cvpair(1,icv)
      !ic = tmp_cvpair(2,icv)
      iv = cvpair(1, w_grp%ncv_start(ipe) + icv - 1)
      ic = cvpair(2, w_grp%ncv_start(ipe) + icv - 1)
      tmp_Cmtrx(1:w_grp%myn_intp_r, icv) = &
        isdf_in%dPsi_intp_loc(1:w_grp%myn_intp_r, ic, 1, 1)* &
        isdf_in%dPsi_intp_loc(1:w_grp%myn_intp_r, iv, 1, 1)
    end do ! icv
    call timacc(75, 2, tsec)
    call timacc(76, 1, tsec)
    call dgemm('t', 'n', ncvdim, blksz, w_grp%myn_intp_r, 1.d0, &
               tmp_Cmtrx(1, 1), w_grp%myn_intp_r, MC_vec(1, 1), w_grp%myn_intp_r, &
               0.d0, tmp_Hvec(1, 1), w_grp%ldncv)
    call timacc(76, 2, tsec)
    call timacc(77, 1, tsec)
    call MPI_REDUCE(tmp_Hvec(1, 1), CMC_vec(1, 1), w_grp%ldncv*blksz, &
                    MPI_DOUBLE, MPI_SUM, ipe, w_grp%comm, err)
    call timacc(77, 2, tsec)
    !if (w_grp%master) print *, w_grp%inode, "tmp_Hvec", tmp_Hvec(1:5)
    !if (w_grp%master) print *, w_grp%inode, "CMC_vec",  CMC_vec(1:5)
  end do ! ipe
  !if (w_grp%master) print *, "Inside apply_matvec_UltraLowMem, CMC_vec ", CMC_vec(1:5)
  !stop
  !if (w_grp%master) print *, "breakpoint 4"
  call timacc(78, 1, tsec)
  do ii = 1, blksz
! Hvec <- R^3/2 * Hvec
    Hvec(1:ncv_loc, ii) = sqrtR(1:ncv_loc)**3*Hvec(1:ncv_loc, ii)
    !print *, "4 Hvec", Hvec(1:5)

! Hvec <- Hvec + R^1/2 * CMC_vec
    Hvec(1:ncv_loc, ii) = Hvec(1:ncv_loc, ii) + &
                          4.d0*sqrtR(1:ncv_loc)*CMC_vec(1:ncv_loc, ii)
    !print *, "5 Hvec", Hvec(1:5)
  end do ! blksz
  call timacc(78, 2, tsec)

end subroutine matvec_isdf_BLOCK

subroutine polynomial_matvec_isdf_BLOCK(vec, Hvec, ncv_loc, &
                                        polynomial, npoly, isdf_in, wfn, cvpair, ncv, sqrtR, &
                                        MC_vec, CMC_vec, tmp_vec, tmp_Hvec, tmp_cvpair, blksz)

  use typedefs
  use mpi_module
  implicit none

  integer, intent(in) :: ncv_loc, npoly, blksz, ncv
  type(ISDF), intent(in) :: isdf_in
  real(dp), intent(in) :: polynomial(npoly + 1), &
                          sqrtR(ncv_loc)
  real(dp), intent(in) :: vec(ncv_loc, blksz)
  type(wavefunction), intent(in) :: wfn
  integer, intent(in) :: cvpair(4, ncv)
  ! polynomial is an array that stores the coefficients from higher order to lower order
  ! polynomial = [a_n, a_{n-1}, a_{n-2}, ..., a_0]
  ! polynomial(1) = a_n
  ! polynomial(2) = a_{n-1}
  ! ...
  ! polynomial(npoly+1) = a_0
  real(dp), intent(out) :: Hvec(ncv_loc, blksz)
  real(dp), intent(inout) :: MC_vec(w_grp%myn_intp_r, blksz), &
                             CMC_vec(w_grp%ldncv, blksz), &
                             tmp_vec(isdf_in%n_intp_r, blksz), tmp_Hvec(w_grp%ldncv, blksz)
  integer, intent(inout) :: tmp_cvpair(2, w_grp%ldncv)
  real(dp) :: Hvec_old(ncv_loc, blksz), tsec(2)
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
  call timacc(58, 1, tsec)
  do k = 1, npoly
    !if (peinf%master) print *, "k =", k, polynomial(k)
    if (k == 1) then
      call matvec_isdf_BLOCK(vec, Hvec, ncv_loc, isdf_in, wfn, cvpair, ncv, &
                             sqrtR, MC_vec, CMC_vec, tmp_vec, tmp_Hvec, tmp_cvpair, blksz)
      Hvec = polynomial(k)*Hvec + polynomial(k + 1)*vec
      !if (peinf%master) print *, "In polynomial_matvec_isdf: Hvec", Hvec(1:10)
    else
      call matvec_isdf_BLOCK(Hvec_old, Hvec, ncv_loc, isdf_in, wfn, cvpair, ncv, &
                             sqrtR, MC_vec, CMC_vec, tmp_vec, tmp_Hvec, tmp_cvpair, blksz)
      Hvec = Hvec + polynomial(k + 1)*vec
      !if (peinf%master) print *, "In polynomial_matvec_isdf: Hvec", Hvec(1:10)
      !stop
    end if
    Hvec_old = Hvec
  end do ! k
  call timacc(58, 2, tsec)

end subroutine polynomial_matvec_isdf_BLOCK

