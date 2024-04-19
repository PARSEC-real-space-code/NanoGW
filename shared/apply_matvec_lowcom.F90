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
subroutine matvec_isdf_lowcom( vec, Hvec, ncv_loc, isdf_in, ncv, sqrtR, &
#ifdef DCU
   !MC_vec, CMC_vec, tmp_vec, &
   blksz, d_PsiV, d_PsiC, d_Cmtrx, d_pvec, d_cvec, &
   hipblasHandle &
#else
   !MC_vec, CMC_vec, tmp_vec, &
   blksz &
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
#endif
  implicit none
#ifdef MPI
  include 'mpif.h'
#endif
  
  integer,  intent(in) :: ncv_loc, blksz, ncv
  real(dp), intent(in) :: vec(ncv_loc, blksz) 
  type(ISDF), intent(in) :: isdf_in
  real(dp), intent(out) :: Hvec(ncv_loc, blksz)
  real(dp), intent(in) :: sqrtR(ncv_loc)
  !real(dp), target, intent(inout) :: MC_vec(w_grp%myn_intp_r, blksz), &
  !  CMC_vec(ncv_loc, blksz), &
  !  tmp_vec(isdf_in%n_intp_r, blksz)
  real(dp), target :: MC_vec(w_grp%myn_intp_r, blksz), &
    CMC_vec(ncv_loc, blksz), &
    tmp_vec(isdf_in%n_intp_r, blksz)
#ifdef DCU
  type(c_ptr), intent(inout) :: d_PsiV, d_PsiC, d_Cmtrx, d_cvec, d_pvec
  type(c_ptr), intent(in) :: hipblasHandle
  real(dp), target :: tmparray(isdf_in%mync, blksz)
  real(dp) :: beta
  real(dp), parameter :: d_one = 1.d0, d_zero = 0.d0
  integer(kind(HIPBLAS_OP_N)), parameter :: &
    transt = HIPBLAS_OP_T, &
    transn = HIPBLAS_OP_N
  integer(c_size_t), parameter :: Bytes_double = 8
  integer(c_size_t) :: Nbytes
#else
  real(dp) :: tmp_Cmtrx(isdf_in%n_intp_r, isdf_in%mync)
#endif
  real(dp) :: tsec(2)
  integer :: iv, ic, icv, incr, intp_end, intp_start, &
    res, err, n_intp_r, ncvdim, ipe, ncvpair, ii, mync, mynv
  ! external functions
  real(dp), external :: ddot

! apply vec <-- R^1/2 vec
  !if (w_grp%master) print *, "1 vec", vec(1:10,1)
  call timacc(68,1,tsec)
  do ii = 1, blksz
  Hvec(1:ncv_loc, ii) = sqrtR(1:ncv_loc) * vec(1:ncv_loc, ii)
  enddo
  call timacc(68,2,tsec)
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
  n_intp_r   = isdf_in%n_intp_r
  intp_start = w_grp%n_intp_start(w_grp%inode)
  intp_end   = w_grp%n_intp_end  (w_grp%inode)
  mync       = isdf_in%mync
  mynv       = isdf_in%mynv
  tmp_vec = 0.d0
  MC_vec  = 0.d0
  CMC_vec = 0.d0
  !if (w_grp%master) print *, "breakpoint 1"
  call timacc(69,1,tsec)
  ! d_cvec ( n_intp_r, jnblock )
  ! d_pvec ( mync, jnblock )
  do iv = 1, mynv
#ifdef DCU
    Nbytes = mync*blksz*Bytes_double
    tmparray(1:mync, 1:blksz) = Hvec( ((iv-1)*mync+1): (iv*mync), 1:blksz )
    call hipCheck(hipMemcpy(d_pvec,  c_loc(tmparray), &
      Nbytes, hipMemcpyHostToDevice))
    ! d_Cmtrx ( n_intp_r, mync )
    call hadamard_prod_mat( d_PsiV, d_PsiC, d_Cmtrx, &
      isdf_in%n_intp_r, iv-1, isdf_in%mynv, isdf_in%mync )
    if (iv .eq. 1) then
      beta = d_zero
    else
      beta = d_one
    endif
    call hipblasCheck(hipblasDgemm(hipblasHandle, transn, transn, n_intp_r, &
      blksz, mync, d_one, d_Cmtrx, n_intp_r, &
      d_pvec, mync, beta, d_cvec, n_intp_r))
#else
    do ic = 1, mync
      tmp_Cmtrx(1:n_intp_r, ic) = isdf_in%PsiC_intp_bl(1:n_intp_r, ic, 1, 1) * &
               isdf_in%PsiV_intp_bl(1:n_intp_r, iv, 1, 1) 
    enddo
    call dgemm('N', 'N', n_intp_r, blksz, mync, 1.d0, tmp_Cmtrx, n_intp_r, &
      Hvec((iv-1)*mync+1, 1), ncv_loc, 1.d0, tmp_vec(1,1), n_intp_r)
#endif
  enddo
#ifdef DCU
  call hipCheck(hipDeviceSynchronize())
  Nbytes = n_intp_r*blksz*Bytes_double
  call hipCheck(hipMemcpy(c_loc(tmp_vec), d_cvec, &
    Nbytes, hipMemcpyDeviceToHost))
#endif
  call timacc(69,2,tsec)
  !if (w_grp%master) print *, "breakpoint 2"
  call timacc(70,1,tsec)
  if ( w_grp%npes .gt. 1) then
    call MPI_allreduce( MPI_IN_PLACE, tmp_vec, isdf_in%n_intp_r*blksz, &
      MPI_DOUBLE, MPI_SUM, w_grp%comm, err )
  endif
  call timacc(70,2,tsec)
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
  call timacc(71,1,tsec)
  !if (w_grp%master) print *, "breakpoint 3"
  call dgemm('N', 'N', w_grp%myn_intp_r, blksz, isdf_in%n_intp_r, 1.d0, &
    isdf_in%Mmtrx_loc(1,1,1,1,1,1,1), &
    w_grp%myn_intp_r, tmp_vec, isdf_in%n_intp_r, 0.d0, &
    MC_vec, w_grp%myn_intp_r)
  call timacc(71,2,tsec)
  !print *, "Mmtrx(1,1:n_intp_r)", isdf_in%Mmtrx(1,1:isdf_in%n_intp_r,1,1,1,1,1)
  !print *, "tmp_vec(1:n_intp_r)", tmp_vec(1:isdf_in%n_intp_r)
  !print *, "Mmtrx(1,1:n_intp_r)*tmp_vec(1:n_intp_r)", &
  ! sum(isdf_in%Mmtrx(1,1:isdf_in%n_intp_r,1,1,1,1,1)*tmp_vec(1:isdf_in%n_intp_r))
  call timacc(70,1,tsec)
  tmp_vec = 0.d0
  tmp_vec(intp_start:intp_end, 1:blksz) = MC_vec(1:w_grp%myn_intp_r, 1:blksz)
  call MPI_allreduce( MPI_IN_PLACE, tmp_vec, isdf_in%n_intp_r*blksz, MPI_DOUBLE, MPI_SUM, &
    w_grp%comm, err )
  call timacc(70,2,tsec)
  !if (w_grp%master) print *, "Mmtrx Cmtrx R^1/2 vec", MC_vec(1:10)
! apply CMC_vec = Cmtrx^T MC_vec
  call timacc(69,1,tsec)
  ! d_cvec ( n_intp_r, jnblock )
  ! d_pvec ( mync, jnblock )
#ifdef DCU
  Nbytes = n_intp_r*blksz*Bytes_double
  call hipCheck(hipMemcpy(d_cvec, c_loc(tmp_vec), &
         Nbytes, hipMemcpyHostToDevice))
#endif
  do iv = 1, mynv
#ifdef DCU
    ! d_Cmtrx ( n_intp_r, mync )
    call hadamard_prod_mat( d_PsiV, d_PsiC, d_Cmtrx, &
      n_intp_r, iv-1, mynv, mync )
    call hipblasCheck(hipblasDgemm(hipblasHandle, transt, transn, mync, &
      blksz, n_intp_r, d_one, d_Cmtrx, n_intp_r, &
      d_cvec, n_intp_r, d_zero, d_pvec, mync))
    call hipCheck(hipDeviceSynchronize())
    Nbytes = mync*blksz*Bytes_double
    call hipCheck(hipMemcpy(c_loc(tmparray), d_pvec, &
      Nbytes, hipMemcpyDeviceToHost))
    CMC_vec( ((iv-1)*mync+1):(iv*mync), 1:blksz ) = &
      tmparray( 1:mync, 1:blksz )
#else
    do ic = 1, mync
      tmp_Cmtrx(1:n_intp_r, ic) = isdf_in%PsiC_intp_bl(1:n_intp_r, ic, 1, 1) * &
         isdf_in%PsiV_intp_bl(1:n_intp_r, iv, 1, 1) 
    enddo
    call dgemm('T', 'N', mync, blksz, n_intp_r, 1.d0, tmp_Cmtrx, n_intp_r, &
       tmp_vec, n_intp_r, 0.d0, CMC_vec((iv-1)*mync+1,1), ncv_loc) 
#endif
  enddo
  call timacc(69,2,tsec)
  !if (w_grp%master) print *, "Inside apply_matvec_UltraLowMem, CMC_vec ", CMC_vec(1:5)
  !stop
  !if (w_grp%master) print *, "breakpoint 4"
  call timacc(72,1,tsec)
  do ii = 1, blksz
!   Hvec <- R^3/2 * Hvec 
    Hvec(1:ncv_loc, ii) = sqrtR(1:ncv_loc)**3 * Hvec(1:ncv_loc, ii)
    !print *, "4 Hvec", Hvec(1:5)
    
!   Hvec <- Hvec + R^1/2 * CMC_vec
    Hvec(1:ncv_loc, ii) = Hvec(1:ncv_loc, ii) + &
      4.d0 * sqrtR(1:ncv_loc) * CMC_vec(1:ncv_loc, ii)
    !print *, "5 Hvec", Hvec(1:5)
  enddo ! blksz
  call timacc(72,2,tsec)

end subroutine matvec_isdf_lowcom

subroutine polynomial_matvec_isdf_lowcom(vec, Hvec, ncv_loc, &
  polynomial, npoly, isdf_in, ncv, sqrtR, &
#ifdef DCU
  !MC_vec, CMC_vec, tmp_vec, &
  blksz, d_PsiV, d_PsiC, d_Cmtrx, d_pvec, d_cvec, &
  hipblasHandle &
#else
  !MC_vec, CMC_vec, tmp_vec, &
  blksz &
#endif
  )

  use typedefs
  use mpi_module
#ifdef DCU
  use hipfort_types
#endif
  implicit none
  
  integer, intent(in) :: ncv_loc, npoly, blksz, ncv
  type(ISDF), intent(in) :: isdf_in
  real(dp), intent(in) :: polynomial(npoly+1), &
   sqrtR(ncv_loc)
  real(dp), intent(in) :: vec(ncv_loc, blksz)
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
  type(c_ptr), intent(inout) :: d_PsiV, d_PsiC, d_Cmtrx, d_cvec, d_pvec
  type(c_ptr), intent(in) :: hipblasHandle
#endif
  real(dp) :: Hvec_old(ncv_loc, blksz), tsec(2)
  integer  :: k

! isdf_in%n_intp_r == w_grp%n_intp_r
! allocate(MC_vec(w_grp%myn_intp_r))
! allocate(CMC_vec(w_grp%myncv))
! allocate(tmp_vec(isdf_in%n_intp_r))
! allocate(tmp_Hvec(w_grp%ldncv))
  !   (a_n H^n + a_{n-1} H^{n-1} + ... + a_1 H + a_0 ) vec
  ! = (H*...(H*(H*(a_n H*vec + a_{n-1} vec) + a_{n-2} vec) + a_{n-3} vec)...a_1 vec) + a_0 vec
  ! Recursively:
  !   Hvec <- a_k H*vec + a_{k-1}
  call timacc(58, 1, tsec)
  do k = 1, npoly
    !if (peinf%master) print *, "k =", k, polynomial(k)
    if (k .eq. 1) then
      call matvec_isdf_lowcom(vec, Hvec, ncv_loc, isdf_in, ncv, sqrtR, &
#ifdef DCU
       !MC_vec, CMC_vec, tmp_vec, &
       blksz, d_PsiV, d_PsiC, d_Cmtrx, d_pvec, d_cvec, &
       hipblasHandle &
#else
       !MC_vec, CMC_vec, tmp_vec, &
       blksz &
#endif
      )
      Hvec = polynomial(k) * Hvec  + polynomial(k+1) * vec
      !if (peinf%master) print *, "In polynomial_matvec_isdf: Hvec", Hvec(1:10)
    else
      call matvec_isdf_lowcom(Hvec_old, Hvec, ncv_loc, isdf_in, ncv, sqrtR, &
#ifdef DCU
       !MC_vec, CMC_vec, tmp_vec, &
       blksz, d_PsiV, d_PsiC, d_Cmtrx, d_pvec, d_cvec, &
       hipblasHandle &
#else
       !MC_vec, CMC_vec, tmp_vec, &
       blksz &
#endif
      )
      Hvec = Hvec + polynomial(k+1) * vec
      !if (peinf%master) print *, "In polynomial_matvec_isdf: Hvec", Hvec(1:10)
      !stop
    endif
    Hvec_old = Hvec
  enddo ! k 
  call timacc(58, 2, tsec)

end subroutine polynomial_matvec_isdf_lowcom

