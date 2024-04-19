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
subroutine matvec_isdf( vec, Hvec, ncv_loc, isdf_in, wfn, tmpCmtrx, sqrtR)

  use typedefs
  use mpi_module
  implicit none
  include 'mpif.h'
  
  integer, intent(in) :: ncv_loc
  real(dp), intent(in) :: vec(ncv_loc) 
  type(ISDF), intent(in) :: isdf_in
  real(dp), intent(in) :: tmpCmtrx(isdf_in%n_intp_r,ncv_loc), &
    sqrtR(ncv_loc)
  real(dp), intent(out) :: Hvec(ncv_loc)
  real(dp), allocatable :: C_vec(:), MC_vec(:), &
    tmp_vec(:), CMC_vec(:)
  type(wavefunction), intent(in) :: wfn
  integer :: cvpair(2,ncv_loc)
  integer :: iv, ic, icv, incr, intp_end, intp_start, &
    res, err, n_intp_loc

! apply vec <-- R^1/2 vec
  !do icv = 1, ncv_loc 
  !  iv = cvpair(1,icv)
  !  ic = cvpair(2,icv)
  !  sqrtR(icv) = (wfn%e1(ic) - wfn%e1(iv))**0.5 
  !enddo
  if (peinf%master) print *, "1 vec", vec(1:10)
  Hvec = sqrtR * vec
  if (peinf%master) print *, "2 Hvec", Hvec(1:10)

! determine the local (n_intp_r) points
  incr = isdf_in%n_intp_r/peinf%npes
  res  = mod(isdf_in%n_intp_r, peinf%npes)
  if (peinf%inode .lt. res) then
      intp_start = peinf%inode*(incr+1)+1
      intp_end   = peinf%inode*(incr+1)+incr+1
  else
      intp_start = res*(incr+1)+(peinf%inode-res)*incr + 1
      intp_end   = res*(incr+1)+(peinf%inode-res)*incr + incr
  endif
  n_intp_loc = intp_end - intp_start + 1
  !print *, intp_start, intp_end
!  Cmtrx is distributed in rows (n_intp_r_loc, ncv_loc)
!  n_intp_r_loc = n_intp_r/nprocs
! Construct Cmtrx 
! Todo: implement dependence on spin and kpt index
  !allocate(tmpCmtrx(isdf_in%n_intp_r, ncv_loc))
  !do icv = 1, ncv_loc
  !  iv = cvpair(1,icv)
  !  ic = cvpair(2,icv)
  !  tmpCmtrx(1:isdf_in%n_intp_r, icv) = &
  !      isdf_in%Psi_intp(1:isdf_in%n_intp_r, iv, 1, 1) &
  !    * isdf_in%Psi_intp(1:isdf_in%n_intp_r, ic, 1, 1)
  !enddo ! icv loop
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
  allocate(tmp_vec(isdf_in%n_intp_r))
! compute C_vec = Cmtrx @ vec 
  call dgemv('N', isdf_in%n_intp_r, ncv_loc, 1.0d0, &
    tmpCmtrx, isdf_in%n_intp_r, Hvec, 1, 0.0d0, tmp_vec, 1)
  if (peinf%master) print *, "Inside matvec_isdf, tmp_vec", tmp_vec(1:5)
  !stop
! compute MC_vec = Mmtrx @ tmp_vec
!  Mtrx is distributed in columns (n_intp_r, n_intp_loc) 
!  Mtrx shape:
!   allocate(isdf_in%Mmtrx( n_intp_r, n_intp_r, nspin, nspin, kpt%nk, 2, gvec%syms%ntrans ))
!  Calculate Mmtrx(1) = <zeta(i)| V_coul(r,r') | zeta(j)> and 
!            Mmtrx(2) = <zeta(i)| f^LDA(r) | zeta(j)>                       
  call MPI_allreduce( MPI_IN_PLACE, tmp_vec, isdf_in%n_intp_r, MPI_DOUBLE, MPI_SUM, &
    peinf%comm, err )
  if (peinf%master) print *, "Cmtrx R^1/2 vec", tmp_vec(1:10)

  allocate(MC_vec(n_intp_loc))
  MC_vec  = 0.d0
  call dgemv('N', n_intp_loc, isdf_in%n_intp_r, 1.0d0, &
    isdf_in%Mmtrx(intp_start,1,1,1,1,1,1), &
    isdf_in%n_intp_r, tmp_vec, 1, 0.0d0, MC_vec, 1)
  !print *, "Mmtrx(1,1:n_intp_r)", isdf_in%Mmtrx(1,1:isdf_in%n_intp_r,1,1,1,1,1)
  !print *, "tmp_vec(1:n_intp_r)", tmp_vec(1:isdf_in%n_intp_r)
  !print *, "Mmtrx(1,1:n_intp_r)*tmp_vec(1:n_intp_r)", &
  ! sum(isdf_in%Mmtrx(1,1:isdf_in%n_intp_r,1,1,1,1,1)*tmp_vec(1:isdf_in%n_intp_r))
  tmp_vec = 0.d0
  tmp_vec(intp_start:intp_end) = MC_vec(1:n_intp_loc)
  call MPI_allreduce( MPI_IN_PLACE, tmp_vec, isdf_in%n_intp_r, MPI_DOUBLE, MPI_SUM, &
    peinf%comm, err )
  if (peinf%master) print *, "Mmtrx Cmtrx R^1/2 vec", tmp_vec(1:10)
  deallocate(MC_vec)
    
! apply CMC_vec = Cmtrx^T MC_vec
  allocate(CMC_vec(ncv_loc))
  call dgemv('T', isdf_in%n_intp_r, ncv_loc, 1.0d0, tmpCmtrx, isdf_in%n_intp_r, tmp_vec, 1, 0.0d0, CMC_vec, 1)
  if (peinf%master) print *, "Cmtrx^T Mmtrx Cmtrx R^1/2 vec", CMC_vec(1:5)
  !print *, "Inside apply_matvec, CMC_vec ", CMC_vec(1:5)

! Hvec <- R^3/2 * Hvec 
  Hvec = sqrtR**3 * Hvec
  !print *, "4 Hvec", Hvec(1:5)

! Hvec <- Hvec + R^1/2 * CMC_vec
  Hvec = Hvec + 4.d0 * sqrtR * CMC_vec
  !print *, "5 Hvec", Hvec(1:5)
  deallocate(CMC_vec, tmp_vec)

end subroutine matvec_isdf

subroutine polynomial_matvec_isdf(vec, Hvec, ncv_loc, polynomial, npoly, isdf_in, wfn, tmp_Cmtrx, sqrtR)

  use typedefs
  use mpi_module
  implicit none
  
  integer, intent(in) :: ncv_loc, npoly
  type(ISDF), intent(in) :: isdf_in
  real(dp), intent(in) :: polynomial(npoly+1), tmp_Cmtrx(isdf_in%n_intp_r), &
   sqrtR(ncv_loc)
  real(dp), intent(in) :: vec(ncv_loc)
  type(wavefunction), intent(in) :: wfn
  integer :: cvpair(2,ncv_loc)
  ! polynomial is an array that stores the coefficients from higher order to lower order
  ! polynomial = [a_n, a_{n-1}, a_{n-2}, ..., a_0]
  ! polynomial(1) = a_n
  ! polynomial(2) = a_{n-1}
  ! ...
  ! polynomial(npoly+1) = a_0
  real(dp), intent(out) :: Hvec(ncv_loc)
  real(dp) :: Hvec_old(ncv_loc)
  integer  :: k

  !   (a_n H^n + a_{n-1} H^{n-1} + ... + a_1 H + a_0 ) vec
  ! = (H*...(H*(H*(a_n H*vec + a_{n-1} vec) + a_{n-2} vec) + a_{n-3} vec)...a_1 vec) + a_0 vec
  ! Recursively:
  !   Hvec <- a_k H*vec + a_{k-1}
  do k = 1, npoly
    !if (peinf%master) print *, "k =", k, polynomial(k)
    if (k .eq. 1) then
        call matvec_isdf(vec, Hvec, ncv_loc, isdf_in, wfn, tmp_Cmtrx, &
         sqrtR)
        Hvec = polynomial(k) * Hvec  + polynomial(k+1) * vec
        !if (peinf%master) print *, "In polynomial_matvec_isdf: Hvec", Hvec(1:10)
    else
        call matvec_isdf(Hvec_old, Hvec, ncv_loc, isdf_in, wfn, tmp_Cmtrx, &
         sqrtR)
        Hvec = Hvec + polynomial(k+1) * vec
        !if (peinf%master) print *, "In polynomial_matvec_isdf: Hvec", Hvec(1:10)
        !stop
    endif
    Hvec_old = Hvec
  enddo ! k 

end subroutine polynomial_matvec_isdf

