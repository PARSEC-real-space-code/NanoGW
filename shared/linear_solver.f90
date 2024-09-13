subroutine dlinear_solver(n, m, A, B, X, inode, verbose, errinfo)

  use myconstants
  implicit none

  integer, intent (in) :: n, m, inode
  real (dp), intent (in) :: A(n,n)
  real (dp), intent (in) :: B(n,m)
  real (dp), intent (out) :: X(n,m)
  integer, intent (in) :: errinfo ! issue: why intent(in)?

  logical :: verbose
  ! Temporary local variables
  real (dp) :: Af(n,n), r(n), c(n), ferr(m), berr(m), work(4*n)
  integer :: ipiv(n), iwork(n)
  character :: equed, fact, trans
  integer :: lda, ldaf
  real (dp) :: rcond
  
  fact = 'E'
  trans = 'N'
  ! syntax: call dgesvx( fact, trans, n, nrhs, a, lda, 
  !           af, ldaf, ipiv, equed, r, c,
  !           b, ldb, x, ldx, rcond, ferr,
  !           berr, work, iwork, info )
  if (verbose) write (*, '(a,i6,a)') "Node", inode, &
      " call dgesvx to solve A*X = B"
  call dgesvx(fact, trans, n, m, A, n, Af, n, ipiv, equed, r, c, B, n, X, n, &
      rcond, ferr, berr, work, iwork, errinfo)
  if (errinfo == 0) then
    if (verbose) print *," dgesvx is successfully excecuted. "
  else if (errinfo < 0) then
    write (*, '(a,i6,a,i2,a)') " Node", inode, " The ", -errinfo, &
         "-th parameter of dgesvx has an illegal value."
  else if (errinfo > 0 .and. errinfo <= n) then
    write (*, '(a,i6,a)') " Node", inode, &
        " The U matrix is singular. LU factorization cannot be completed."

  else ! errinfo = n+1
    if (verbose) then 
      write (*, '(a,i6,a,a)') " Node", inode, &
          " The U matrix is nonsingular, but rcond is less than ", &
          "machine precision."
      write (*, '(a,a,f12.6)') " rcond is the reciprocal of condition ", &
          "number: = ", rcond
    end if
  end if

  return
end subroutine dlinear_solver
