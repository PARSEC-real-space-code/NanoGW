subroutine polyfitcheb_sqrt(xmin, xmax, npoly, nx, inode, beta, fit_error, verbosity)
  
  use typedefs
  implicit none
  !
  ! variable types
  !
  real(dp), intent(inout) :: xmin, xmax
  integer, intent(in) :: npoly, inode
  real(dp), intent(out) :: beta(npoly+1)
  integer, intent(in) :: nx
  real(dp), intent(out) :: fit_error
  integer, parameter :: mx = 2000
  real(dp) :: dx, tmpx, xx(nx), yy(nx), pxx(nx, npoly+1), &
    M(npoly+1, npoly+1), beta2(npoly+1,1), tmp_xx, tmp_pxx, &
    y(npoly+1), y2(npoly+1,1), fit_yy, emid, scal, xmin_new, xmax_new
  integer :: errinfo, ip, ix, ii
  logical :: verbosity

  ! print 
  if (xmin .ge. xmax) then
      ! swap xmin and xmax
      tmpx = xmin
      xmin = xmax
      xmax = tmpx
  endif
  if (xmin .le. 0.d0) then
     print *, " emin ", xmin, " eV. "
     call die("emin is smaller than zero!")
  endif
  xmin_new = xmin
  xmax_new = xmax
  if (xmin <= 1e-4) then
      write(6, *) "xmin = ", xmin
      call die("Polyfit_sqrt xmin is too small")
  endif
  emid = 0.5*(xmin_new + xmax_new)
  scal = 0.5*(xmax_new - xmin_new)
  !
  dx = pi/(nx-1)
  pxx(1:nx, 1) = 1.d0

  do ix = 1, nx
    xx(ix) = emid - scal * cos(dx*(ix-1)) ! Use chebyshev nodes for better fitting
    yy(ix) = sqrt(xx(ix))
    do ip = 2, npoly+1
      pxx(ix, ip) = xx(ix)**(ip-1)
    enddo
  enddo
  ! Do the following computations
  ! M = pxx^T pxx
  ! y = pxx^T yy
  ! solve M beta = y
  call dgemm('t', 'n', npoly+1, npoly+1, nx, 1.d0, pxx, nx, pxx, nx, &
  0.d0, M, npoly+1)
  !print *, "M = "
  !do ii = 1, npoly+1
  !  print *, M(ii, 1:npoly+1)
  !enddo
  call dgemv('t', nx, npoly+1, 1.d0, pxx, nx, yy, 1, 0.d0, y, 1 )
  !print *, "y = "
  !print *, y(1:npoly+1)
  y2(1:npoly+1, 1) = y(1:npoly+1)
  call dlinear_solver(npoly+1, 1, M, y, beta2, inode, .False., errinfo )
  
  beta(1:npoly+1) = beta2(1:npoly+1, 1)
  if (verbosity) then
    print *, "beta = "
    print *, beta2(1:npoly+1, 1)
  
    print *, "Compare"
    ! check the correctness
    print *, "x sqrt(x) fitted_sqrt rel. err."
    do ix = 1, nx
      fit_yy = 0.0
      do ip = 1, npoly+1
        fit_yy = fit_yy + pxx(ix,ip)*beta(ip)
      enddo
      write(6, '( f10.5  f10.5  f10.5   f6.1 a )') xx(ix), yy(ix), fit_yy, &
        (fit_yy-yy(ix))/yy(ix)*100, "%" 
    enddo
  endif

  fit_error = 0.d0
  do ix = 1, mx
    tmp_xx = xmin + (ix-1) * (xmax-xmin)/(mx-1)
    tmp_pxx = 1.d0
    fit_yy = 0.d0
    do ip = 1, npoly+1
      fit_yy = fit_yy + tmp_pxx * beta(ip)
      if (ip < npoly+1) tmp_pxx = tmp_pxx * tmp_xx
    enddo
    fit_error = fit_error + abs( (sqrt(tmp_xx) - fit_yy)/sqrt(tmp_xx) ) / dble(mx)
  enddo

end subroutine
