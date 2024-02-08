!===================================================================
!
! Ad hoc implementation of zdotu (BLAS function). Gfortran compiler
! seems to crash with the original BLAS function.
!
! Copyright (C) 2009 Murilo L. Tiago, http://users.ices.utexas.edu/~mtiago
! This file is part of RGWBS. It is distributed under the GPL v1.
!
!-------------------------------------------------------------------
function zdot_u(ndim,xx,incx,yy,incy) result(dot_prod)

  use myconstants
  implicit none

  integer, intent(in) :: ndim, incx, incy
  complex(dpc), intent(in) :: xx(ndim*incx-incx+1), yy(ndim*incy-incy+1)
  complex(dpc) :: dot_prod
  integer :: ii

  dot_prod = zzero
  do ii = 1, ndim
     dot_prod = dot_prod + xx((ii-1)*incx + 1) * yy((ii-1)*incy + 1)
  enddo

end function zdot_u
!===================================================================
!
! Ad hoc implementation of zdotc (BLAS function). Gfortran compiler
! seems to crash with the original BLAS function.
!
!-------------------------------------------------------------------
function zdot_c(ndim,xx,incx,yy,incy) result(dot_prod)

  use myconstants
  implicit none

  integer, intent(in) :: ndim, incx, incy
  complex(dpc), intent(in) :: xx(ndim*incx-incx+1), yy(ndim*incy-incy+1)
  complex(dpc) :: dot_prod
  integer :: ii

  dot_prod = zzero
  do ii = 1, ndim
     dot_prod = dot_prod + conjg( xx((ii-1)*incx + 1) ) * yy((ii-1)*incy + 1)
  enddo

end function zdot_c
!===================================================================
