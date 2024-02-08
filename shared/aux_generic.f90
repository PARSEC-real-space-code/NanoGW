!===================================================================
!
! Wrapper to zdotu (BLAS function).
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
  complex(dpc), external :: zdotu

  dot_prod = zdotu(ndim,xx,incx,yy,incy)

end function zdot_u
!===================================================================
!
! Wrapper to zdotc (BLAS function).
!
!-------------------------------------------------------------------
function zdot_c(ndim,xx,incx,yy,incy) result(dot_prod)

  use myconstants
  implicit none

  integer, intent(in) :: ndim, incx, incy
  complex(dpc), intent(in) :: xx(ndim*incx-incx+1), yy(ndim*incy-incy+1)
  complex(dpc) :: dot_prod
  complex(dpc), external :: zdotc

  dot_prod = zdotc(ndim,xx,incx,yy,incy)

end function zdot_c
!===================================================================
