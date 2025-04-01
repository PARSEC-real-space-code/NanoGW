program test_polyfit

  implicit none
  !
  ! variable types
  !
  integer, parameter :: dp = kind(1.0d0)
  integer, parameter :: dpc = kind((1.0d0, 1.0d0))

  real(8) :: xmin, xmax
  integer, parameter :: npoly = 8
  real(8) :: a(npoly + 1)

  print *, dp, dpc
  xmin = 0.05d0
  xmax = 8.d0
  ! print *, xmin, xmax

  call polyfit_sqrt(xmin, xmax, npoly, 0, a, .true.)

  print *, a(1:npoly + 1)

end program test_polyfit
