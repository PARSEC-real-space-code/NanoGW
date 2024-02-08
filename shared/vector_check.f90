!===============================================================
!
! Compares two input real arrays of length nn. If they are
! different, return an error message and lcheck = .false.
! If they are identical, return lcheck = .true.
!
! Copyright (C) 2009 Murilo L. Tiago, http://users.ices.utexas.edu/~mtiago
! This file is part of RGWBS. It is distributed under the GPL v1.
!
!---------------------------------------------------------------
subroutine vector_check(filename,label,nn,vec_1,vec_2,lcheck)

  use myconstants
  implicit none

  ! name of checkpoint file, name of field in parent routine
  character (len=*), intent(in) :: filename, label
  ! length of arrays
  integer, intent(in) :: nn
  ! arrays to be compared
  real(dp), dimension(nn), intent(in) :: vec_1, vec_2
  ! true if arrays are identical, false otherwise
  logical, intent(out) :: lcheck

  integer :: ii
  real(dp), parameter :: TOL_D = 1.d-8

  lcheck = .true.
  do ii = 1, nn
     if (abs(vec_1(ii) - vec_2(ii)) > TOL_D) then
        lcheck = .false.
        exit
     endif
  enddo

  if (.not. lcheck) then
     write(6,*) ' WARNING: parameter mismatch in checkpoint file'
     write(6,*) ' Value of ', label, ' in memory is ', vec_1
     write(6,*) ' Value of ', label, ' from ', filename, ' is ', vec_2
     write(6,*) ' Ignoring checkpoint file.'
  endif

end subroutine vector_check
!===============================================================
