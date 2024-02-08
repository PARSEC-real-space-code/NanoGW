!===============================================================
!
! Compares two input integer arrays of length nn. If they are
! different, return an error message and lcheck = .false.
! If they are identical, return lcheck = .true.
!
! Copyright (C) 2009 Murilo L. Tiago, http://users.ices.utexas.edu/~mtiago
! This file is part of RGWBS. It is distributed under the GPL v1.
!
!---------------------------------------------------------------
subroutine map_check(filename,label,nn,map_1,map_2,lcheck)

  implicit none

  ! name of checkpoint file, name of field in parent routine
  character (len=*), intent(in) :: filename, label
  ! length of arrays
  integer, intent(in) :: nn
  ! arrays to be compared
  integer, intent(in), dimension(nn) :: map_1, map_2
  ! true if arrays are identical, false otherwise
  logical, intent(out) :: lcheck

  integer :: ii

  lcheck = .true.
  do ii = 1, nn
     if (map_1(ii) /= map_2(ii)) then
        lcheck = .false.
        exit
     endif
  enddo

  if (.not. lcheck) then
     write(6,*) ' WARNING: parameter mismatch in checkpoint file '
     write(6,*) ' Value of ', label, ' in memory is ', map_1
     write(6,*) ' Value of ', label, ' from ', filename, ' is ', map_2
     write(6,*) ' Ignoring checkpoint file.'
  endif

end subroutine map_check
!===============================================================
