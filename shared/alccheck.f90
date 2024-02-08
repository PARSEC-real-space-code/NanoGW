!===============================================================
!
! Checks allocation of array with size size, name wst. If
! allocation fails (istat /= 0), an error message is writen in
! standard output and the calculation is aborted in an abrupt and
! possibly nasty way.
!
! All "big arrays" (i.e., with size proportional to grid size or
! square of number of eigenvalues) should have allocation checked.
!
! Copyright (C) 2009 Murilo L. Tiago, http://users.ices.utexas.edu/~mtiago
! This file is part of RGWBS. It is distributed under the GPL v1.
!
!---------------------------------------------------------------
subroutine alccheck(wst,wloc,isize,istat)
  implicit none

  ! name of array, name of parent subroutine
  character (len=*), intent(in) :: wst, wloc
  ! size of array, allocation status (see above)
  integer, intent(in) :: isize, istat

  character (len=600) :: lastwords

  if(istat /= 0) then
     write(lastwords,'(a,a,a,a,a,i12)') 'Cannot allocate ', wst, &
          ' in ', wloc, ' , ARRAY SIZE= ', isize
     call die(lastwords)
  endif

end subroutine alccheck
!===============================================================
