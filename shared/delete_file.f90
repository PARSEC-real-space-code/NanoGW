!===============================================================
!
! Delete specified file.
!
! Copyright (C) 2009 Murilo L. Tiago, http://users.ices.utexas.edu/~mtiago
! This file is part of RGWBS. It is distributed under the GPL v1.
!
!-------------------------------------------------------------------
subroutine delete_file(itape,filnam)

  implicit none
  ! arguments
  ! number of file (must be an unused valid unit number)
  integer, intent(in) :: itape
  ! name of file to be deleted
  character(len=*), intent(in) :: filnam

  open(itape,file=filnam(1:len_trim(filnam)))
  close(itape,status='DELETE')

end subroutine delete_file
!===============================================================
