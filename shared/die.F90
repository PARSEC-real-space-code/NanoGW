!===================================================================
!
! Print text in "lastwords" and abort the calculation. All
! crashes during execution should go through this subroutine.
!
! Copyright (C) 2009 Murilo L. Tiago, http://users.ices.utexas.edu/~mtiago
! This file is part of RGWBS. It is distributed under the GPL v1.
!
!-------------------------------------------------------------------
subroutine die(lastwords)

  use mpi_module
  implicit none
#ifdef MPI
  include 'mpif.h'
#endif

  ! arguments
  character (len=*), intent(in) :: lastwords

  ! local variables
  character (len=26) :: datelabel
#ifdef MPI
  integer :: info
#endif

  write(6,'(1x,a)') lastwords(1:len(trim(lastwords)))
  call get_date(datelabel)
  write(6,'(a,i5,3a)') ' PE # ', peinf%inode, ' stops. ', datelabel, ' UTC'
  call flush(6)
#ifdef MPI
  call MPI_ABORT(peinf%comm,info)
#endif
  stop

end subroutine die
!===================================================================
