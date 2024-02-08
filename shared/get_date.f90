!===================================================================
!
! Returns the date (YYYY MMM DD) and time (hh:mm:ss tzone) in z 26
! character string.
!
! Copyright (C) 2009 Murilo L. Tiago, http://users.ices.utexas.edu/~mtiago
! This file is part of RGWBS. It is distributed under the GPL v1.
!
!-------------------------------------------------------------------
subroutine get_date(datelabel)

  implicit none

  ! arguments
  integer :: values(8)
  character (len=26), intent(out) :: datelabel

  ! local variables
  character (len=3) :: month(12) = (/'Jan','Feb','Mar','Apr','May', &
       'Jun','Jul','Aug','Sep','Oct','Nov','Dec'/)

  call date_and_time(VALUES=values)

  write(datelabel,'(i4,1x,a,1x,i2,1x,i2,a,i2,a,i2,1x,i5)') &
       values(1), month(values(2)), values(3), &
       values(5), ':', values(6), ':', values(7), values(4)

  return
end subroutine get_date
!===================================================================
