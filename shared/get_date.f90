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
  character(len=26), intent(out) :: datelabel

  integer :: tz_hours, tz_minutes
  character(len=1) :: tz_sign

  call date_and_time(VALUES=values)

  ! time zone information
  if (values(4) >= 0) then
    tz_sign = "+"
  else
    tz_sign = "-"
  end if
  tz_hours = abs(values(4)) / 60
  tz_minutes = mod(abs(values(4)), 60)

  write(datelabel, '(I4.4,"-",I2.2,"-",I2.2," ",I2.2,":",I2.2,":",I2.2,1X,A,I2.2,I2.2)') &
    values(1), values(2), values(3), values(5), values(6), values(7), tz_sign, tz_hours, tz_minutes

end subroutine get_date
!===================================================================
