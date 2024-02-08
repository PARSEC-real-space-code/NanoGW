!===================================================================
!
! If verbose is true, print time elapsed since beginning of
! calculation, in seconds.
! Otherwise, do not print anything out.
! The time counter is initialized in the first call.
!
! Copyright (C) 2009 Murilo L. Tiago, http://users.ices.utexas.edu/~mtiago
! This file is part of RGWBS. It is distributed under the GPL v1.
!
!---------------------------------------------------------------
subroutine stopwatch(verbose,label)

  use myconstants
  implicit none

  ! arguments
  logical, intent(in) :: verbose
  character (len=*), intent(in) :: label

  ! local variables
  real(dp), save :: start_time = zero
  real(dp) :: this_time

  if (start_time == zero) call cpu_time(start_time)

  if (verbose) then
     call cpu_time(this_time)
     write(6,'(a,f11.2,a,a)') ' ELAPSED TIME = ', &
          this_time - start_time, ' s. ', label(1:len_trim(label))
     call flush(6)
  endif

end subroutine stopwatch
!===================================================================
