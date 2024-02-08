!===================================================================
!
! Timing subroutine. Return elapsed cpu and wall clock times in sec.
! Also return the number of times the counter has been called.
!
! Depending on value of "option" routine will:
! (1) Start with new incremental time slice for accumulator n
!   also increase by one the counter for this accumulator.
! (2) Stop time slice; add time to accumlator n.
! (3) Report accumulated time for accumulator n
!   and number of time that the routine has been called.
!
! INPUT:
!  n = index of accumulator (distinguish what is being timed);
!  option = see comment above
!
! OUTPUT:
!  on option = 3 : 
!    tottim(1) = accumulated CPU time for accumulator n.
!    tottim(2) = accumulated wallclock time for accumulator n.
!    option = number of times that accumulator n has been incremented.
!  There is no output or any other input value of option 3.
!
! Initial version from the Berkeley GW package (G-M Rignanese, X. Blase, 1998 ?).
! Revised by M. Tiago (2004).
!
! Copyright (C) 2009 Murilo L. Tiago, http://users.ices.utexas.edu/~mtiago
! This file is part of RGWBS. It is distributed under the GPL v1.
!
!-------------------------------------------------------------------
subroutine timacc(n,option,tottim)    

  use myconstants
  implicit none

  ! arguments
  integer, intent(in) :: n
  integer, intent(inout) :: option
  real(dp),intent(out) :: tottim(2)

  ! local variables
  character (len=600) :: lastwords
  integer, parameter :: MTIM = 100 ! Maximum number of "timing slots" available
  integer, save :: ncount(MTIM) = 0
  real(dp), dimension(2,MTIM), save :: acctim = zero, tzero = zero
  real(dp) :: cpu, wall
  integer :: values(8)

  !-------------------------------------------------------------------
  ! Does n have valid value ?
  !
  if (n < 1 .or. n > MTIM) then
     write(lastwords,'(a,i6,a,i8)') ' timacc: dim mtim=', MTIM, &
          ' but input n= ', n
     call die(lastwords)
  end if

  !-------------------------------------------------------------------
  ! What time is it?
  !
  call cpu_time(cpu)
  call date_and_time(VALUES=values)
  wall = ((values(3)*24.0d0+values(5))*60.0d0 + values(6))*60.0d0 &  
       + values(7) + values(8)*0.001d0

  select case (option)
  case(1)
     ! Start stopwatch for accumulator n.
     tzero(1,n) = cpu
     tzero(2,n) = wall
  case(2)
     ! Stop stopwatch and accumulate elapsed time for n.
     acctim(1,n) = acctim(1,n) + cpu - tzero(1,n)
     acctim(2,n) = acctim(2,n) + wall - tzero(2,n)
     ncount(n) = ncount(n) + 1
  case(3)
     ! Return accumulated time for n.
     tottim(1) = acctim(1,n)
     tottim(2) = acctim(2,n)
     option = ncount(n)
  case default
     write(lastwords,'(a,i10)' ) ' timacc: input option not valid, = ' , option
     call die(lastwords)
  end select

end subroutine timacc
!===================================================================
