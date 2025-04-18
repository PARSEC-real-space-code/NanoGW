!===================================================================
!
! Prints out timing report.
! If verbose = .true., prints out to unit 6. Otherwise, does some
! time accounting but produces no output.
! Labels for counters are stored in routnam.
! Counter 1 is for total time.
! Counters 2 to count1+1 are for subroutines called by main program.
! Counters count1+2 to count1+count2+1 are for subroutines called by
! other subroutines.
! INPUT:
!    verbose : output flag, see above
!    comm : MPI communicator
!    count1 : number of outer subroutines, see above
!    count2 : number of inner subroutines, see above
!    routnam : name of all subroutines.
!
! Copyright (C) 2009 Murilo L. Tiago, http://users.ices.utexas.edu/~mtiago
! This file is part of RGWBS. It is distributed under the GPL v1.
!
!-------------------------------------------------------------------
subroutine finalize(verbose, comm, count1, count2, routnam, timerlist)

  use myconstants
  implicit none

#ifdef MPI
  include 'mpif.h'
#endif
  logical, intent(in) :: verbose
  integer, intent(in) :: comm, count1, count2, timerlist(count1 + count2)
  character(len=40), intent(in) :: routnam(count1 + count2)

  character(len=26) :: datelabel
  integer :: ii, jj
  real(dp) :: tsec(2)
#ifdef MPI
  integer :: info
  real(dp) :: tmin(2), tmax(2)
#endif

!-------------------------------------------------------------------

  if (verbose) then
    write (6, '(/,a,/)') repeat('-', 96)
#ifdef MPI
    write (6, '(33x,a10,22x,a10,19x,a1)') 'CPU [s]', 'WALL [s]', '#'
    write (6, '(23x,3a10,2x,3a10)') 'min.', 'master', 'max.', 'min.', 'master', 'max.'
#else
    write (6, '(22x,a10,2x,a10,7x,a1)') 'CPU [s]', 'WALL [s]', '#'
#endif
  end if
  do ii = 1, count1 + count2
    if (ii == count1 + 1 .and. verbose) write (6, *)
    jj = 3
    call timacc(timerlist(ii), jj, tsec)
#ifdef MPI
    call MPI_ALLREDUCE(tsec, tmin, 2, MPI_DOUBLE_PRECISION, MPI_MIN, comm, info)
    call MPI_ALLREDUCE(tsec, tmax, 2, MPI_DOUBLE_PRECISION, MPI_MAX, comm, info)
#endif
    if (verbose .and. jj > 0) then
#ifdef MPI
      write (6, '(1x,a22,3f10.2,2x,3f10.2,i10)') &
        routnam(ii), tmin(1), tsec(1), tmax(1), tmin(2), tsec(2), tmax(2), jj
#else
      write (6, '(1x,a22,f10.2,2x,f10.2,i10)') routnam(ii), tsec(1), tsec(2), jj
#endif
    end if
  end do ! ii
  jj = 3
  call timacc(1, jj, tsec)
#ifdef MPI
  call MPI_ALLREDUCE(tsec, tmin, 2, MPI_DOUBLE_PRECISION, MPI_MIN, comm, info)
  call MPI_ALLREDUCE(tsec, tmax, 2, MPI_DOUBLE_PRECISION, MPI_MAX, comm, info)
#endif

  if (verbose) then
#ifdef MPI
    write (6, '(/,1x,a6,16x,3f10.2,2x,3f10.2,/)') &
      'TOTAL:', tmin(1), tsec(1), tmax(1), tmin(2), tsec(2), tmax(2)
#else
    write (6, '(/,1x,a6,16x,f10.2,2x,f10.2,/)') 'TOTAL:', tsec(1), tsec(2)
#endif
    call get_date(datelabel)
    write (6, '(3a,/,a)') ' Finished ', datelabel, ' UTC', repeat('-', 96)
  end if

#ifdef MPI
  call MPI_FINALIZE(info)
#endif

  return
end subroutine finalize
!===================================================================
