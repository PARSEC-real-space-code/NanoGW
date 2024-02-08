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
subroutine finalize(verbose,comm,count1,count2,routnam,timerlist)

  use myconstants
  implicit none

#ifdef MPI
  include 'mpif.h'
#endif
  logical, intent(in) :: verbose
  integer, intent(in) :: comm, count1, count2, timerlist(count1+count2)
  character (len=40), intent(in) :: routnam(count1+count2)

  character (len=26) :: datelabel
  integer :: ii, jj
  real(dp) :: tsec(2)
#ifdef MPI
  integer :: info
  real(dp) :: tmin(2), tmax(2)
#endif

!-------------------------------------------------------------------

  if (verbose) write(6,'(/,a,/,48x,a10,3x,a10,4x,a8,/)') &
       repeat('-',65), 'CPU [s]', 'WALL [s]', '#'
  do ii = 1, count1 + count2
     if (ii == count1 .and. verbose) write(6,*)
     jj = 3
     call timacc(timerlist(ii),jj,tsec)
#ifdef MPI
     call MPI_ALLREDUCE(tsec,tmin,2,MPI_DOUBLE_PRECISION,MPI_MIN,comm,info)
     call MPI_ALLREDUCE(tsec,tmax,2,MPI_DOUBLE_PRECISION,MPI_MAX,comm,info)
#endif
     if (verbose) then
#ifdef MPI
        write(6,'(1x,a40,a,f10.3,3x,f10.3)') &
             routnam(ii), '( min. )', tmin(1), tmin(2)
        write(6,'(41x,a,f10.3,3x,f10.3,3x,i8)') &
             '(master)', tsec(1), tsec(2), jj
        write(6,'(41x,a,f10.3,3x,f10.3)') '( max. )', tmax(1), tmax(2)
#else
        write(6,'(1x,a40,8x,f10.3,3x,f10.3,3x,i8)') &
             routnam(ii), tsec(1), tsec(2), jj
#endif
     endif
  enddo
  jj = 3
  call timacc(1,jj,tsec)
#ifdef MPI
  call MPI_ALLREDUCE(tsec,tmin,2,MPI_DOUBLE_PRECISION,MPI_MIN,comm,info)
  call MPI_ALLREDUCE(tsec,tmax,2,MPI_DOUBLE_PRECISION,MPI_MAX,comm,info)
#endif

  if (verbose) then
#ifdef MPI
     write(6,'(/,1x,a40,a,f10.3,3x,f10.3)') &
          'TOTAL', '( min. )', tmin(1), tmin(2)
     write(6,'(41x,a,f10.3,3x,f10.3)') '(master)', tsec(1), tsec(2)
     write(6,'(41x,a,f10.3,3x,f10.3,/)') '( max. )', tmax(1), tmax(2)
#else
     write(6,'(/,1x,a40,8x,f10.3,3x,f10.3,/)') &
          'TOTAL', tsec(1), tsec(2)
#endif
     call get_date(datelabel)
     write(6,'(3a,/,a)') ' Finished ', datelabel, ' UTC', repeat('-',65)
  endif

#ifdef MPI
  call MPI_FINALIZE(info)
#endif

  return
  end subroutine finalize
!===================================================================
