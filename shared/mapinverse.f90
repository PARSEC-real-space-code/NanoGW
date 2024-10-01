#include "../shared/mycomplex.h"
!===================================================================
!
! For input integer sequences map1(1:n1,m) and map2(1:n2,m),
! determines an associated sequence map1to2 so that
!  map1to2(i1) = i2 for all indices (i1,i2) such that map1(i1,:) = map2(i2,:)
!  map1to2(i1) = 0 for all values of i1 such that map1(i1,:) has no
!                  correspondent in map2
!
! Copyright (C) 2009 Murilo L. Tiago, http://users.ices.utexas.edu/~mtiago
! This file is part of RGWBS. It is distributed under the GPL v1.
!
!-------------------------------------------------------------------
subroutine mapinverse(m, n1, map1, n2, map2, map1to2)

  use mpi_module
  implicit none
  include 'mpif.h'

  ! arguments
  integer, intent(in) :: m, n1, n2
  integer, intent(in) :: map1(m, n1), map2(m, n2)
  integer, intent(inout) :: map1to2(n1)

  ! local variables
  integer :: i1, i2, mm, mapdiff, n_start, n_end, &
             err, incr, res
  logical :: found

  ! Weiwei Comment: the old code becomes extremely slow
  ! when n1 and n2 are very large
  !--old code start here--
!  map1to2 = 0
!  do i1 = 1,n1
!     do i2 = 1,n2
!        found = .true.
!        do mm = 1, m
!           if (map1(mm,i1) /= map2(mm,i2)) then
!              found = .false.
!           endif
!        enddo
!        if (found) then
!           map1to2(i1) = i2
!           exit
!        endif
!     enddo
!  enddo
  !--old code end here--

  !--new code start here--
  map1to2 = 0
  incr = n1/peinf%npes
  res = mod(n1, peinf%npes)
  if (peinf%inode < res) then
    n_start = peinf%inode*(incr + 1) + 1
    n_end = peinf%inode*(incr + 1) + incr + 1
  else
    n_start = res*(incr + 1) + (peinf%inode - res)*incr + 1
    n_end = res*(incr + 1) + (peinf%inode - res)*incr + incr
  end if
  ! print *, "inode ", peinf%inode, n_start, n_end
  do i1 = n_start, n_end
    do i2 = 1, n2
      mapdiff = sum(abs(map2(1:m, i2) - map1(1:m, i1)))
      if (mapdiff == 0) then
        map1to2(i1) = i2
        exit
      end if
    end do
  end do
  call MPI_BARRIER(peinf%comm, err)
  call mpi_allreduce(MPI_IN_PLACE, map1to2, n1, MPI_INTEGER, MPI_SUM, &
                     peinf%comm, err)
! call int_psum( n1, peinf%npes, peinf%comm, map1to2 )
  !--new code end here--

end subroutine mapinverse
!===================================================================
