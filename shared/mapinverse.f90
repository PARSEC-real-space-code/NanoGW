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
subroutine mapinverse(m,n1,map1,n2,map2,map1to2)

  implicit none

  ! arguments
  integer, intent(in) :: m, n1, n2
  integer, intent(in) :: map1(m,n1), map2(m,n2)
  integer, intent(inout) :: map1to2(n1)

  ! local variables
  integer :: i1, i2, mm
  logical :: found

  map1to2 = 0
  do i1 = 1,n1
     do i2 = 1,n2
        found = .true.
        do mm = 1, m
           if (map1(mm,i1) /= map2(mm,i2)) then
              found = .false.
           endif
        enddo
        if (found) then
           map1to2(i1) = i2
           exit
        endif
     enddo
  enddo

end subroutine mapinverse
!===================================================================
