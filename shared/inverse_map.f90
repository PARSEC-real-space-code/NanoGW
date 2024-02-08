!===================================================================
!
! Given integer sequences map(1:nn) and invm(1:mm), return a new
! sequence stored in map such that:
!
!    map_new(i) = invm(map_old(i))
!
! where invm(map_old(i)) is assumed non-zero for all i =1, ..., n
!
! Copyright (C) 2009 Murilo L. Tiago, http://users.ices.utexas.edu/~mtiago
! This file is part of RGWBS. It is distributed under the GPL v1.
!
!-------------------------------------------------------------------
subroutine inverse_map(mm,invm,nn,map)

  implicit none

  ! arguments
  integer, intent(in) :: mm, nn, invm(mm)
  integer, intent(inout) :: map(nn)

  ! local variables
  integer :: tmpmap(nn), ii

  if (nn == 0) return

  tmpmap = map
  do ii = 1, nn
     map(ii) = invm( tmpmap(ii) )
  enddo

end subroutine inverse_map
!===================================================================
