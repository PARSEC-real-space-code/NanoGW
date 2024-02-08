!===================================================================
!
!  Non-recursive implementation of the quick sort method.
!
! INPUT:
!    nn : length of array to be sorted
!    arrin : (double precision) array to be sorted
!
! OUTPUT:
!    indx : sorting index. Lowest component in input array is arrin(indx(1)),
!           second lowest component is arrin(indx(2)) and so on.
!
!  Author: Yunkai Zhou, University of Minnesota, 2005
!
! Copyright (C) 2005 Finite Difference Research Group
! This file is part of parsec, http://www.ices.utexas.edu/parsec/
!
!-------------------------------------------------------------------
subroutine quicksort(nn, arrin, indx)

  use myconstants
  implicit none

  ! arguments
  integer, intent(in) :: nn                ! length of input array
  real(dp), intent(in) :: arrin(0:nn-1)    ! the array to be sorted
  integer, intent(out) :: indx(0:nn-1)     ! return the index of sorted arrin

  ! local variables
  integer :: ii, igap, jj, itemp
  real(dp) :: temp, xtmp(0:nn-1)

  !-------------------------------------------------------------------

  do ii = 0, nn - 1
     indx(ii) = ii + 1
  end do

  if (nn <= 1) return

  xtmp = arrin
  igap = nn / 2

  do while (igap > 0) 
              
     do ii = igap, nn - 1
        jj = ii - igap

        do while (jj >= 0) 
           ! Always sort the array in non-decreasing order 
           ! according to xtmp. Need to reorder xtmp since the
           ! next comparisons are based on the changing xtmp.
           if (xtmp(jj) > xtmp(jj+igap)) then
              temp = xtmp(jj)
              xtmp(jj) = xtmp(jj+igap)
              xtmp(jj+igap) = temp

              itemp = indx(jj)
              indx(jj) = indx(jj+igap)
              indx(jj+igap) = itemp
           else
              exit
           end if
           jj = jj - igap
        end do
     end do
     igap = igap / 2
  end do

end subroutine quicksort
!===================================================================
