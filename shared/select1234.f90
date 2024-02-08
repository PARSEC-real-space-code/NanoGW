!===================================================================
!
! For a quadruplet m1,m2,m3,m4 (absolute indices), look for
! equivalent quadruplets and figure out if they are also being
! computed. If other quadruplets are found, compute this one only
! if m3 is less than or equal to m4, and m1 is less or equal to m2.
! In this case, outflag returns .true.
!
! Copyright (C) 2009 Murilo L. Tiago, http://users.ices.utexas.edu/~mtiago
! This file is part of RGWBS. It is distributed under the GPL v1.
!
!-------------------------------------------------------------------
subroutine select1234(mat,syms,m1,m2,m3,m4,csp,rep_1,rep_2,rep_3,rep_4, &
     sympair,outflag)

  use typedefs
  implicit none

  ! arguments
  ! kernel structure
  type (kernelinfo), intent(in) :: mat
  ! Abelian group
  type (symmetries), intent(in) :: syms
  ! order of orbitals 1, 2, 3, 4, and spin channel
  integer, intent(in) :: m1, m2, m3, m4, csp
  ! representations of orbitals 1, 2, 3, 4
  integer, intent(in), dimension(0:syms%ngr) :: rep_1, rep_2, rep_3, rep_4
  ! true if symmetric pair (3, 4, 1, 2) will be calculated
  logical, intent(in) :: sympair
  ! true if pair (1, 2, 3, 4) needs to be calculated
  logical, intent(out) :: outflag

  ! local variables
  integer :: jj, jstart, jstop, igr, rep_prod

  outflag = .true.
  !-------------------------------------------------------------------
  ! Make sure that spin bounds are correct. If sympair is true, then
  ! the spin of (m1,m2) is equal to the spin of (m3,m4) by construction.
  ! Now, must make sure that the search of equivalent pairs have matching
  ! spins on the columns.
  if (csp == 1) then
     jstart = 1
     jstop = mat%ncol_up
  else
     jstart = mat%ncol_up + 1
     jstop = mat%ncol
  endif
  !-------------------------------------------------------------------
  ! Is the pair (m4,m3) going to be calculated as well?
  ! If so, do (m3,m4) only if m4>=m3. Otherwise, return with outflag=.false.
  !
  do jj = jstart, jstop
     if ((m3 == mat%col(2,jj)).and.(m4 == mat%col(1,jj)).and.(m4 < m3)) then
        outflag = .false.
        return
     endif
  enddo
  !-------------------------------------------------------------------
  ! Must update the spin bounds in order to search for equivalent pairs
  ! along rows.
  if (csp == 1) then
     jstart = 1
     jstop = mat%nrow_up
  else
     jstart = mat%nrow_up + 1
     jstop = mat%nrow
  endif
  !-------------------------------------------------------------------
  ! Is matrix element for the quadruplet (m1,m2,m3,m4) going to be
  ! calculated with (m3,m4,m1,m2)?
  ! if so, do (m3,m4,m1,m2) if m3>m1 or m3=m1, m4>=m2. Otherwise, return 
  ! with outflag=.false.
  ! 
  do jj = jstart, jstop
     if ((m4 == mat%row(2,jj)) .and. (m3 == mat%row(1,jj)) .and. sympair) then
        if (m1 > m3) then
           outflag = .false.
           return
        endif
        if (( m3 == m1) .and. (m4 < m2)) then
           outflag = .false.
           return
        endif
     endif
  enddo

  !-------------------------------------------------------------------
  ! If there are non-Abelian groups, calculate the products of representations.
  ! If there is no match, skip this matrix element.
  !
  do igr = 1, syms%ngr
     rep_prod = 0
     do jj = 1, syms%nrep(igr)
        rep_prod = rep_prod + syms%g_prod(rep_1(igr),rep_2(igr),jj,igr) * &
             syms%g_prod(rep_3(igr),rep_4(igr),jj,igr)
     enddo
     if (rep_prod == 0) outflag = .false.
  enddo

end subroutine select1234
!===================================================================

