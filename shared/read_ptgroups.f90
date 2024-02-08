!===================================================================
!
! Read parameters for symmetry groups from rgwbs.in file.
!
! INPUT:
!    none
!
! OUTPUT:
!    syms%ngr : number of non-Abelian groups
!    syms%grfilename : names of input files where non-Abelian groups
!                      are specified.
!
! Copyright (C) 2009 Murilo L. Tiago, http://users.ices.utexas.edu/~mtiago
! This file is part of RGWBS. It is distributed under the GPL v1.
!
!-------------------------------------------------------------------
subroutine read_ptgroups(syms)
  !
  ! Read definition of point groups.
  !
  use typedefs
  use esdf
  implicit none

  ! arguments
  type (symmetries), intent(inout) :: syms

  ! local variables
  integer :: ii

  call esdf_init('rgwbs.in')

  if (esdf_block ('point_group_tables',syms%ngr )) then
     allocate(syms%grfilename(syms%ngr))
     do ii = 1, syms%ngr
        read(block_data(ii),*) syms%grfilename(ii)
     enddo
  else
     syms%ngr = 0
  endif

  call esdf_close

end subroutine read_ptgroups
!===================================================================
