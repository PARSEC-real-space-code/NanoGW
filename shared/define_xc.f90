!===================================================================
!
! Define exchange-correlation functional in XC module.
!
! Copyright (C) 2009 Murilo L. Tiago, http://users.ices.utexas.edu/~mtiago
! This file is part of RGWBS. It is distributed under the GPL v1.
!
!-------------------------------------------------------------------
subroutine define_xc(nspin,x_id,c_id,init_gr)

  use myconstants
  use xc_functionals
  implicit none

  ! arguments
  ! number of spins, exchange flag, correlation flag
  integer, intent(in) :: nspin, x_id, c_id
  ! init_gr : true if functional includes spacial derivatives
  ! (e.g. GGA), left with input value otherwise
  logical, intent(inout) :: init_gr

  call xc_init(nspin,x_id,c_id,0,zero,one,.false.,xc_dft)

  if (.not. init_gr) init_gr = xc_dft%has_grad

end subroutine define_xc
!===================================================================
