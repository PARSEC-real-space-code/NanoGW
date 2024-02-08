!===================================================================
!
! Read type of exchange-correlation functional from input.
!
! Copyright (C) 2009 Murilo L. Tiago, http://users.ices.utexas.edu/~mtiago
! This file is part of RGWBS. It is distributed under the GPL v1.
!
!-------------------------------------------------------------------
subroutine paratec_xc(nspin,x_id,c_id)

  use myconstants
  use xc_functionals
  use esdf
  implicit none

  ! arguments
  ! number of spin channels
  integer, intent(in) :: nspin
  ! exchange flag, correlation flag
  integer, intent(out) :: x_id, c_id

  ! local variables
  character (len=80) :: keyword

  call esdf_init('input')

  keyword = esdf_reduce(esdf_string ('exchange_correlation', 'ceperley_alder'))

  select case (trim(keyword))
  case ('ceperley_alder')
     x_id = XC_LDA_X
     c_id = XC_LDA_C_PZ
  case ('perdew_burke_ernzerhof')
     x_id = XC_GGA_X_PBE
     c_id = XC_GGA_C_PBE
  case default
     write(6,*) 'ERROR: unknown correlation type:', trim(keyword)
     call die('STOP.')
  end select

  call esdf_close

end subroutine paratec_xc
!===================================================================
