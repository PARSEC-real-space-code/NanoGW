!===================================================================
!
! Read type of exchange-correlation functional from parsec.in.
!
! Copyright (C) 2009 Murilo L. Tiago, http://users.ices.utexas.edu/~mtiago
! This file is part of RGWBS. It is distributed under the GPL v1.
!
!-------------------------------------------------------------------
subroutine parsec_xc(nspin,x_id,c_id)

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

  call esdf_init('parsec.in')

  if (nspin == 1) then
     keyword = esdf_reduce(esdf_string ('correlation_type', 'ca'))
  else
     keyword = esdf_reduce(esdf_string ('correlation_type', 'pl'))
  endif

  select case (trim(keyword))
  case ('xa')
     x_id = XC_LDA_X
     c_id = XC_LDA_C_XALPHA
  case ('wi')
     x_id = XC_LDA_X
     c_id = XC_LDA_C_WIGNER
  case ('hl')
     x_id = XC_LDA_X
     c_id = XC_LDA_C_HL
  case ('ca','pz','lda')
     x_id = XC_LDA_X
     c_id = XC_LDA_C_PZ
  case ('pl','pw92','pwlda')
     x_id = XC_LDA_X
     c_id = XC_LDA_C_PW
  case ('pb','pbe')
     x_id = XC_GGA_X_PBE
     c_id = XC_GGA_C_PBE
  case ('blyp')
     x_id = XC_GGA_X_B88
     c_id = XC_GGA_C_LYP
  case default
     write(6,*) 'ERROR: unknown correlation type:', trim(keyword)
     call die('STOP.')
  end select

  call esdf_close

end subroutine parsec_xc
!===================================================================
