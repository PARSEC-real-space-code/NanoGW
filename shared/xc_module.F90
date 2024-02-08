!===================================================================
!
! Interface module for libraries libxc and libstrig_f. This module
! contains functions and subroutines that define and compute
! exchange-correlation functionals of several types (LDA, GGA, meta-GGA,
! hybrids). Currently, the module recognizes only some of the LDA and
! GGA functionals. Better support should be added if necessary.
! Libraries libxc and libstring_f are developed by the Octopus development
! group (see http://www.tddft.org/programs/octopus for documentation,
! mailing list, download etc.).
!
! This interface is based on the XC interface distributed in the APE code
! (see http://www.tddft.org/programs/APE, by Micael Oliveira and Fernando
! Nogueira).
!
! Warnings: 
!  1. The input density is assumed in a.u.^-3 (Bohr^-3) and the output
! energies are in Ryd.
!
!  2. There should be a call to xc_end at the end of calculation.
!
! Copyright (C) 2009 Murilo L. Tiago, http://users.ices.utexas.edu/~mtiago
! This file is part of RGWBS. It is distributed under the GPL v1.
!
!-------------------------------------------------------------------
module lib_xc_m
  implicit none

  private
  public ::                         &
    xc_f90_info_number,             &
    xc_f90_info_kind,               &
    xc_f90_info_name,               &
    xc_f90_info_family,             &
    xc_f90_info_refs,               &
    xc_f90_family_from_id,          &
    xc_f90_lda_init,                &
    xc_f90_lda_vxc,                 &
    xc_f90_lda_fxc,                 &
    xc_f90_lda_kxc,                 &
    xc_f90_lda_end,                 &
    xc_f90_lca_init,                &
    xc_f90_lca,                     &
    xc_f90_gga_init,                &
    xc_f90_gga,                     &
    xc_f90_gga_end,                 &
    xc_f90_gga_lb,                  &
    xc_f90_mgga_init,               &
    xc_f90_mgga,                    &
    xc_f90_mgga_end

  ! Families of xc functionals
  integer, public, parameter ::     &
    XC_FAMILY_UNKNOWN       =  -1,  &
    XC_FAMILY_LDA           =   1,  &
    XC_FAMILY_GGA           =   2,  &
    XC_FAMILY_MGGA          =   4,  &
    XC_FAMILY_LCA           =   8,  &
    XC_FAMILY_OEP           =  16

  integer, public, parameter ::     &
    XC_UNPOLARIZED          =   1,  &  ! Spin unpolarized
    XC_POLARIZED            =   2      ! Spin polarized

  integer, public, parameter ::     &
    XC_NON_RELATIVISTIC     =   0,  &  ! Functional includes or not realtivistic
    XC_RELATIVISTIC         =   1      ! corrections. Only available in some functionals.

  ! Kinds
  integer, public, parameter ::     &
    XC_EXCHANGE             =   0,  &
    XC_CORRELATION          =   1,  &
    XC_EXCHANGE_CORRELATION =   2

  ! the LDAs
  integer, public, parameter ::     &
    XC_LDA_X                =   1,  &  ! Exchange
    XC_LDA_C_WIGNER         =   2,  &  ! Wigner parametrization
    XC_LDA_C_RPA            =   3,  &  ! Random Phase Approximation
    XC_LDA_C_HL             =   4,  &  ! Hedin & Lundqvist
    XC_LDA_C_GL             =   5,  &  ! Gunnarson & Lundqvist
    XC_LDA_C_XALPHA         =   6,  &  ! Slaters Xalpha
    XC_LDA_C_VWN            =   7,  &  ! Vosko, Wilk, & Nussair
    XC_LDA_C_VWN_RPA        =   8,  &  ! Vosko, Wilk, & Nussair (RPA)
    XC_LDA_C_PZ             =   9,  &  ! Perdew & Zunger
    XC_LDA_C_PZ_MOD         =  10,  &  ! Perdew & Zunger (Modified)
    XC_LDA_C_OB_PZ          =  11,  &  ! Ortiz & Ballone (PZ)
    XC_LDA_C_PW             =  12,  &  ! Perdew & Wang
    XC_LDA_C_PW_MOD         =  13,  &  ! Perdew & Wang (Modified)
    XC_LDA_C_OB_PW          =  14,  &  ! Ortiz & Ballone (PW)
    XC_LDA_C_AMGB           =  15      ! Attacalite et al

  ! the GGAs
  integer, public, parameter ::     &
    XC_GGA_X_PBE            = 101,  &  ! Perdew, Burke & Ernzerhof exchange
    XC_GGA_X_PBE_R          = 102,  &  ! Perdew, Burke & Ernzerhof exchange (revised)
    XC_GGA_X_B86            = 103,  &  ! Becke 86 Xalpha,beta,gamma
    XC_GGA_X_B86_R          = 104,  &  ! Becke 86 Xalpha,beta,gamma reoptimized
    XC_GGA_X_B86_MGC        = 105,  &  ! Becke 88 Xalfa,beta,gamma (with mod. grad. correction)
    XC_GGA_X_B88            = 106,  &  ! Becke 88
    XC_GGA_X_G96            = 107,  &  ! Gill 96
    XC_GGA_X_PW86           = 108,  &  ! Perdew & Wang 86
    XC_GGA_X_PW91           = 109,  &  ! Perdew & Wang 91
    XC_GGA_X_OPTX           = 110,  &  ! Handy & Cohen OPTX 01
    XC_GGA_X_DK87_R1        = 111,  &  ! dePristo & Kress 87 version R1
    XC_GGA_X_DK87_R2        = 112,  &  ! dePristo & Kress 87 version R1
    XC_GGA_X_LG93           = 113,  &  ! Lacks & Gordon 93
    XC_GGA_X_FT97_A         = 114,  &  ! Filatov & Thiel 97 (version A)
    XC_GGA_X_FT97_B         = 115,  &  ! Filatov & Thiel 97 (version B)
    XC_GGA_C_PBE            = 130,  &  ! Perdew, Burke & Ernzerhof correlation
    XC_GGA_C_LYP            = 131,  &  ! Lee, Yang & Parr
    XC_GGA_C_P86            = 132,  &  ! Perdew 86
    XC_GGA_XC_LB            = 160      ! van Leeuwen & Baerends


  ! the meta-GGAs
  integer, public, parameter ::     &
    XC_MGGA_X_TPSS          = 201,  &  ! Perdew, Tao, Staroverov & Scuseria exchange
    XC_MGGA_C_TPSS          = 202      ! Perdew, Tao, Staroverov & Scuseria correlation

  ! the LCAs
  integer, public, parameter ::     &
    XC_LCA_OMC              = 301,  &  ! Orestes, Marcasso & Capelle
    XC_LCA_LCH              = 302      ! Lee, Colwell & Handy

  ! info
  interface
    integer function xc_f90_info_number(info)
      integer (kind=8), intent(in) :: info
    end function xc_f90_info_number

    integer function xc_f90_info_kind(info)
      integer (kind=8), intent(in) :: info
    end function xc_f90_info_kind

    subroutine xc_f90_info_name(info, s)
      integer (kind=8), intent(in)  :: info
      character(len=*), intent(out) :: s
    end subroutine xc_f90_info_name

    integer function xc_f90_info_family(info)
      integer (kind=8), intent(in)  :: info
    end function xc_f90_info_family

    subroutine xc_f90_info_refs(info, n, s)
      integer (kind=8), intent(in)    :: info
      integer,          intent(inout) :: n
      character(len=*), intent(out)   :: s
    end subroutine xc_f90_info_refs
  end interface


  ! functionals
  interface
    integer function xc_f90_family_from_id(id)
      integer, intent(in) :: id
    end function xc_f90_family_from_id
  end interface


  ! We will use the same public interface (xc_lda_init) for the three C procedures
  interface xc_f90_lda_init
    subroutine xc_f90_lda_init_(p, info, functional, nspin)
      integer (kind=8), intent(out) :: p
      integer (kind=8), intent(out) :: info
      integer,          intent(in)  :: functional
      integer,          intent(in)  :: nspin
    end subroutine xc_f90_lda_init_

    subroutine xc_f90_lda_x_init(p, info, functional, nspin, dim, irel)
      integer (kind=8), intent(out) :: p
      integer (kind=8), intent(out) :: info
      integer,   intent(in)  :: functional
      integer,   intent(in)  :: nspin  ! XC_UNPOLARIZED or XC_POLARIZED
      integer,   intent(in)  :: dim    ! 2 or 3 dimensions
      integer,   intent(in)  :: irel   ! XC_NON_RELATIVISTIC or XC_RELATIVISTIC
    end subroutine xc_f90_lda_x_init

    subroutine xc_f90_lda_c_xalpha_init(p, info, functional, nspin, dim, alpha)
      use myconstants
      integer (kind=8), intent(out) :: p
      integer (kind=8), intent(out) :: info
      integer,   intent(in)  :: functional
      integer,   intent(in)  :: nspin  ! XC_UNPOLARIZED or XC_POLARIZED
      integer,   intent(in)  :: dim    ! 2 or 3 dimensions
      real(dp),   intent(in)  :: alpha  ! Ec = alpha Ex
    end subroutine xc_f90_lda_c_xalpha_init
  end interface

  interface
    subroutine xc_f90_lda_end(p)
      integer (kind=8), intent(inout) :: p
    end subroutine xc_f90_lda_end

    subroutine xc_f90_lda_vxc(p, rho, e, v)
      use myconstants
      integer (kind=8), intent(in)  :: p
      real(dp),     intent(in)  :: rho   ! rho(nspin) the density
      real(dp),     intent(out) :: e     ! the energy per unit particle
      real(dp),     intent(out) :: v     ! v(nspin) the potential
    end subroutine xc_f90_lda_vxc

    subroutine xc_f90_lda_fxc(p, rho, fxc)
      use myconstants
      integer (kind=8), intent(in)  :: p
      real(dp),     intent(in)  :: rho   ! rho(nspin) the density
      real(dp),     intent(out) :: fxc   ! v(nspin,nspin) the xc kernel
    end subroutine xc_f90_lda_fxc

    subroutine xc_f90_lda_kxc(p, rho, kxc)
      use myconstants
      integer (kind=8), intent(in)  :: p
      real(dp),     intent(in)  :: rho   ! rho(nspin) the density
      real(dp),     intent(out) :: kxc
    end subroutine xc_f90_lda_kxc

  end interface


  ! We will use the same public procedure for the two C procedures.
  interface xc_f90_gga_init
    subroutine xc_f90_gga_init_(p, info, functional, nspin)
      integer (kind=8), intent(out) :: p
      integer (kind=8), intent(out) :: info
      integer,   intent(in)  :: functional
      integer,   intent(in)  :: nspin
    end subroutine xc_f90_gga_init_

    subroutine xc_f90_gga_lb_init(p, info, functional, nspin, modified, threshold)
      use myconstants
      integer (kind=8), intent(out) :: p
      integer (kind=8), intent(out) :: info
      integer,   intent(in)  :: functional
      integer,   intent(in)  :: nspin
      integer,   intent(in)  :: modified
      real(dp),  intent(in)  :: threshold
    end subroutine xc_f90_gga_lb_init
 end interface

  interface
    subroutine xc_f90_gga_end(p)
      integer (kind=8), intent(inout) :: p
    end subroutine xc_f90_gga_end

    subroutine xc_f90_gga(p, rho, grho, e, dedd, dedgd)
      use myconstants
      integer (kind=8), intent(in)  :: p
      real(dp),     intent(in)  :: rho   ! rho(nspin) the density
      real(dp),     intent(in)  :: grho  ! grho(3,nspin) the gradient of the density
      real(dp),     intent(out) :: e     ! the energy per unit particle
      real(dp),     intent(out) :: dedd  ! dedd(nspin) the derivative of the energy
      ! in terms of the density
      real(dp),     intent(out) :: dedgd ! and in terms of the gradient of the density
    end subroutine xc_f90_gga

    subroutine xc_f90_gga_lb(p, rho, grho, r, ip, qtot, dedd)
      use myconstants
      integer (kind=8), intent(in)  :: p
      real(dp),     intent(in)  :: rho   ! rho(nspin) the density
      real(dp),     intent(in)  :: grho  ! grho(3,nspin) the gradient of the density
      real(dp),     intent(in)  :: r     ! distance from center of finite system
      real(dp),     intent(in)  :: ip    ! ionization potential
      real(dp),     intent(in)  :: qtot  ! total charge
      real(dp),     intent(out) :: dedd
    end subroutine xc_f90_gga_lb
  end interface

  ! the meta-GGAs
  interface
    subroutine xc_f90_mgga_init(p, info, functional, nspin)
      integer (kind=8), intent(out) :: p
      integer (kind=8), intent(out) :: info
      integer,   intent(in)  :: functional
      integer,   intent(in)  :: nspin
    end subroutine xc_f90_mgga_init

    subroutine xc_f90_mgga_end(p)
      integer (kind=8), intent(inout) :: p
    end subroutine xc_f90_mgga_end

    subroutine xc_f90_mgga(p, rho, grho, tau, e, dedd, dedgd, dedtau)
      use myconstants
      integer (kind=8), intent(in)  :: p
      real(dp),     intent(in)  :: rho   ! rho(nspin) the density
      real(dp),     intent(in)  :: grho  ! grho(3,nspin) the gradient of the density
      real(dp),     intent(in)  :: tau   ! tau(nspin) the kinetic energy density
      real(dp),     intent(out) :: e     ! the energy per unit particle
      real(dp),     intent(out) :: dedd  ! dedd(nspin) the derivative of the energy
      ! in terms of the density
      real(dp),     intent(out) :: dedgd ! in terms of the gradient of the density
      real(dp),     intent(out) :: dedtau! and in terms of tau
    end subroutine xc_f90_mgga
  end interface

  ! the LCAs
  interface
    subroutine xc_f90_lca_init(p, info, functional, nspin)
      integer (kind=8), intent(out) :: p
      integer (kind=8), intent(out) :: info
      integer,   intent(in)  :: functional
      integer,   intent(in)  :: nspin
    end subroutine xc_f90_lca_init

    subroutine xc_f90_lca(p, rho, v, e, dedd, dedv)
      use myconstants
      integer (kind=8), intent(in)  :: p
      real(dp),     intent(in)  :: rho   ! rho(nspin) the density
      real(dp),     intent(in)  :: v     ! v(3,nspin) the vorticity
      real(dp),     intent(out) :: e     ! the energy per unit particle
      real(dp),     intent(out) :: dedd  ! dedd(nspin) the derivative of the energy
      ! in terms of the density
      real(dp),     intent(out) :: dedv  ! and in terms of the vorticity
    end subroutine xc_f90_lca
  end interface

end module lib_xc_m

module xc_functionals

  use myconstants
  use lib_xc_m
  implicit none

  type xc_type
     logical :: has_grad
     integer :: x_id
     integer :: c_id
     integer :: x_fam
     integer :: c_fam

     integer :: x_irel
     integer :: nspin

     integer :: lb94_j
     real(dp) :: lb94_alpha
     real(dp) :: xalpha

     integer (kind=8) :: x_func ! the pointer used to call the library
     integer (kind=8) :: x_info ! information about the functional
     integer (kind=8) :: c_func ! the pointer used to call the library
     integer (kind=8) :: c_info ! information about the functional

  end type xc_type

  type(xc_type) :: xc_dft

contains

  subroutine xc_init(nspin, x_id, c_id, lb94_j, lb94_alpha, xalpha, x_rel, xc_model)
    !-----------------------------------------------------------------------!
    ! Initializes exchange-correlation model data. The exchange-correlation !
    ! models are read from the input file.                                  !
    !-----------------------------------------------------------------------!
    use lib_xc_m
    implicit none
    type(xc_type), intent(out) :: xc_model
    integer, intent(in) :: nspin, x_id, c_id, lb94_j
    real(dp), intent(in) :: lb94_alpha, xalpha
    logical, intent(in) :: x_rel

    character (len=800) :: lastwords

    xc_model%nspin = nspin
    xc_model%x_irel = XC_NON_RELATIVISTIC
    xc_model%x_id = x_id
    xc_model%c_id = c_id
    xc_model%has_grad = .false.

    if(xc_model%x_id /= 0) then
      ! get the family of the functional
      xc_model%x_fam = xc_f90_family_from_id(xc_model%x_id)

      if(xc_model%x_fam == XC_FAMILY_UNKNOWN .or. xc_model%x_fam == XC_FAMILY_MGGA) then
        write(lastwords, '(a,i3,a,a)') "'", xc_model%x_id, &
             "' is not a known exchange functional!", &
             "Please check the manual for a list of possible values."
        call die(0,lastwords)
      end if
      if (xc_model%x_fam == XC_FAMILY_GGA) xc_model%has_grad = .true.
    end if

   if(xc_model%c_id /= 0) then
      ! get the family of the functional
      xc_model%c_fam = xc_f90_family_from_id(xc_model%c_id)

      if(xc_model%c_fam == XC_FAMILY_UNKNOWN .or. xc_model%c_fam == XC_FAMILY_MGGA .or. xc_model%c_id == XC_LDA_C_AMGB) then
        write(lastwords, '(a,i3,a)') "'", xc_model%c_id, &
             "' is not a known correlation functional!", &
             "Please check the manual for a list of possible values."
        call die(0,lastwords)
     end if
      if (xc_model%c_fam == XC_FAMILY_GGA) xc_model%has_grad = .true.
    end if

    !Extra variables
    if (xc_model%x_fam == XC_FAMILY_LDA) then
       if (x_rel) xc_model%x_irel = XC_RELATIVISTIC
    end if
    if (xc_model%x_fam == XC_FAMILY_GGA .and. xc_model%x_id == XC_GGA_XC_LB) then
       xc_model%lb94_j = lb94_j
       xc_model%lb94_alpha = lb94_alpha
    end if
    if (xc_model%c_fam == XC_FAMILY_LDA .and. xc_model%c_id == XC_LDA_C_XALPHA) then
       xc_model%xalpha = xalpha
    end if

    !Initialize
    call x_init(xc_model)
    call c_init(xc_model)

  end subroutine xc_init

  subroutine x_init(xc_model)
    !-----------------------------------------------------------------------!
    ! Initialize the libxc objects of the exchange part of xc_model.        !
    !-----------------------------------------------------------------------!
    use lib_xc_m
    implicit none
    type(xc_type), intent(inout) :: xc_model

    select case(xc_model%x_fam)
    case(XC_FAMILY_LDA)
      call xc_f90_lda_init(xc_model%x_func, xc_model%x_info, XC_LDA_X, xc_model%nspin, 3, xc_model%x_irel)

    case(XC_FAMILY_GGA)
      if(xc_model%x_id == XC_GGA_XC_LB) then
        call xc_f90_gga_init(xc_model%x_func, xc_model%x_info, xc_model%x_id, xc_model%nspin, xc_model%lb94_j, xc_model%lb94_alpha)
      else
        call xc_f90_gga_init(xc_model%x_func, xc_model%x_info, xc_model%x_id, xc_model%nspin)
      end if

    end select

  end subroutine x_init

  subroutine c_init(xc_model)
    !-----------------------------------------------------------------------!
    ! Initialize the libxc objects of the correlation part of xc_model.     !
    !-----------------------------------------------------------------------!
    type(xc_type), intent(inout) :: xc_model

    select case(xc_model%c_fam)
    case(XC_FAMILY_LDA)

       if(xc_model%c_id /= XC_LDA_C_XALPHA) then
          call xc_f90_lda_init(xc_model%c_func, xc_model%c_info, xc_model%c_id, xc_model%nspin)
       else
          call xc_f90_lda_init(xc_model%c_func, xc_model%c_info, XC_LDA_C_XALPHA, xc_model%nspin, 3, xc_model%xalpha)
       end if

    case(XC_FAMILY_GGA)
       call xc_f90_gga_init(xc_model%c_func, xc_model%c_info, xc_model%c_id, xc_model%nspin)

    end select

  end subroutine c_init

  subroutine xc_end(xc_model)
    !-----------------------------------------------------------------------!
    ! Frees all memory associated to the xc_model.                          !
    !-----------------------------------------------------------------------!
    type(xc_type), intent(inout) :: xc_model

    select case (xc_model%x_fam)
    case (XC_FAMILY_LDA); call xc_f90_lda_end(xc_model%x_func)
    case (XC_FAMILY_GGA); call xc_f90_gga_end(xc_model%x_func)
    end select

    select case (xc_model%c_fam)
    case (XC_FAMILY_LDA); call xc_f90_lda_end(xc_model%c_func)
    case (XC_FAMILY_GGA); call xc_f90_gga_end(xc_model%c_func)
    end select

    xc_model%nspin  = 0
    xc_model%x_id   = 0; xc_model%c_id = 0
    xc_model%x_fam  = 0; xc_model%c_fam = 0
    xc_model%x_irel = 0
    xc_model%lb94_j = 0
    xc_model%lb94_alpha = zero
    xc_model%xalpha = zero

  end subroutine xc_end

  function x_model_name(xc_model)
    !-----------------------------------------------------------------------!
    ! Returns the name of the exchange functional.                          !
    !-----------------------------------------------------------------------!
    type(xc_type), intent(in) :: xc_model
    character(120) :: x_model_name

    if (xc_model%x_id == 0) then
      x_model_name = "None"
    else
      call xc_f90_info_name(xc_model%x_info, x_model_name)
    end if

  end function x_model_name

  function c_model_name(xc_model)
    !-----------------------------------------------------------------------!
    ! Returns the name of the correlation functional.                       !
    !-----------------------------------------------------------------------!
    type(xc_type), intent(in) :: xc_model
    character(120) :: c_model_name

    if (xc_model%c_id == 0) then
       c_model_name = "None"
    else
       call xc_f90_info_name(xc_model%c_info, c_model_name)
    end if

  end function c_model_name
  !===================================================================
  !
  ! Calculates the exchange potential and exchange energy for a given
  ! electron density. Both input electron density and output potential
  ! are calculated on real space in the irreducible wedge.
  !
  ! INPUT:
  !    vxc(:,s) : electron density for spin channel s, in units of (bohr)^-3
  !
  ! OUTPUT:
  !    vxc(:,s) : exchange potential for spin channel s, in units of Ry.
  !    exc      : exchange energy, in units of Ry.
  !
  !-------------------------------------------------------------------
  subroutine vx_get(gvec,xc_model,nspin,ntrans,vxc,exc)

    use myconstants
    use typedefs
    use lib_xc_m
    implicit none

    ! arguments
    ! real-space grid
    type (gspace), intent(inout) :: gvec
    ! exchange-correlation functional
    type(xc_type), intent(in) :: xc_model
    ! number of spin channels and number of irreducible representations
    integer, intent(in) :: nspin, ntrans
    real(dp), intent(inout) :: vxc(gvec%nr,nspin)
    real(dp), intent(out) :: exc

    ! local variables
    ! gradient of electron density.
    real(dp), dimension (:,:,:), allocatable :: grad
    ! minus laplacian of electron density.
    real(dp), dimension (:,:), allocatable :: mlapl
    ! energy derivative, as defined in libxc
    real(dp), dimension (:,:), allocatable :: vsigma
    ! gradient of energy derivative, Grad[ vsigma ]
    real(dp), dimension (:,:,:), allocatable :: gradvs
    ! allocation check.
    integer :: alcstat
    ! local number of grid points.
    integer :: ndim
    ! counters.
    integer :: ii, isp
    ! temporary variables.
    real(dp) :: density(nspin), sigma(2*nspin-1), vsigma_l(2*nspin-1), &
         ex_tmp, vxc_tmp(nspin)
    ! spin constants:
    integer, parameter :: UP = 1, DOWN = 2

    !---------------------------------------------------------------
    ! Initialize internal variables.
    ndim = gvec%nr
    exc = zero

    !---------------------------------------------------------------
    ! With GGA, we must calculate the gradient and laplacian of the
    ! electron density (rho==exc) at each grid point. Notice that
    ! get_lap calculates the negative of the laplacian.
    !
    if (xc_model%has_grad) then
       allocate(mlapl(ndim,nspin),stat=alcstat)
       call alccheck('mlapl','get_xc',ndim*nspin,alcstat)
       allocate(grad(3,ntrans*ndim,nspin),stat=alcstat)
       call alccheck('grad','get_xc',3*ntrans*ndim*nspin,alcstat)
       do isp = 1, nspin
          call get_lap_sym(gvec,vxc(1,isp),mlapl(1,isp))
          call get_grad_sym(gvec,vxc(1,isp),grad(1,1,isp))
       enddo
       ! Allocate auxiliary arrays.
       allocate(vsigma(ndim,2*nspin-1),stat=alcstat)
       call alccheck('vsigma','get_xc',ndim*(2*nspin-1),alcstat)
       vsigma = zero
       allocate(gradvs(3,ntrans*ndim,2*nspin-1),stat=alcstat)
       call alccheck('gradvs','get_xc',3*ntrans*ndim*(2*nspin-1),alcstat)
    endif

    !---------------------------------------------------------------
    ! Calculate density part of exchange.
    !
    do ii = 1, ndim
       density(1:nspin) = vxc(ii,1:nspin)
       if (xc_model%has_grad) then
          sigma = zero
          sigma(1) = dot_product(grad(:,ii,UP),grad(:,ii,UP))
          if (nspin > 1) then
             sigma(2) = dot_product(grad(:,ii,UP),grad(:,ii,DOWN))
             sigma(3) = dot_product(grad(:,ii,DOWN),grad(:,ii,DOWN))
          endif
       endif

       select case(xc_model%x_fam)
       case(XC_FAMILY_LDA)
          call xc_f90_lda_vxc(xc_model%x_func, density(1), ex_tmp, vxc_tmp(1))
       case(XC_FAMILY_GGA)
          if(xc_model%x_id == XC_GGA_XC_LB) then
             ! Not yet implemented.
          else
             call xc_f90_gga(xc_model%x_func, density(1), sigma(1), ex_tmp, &
                  vxc_tmp(1), vsigma_l(1))
          end if
       case default
          ex_tmp = zero
          vxc_tmp = zero
          if (xc_model%has_grad) vsigma = zero
       end select
       exc = exc + ex_tmp
       vxc(ii,1:nspin) = vxc_tmp(1:nspin)
       if (xc_model%has_grad) then
          vsigma(ii,1) = vsigma_l(1)
          if (nspin > 1) vsigma(ii,2:3) = vsigma_l(2:3)
       endif
    enddo

    !---------------------------------------------------------------
    ! Calculate gradient part of exchange.
    !
    if (xc_model%has_grad) then
       do isp = 1, 2*nspin - 1
          call get_grad_sym(gvec,vsigma(1,isp),gradvs(1,1,isp))
       enddo
       do ii = 1, ndim
             vxc(ii,1) = vxc(ii,1) + two*vsigma(ii,1)*mlapl(ii,1) - &
                  two*dot_product(gradvs(:,ii,1),grad(:,ii,1))
          if (nspin > 1) then
             vxc(ii,1) = vxc(ii,1) + vsigma(ii,2)*mlapl(ii,2) - &
                  dot_product(gradvs(:,ii,2),grad(:,ii,2))
             vxc(ii,2) = vxc(ii,2) + two*vsigma(ii,3)*mlapl(ii,2) - &
                  two*dot_product(gradvs(:,ii,3),grad(:,ii,2)) + &
                  vsigma(ii,2)*mlapl(ii,1) - &
                  dot_product(gradvs(:,ii,2),grad(:,ii,1))
          endif
       enddo
    endif
    ! Convert units hartree -> rydbergs.
    vxc = vxc * two
    exc = exc * two * real(ntrans,dp)

    ! Free up memory.
    if (xc_model%has_grad) then
       deallocate(grad)
       deallocate(mlapl)
       deallocate(gradvs)
       deallocate(vsigma)
    endif

  end subroutine vx_get
  !===================================================================
  !
  ! Calculates the exchange-correlation potential and exchange-correlation
  ! energy for a given electron density. Both input electron density and
  ! output potential are calculated on real space in the irreducible wedge.
  !
  ! INPUT:
  !    vxc(:,s) : electron density for spin channel s, in units of (bohr)^-3
  !
  ! OUTPUT:
  !    vxc(:,s) : exchange-correlation potential for spin channel s, in units of Ry.
  !    exc      : exchange-correlation energy, in units of Ry.
  !
  !-------------------------------------------------------------------
  subroutine vxc_get(gvec,xc_model,nspin,ntrans,vxc,exc)

    use myconstants
    use typedefs
    use lib_xc_m
    implicit none

    ! arguments
    ! real-space grid
    type (gspace), intent(inout) :: gvec
    ! exchange-correlation functional
    type(xc_type), intent(in) :: xc_model
    ! number of spin channels and number of irreducible representations
    integer, intent(in) :: nspin, ntrans
    real(dp), intent(inout) :: vxc(gvec%nr,nspin)
    real(dp), intent(out) :: exc

    ! local variables
    ! gradient of electron density.
    real(dp), dimension (:,:,:), allocatable :: grad
    ! minus laplacian of electron density.
    real(dp), dimension (:,:), allocatable :: mlapl
    ! energy derivative, as defined in libxc
    real(dp), dimension (:,:), allocatable :: vsigma
    ! gradient of energy derivative, Grad[ vsigma ]
    real(dp), dimension (:,:,:), allocatable :: gradvs
    ! allocation check.
    integer :: alcstat
    ! local number of grid points.
    integer :: ndim
    ! counters.
    integer :: ii, isp
    ! temporary variables.
    real(dp) :: density(nspin), sigma(2*nspin-1), vsigma_l(2*nspin-1), &
         ex_tmp, vxc_tmp(nspin)
    ! spin constants:
    integer, parameter :: UP = 1,DOWN = 2

    !---------------------------------------------------------------
    ! Initialize internal variables.
    ndim = gvec%nr
    exc = zero

    !---------------------------------------------------------------
    ! With GGA, we must calculate the gradient of the electron
    ! density (rho==exc) at each grid point.
    !
    if (xc_model%has_grad) then
       allocate(mlapl(ndim,nspin),stat=alcstat)
       call alccheck('mlapl','get_xc',ndim*nspin,alcstat)
       allocate(grad(3,ntrans*ndim,nspin),stat=alcstat)
       call alccheck('grad','get_xc',3*ntrans*ndim*nspin,alcstat)
       do isp = 1, nspin
          call get_lap_sym(gvec,vxc(1,isp),mlapl(1,isp))
          call get_grad_sym(gvec,vxc(1,isp),grad(1,1,isp))
       enddo
       ! Allocate auxiliary arrays.
       allocate(vsigma(ndim,2*nspin-1),stat=alcstat)
       call alccheck('vsigma','get_xc',ndim*(2*nspin-1),alcstat)
       vsigma = zero
       allocate(gradvs(3,ntrans*ndim,2*nspin-1),stat=alcstat)
       call alccheck('gradvs','get_xc',3*ntrans*ndim*(2*nspin-1),alcstat)
    endif

    !---------------------------------------------------------------
    ! Calculate density part of exchange-correlation.
    !
    do ii = 1, ndim
       density(1:nspin) = vxc(ii,1:nspin)
       if (xc_model%has_grad) then
          sigma = zero
          sigma(1) = dot_product(grad(:,ii,UP),grad(:,ii,UP))
          if (nspin > 1) then
             sigma(2) = dot_product(grad(:,ii,UP),grad(:,ii,DOWN))
             sigma(3) = dot_product(grad(:,ii,DOWN),grad(:,ii,DOWN))
          endif
       endif

       ! Exchange.
       select case(xc_model%x_fam)
       case(XC_FAMILY_LDA)
          call xc_f90_lda_vxc(xc_model%x_func, density(1), ex_tmp, vxc_tmp(1))
       case(XC_FAMILY_GGA)
          if(xc_model%x_id == XC_GGA_XC_LB) then
             ! Not yet implemented.
          else
             call xc_f90_gga(xc_model%x_func, density(1), sigma(1), ex_tmp, &
                  vxc_tmp(1), vsigma_l(1))
          end if
       case default
          ex_tmp = zero
          vxc_tmp = zero
          if (xc_model%has_grad) vsigma = zero
       end select
       exc = exc + ex_tmp * sum(density)
       vxc(ii,1:nspin) = vxc_tmp(1:nspin)
       if (xc_model%has_grad) then
          vsigma(ii,1) = vsigma_l(1)
          if (nspin > 1) vsigma(ii,2:3) = vsigma_l(2:3)
       endif

       ! Correlation.
       select case(xc_model%c_fam)
       case(XC_FAMILY_LDA)
          call xc_f90_lda_vxc(xc_model%c_func, density(1), ex_tmp, vxc_tmp(1))
       case(XC_FAMILY_GGA)
          if(xc_model%c_id == XC_GGA_XC_LB) then
             ! Not yet implemented.
          else
             call xc_f90_gga(xc_model%c_func, density(1), sigma(1), ex_tmp, &
                  vxc_tmp(1), vsigma_l(1))
          end if

       case default
          ex_tmp = zero
          vxc_tmp = zero
          if (xc_model%has_grad) vsigma = zero
       end select
       exc = exc + ex_tmp * sum(density)
       vxc(ii,1:nspin) = vxc(ii,1:nspin) + vxc_tmp(1:nspin)
       if (xc_model%has_grad) then
          vsigma(ii,1) = vsigma(ii,1) + vsigma_l(1)
          if (nspin > 1) vsigma(ii,2:3) = vsigma(ii,2:3) + vsigma_l(2:3)
       endif
    enddo
    !---------------------------------------------------------------
    ! Calculate gradient part of exchange-correlation.
    !
    if (xc_model%has_grad) then
       do isp = 1, 2*nspin - 1
          call get_grad_sym(gvec,vsigma(1,isp),gradvs(1,1,isp))
       enddo
       do ii = 1, ndim
             vxc(ii,1) = vxc(ii,1) + two*vsigma(ii,1)*mlapl(ii,1) - &
                  two*dot_product(gradvs(:,ii,1),grad(:,ii,1))
          if (nspin > 1) then
             vxc(ii,1) = vxc(ii,1) + vsigma(ii,2)*mlapl(ii,2) - &
                  dot_product(gradvs(:,ii,2),grad(:,ii,2))
             vxc(ii,2) = vxc(ii,2) + two*vsigma(ii,3)*mlapl(ii,2) - &
                  two*dot_product(gradvs(:,ii,3),grad(:,ii,2)) + &
                  vsigma(ii,2)*mlapl(ii,1) - &
                  dot_product(gradvs(:,ii,2),grad(:,ii,1))
          endif
       enddo
    endif
    ! Convert units hartree -> rydbergs.
    vxc = vxc * two
    exc = exc * two * real(ntrans,dp)

    ! Free up memory.
    if (xc_model%has_grad) then
       deallocate(grad)
       deallocate(mlapl)
       deallocate(gradvs)
       deallocate(vsigma)
    endif

  end subroutine vxc_get
  !===================================================================
  !
  ! Calculates the exchange-correlation kernel (functional derivative of
  ! XC potential) for a given electron density. Both input electron density
  ! and output kernel are calculated on real space in the irreducible wedge.
  ! The functional is assumed to be LDA, or LSDA.
  !
  ! INPUT:
  !    fxc(:,s,1) : electron density for spin channel s, in units of (bohr)^-3
  !    fxc(:,:,2) is not referenced if nspin > 1
  !    kflag : kernel flag
  !           kflag = 3  -> calculate triplet kernel, as defined in Eq. 50
  !                         of Vasiliev, Ogut, and Chelikowsky, Phys. Rev. B
  !                         65, 115416 (2002). It implies nspin = 1.
  !           otherwise  -> calculate singlet kernel or spin polarized kernel.
  !                         see Eq. 5.13 of Onida, Reining, and Rubio, Rev.
  !                         Mod. Phys. 74, 601 (2002).
  !
  ! OUTPUT:
  !    fxc(:,s,s') : exchange-correlation kernel for spin channel s and s',
  !                  in units of Ry* (bohr)^-3.
  !
  !-------------------------------------------------------------------
  subroutine fxc_get(xc_model,nspin,ndim,kflag,fxc)

    use myconstants
    use lib_xc_m
    implicit none

    ! arguments
    ! exchange-correlation functional
    type(xc_type), intent(in) :: xc_model
    integer, intent(in) :: &
         nspin, &    ! number of spin channels
         ndim, &     ! size of real-space grid
         kflag       ! kernel flag (see above)
    real(dp), intent(inout) :: fxc(ndim,nspin,nspin)

    ! local variables
    ! counters.
    integer :: ii
    ! temporary variables.
    real(dp) :: fx(2*nspin-1), fc(2*nspin-1), rho(nspin)

    !---------------------------------------------------------------

    do ii = 1, ndim
       if (kflag == 3) then
          fxc(ii,1,1) = vspin_LDA(fxc(ii,1,1))
       elseif (kflag == 1 .or. kflag == 2) then
          rho(:) = fxc(ii,:,1)
          call xc_f90_lda_fxc(xc_model%x_func,rho(1),fx(1))
          call xc_f90_lda_fxc(xc_model%c_func,rho(1),fc(1))
          fxc(ii,1,1) = fx(1) + fc(1)
          if (nspin > 1) then
             fxc(ii,2,2) = fx(3) + fc(3)
             fxc(ii,1,2) = fx(2) + fc(2)
             fxc(ii,2,1) = fxc(ii,1,2)
          endif
       endif
    enddo
    fxc = fxc * two

  end subroutine fxc_get
  !===================================================================
  !
  ! vspin_LDA(n) = TDLDA potential for triplet states
  ! (see Vasiliev et al, PRB 2002).
  !
  !---------------------------------------------------------------
  function vspin_LDA(n) result(vspin)

    use myconstants
    implicit none

    ! arguments
    real(dp), intent(in) :: n
    real(dp) :: vspin

    ! local variables
    real(dp) :: rs, lrs, exc, excs

    ! Exchange constant, X1 = - ( 9*pi/4 )^(1/3) * 3/( 2*pi )
    real(dp), parameter :: ca_x1 = -0.91633058656629d0
    ! Xf = -( 2 / (9*pi) )^(1/3) * 2
    real(dp), parameter :: ca_xf = -0.827133987866d0
    ! High density (rs<1) constants, LSDA
    real(dp), parameter :: ca_c1 = 0.0622d0
    real(dp), parameter :: ca_c2 = 0.0960d0
    real(dp), parameter :: ca_c3 = 0.0040d0
    real(dp), parameter :: ca_c4 = 0.0232d0
    real(dp), parameter :: ca_c5 = 0.0192d0

    real(dp), parameter :: ca_as =  0.0311d0
    real(dp), parameter :: ca_bs = -0.0538d0
    real(dp), parameter :: ca_cs =  0.0014d0
    real(dp), parameter :: ca_ds = -0.0096d0
    ! Low density (rs>1) constants, LSDA
    real(dp), parameter :: ca_g = -0.2846d0
    real(dp), parameter :: ca_b1 = 1.0529d0
    real(dp), parameter :: ca_b2 = 0.3334d0
    real(dp), parameter :: ca_gs = -0.1686d0
    real(dp), parameter :: ca_b1s = 1.3981d0
    real(dp), parameter :: ca_b2s = 0.2611d0

    vspin = zero

    ! Exit if the density is "zero"
    if (real(n,dp) <= zero) return

    ! rs
    rs = (three/(four*pi*real(n,dp)))**third
    ! exchange
    vspin = four*four/nine/three*ca_x1*(rs**2)*pi
    ! correlation: high density
    if (rs < one) then
       lrs = log(rs)
       exc = ca_c1*lrs - ca_c2 + ca_c3*rs*lrs - ca_c4*rs
       excs = ca_as*lrs + ca_bs + ca_cs*rs*lrs + ca_ds*rs
       ! correlation: low density
    else
       exc = ca_g/(one+ca_b1*sqrt(rs)+ca_b2*rs)
       excs = ca_gs/(one+ca_b1s*sqrt(rs)+ca_b2s*rs)
    endif
    vspin = vspin + four*(excs-exc)/(nine*n*0.259921049890)

  end function vspin_LDA

end module xc_functionals
!===================================================================
