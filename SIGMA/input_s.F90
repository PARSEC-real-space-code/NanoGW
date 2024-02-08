!===================================================================
!
! Read input parameters from file rgwbs.in
! Only search for parameters specific to BSESOLV.
!
! Copyright (C) 2009 Murilo L. Tiago, http://users.ices.utexas.edu/~mtiago
! This file is part of RGWBS. It is distributed under the GPL v1.
!
!-------------------------------------------------------------------
subroutine input_s(sig,kpt_sig,snorm,writeqp,readvxc,readocc,cohsex, &
     nooffd, hqp_sym,n_it,chkpt,static_type,sig_en, &
     max_conv,xbuff,ecuts,qpmix,sig_cut,verbose)

  use typedefs
  use esdf
  use mpi_module
  implicit none

  ! arguments
  ! self-energy data
  type(siginfo), intent(out) :: sig
  ! k-point type for self-energy
  type(kptinfo), intent(out) :: kpt_sig
  logical, intent(out) :: &
       snorm, &      ! true if sum rule is renormalized
       writeqp, &    ! true if QP eigenvectors are printed in file parsec_qp.dat
       readvxc, &    ! true if Vxc matrix elements are read from file sigma_mtxel.dat
       readocc, &    ! true if orbital occupancies should be read
       cohsex, &     ! true if COHSEX (=static) approximation is used
       nooffd, &     ! true if only diagonal of self energy matrix to be computed
       hqp_sym       ! true if H_qp is symmetrized
  logical, intent(in)  :: &
       verbose       ! true for printing out extra info
  integer, intent(out) :: &
       n_it, &       ! number of self-consistent iterations
       chkpt, &      ! checkpoint flag
       sig_en, &     ! number of energy points where self-energy is calculated
       static_type   ! type of static correction
  real(dp), intent(out) :: &
       max_conv, &   ! tolerance in convergence test, during self-consistency 
       xbuff, &      ! size of buffer arrays in wpol0
       ecuts, &      ! resolution in energy denominators
       qpmix, &      ! mixing factor used in self-consistency
       sig_cut       ! cut-off factor in the update of self-energy

  ! local variables
  logical :: ifound
  character (len=800) :: tstring
  integer :: ii, i1, i2, jj, nlines, min_s, max_s
  integer :: diag(maxdata), offdiag(maxdata,2)
  integer, parameter :: MIN_INIT = 88888888
  real(dp) :: dtmp
  logical, dimension(:,:), allocatable :: inv_offd

  !-----------------------------------------------------------------------
  ! Initialize mandatory input info.
  !
  sig%ndiag = -1
  sig%ndiag_s = -1
  sig%noffd = -1
  sig%noffd_s = -1
  sig%nmap = -1

  !-----------------------------------------------------------------------
  ! Start searching for keywords in input file.
  !
  call esdf_init('rgwbs.in')

  sig%nmax_c = esdf_integer('max_number_states',-1)

  min_s = esdf_integer('qp_min',MIN_INIT)
  max_s = esdf_integer('qp_max',-1)

  sig%deltae = esdf_physical('energy_range',20.d0,'eV')
  sig%deltae = sig%deltae / ryd

  sig%nen = esdf_integer('energy_data_points',21)

  snorm = esdf_defined('renormalize_sumrule')

  max_conv = esdf_physical('convergence_potential',1.d-3,'eV')

  n_it = esdf_integer('number_iterations',0)

  chkpt = -1
  ifound = esdf_defined('read_checkpoint')
  if (n_it == 0) ifound = .true.
  if (ifound) chkpt = 0

  writeqp = esdf_defined('write_qp_wavefunctions')

  readvxc = esdf_defined('read_vxc_matrix_elements')

  readocc = esdf_defined('read_orbital_occupancies')

  cohsex = esdf_defined('cohsex_approximation')

  nooffd = esdf_defined('only_diagonal')

  tstring = esdf_reduce(esdf_string('exchange_correlation','gw'))
  select case (trim(tstring))
  case ('gw')
     sig%xc = XC_GW
  case ('hartree_fock')
     sig%xc = XC_HF
  case ('b3lyp')
     sig%xc = XC_B3LYP
  case ('blyp')
     sig%xc = XC_GGA_BLYP
  case ('lda_ca')
     sig%xc = XC_LDA_CA
  case ('gga_pbe')
     sig%xc = XC_GGA_PBE
  case default
     write(6,*) 'ERROR: unknown exchange-correlation type:', trim(tstring)
     call die('STOP.')
  end select

  if (sig%xc /= XC_GW) chkpt = -1

  xbuff = esdf_physical('scratch_disk_size',zero,'MB')

  ecuts = esdf_physical('dynamic_energy_resolution',0.005d0,'Ry')

  qpmix = esdf_double('qp_mixing_param',one)

  sig_cut = esdf_physical('self_energy_cutoff',zero,'eV')

  hqp_sym = (.not. esdf_defined('no_hqp_symmetrize'))

  static_type = esdf_integer('static_type',2)
  if ( static_type .lt. 1 .or. static_type .gt. 2) then
     write(6, *) "WARNINING: static_type is not right."
     write(6, *) "   Please choose static_type = 1 or 2."
     static_type = 2
     write(6, *) "   Program will use the default value 2."
  endif

  tstring = esdf_reduce(esdf_string('sigma_energy','left'))
  select case (trim(tstring))
  case ('left')
     sig_en = SIG_LEFT
  case ('right')
     sig_en = SIG_RIGHT
  case ('average')
     sig_en = SIG_AV
  case default
     write(6,*) 'ERROR: unknown choice of sigma_energy :', trim(tstring)
     call die('STOP.')
  end select

  if(esdf_block('sigma_kpoints',kpt_sig%nk))then
     allocate(kpt_sig%fk(3,kpt_sig%nk))
     do ii = 1, kpt_sig%nk
        read(block_data(ii),*) (kpt_sig%fk(jj,ii),jj=1,3), dtmp
        kpt_sig%fk(:,ii) = kpt_sig%fk(:,ii) / dtmp
     enddo
  else
     kpt_sig%nk = 1
     allocate(kpt_sig%fk(3,kpt_sig%nk))
     kpt_sig%fk = zero
  endif

  if (esdf_block('diag',nlines)) then
     if (sig%ndiag < 0) sig%ndiag = 0
     call esdf_parse_block('diag',nlines,sig%ndiag,maxdata,diag)
  endif

  if (esdf_block('offdiag',nlines)) then
     if (sig%noffd < 0) sig%noffd = 0
     call esdf_parse_block_2('offdiag',nlines,sig%noffd,maxdata,offdiag)
  endif

  call esdf_close

  !-----------------------------------------------------------------------
  ! Check if indices of diagonal/off-diagonal matrix elements of sigma
  ! are defined.
  !

  if (nooffd) then
     offdiag = -1
     sig%noffd = 0
  endif

  if (sig%ndiag > 0) min_s = min(min_s,minval(diag(1:sig%ndiag)))
  if (sig%noffd > 0) min_s = min(min_s,minval(offdiag(1:sig%noffd,:)))
  if (sig%ndiag > 0) max_s = max(max_s,maxval(diag(1:sig%ndiag)))
  if (sig%noffd > 0) max_s = max(max_s,maxval(offdiag(1:sig%noffd,:)))

  if (min_s == MIN_INIT) min_s = 1

  if (max_s > maxdata) call die('ERROR: Too many matrix elements '// &
       'Increase value of internal parameter maxdata in typedefs.F90.')
  if (max_s > 0) then
     sig%nmap = max_s - min_s + 1
     allocate(sig%map(sig%nmap))
     do ii = 1, sig%nmap
        sig%map(ii) = ii + min_s - 1
     enddo

     sig%ndiag_s = max(sig%ndiag,0)
     do ii = min_s, max_s
        ifound = .false.
        do jj = 1, sig%ndiag
           if (diag(jj) == ii) then
              ifound = .true.
              exit
           endif
        enddo
        if (.not. ifound) then
           sig%ndiag_s = sig%ndiag_s + 1
           diag(sig%ndiag_s) = ii
        endif
     enddo

     if (nooffd) then
        sig%noffd_s = 0
     else
        sig%noffd_s = max(sig%noffd,0)
        allocate(inv_offd(min_s:max_s,min_s:max_s))
        inv_offd = .true.
        do jj = 1, sig%noffd
           inv_offd(offdiag(jj,1),offdiag(jj,2)) = .false.
           inv_offd(offdiag(jj,2),offdiag(jj,1)) = .false.
        enddo
        do i1 = min_s, max_s
           do i2 = i1 + 1, max_s
              if (inv_offd(i1,i2)) then
                 sig%noffd_s = sig%noffd_s + 1
                 offdiag(sig%noffd_s,1) = i1
                 offdiag(sig%noffd_s,2) = i2
              endif
           enddo
        enddo
        deallocate(inv_offd)
     endif
  endif

  if(sig%ndiag_s > 0) then
     allocate(sig%diag(sig%ndiag_s))
     do ii = 1, sig%ndiag_s
        jj = 0
        do jj = 1, sig%nmap
           if (sig%map(jj) == diag(ii)) exit
        enddo
        sig%diag(ii) = jj
     enddo
  endif

  if(sig%noffd_s > 0) then
     allocate(sig%off1(sig%noffd_s))
     allocate(sig%off2(sig%noffd_s))
     do ii=1,sig%noffd_s
        jj = 0
        do jj = 1, sig%nmap
           if (sig%map(jj) == offdiag(ii,1)) exit
        enddo
        sig%off1(ii) = jj
        jj = 0
        do jj = 1, sig%nmap
           if (sig%map(jj) == offdiag(ii,2)) exit
        enddo
        sig%off2(ii) = jj
     enddo
  endif

  if(peinf%master .and. verbose) then
     write(6,'(" in input_s() ")')
     write(6,'("    ii    sig%map")')
     do ii = 1, sig%nmap
        write(6,'(i8,i8)') ii,sig%map(ii)
     enddo
   
     write(6,'("    ii    sig%diag")')
     do ii = 1, sig%ndiag_s
        write(6,'(i8,i8)') ii,sig%diag(ii)
     enddo
   
     write(6,'("    ii    sig%offdiag")')
     do ii = 1, sig%noffd_s
        write(6,'(i8,i8,i8)') ii,sig%off1(ii),sig%off2(ii)
     enddo
  endif
  call MPI_Barrier(peinf%comm,ii)

end subroutine input_s
!===================================================================
