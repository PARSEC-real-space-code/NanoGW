!===================================================================
!
! Read input parameters from file rgwbs.in
! Only search for parameters specific to BSESOLV.
!
! Copyright (C) 2009 Murilo L. Tiago, http://users.ices.utexas.edu/~mtiago
! This file is part of RGWBS. It is distributed under the GPL v1.
!
!-------------------------------------------------------------------
subroutine input_b(bsepol_in,q_bse,writeig,trip_flag,trunc_c,mix_flag, &
     noxchange,snorm,readocc,hf,bsecut,eref,ecutb2)

  use typedefs
  use esdf
  implicit none

  ! arguments
  ! input polarizability data
  type(polinfo), intent(out), dimension(2) :: bsepol_in
  ! input q-vector
  type(qptinfo), intent(out) :: q_bse
  logical, intent(out) :: &
       writeig, &    ! true if BSE eigenvectors are printed in file bse_diag.dat
       trip_flag, &  ! true if BSE is solved for triplet excitations only
       trunc_c, &    ! true if the long-wavelength part of the Coulomb
                     ! interaction is removed
       mix_flag, &   ! true if Tamm-Dancof approximation is not used
       noxchange, &  ! true if exchange part of kernel is ignored
       snorm, &      ! true if sum rule is renormalized
       readocc, &    ! true if orbital occupancies should be read
       hf            ! true if BSE is solved within the Hartree-Fock approximation
  real(dp), intent(out) :: &
       bsecut, &     ! energy cutoff in BSE
       eref, &       ! energy reference at which the direct kernel is calculated
       ecutb2        ! resolution in energy denominators

  ! local variables
  integer :: ii, jj, nlines
  real(dp) :: dtmp
  integer, dimension(maxdata,2) :: vmap_b, cmap_b

  !-----------------------------------------------------------------------
  ! Initialize mandatory input info.
  !
  bsepol_in(:)%nval = -1
  bsepol_in(:)%ncond = -1

  !-----------------------------------------------------------------------
  ! Start searching for keywords in input file.
  !
  call esdf_init('rgwbs.in')

  eref = esdf_physical('energy_reference',-one,'eV')

  bsecut = esdf_physical('bse_cutoff',-one,'eV')
  bsecut = bsecut/ryd

  writeig = (.not. esdf_defined('no_eigenvectors'))

  trip_flag = esdf_defined('bse_triplet_kernel')

  trunc_c =  (esdf_defined('truncate_coulomb'))

  mix_flag = esdf_defined('use_mixing')

  noxchange = esdf_defined('no_exchange')

  snorm = esdf_defined('renormalize_sumrule')

  readocc = esdf_defined('read_orbital_occupancies')

  if ( esdf_reduce(esdf_string('exchange_correlation','gw')) == &
       'hartree_fock' ) then
     hf = .true.
  else
     hf = .false.
  endif

  ecutb2 = esdf_physical('bse_energy_resolution',zero,'eV')

  if (esdf_block('bse_valence',nlines)) then
     if (bsepol_in(1)%nval < 0) bsepol_in(1)%nval = 0
     call esdf_parse_block('bse_valence', &
          nlines,bsepol_in(1)%nval,maxdata,vmap_b(1,1))
     bsepol_in(2)%nval = bsepol_in(1)%nval
     vmap_b(1:bsepol_in(1)%nval,2) = vmap_b(1:bsepol_in(1)%nval,1)
  endif

  if (esdf_block('bse_valence_up',nlines)) then
     if (bsepol_in(1)%nval < 0) bsepol_in(1)%nval = 0
     call esdf_parse_block('bse_valence_up',nlines, &
          bsepol_in(1)%nval,maxdata,vmap_b(1,1))
  endif

  if (esdf_block('bse_valence_down',nlines)) then
     if (bsepol_in(2)%nval < 0) bsepol_in(2)%nval = 0
     call esdf_parse_block('bse_valence_down',nlines, &
          bsepol_in(2)%nval,maxdata,vmap_b(1,2))
  endif

  if (esdf_block('bse_conduction',nlines)) then
     if (bsepol_in(1)%ncond < 0) bsepol_in(1)%ncond = 0
     call esdf_parse_block('bse_conduction',nlines, &
          bsepol_in(1)%ncond,maxdata,cmap_b(1,1))
     bsepol_in(2)%ncond = bsepol_in(1)%ncond
     cmap_b(1:bsepol_in(1)%ncond,2) = cmap_b(1:bsepol_in(1)%ncond,1)
  endif

  if (esdf_block('bse_conduction_up',nlines)) then
     if (bsepol_in(1)%ncond < 0) bsepol_in(1)%ncond = 0
     call esdf_parse_block('bse_conduction_up',nlines, &
          bsepol_in(1)%ncond,maxdata,cmap_b(1,1))
  endif

  if (esdf_block('bse_conduction_down',nlines)) then
     if (bsepol_in(2)%ncond < 0) bsepol_in(2)%ncond = 0
     call esdf_parse_block('bse_conduction_down',nlines, &
          bsepol_in(2)%ncond,maxdata,cmap_b(1,2))
  endif

  if(esdf_block('qpoints_bse',q_bse%nk)) then
     if (q_bse%nk > 1) call die( &
          'ERROR: number of q-points in BSE must not be greater than one. Stop.')
     allocate(q_bse%fk(3,q_bse%nk))
     allocate(q_bse%zerok(q_bse%nk))
     do ii = 1, q_bse%nk
        read(block_data(ii),*) (q_bse%fk(jj,ii),jj=1,3), dtmp, jj
        q_bse%fk(:,ii) = q_bse%fk(:,ii) / dtmp
        if (jj > 0) then
           q_bse%zerok(ii) = .true.
        else
           q_bse%zerok(ii) = .false.
        endif
     enddo
  else
     q_bse%nk = 1
     allocate(q_bse%fk(3,q_bse%nk))
     q_bse%fk = zero
     allocate(q_bse%zerok(q_bse%nk))
     q_bse%zerok = .true.
  endif

  call esdf_close

  !-----------------------------------------------------------------------
  ! Check if indices of valence and conduction states are defined.
  !
  do ii = 1, 2
     if (bsepol_in(ii)%nval /= 0) then
        allocate(bsepol_in(ii)%vmap(bsepol_in(ii)%nval))
        bsepol_in(ii)%vmap(1:bsepol_in(ii)%nval) = &
             vmap_b(1:bsepol_in(ii)%nval,ii)
     endif
     if (bsepol_in(ii)%ncond /= 0) then
        allocate(bsepol_in(ii)%cmap(bsepol_in(ii)%ncond))
        bsepol_in(ii)%cmap(1:bsepol_in(ii)%ncond) = &
             cmap_b(1:bsepol_in(ii)%ncond,ii)
     endif
  enddo

  return

end subroutine input_b
!===================================================================
