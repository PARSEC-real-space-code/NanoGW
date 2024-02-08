!===================================================================
!
! Read input parameters from file rgwbs.in
! This is a parsing subroutine that searches for meaningful keywords
! in input file and initializes parameters according to input.
!
! If no file is found, use default values and quietly fold its tent...
!
! Copyright (C) 2009 Murilo L. Tiago, http://users.ices.utexas.edu/~mtiago
! This file is part of RGWBS. It is distributed under the GPL v1.
!
!-------------------------------------------------------------------
subroutine input_g(pol_in,qpt,tdldacut,nbuff,lcache,wgr_npes, &
      nolda,tamm_d,rgr_num,dft_code,doisdf,n_intp,intp_type,isdf_type, &
      lessmemory, fastselect, verbose)

  use typedefs
  use esdf
  use mpi_module
  implicit none

  ! arguments
  ! polarizability type, stores ranges of orbitals and number of orbitals
  type(polinfo), intent(out), dimension(2) :: pol_in
  ! number and coordinates of q-points
  type (qptinfo), intent(out) :: qpt
  integer, intent(out) :: &
       nbuff, &      ! length of buffer arrays in kernel
       lcache, &     ! length of cache arrays in integrates
       wgr_npes, &   ! number of PEs in wfn groups
       rgr_num, &    ! number of groups where representations/q-points
                     ! are distributed
       dft_code      ! specification of DFT code

  real(dp), intent(out) :: tdldacut  ! energy cutoff in TDLDA
  logical, intent(out) :: &
       nolda, &      ! true if LDA kernel is not used
       tamm_d        ! true if Tamm-Dancof approximation is used
  logical, intent(in) :: verbose ! if true, prinout extra debug info

  ! local variables
  character (len=80) :: strflag
  integer :: ii, jj, nlines
  real(dp) :: sbuff, dtmp
  integer, dimension(maxdata,2) :: vmap, cmap

  ! variables for ISDF method
  logical, intent(out) :: doisdf, lessmemory, fastselect
  integer, intent(out) :: n_intp
  integer, intent(out) :: intp_type, isdf_type

  !-----------------------------------------------------------------------
  ! Initialize input info.
  !
  pol_in(:)%nval = -1
  pol_in(:)%ncond = -1

  !-----------------------------------------------------------------------
  ! Start reading input parameters and search for ESDF keywords.
  !
  call esdf_init('rgwbs.in')

  dtmp = -one
  tdldacut = esdf_physical('tdlda_cutoff',dtmp,'eV')
  tdldacut = tdldacut/ryd

  sbuff = esdf_physical('buffer_size',zero,'MB')
  if (sbuff == zero) then
     nbuff = 50000
  else
     nbuff = nint(1024.0 * 128.0 * sbuff)
  endif

  lcache = esdf_integer('cache_size',4000)

  nolda = esdf_defined('no_lda_kernel')

  tamm_d = esdf_defined('tamm_dancoff')

  ! read variables for ISDF method
  doisdf = (esdf_defined('doisdf'))

  lessmemory = (esdf_defined('lessmemory'))

  fastselect = (esdf_defined('fastselect'))

  n_intp = esdf_integer('num_isdf_points',-1)

  intp_type = esdf_integer('intp_type',1)

  isdf_type = esdf_integer('isdf_type',1)

  ! ------

  ii = 1
  rgr_num = esdf_integer('distribute_representations',ii)

  strflag = esdf_reduce(esdf_string('dft_program', 'parsec'))
  select case (trim(strflag))
  case ('parsec','PARSEC')
     dft_code = PARSEC
  case ('paratec','PARATEC')
     dft_code = PARATEC
  case default
     write(6,*) 'ERROR: unknown DFT program : ',trim(strflag)
  end select

  ii = 1
  wgr_npes = esdf_integer('distribute_wavefunctions',ii)

  if (esdf_block('tdlda_valence',nlines)) then
     if (pol_in(1)%nval < 0) pol_in(1)%nval = 0
     call esdf_parse_block('tdlda_valence', &
          nlines,pol_in(1)%nval,maxdata,vmap(1,1))
     pol_in(2)%nval = pol_in(1)%nval
     vmap(1:pol_in(1)%nval,2) = vmap(1:pol_in(1)%nval,1)
  endif

  if (esdf_block('tdlda_valence_up',nlines)) then
     if (pol_in(1)%nval < 0) pol_in(1)%nval = 0
     call esdf_parse_block('tdlda_valence_up',nlines, &
          pol_in(1)%nval,maxdata,vmap(1,1))
  endif

  if (esdf_block('tdlda_valence_down',nlines)) then
     if (pol_in(2)%nval < 0) pol_in(2)%nval = 0
     call esdf_parse_block('tdlda_valence_down',nlines, &
          pol_in(2)%nval,maxdata,vmap(1,2))
  endif

  if (esdf_block('tdlda_conduction',nlines)) then
     if (pol_in(1)%ncond < 0) pol_in(1)%ncond = 0
     call esdf_parse_block('tdlda_conduction',nlines, &
          pol_in(1)%ncond,maxdata,cmap(1,1))
     pol_in(2)%ncond = pol_in(1)%ncond
     cmap(1:pol_in(1)%ncond,2) = cmap(1:pol_in(1)%ncond,1)
  endif

  if (esdf_block('tdlda_conduction_up',nlines)) then
     if (pol_in(1)%ncond < 0) pol_in(1)%ncond = 0
     call esdf_parse_block('tdlda_conduction_up',nlines, &
          pol_in(1)%ncond,maxdata,cmap(1,1))
  endif

  if (esdf_block('tdlda_conduction_down',nlines)) then
     if (pol_in(2)%ncond < 0) pol_in(2)%ncond = 0
     call esdf_parse_block('tdlda_conduction_down',nlines, &
          pol_in(2)%ncond,maxdata,cmap(1,2))
  endif

  if (esdf_block('qpoints',qpt%nk))then
     allocate(qpt%fk(3,qpt%nk))
     allocate(qpt%zerok(qpt%nk))
     do ii = 1, qpt%nk
        read(block_data(ii),*) (qpt%fk(jj,ii),jj=1,3), dtmp, jj
        qpt%fk(:,ii) = qpt%fk(:,ii) / dtmp
        if (jj > 0) then
           qpt%zerok(ii) = .true.
        else
           qpt%zerok(ii) = .false.
        endif
     enddo
  else
     qpt%nk = 1
     allocate(qpt%fk(3,qpt%nk))
     qpt%fk = zero
     allocate(qpt%zerok(qpt%nk))
     qpt%zerok = .true.
  endif

  call esdf_close

  !-----------------------------------------------------------------------
  ! Check if indices of valence and conduction states are defined.
  !
  do ii = 1, 2
     if (pol_in(ii)%nval /= 0) then
        allocate(pol_in(ii)%vmap(pol_in(ii)%nval))
        pol_in(ii)%vmap(1:pol_in(ii)%nval) = vmap(1:pol_in(ii)%nval,ii)
     endif
     if (pol_in(ii)%ncond /= 0) then
        allocate(pol_in(ii)%cmap(pol_in(ii)%ncond))
        pol_in(ii)%cmap(1:pol_in(ii)%ncond) = cmap(1:pol_in(ii)%ncond,ii)
     endif
  enddo
  
  ! W Gao dbg 
  if(peinf%master .and. verbose) then
     do ii = 1, 2
        if (pol_in(ii)%nval /= 0) then
           write(6,'("pol_in(",i2,")%nval = ",i8)') ii,pol_in(ii)%nval
           write(6,'("     j      vmap(j) ")')
           do jj = 1, pol_in(ii)%nval
              write(6,'(i8,i8)') jj, pol_in(ii)%vmap(jj)
           enddo
        endif
        write(6,'()')
        if (pol_in(ii)%ncond /= 0) then
           write(6,'("pol_in(",i2,")%ncond = ",i8)') ii,pol_in(ii)%ncond
           write(6,'("     i      cmap(i) ")')
           do jj = 1, pol_in(ii)%ncond
              write(6,'(i8,i8)') jj, pol_in(ii)%cmap(jj)
           enddo
        endif
        write(6,'()')
     enddo
  endif
        
end subroutine input_g
!===================================================================
