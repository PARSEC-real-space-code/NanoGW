!===================================================================
!
! Read in user parameters from parsec.in and define cluster structure.
!
! OUTPUT:
!    psp : pseudopotentials for each atom type
!    norder : order of finite differences
!    ierr : error flag (0 for successful exit, positive otherwise)
!
! Copyright (C) 2009 Murilo L. Tiago, http://users.ices.utexas.edu/~mtiago
! This file is part of RGWBS. It is distributed under the GPL v1.
!
!-------------------------------------------------------------------
subroutine parsec_atoms(gvec,verbose,norder,ierr)

  use typedefs
  use esdf
  use psp_module
  implicit none

  ! arguments
  type (gspace), intent(inout):: gvec
  ! true for master processor
  logical, intent(in) :: verbose
  ! finite difference order
  integer, intent(out) :: norder
  ! error flag
  integer, intent(out) :: ierr

  ! local variables
  ! default values for the different types
  integer, parameter :: integer_default = 0
  real(dp), parameter :: double_default  = zero
  character, parameter :: string_default = ' '

  ! scale used for atomic coordinates
  real(dp) :: cscale(3,3)

  ! string flags
  character (len=60) :: strflag
  character (len=1) :: lflag

  ! counters, temporary variables
  integer :: ii, jj, ity, itmp, iat
  real(dp) :: temp0, vtmp(3)

  ! pspformat = MARTINS_NEW  --> read *_POTRE.DAT files (with header info)
  ! pspformat = MARTINS_WANG --> read *_POTRW.DAT files (with header info)
  ! pspformat = FHIPP        --> read *_FHIPP.DAT files (with header info)
  integer :: pspformat

  ! extension that turns atom name into the pseudopotential file name
  character(len=90) :: fname

  !-------------------------------------------------------------------

  ! Initialize error flag.
  ierr = 0

  ! User parameter input from file 'parsec_input'	
  ! Progress output written into 'parsec_output'
  call esdf_init ('parsec.in')

  ! Read in order of finite difference expansion.
  itmp = 12
  norder = esdf_integer ('Expansion_Order',itmp)

  norder = norder / 2

  if (verbose) write(6,'(a,/,a,/)') ' Atom data :', ' -----------'
  ! Read in and write number of atom types.
  type_num = esdf_integer('Atom_Types_Num',integer_default)
  if (type_num == integer_default) then
     if (verbose) then
        write(6,*) 'ERROR: unknown number of atom types. '
        write(6,*) ' Input line "Atom_Types_Num" not found!'
        write(6,*) 'STOP in read_cluster'
     endif
     ierr = 1
     return
  endif

  strflag = esdf_reduce(esdf_string ('Coordinate_Unit','cartesian_bohr'))
  ! cscale : for periodic systems, hold unit lattice vectors in column-wise form
  cscale = zero

  select case (trim(strflag))
  case ('cartesian_bohr')
     if (verbose) write(6,*) 'Atom coordinates input in Bohr radii.'
     do ii = 1, 3
        cscale(ii,ii) = one
     enddo
  case ('cartesian_ang')
     if (verbose) write(6,*) 'Atom coordinates input in angstroms.'
     do ii = 1, 3
        cscale(ii,ii) = one/angs
     enddo
  case ('lattice_vectors')
     if (verbose) write(6,*) 'Atom coordinates input in lattice vector units.'
     cscale = -one
  end select
  temp0 = esdf_physical('Coordinate_Scale',one,'bohr')
  if (temp0 /= one) then
     if (verbose) then
        write(6,*) 'Rescaling atomic coordinates by ',temp0,' a.u.'
        if (trim(strflag) /= 'cartesian_bohr') write(6,*) &
             'WARNING: coordinate_unit input is superseded by coordinate_scale!'
     endif
     cscale(:,:) = zero
     do ii = 1, 3
        cscale(ii,ii) = temp0
     enddo
  endif

  allocate(psp(type_num))

  ! Read name and atom coordinates.
  do ity = 1, type_num
     psp(ity)%name = esdf_string ('Atom_Type', string_default)
     if (esdf_block ('Atom_Coord',psp(ity)%natmi )) then   
        allocate (psp(ity)%ratm (3,psp(ity)%natmi))
        do iat = 1, psp(ity)%natmi
           read(block_data(iat),*) (psp(ity)%ratm(jj,iat),jj=1,3)
        enddo
     else
        if (verbose) then
           write(6,*) 'ERROR: no atom coordinates were found'
           write(6,*) 'for atom type ', psp(ity)%name
           write(6,*) 'STOP in read_cluster.'
        endif
        ierr = 1
        return
     endif
  enddo

  ! Rescale atomic coordinates if needed.
  if (maxval(abs(cscale)) > zero) then
     if (cscale(1,1) < zero) cscale = gvec%avec
     do ity = 1, type_num
        do iat = 1, psp(ity)%natmi
           call dmatvec3('N',cscale,psp(ity)%ratm(1,iat),vtmp)
           psp(ity)%ratm(:,iat) = vtmp
        enddo
     enddo
  endif

  if (verbose) then
     write(6,'(/,a,/,1x,3(10x,a),10x,a)') &
          '  Atomic coordinates (all atoms) are :',  &
          'x [a.u.]', 'y [a.u.]', 'z [a.u.]', ' Type'
     do ity = 1, type_num
        do iat = 1, psp(ity)%natmi
           write(6,'(i4,3(1x,f17.9),9x,a)') iat, &
                psp(ity)%ratm(:,iat), psp(ity)%name
        enddo
     enddo
     write(6,'(/,a,i4)') ' Total number of atoms = ', sum(psp(:)%natmi)
     write(6,'(a,i4)') ' Total number of atom types = ', type_num
  endif

  do ity = 1, type_num
     if (verbose) &
          write(6,'(/,a,a,a,/,a,/)') ' Pseudopotential messages, element ', &
          trim(psp(ity)%name),' :', ' --------------------------------------'
     strflag = esdf_reduce(esdf_string ('Pseudopotential_Format','martins_new'))
     select case(trim(strflag))
     case ('martins_new')
        pspformat = MARTINS_NEW
        fname = trim(psp(ity)%name)//'_POTRE.DAT'
        if (verbose) write(6,*) 'Pseudopotential format : new Martins'
     case ('martins_wang')
        pspformat = MARTINS_WANG
        fname = trim(psp(ity)%name)//'_POTRW.DAT'
        if (verbose) write(6,*) 'Pseudopotential format : Martins-Wang'
     case ('fhipp')
        pspformat = FHIPP
        fname = trim(psp(ity)%name)//'_FHIPP.DAT'
        if (verbose) write(6,*) 'Pseudopotential format : FHIPP'
     case default
        if (verbose) write(6,'(a,a,/,a)') &
             ' ERROR. This pseudopotential format is not recognized: ', &
             trim(strflag), ' STOP.'
        ierr = 1
        return
     end select
     !
     ! Read other pseudopotential parameters, if available
     ! WARNING: if pseudopotential files have header info, the number
     ! of angular momentum channels must be consistent with what is
     ! inside the psp files.
     !
     if ( pspformat /= MARTINS_WANG ) then
        lflag = esdf_reduce(esdf_string('Local_Component', string_default))
        select case (trim(lflag))
        case ('s')
           psp(ity)%loc = 1
        case ('p')
           psp(ity)%loc = 2
        case ('d')
           psp(ity)%loc = 3
        case ('f')
           psp(ity)%loc = 4
        case default
           if (verbose) then
              write(6,'(/,a)') ' ERROR: problem with this atom type'
              write(6,*) 'Impossible local pseudopotential : ', lflag
              write(6,*) 'Must be s, p, d, or f'
              write(6,*) 'STOP in read_cluster'
           endif
           ierr = 1
           return
        end select
     endif

     if (pspformat /= MARTINS_NEW ) psp(ity)%rcore = &
          esdf_physical ('Core_Cutoff_Radius',double_default,'bohr')
     !
     ! Read pseudopotentials.
     !
     call read_psp(psp(ity),verbose,fname,pspformat,gvec%hcub,ierr)
     if (ierr /= 0) return
  enddo

  call esdf_close

end subroutine parsec_atoms
!===================================================================
