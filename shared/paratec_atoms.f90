!===================================================================
!
! Read in user parameters from input and define structure cluster.
!
! OUTPUT:
!    psp : pseudopotentials for each atom type
!    gvec%avec, bvec, bdot : unit lattice vectors in real space and reciprocal space
!    ierr : error flag (0 for successful exit, positive otherwise)
!
! Copyright (C) 2009 Murilo L. Tiago, http://users.ices.utexas.edu/~mtiago
! This file is part of RGWBS. It is distributed under the GPL v1.
!
!-------------------------------------------------------------------
subroutine paratec_atoms(gvec,verbose,ierr)

  use typedefs
  use esdf
  use psp_module
  implicit none

  ! arguments
  type (gspace), intent(inout):: gvec
  ! true for master processor
  logical, intent(in) :: verbose
  ! error flag
  integer, intent(out) :: ierr

  ! local variables
  character (len=800) :: lastwords
  ! extension that turns atom name into the pseudopotential file name
  character(len=90) :: fname
  ! counters, temporary variables
  logical :: scaleflag
  integer :: ii, jj, ity, iat, pspformat
  real(dp) :: temp0, vtmp(3), avec(3,3)
  ! string flag
  character (len=60) :: strflag, ldata
  ! default values for the different types
  integer, parameter :: integer_default = 0
  real(dp), parameter :: double_default  = zero
  character, parameter :: string_default = ' '

  !-------------------------------------------------------------------

  ! Initialize error flag.
  ierr = 0

  ! User parameter input from file 'input'
  call esdf_init ('input')

  if (esdf_block('latticevecs',jj)) then
     do ii = 1, 3
        read (block_data(ii),*) strflag, gvec%avec(:,ii)
     enddo
     if (jj == 4) then
        read(block_data(jj),*) strflag, vtmp(1)
        avec = gvec%avec
        call mtrxin(avec,temp0,ii)
        if (ii /= 0) call die('ERROR in mtrxin. Stop.')
        temp0 = exp( log(vtmp(1)/abs(temp0)) * third )
        gvec%avec = gvec%avec * temp0
        if ( abs(vtmp(1) - gvec%celvol) > 1.d-6 ) then
           write(6,*) 'ERROR: volume mismatch: ', vtmp(1), gvec%celvol
           ierr = 1
           return
        endif
     endif
  endif
  !
  ! Calculate unit reciprocal lattice vectors.
  !
  gvec%bvec = gvec%avec
  call mtrxin(gvec%bvec,temp0,ii)
  if (ii /= 0) call die('ERROR in mtrxin. Stop.')
  gvec%bvec = two * pi * gvec%bvec
  avec = matmul(gvec%bvec,transpose(gvec%bvec))
  temp0 = maxval(abs(avec - gvec%bdot))
  if ( temp0 > 1.d-6 ) then
     write(lastwords,*) 'ERROR: metric mismatch. ', temp0, avec, gvec%bdot
     ierr = 1
     return
  endif

  if (verbose) write(6,'(a,/,a,/)') ' Atom data :', ' -----------'
  ! Determine number of atoms and number of chemical species.
  type_num = 0
  if (esdf_block('coordinates',jj)) then
    do ii = 1, jj
       read(block_data(ii),*) strflag
       if (trim(strflag) == 'newtype') type_num = type_num + 1
    enddo
    allocate(psp(type_num))
    ity = 0
    psp(:)%natmi = 0
    do ii = 1, jj
       read(block_data(ii),*) strflag
       if (trim(strflag) == 'newtype') then
          ity = ity + 1
       elseif (trim(strflag) == 'coord') then
          psp(ity)%natmi = psp(ity)%natmi + 1
       endif
    enddo
    ity = 0
    do ii = 1, jj
       read(block_data(ii),*) strflag, ldata
       if (trim(strflag) == 'newtype') then
          ity = ity + 1
          iat = 0
          psp(ity)%name = trim(ldata)
          allocate(psp(ity)%ratm(3,psp(ity)%natmi))
       elseif (trim(strflag) == 'coord') then
          iat = iat + 1
          read(block_data(ii),*) ldata, psp(ity)%ratm(:,iat)
       endif
    enddo
  endif

  if (type_num == 0) then
     if (verbose) then
        write(6,*) 'ERROR: unknown number of atom types. '
        write(6,*) ' Input line "coordinates" not found!'
        write(6,*) 'STOP in usrinput'
     endif
     ierr = 1
     return
  endif

  !  Rescale atomic coordinates if needed.
  scaleflag = esdf_defined('coordinates_absolute')
  if (scaleflag) then
     if (verbose) write(6,*) 'Atom coordinates input in Bohr radii.'
  else
     if (verbose) write(6,*) &
          'Atom coordinates input in lattice vector units.'
     do ity = 1, type_num
        do iat = 1, psp(ity)%natmi
           call dmatvec3('N',gvec%avec,psp(ity)%ratm(1,iat),vtmp)
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

  pspformat = esdf_integer ('pp_format',MARTINS_KB)
  if (pspformat == MARTINS_KB) then
     if (verbose) then
        write(6,*) 'ERROR: Original PSP format in reciprocal space'
        write(6,*) 'is not supported.'
        write(6,*) 'Supported formats are L-W Wang, FHIPP only.'
        write(6,*) 'STOP in paratec_atoms.'
     endif
     ierr = 1
     return
  endif

  if (pspformat == FHIPP) then
     if (esdf_block('pseudopotential',jj)) then
        do ity = 1, type_num
           read (block_data(ity),*) strflag, psp(ity)%loc
        enddo
     endif
  endif

  do ity = 1, type_num

     if (verbose) &
          write(6,'(/,a,a,a,/,a,/)') ' Pseudopotential messages, element ', &
          trim(psp(ity)%name),' :', ' ', repeat('-',38)

     fname = trim(psp(ity)%name)//'_POT.DAT'
     select case(pspformat)
     case (MARTINS_WANG)
        if (verbose) write(6,*) 'Pseudopotential format : Martins-Wang'
     case (FHIPP)
        if (verbose) write(6,*) 'Pseudopotential format : FHIPP'
     case default
        if (verbose) write(6,'(a,i5,/,a)') &
             ' ERROR. This pseudopotential format is not recognized: ', &
             pspformat, ' STOP.'
        ierr = 1
        return
     end select

     psp(ity)%rcore = zero

     call read_psp(psp(ity),verbose,fname,pspformat,gvec%hcub,ierr)
     if (ierr /= 0) return
  enddo

  call esdf_close

end subroutine paratec_atoms
!===================================================================
