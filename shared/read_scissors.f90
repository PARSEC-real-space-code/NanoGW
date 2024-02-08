!===================================================================
!
! Aply a scissors operator with constant and slope. Operator is read
! from file scissors.in. Each line on that file corresponds to one
! operator, with format: jsp, imin, imax, eref, const, slope
!
! jsp   : spin channel (1 or 2)
! imin  : lowest orbital (integer)
! imax  : highest orbital (integer)
! eref  : reference energy, in eV (real)
! const : constant shift, in eV (real)
! slope : slope (real)
!
! Eigenvalues of all orbitals between imin and imax and with spin "jsp"
! are modified according to this rule:
!
!            E_new = E_old + ( E_old - E_ref ) * slope + const
!
! If no scissors operator is found in rgwbs.in, then E_new = E_old.
!
! Copyright (C) 2009 Murilo L. Tiago, http://users.ices.utexas.edu/~mtiago
! This file is part of RGWBS. It is distributed under the GPL v1.
!
!-------------------------------------------------------------------
subroutine read_scissors(verbose,isp,nstate,e_old,e_new)

  use myconstants
  use esdf
  implicit none

  ! arguments
  ! output flag
  logical, intent(in) :: verbose
  ! input spin component, length of eigenvalue arrays
  ! scissors operators in rgwbs.in with spin component different from isp are ignored
  integer, intent(in) :: isp, nstate
  ! input list of eigenvalues
  real(dp), intent(in) :: e_old(nstate)
  ! output list of eigenvalues
  real(dp), intent(out) :: e_new(nstate)

  ! local variables
  integer :: ii, jj, jsp, imin, imax, nscissors, ierr
  real(dp) :: eref, const, slope

  call dcopy(nstate,e_old,1,e_new,1)

  call esdf_init('rgwbs.in')
  if (esdf_block('scissors',nscissors)) then
     if (verbose) write(6,'(/,a,/)') repeat('*',65)
     do jj = 1, nscissors
        read(block_data(jj),*,iostat=ierr) jsp, imin, imax, const, eref, slope
        if (ierr /= 0) call die('ERROR in reading scissors operator')
        if (isp /= jsp) cycle
        if (verbose) then
           write(6,*) 'Warning. Applying scissors operator for spin ',isp
           write(6,*) ' levels ',imin,' to ',imax
           write(6,'(a,g20.10,a)') ' constant shift = ',const,' eV'
           write(6,'(a,g20.10,a)') ' Reference energy ',eref,' eV'
           write(6,'(a,g20.10)') ' slope = ',slope
        endif
        do ii = 1, nstate
           if (ii <= imax .and. ii >= imin) &
              e_new(ii) = e_old(ii) + &
              (e_old(ii) - eref/ryd) * slope + const/ryd
        enddo
     enddo
     if (verbose) write(6,'(/,a,/)') repeat('*',65)
  endif

  call esdf_close

end subroutine read_scissors
!===================================================================
