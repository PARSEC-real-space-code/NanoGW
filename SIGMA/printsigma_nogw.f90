!===================================================================
!
! Write data to disk for Hartree-Fock or DFT calculations.
!
! Copyright (C) 2009 Murilo L. Tiago, http://users.ices.utexas.edu/~mtiago
! This file is part of RGWBS. It is distributed under the GPL v1.
!
!-------------------------------------------------------------------
subroutine printsigma_nogw(nrep,sig,wfn,q_p,iord,qpmix,max_sig,file_name)

  use typedefs
  use mpi_module
  implicit none

  ! arguments
  ! number of irreducible representations
  integer, intent(in) :: nrep
  ! self-energy for the current spin channel and current k-point
  type (siginfo), intent(in) :: sig
  ! electron wavefunctions for the current spin channel and current k-point
  type (wavefunction), intent(in) :: wfn
  ! quasi-particle data for the current k-point, all representations
  type (qpinfo), dimension(nrep), intent(inout) :: q_p
  ! list of electron orbitals ordered according to representation
  integer, intent(in) :: iord(wfn%nstate)
  ! mixing parameter
  real(dp), intent(in) :: qpmix
  ! maximum amplitude of self-energy correction, in eV. The self-energy
  ! correction is defined as Sigma_i - Sigma_(i-1) where Sigma_(i-1) is
  ! the self-energy at previous iteration, or Vxc if i is the first
  ! SCGW iteration.
  real(dp), intent(inout) :: max_sig
  ! name of output file
  character (len=*), intent(in) :: file_name

  ! local variables
  integer, parameter :: jtape = 23
  integer :: ii, jj, i1, i2, j1, j2, isig, irp
  real(dp) :: elda
  complex(dpc) :: esig

  !-------------------------------------------------------------------
  ! Save self-energy, diagonal part into QP structures.
  !
  do isig = 1, sig%ndiag
     ii = sig%map(sig%diag(isig))
     jj = iord(ii)
     irp = wfn%irep(ii)
     q_p(irp)%sigmaqp(jj,jj) = sig%xdiag(isig) + sig%scsdiag(isig)
     q_p(irp)%hqp(jj,jj) = wfn%e0(ii)*ryd + &
          qpmix * (q_p(irp)%sigmaqp(jj,jj) - sig%vxcdiag(isig))
  enddo
  !-------------------------------------------------------------------
  ! Save off-diagonal part into QP structures.
  !
  do isig = 1,sig%noffd
     i1 = sig%map(sig%off1(isig))
     j1 = iord(i1)
     i2 = sig%map(sig%off2(isig))
     j2 = iord(i2)
     irp = wfn%irep(i1)
     q_p(irp)%sigmaqp(j1,j2) = sig%xoffd(isig) + sig%scsoffd(isig)
     q_p(irp)%hqp(j1,j2) = qpmix * (q_p(irp)%sigmaqp(j1,j2) - sig%vxcoffd(isig))
     q_p(irp)%sigmaqp(j2,j1) = conjg( sig%xoffd(isig) + sig%scsoffd(isig) )
     q_p(irp)%hqp(j2,j1) = qpmix * &
          (q_p(irp)%sigmaqp(j2,j1) - conjg( sig%vxcoffd(isig) ))
  enddo
  !-------------------------------------------------------------------
  ! Print self-energy tables.
  !
  if (peinf%master) then
     open(jtape,file=file_name,form='formatted',status='unknown')
     if (sig%ndiag > 0) then
        write(jtape,'(/,2a,/)') ' OUTPUT SELF-ENERGY DIAGONAL ', &
                   'MATRIX ELEMENTS (eV) :'
        write(jtape,'(2a)') '    n  rep     occ       ', &
             'E_0      V_xc   Sigma_x   Sigma_c  Sigma-V_xc    E_qp  '
        do isig = 1, sig%ndiag
           ii = sig%map(sig%diag(isig))
           elda = wfn%e0(ii)*ryd
           esig = sig%xdiag(isig) + sig%scsdiag(isig) - sig%vxcdiag(isig)
           write(jtape,'(2i5,f8.3,6f10.3)') ii,wfn%irep(ii),wfn%occ1(ii), &
                elda,real(sig%vxcdiag(isig),dp),real(sig%xdiag(isig),dp), &
                real(sig%scsdiag(isig),dp), &
                real(esig,dp),elda + real(esig,dp)
           if (max_sig < abs(esig)) max_sig = abs(esig)
        enddo
        write(jtape,'(/,a)') repeat('-',65)
     endif

     if (sig%noffd > 0) then
        write(jtape,'(/,3a,/)') ' OUTPUT SELF-ENERGY OFF-DIAGONAL ', &
             'MATRIX ELEMENTS (eV) :'
        write(jtape,'(2a)') '   n1   n2  rep        V_xc       ', &
             'Sigma_x      Sigma_c   Sigma-V_xc'
        do isig = 1,sig%noffd
           esig = sig%xoffd(isig) + sig%scsoffd(isig) - sig%vxcoffd(isig)
           write(jtape,'(3i5,4f13.3)') sig%map(sig%off1(isig)), &
                sig%map(sig%off2(isig)),wfn%irep(sig%map(sig%off1(isig))), &
                real(sig%vxcoffd(isig),dp),real(sig%xoffd(isig),dp), &
                real(sig%scsoffd(isig),dp),real(esig,dp)
        enddo
        write(jtape,'(/,a)') repeat('-',65)
     endif
     close(jtape)
  endif

end subroutine printsigma_nogw
!===================================================================
