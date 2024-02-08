!===================================================================
!
! Write data to disk and save self-energy data on QP structure.
! At output, these are the quantities modified in this subroutine:
!  q_p(:)%sigmai : imaginary part of self-energy at QP energy
!  q_p(:,:)%sigmaqp : self-energy correction in the basis of input orbitals
!  q_p(:,:)%hqp : QP Hamiltonian in the basis of input orbitals.
!
!  H_qp (i,j) = E_0 * delta(i,j) + Sigma(i,j) - V_xc(i,j)
!
! Notice that H_qp and Sigma are not necessarily Hermitian matrices (they
! are symmetrized later).If Sigma_qp(i,j) is not available because
! orbitals i and/or j are not included in the problem, then
! Sigma_qp(i,j) and H_qp(i,j) are left unchanged.
! Also, H_qp and Sigma are stored in units of eV.
!
! Copyright (C) 2009 Murilo L. Tiago, http://users.ices.utexas.edu/~mtiago
! This file is part of RGWBS. It is distributed under the GPL v1.
!
!-------------------------------------------------------------------
subroutine printsigma(nrep,sig_in,sig,wfn,q_p,itape,isp,ik, &
     iord,qpmix,max_sig,file_name)

  use typedefs
  use mpi_module
  implicit none

  ! arguments
  integer, intent(in) :: &
       nrep, &       ! number of irreducible representations
       itape, &      ! output unit
       isp, &        ! current spin channel
       ik            ! current k-point
  ! self-energy for the current spin channel and current k-point
  type (siginfo), intent(in) :: sig_in, sig
  ! electron wavefunctions for the current spin channel and current k-point
  type (wavefunction), intent(in) :: wfn
  ! quasi-particle data for the current k-point, all representations
  type (qpinfo), dimension(nrep), intent(inout) :: q_p
  ! list of electron orbitals ordered according to representation
  integer, intent(in) :: iord(sig%nmap)
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
  integer :: ii, jj, i1, i2, j1, j2, isig, irp, ierr, jtape
  real(dp) ::  edft, en0, xtmp
  complex(dpc) :: sigmaoff(2), ztmp
  real(dp), dimension(:), allocatable :: zfac, eqp, &
       xtable
  complex(dpc), dimension(:), allocatable :: sigma0, ftable, &
       sexdiag, scdiag, sgdiag, sigmai, sigmai_table
  real(dp), parameter :: EN_ACC = 0.01

  !-------------------------------------------------------------------
  ! Allocate and initialize buffers.
  !
  if (sig%ndiag_s > 0) then
     allocate(sigma0(sig%ndiag_s))
     sigma0 = zzero
     allocate(eqp(sig%ndiag_s))
     eqp = zero
  endif
  if (sig%ndiag > 0) then
     allocate(zfac(sig%ndiag))
     zfac = zero
     allocate(sexdiag(sig%ndiag))
     sexdiag = zzero
     allocate(scdiag(sig%ndiag))
     scdiag = zzero
     allocate(sgdiag(sig%ndiag))
     sgdiag = zzero
     allocate(sigmai(sig%ndiag))
     sigmai = zzero
  endif
  allocate(ftable(sig_in%nen))
  allocate(xtable(sig_in%nen))
  allocate(sigmai_table(sig_in%nen))
  !-------------------------------------------------------------------
  ! Determine the quasi-particle energy E_qp by solving the equation:
  !
  !   E_qp = E_0 + Sigma(E_qp) - Sigma_old         (*)
  !
  ! where Sigma_old = Vxc if this is the first SCGW iteration, or
  ! the self-energy in previous iteration otherwise.
  ! Equation (*) is solved using a Newton-Raphson algorithm and
  ! spline interpolation of Sigma. The set of quasi-particle energies
  ! to be determined is distributed over PEs. At the end, PEs share data.
  !
  do isig = 1 + peinf%inode, sig%ndiag, peinf%npes
     edft =  wfn%e0(sig%map(sig%diag(isig)))*ryd
     en0 =  wfn%e1(sig%map(sig%diag(isig)))*ryd
     do ii = 1, sig_in%nen
        xtable(ii) = en0 - sig_in%deltae*half*ryd + &
             sig_in%deltae*ryd / real(sig_in%nen - 1,dp) * real(ii-1,dp)
        ftable(ii) = sig%xdiag(isig) + sig%scdiag(ii,isig) + &
             sig%sgdiag(ii,isig) - sig%vxcdiag(isig) + edft - xtable(ii)
        sigmai_table(ii) = sig%scdiag(ii,isig + sig%ndiag) + &
             sig%sgdiag(ii,isig + sig%ndiag)
     enddo
     eqp(isig) = en0
     call newtonr(sig_in%nen,xtable,ftable,eqp(isig),EN_ACC,ztmp,ierr)
     if (ierr /= 0) then
        write(6,*) 'ERROR!! Newton-Raphson exceeded maximum ', &
             'number of iterations (',ierr,') at orbital ', &
             sig%map(sig%diag(isig))
        write(6,*) ' Value of ( Sigma[E_qp] - Vxc - E_qp + E_dft ) ', &
             'is ',ztmp,' eV'
        write(6,*) 'Quasiparticle energy may be incorrect!.'
     endif

     call spline_evaluate(sig_in%nen,xtable,ftable,eqp(isig),sigma0(isig),zfac(isig))
     call spline_evaluate(sig_in%nen,xtable,ftable,edft,sigma0(isig),xtmp)
     sigma0(isig) = sigma0(isig) + sig%vxcdiag(isig)
     zfac(isig) = -one/zfac(isig)
     call spline_evaluate(sig_in%nen,xtable,sig%sexdiag(1,isig),edft,sexdiag(isig),xtmp)
     call spline_evaluate(sig_in%nen,xtable,sig%scdiag(1,isig),edft,scdiag(isig),xtmp)
     call spline_evaluate(sig_in%nen,xtable,sig%sgdiag(1,isig),edft,sgdiag(isig),xtmp)
     call spline_evaluate(sig_in%nen,xtable,sigmai_table,eqp(isig),sigmai(isig),xtmp)

  enddo
  do isig = sig%ndiag + 1 + peinf%inode, sig%ndiag_s, peinf%npes
     sigma0(isig) = sig%xdiag(isig) + sig%scsdiag(isig) + sig%sgsdiag(isig)
     eqp(isig) = wfn%e0(sig%map(sig%diag(isig)))*ryd + &
          real(sigma0(isig),dp) - real(sig%vxcdiag(isig),dp)
  enddo
  if (sig%ndiag_s > 0) then
     call zpsum(sig%ndiag_s,peinf%npes,peinf%comm,sigma0)
     call dpsum(sig%ndiag_s,peinf%npes,peinf%comm,eqp)
  endif
  if (sig%ndiag > 0) then
     call dpsum(sig%ndiag,peinf%npes,peinf%comm,zfac)
     call zpsum(sig%ndiag,peinf%npes,peinf%comm,sexdiag)
     call zpsum(sig%ndiag,peinf%npes,peinf%comm,scdiag)
     call zpsum(sig%ndiag,peinf%npes,peinf%comm,sgdiag)
     call zpsum(sig%ndiag,peinf%npes,peinf%comm,sigmai)
  endif

  !-------------------------------------------------------------------
  ! Save diagonal part into QP structures.
  !
  do isig = 1,sig%ndiag_s
     ii = sig%map(sig%diag(isig))
     jj = iord(sig%diag(isig))
     irp = wfn%irep(ii)
     q_p(irp)%hqp(jj,jj) = wfn%e1(ii)*ryd + qpmix * (eqp(isig) - wfn%e1(ii)*ryd)
     q_p(irp)%sigmaqp(jj,jj) = eqp(isig) - wfn%e0(ii)*ryd + sig%vxcdiag(isig)
  enddo
  do isig = 1,sig%ndiag
     jj = iord(sig%diag(isig))
     irp = wfn%irep(sig%map(sig%diag(isig)))
     q_p(irp)%sigmai(jj) = real(sigmai(isig),dp)
  enddo

  !-------------------------------------------------------------------
  ! Master PE prints data.
  !
  if (peinf%master) then
     jtape = itape + 1
     open(jtape,file=file_name,form='formatted',status='unknown')
     if (sig%ndiag > 0) then
        write(jtape,'(a)') repeat('-',65)
        write(jtape,'(/,2a)') ' OUTPUT SELF-ENERGY DIAGONAL MATRIX ', &
             'ELEMENTS (eV)'
        write(jtape,'(3a)') '     n   occ        E_dft     V_xc   ', &
             'Sigma_x  Sigma_sx   Sigma_c   Sigma_g  Sigma-V_xc   ', &
             'zfac  Im[Sigma]      E_0      E_qp  '
        do isig = 1,sig%ndiag
           en0 =  wfn%e1(sig%map(sig%diag(isig)))*ryd
           edft =  wfn%e0(sig%map(sig%diag(isig)))*ryd
           write(jtape,'(i6,f8.3,11(1x,f9.3))') sig%map(sig%diag(isig)), &
                wfn%occ1(sig%map(sig%diag(isig))),edft, &
                real(sig%vxcdiag(isig),dp), real(sig%xdiag(isig),dp),&
                real(sexdiag(isig),dp),real(scdiag(isig),dp), &
                real(sgdiag(isig),dp),real(sigma0(isig)-sig%vxcdiag(isig),dp), &
                zfac(isig),real(sigmai(isig),dp),en0,eqp(isig)
           write(itape,'(2i8,f10.3,2i4,2f10.3)') sig%map(sig%diag(isig)), &
                sig%map(sig%diag(isig)),eqp(isig),isp,ik,edft,en0
           en0 = en0 - eqp(isig)
           if (max_sig < abs(en0)) max_sig = abs(en0)
        enddo
     endif

     if (sig%ndiag_s > sig%ndiag) then
        write(jtape,'(/,a)') repeat('-',65)
        write(jtape,'(/,2a)') ' OUTPUT SELF-ENERGY DIAGONAL MATRIX ', &
             'ELEMENTS (eV), COHSEX APPROXIMATION'
        write(jtape,'(3a)') '     n   occ        E_dft      V_xc   ', &
             'Sigma_x  Sigma_c   Sigma_g  Sigma-V_xc   ', &
             'E_0      E_qp  '
        do isig = sig%ndiag + 1,sig%ndiag_s
           en0 =  wfn%e1(sig%map(sig%diag(isig)))*ryd
           edft =  wfn%e0(sig%map(sig%diag(isig)))*ryd
           write(jtape,'(i6,f8.3,10f10.3)') sig%map(sig%diag(isig)), &
                wfn%occ1(sig%map(sig%diag(isig))),edft, &
                real(sig%vxcdiag(isig),dp),real(sig%xdiag(isig),dp), &
                real(sig%scsdiag(isig),dp),real(sig%sgsdiag(isig),dp), &
                real(sigma0(isig)-sig%vxcdiag(isig),dp),en0,eqp(isig)
           write(itape,'(2i8,f10.3,2i4,2f10.3)') sig%map(sig%diag(isig)), &
                sig%map(sig%diag(isig)),eqp(isig),isp,ik,edft,en0
           en0 = en0 - eqp(isig)
           if (max_sig < abs(en0)) max_sig = abs(en0)
        enddo
     endif

#ifdef DEBUG
     write(jtape,'(/,a)') repeat('-',65)
     write(jtape,'(/,2a)') ' OUTPUT SELF-ENERGY DIAGONAL MATRIX ', &
          'ELEMENTS (eV)'
     write(jtape,'(2a)') '    n     E_dft   Sigma_x      zfac  ', &
          'Sigma_sx     Sigma      V_xc       E_0      E_qp  '
     do isig = 1, sig%ndiag
        en0 =  wfn%e1(sig%map(sig%diag(isig)))*ryd
        edft =  wfn%e0(sig%map(sig%diag(isig)))*ryd
        write(jtape,'(i6,8f10.3)') sig%map(sig%diag(isig)), &
             edft,real(sig%xdiag(isig),dp),zfac(isig),real(sexdiag(isig),dp), &
             real(sigma0(isig),dp),real(sig%vxcdiag(isig),dp),en0,eqp(isig)
     enddo
#endif
  endif

  deallocate(sigmai_table)
  deallocate(xtable)
  deallocate(ftable)
  if (sig%ndiag_s > 0) then
     deallocate(sigma0)
     deallocate(eqp)
  endif
  if (sig%ndiag > 0) then
     deallocate(zfac)
     deallocate(sexdiag)
     deallocate(scdiag)
     deallocate(sgdiag)
     deallocate(sigmai)
  endif

  !-------------------------------------------------------------------
  ! Work now with off-diagonal part of self-energy.
  !
  if (sig%noffd_s > 0) then
     ! Save off-diagonal part into QP structures.
     do isig = 1, sig%noffd
        i1 = sig%map(sig%off1(isig))
        j1 = iord(sig%off1(isig))
        i2 = sig%map(sig%off2(isig))
        j2 = iord(sig%off2(isig))
        irp = wfn%irep(i1)
        q_p(irp)%sigmaqp(j1,j2) = sig%xoffd(isig) + &
             sig%scoffd(1,isig) + sig%sgoffd(1,isig)
        q_p(irp)%hqp(j1,j2) = qpmix * ( q_p(irp)%sigmaqp(j1,j2) - &
             sig%vxcoffd(isig) )
        q_p(irp)%sigmaqp(j2,j1) = conjg( sig%xoffd(isig) ) + &
             sig%scoffd(2,isig) + sig%sgoffd(2,isig)
        q_p(irp)%hqp(j2,j1) = qpmix * ( q_p(irp)%sigmaqp(j2,j1) - &
             conjg( sig%vxcoffd(isig) ) )
     enddo

     do isig = sig%noffd + 1, sig%noffd_s
        i1 = sig%map(sig%off1(isig))
        j1 = iord(sig%off1(isig))
        i2 = sig%map(sig%off2(isig))
        j2 = iord(sig%off2(isig))
        irp = wfn%irep(i1)
        q_p(irp)%sigmaqp(j1,j2) = sig%xoffd(isig) + &
             sig%scsoffd(isig) + sig%sgsoffd(isig)
        q_p(irp)%hqp(j1,j2) = qpmix * ( q_p(irp)%sigmaqp(j1,j2) - &
             sig%vxcoffd(isig) )
        q_p(irp)%sigmaqp(j2,j1) = conjg( sig%xoffd(isig) + &
             sig%scsoffd(isig) + sig%sgsoffd(isig) )
        q_p(irp)%hqp(j2,j1) = qpmix * ( q_p(irp)%sigmaqp(j2,j1) - &
             conjg( sig%vxcoffd(isig) ) )
     enddo

     ! Print data to screen.
     if (peinf%master) then
        if (sig%ndiag > 0) then
           write(jtape,'(a,/)') repeat('-',65)
           write(jtape,'(/,2a)') ' OUTPUT SELF-ENERGY OFF-DIAGONAL ', &
                'MATRIX ELEMENTS (eV)'
           write(jtape,'(2a)') '    n1    n2     V_xc    Sigma_x  Sigma_c1', &
                '  Sigma_c2  Sigma_g1  Sigma_g2  Sigma-V_xc1  Sigma-V_xc2'
           do isig = 1, sig%noffd
              sigmaoff(1) = sig%xoffd(isig) + sig%scoffd(1,isig) + &
                   sig%sgoffd(1,isig) - sig%vxcoffd(isig)
              sigmaoff(2) = conjg( sig%xoffd(isig) ) + sig%scoffd(2,isig) + &
                   sig%sgoffd(2,isig) - conjg( sig%vxcoffd(isig) )
              write(jtape,'(2i6,6f10.3,2f13.3)') sig%map(sig%off1(isig)), &
                   sig%map(sig%off2(isig)),real(sig%vxcoffd(isig),dp), &
                   real(sig%xoffd(isig),dp),real(sig%scoffd(:,isig),dp), &
                   real(sig%sgoffd(:,isig),dp),real(sigmaoff(:),dp)
              write(itape,'(2i8,f10.3,2i4,f10.3)') sig%map(sig%off1(isig)), &
                   sig%map(sig%off2(isig)),real(sigmaoff(1),dp),isp,ik,zero
              write(itape,'(2i8,f10.3,2i4,f10.3)') sig%map(sig%off2(isig)), &
                   sig%map(sig%off1(isig)),real(sigmaoff(2),dp),isp,ik,zero
           enddo
        endif

        if (sig%noffd_s > sig%noffd) then
           write(jtape,'(/,a)') repeat('-',65)
           write(jtape,'(/,2a)') ' OUTPUT SELF-ENERGY OFF-DIAGONAL ', &
                'MATRIX ELEMENTS (eV), COHSEX APPROXIMATION'
           write(jtape,'(2a)') '    n1    n2     V_xc    Sigma_x  ', &
                'Sigma_c   Sigma_g   Sigma-V_xc1  Sigma-V_xc2'
           do isig = sig%noffd + 1, sig%noffd_s
              sigmaoff(1) = sig%xoffd(isig) + sig%scsoffd(isig) + &
                   sig%sgsoffd(isig) - sig%vxcoffd(isig)
              sigmaoff(2) = conjg( sig%xoffd(isig) + sig%scsoffd(isig) + &
                   sig%sgsoffd(isig) - sig%vxcoffd(isig) )
              write(jtape,'(2i6,4f10.3,2f13.3)') sig%map(sig%off1(isig)), &
                   sig%map(sig%off2(isig)),real(sig%vxcoffd(isig),dp), &
                   real(sig%xoffd(isig),dp),real(sig%scsoffd(isig),dp), &
                   real(sig%sgsoffd(isig),dp),real(sigmaoff(:),dp)
              write(itape,'(2i8,f10.3,2i4,f10.3)') sig%map(sig%off1(isig)), &
                   sig%map(sig%off2(isig)),real(sigmaoff(1),dp),isp,ik,zero
              write(itape,'(2i8,f10.3,2i4,f10.3)') sig%map(sig%off2(isig)), &
                   sig%map(sig%off1(isig)),real(sigmaoff(1),dp),isp,ik,zero
           enddo
        endif

#ifdef DEBUG
        write(jtape,'(/,a)') repeat('-',65)
        write(jtape,'(/,3a)') ' OUTPUT SELF-ENERGY OFF-DIAGONAL MATRIX ', &
             'ELEMENTS (eV)'
        write(jtape,'(2a)') '    n1    n2   Sigma_x            ', &
             'Sigma_c     Sigma          V_xc Sigma-V_xc'
        do isig = 1,sig%noffd
           sigmaoff(1) = sig%xoffd(isig) + sig%scoffd(1,isig) + &
                sig%sgoffd(1,isig) - sig%vxcoffd(isig)
           sigmaoff(2) = conjg( sig%xoffd(isig) ) + sig%scoffd(2,isig) + &
                sig%sgoffd(2,isig) - conjg( sig%vxcoffd(isig) )
           write(jtape,'(2i6,4f10.3,2f13.3)') sig%map(sig%off1(isig)), &
                sig%map(sig%off2(isig)),real(sig%xoffd(isig),dp),zero, &
                real(sig%scoffd(1,isig) + sig%sgoffd(1,isig),dp), &
                real(sigmaoff(1) + sig%vxcoffd(isig),dp), &
                real(sig%vxcoffd(isig),dp),real(sigmaoff(1),dp)
           write(jtape,'(2i6,4f10.3,2f13.3)') sig%map(sig%off2(isig)), &
                sig%map(sig%off1(isig)),real(sig%xoffd(isig),dp),zero, &
                real(sig%scoffd(2,isig) + sig%sgoffd(2,isig),dp), &
                real(sigmaoff(2) + sig%vxcoffd(isig),dp), &
                real(sig%vxcoffd(isig),dp),real(sigmaoff(2),dp)
        enddo
#endif
        write(jtape,'(/,a)') repeat('-',65)
     endif
  endif
  if (peinf%master) close(jtape)

end subroutine printsigma
!===================================================================
