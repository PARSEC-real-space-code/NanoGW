!===================================================================
!
! Write checkpoint file sigma.chkpt.dat.
! All output is written to disk. All parameters are returned
! with input values.
!
! Copyright (C) 2009 Murilo L. Tiago, http://users.ices.utexas.edu/~mtiago
! This file is part of RGWBS. It is distributed under the GPL v1.
!
!-------------------------------------------------------------------
subroutine write_sigma(nkpt,nspin,kpt,sig_in,sig,chkpt)

  use typedefs
  implicit none

  ! arguments
  integer, intent(in) :: nkpt, nspin, chkpt
  type (kptinfo), intent(in) :: kpt
  type (siginfo), intent(in) :: sig_in
  type (siginfo), dimension(nspin,nkpt), intent(in) :: sig

  ! local variables
  integer :: isp, ik, ii, jk
  ! number of output unit
  integer, parameter :: iunit = 60

  open(iunit,file='sigma.chkpt.dat',form='unformatted')
  write(iunit) chkpt, nspin, nkpt, sig_in%xc

  do ik = 1, nkpt
     do isp = 1, nspin
        jk = sig(isp,ik)%indxk
        write(iunit) sig(isp,ik)%nmap, sig(isp,ik)%indxk
        write(iunit) sig(isp,ik)%ndiag_s, sig(isp,ik)%ndiag
        write(iunit) sig(isp,ik)%noffd_s, sig(isp,ik)%noffd
        write(iunit) sig_in%nen, sig_in%deltae
        if (sig(isp,ik)%nmap > 0) &
             write(iunit) (sig(isp,ik)%map(ii),ii=1,sig(isp,ik)%nmap)
        if (sig(isp,ik)%ndiag_s > 0) then
           write(iunit) (sig(isp,ik)%diag(ii),ii=1,sig(isp,ik)%ndiag_s)
           write(iunit) (kpt%wfn(isp,jk)%e0(sig(isp,ik)%map( &
                sig(isp,ik)%diag(ii))),ii=1,sig(isp,ik)%ndiag_s)
           write(iunit) (kpt%wfn(isp,jk)%e1(sig(isp,ik)%map( &
                sig(isp,ik)%diag(ii))),ii=1,sig(isp,ik)%ndiag_s)
           write(iunit) (sig(isp,ik)%vxcdiag(ii),ii=1,sig(isp,ik)%ndiag_s)
           write(iunit) (sig(isp,ik)%xdiag(ii),ii=1,sig(isp,ik)%ndiag_s)
           write(iunit) (sig(isp,ik)%scsdiag(ii),ii=1,sig(isp,ik)%ndiag_s)
           if (sig_in%xc == XC_GW) then
              write(iunit) (sig(isp,ik)%sgsdiag(ii),ii=1,sig(isp,ik)%ndiag_s)
           endif
        endif
        if (sig_in%xc == XC_GW .and. sig(isp,ik)%ndiag > 0) then
           write(iunit) (sig(isp,ik)%scdiag(:,ii),ii=1,2*sig(isp,ik)%ndiag)
           write(iunit) (sig(isp,ik)%sexdiag(:,ii),ii=1,2*sig(isp,ik)%ndiag)
           write(iunit) (sig(isp,ik)%sgdiag(:,ii),ii=1,2*sig(isp,ik)%ndiag)
        endif
        if (sig(isp,ik)%noffd_s > 0) then
           write(iunit) (sig(isp,ik)%off1(ii),ii=1,sig(isp,ik)%noffd_s)
           write(iunit) (sig(isp,ik)%off2(ii),ii=1,sig(isp,ik)%noffd_s)
           write(iunit) (sig(isp,ik)%vxcoffd(ii),ii=1,sig(isp,ik)%noffd_s)
           write(iunit) (sig(isp,ik)%xoffd(ii),ii=1,sig(isp,ik)%noffd_s)
           write(iunit) (sig(isp,ik)%scsoffd(ii),ii=1,sig(isp,ik)%noffd_s)
           if (sig_in%xc == XC_GW) then
              write(iunit) (sig(isp,ik)%sgsoffd(ii),ii=1,sig(isp,ik)%noffd_s)
           endif
        endif
        if (sig_in%xc == XC_GW .and. sig(isp,ik)%noffd > 0) then
           write(iunit)  (sig(isp,ik)%scoffd(:,ii),ii=1,2*sig(isp,ik)%noffd)
           write(iunit) (sig(isp,ik)%sgoffd(:,ii),ii=1,2*sig(isp,ik)%noffd)
        endif
     enddo
  enddo

  close(iunit)
  write(6,'(/,a,i5,/)') ' Checkpoint control = ',chkpt

end subroutine write_sigma
!===================================================================
