!===================================================================
!
! Write matrix terms of the self-energy in standard output (unit 6).
! Useful during debug only.
!
! Copyright (C) 2009 Murilo L. Tiago, http://users.ices.utexas.edu/~mtiago
! This file is part of RGWBS. It is distributed under the GPL v1.
!
!-------------------------------------------------------------------
subroutine save_sigma(nkpt,nspin,sig_in,sig)

  use typedefs
  implicit none

  type (siginfo), intent(in) :: sig_in
  type (siginfo),dimension(nspin,nkpt), intent(in) :: sig
  ! number of k-points and of spin components
  integer, intent(in) :: nkpt, nspin

  ! index of current spin component and of current k-point
  integer :: isp, ik

  do ik = 1, nkpt
     do isp = 1, nspin
#ifdef DEBUG
        if (sig(isp,ik)%ndiag_s > 0) then
           write(6,10) 'Exchange', isp, ik
           write(6,'(5g14.6)') real(sig(isp,ik)%xdiag,dp)*ryd, &
                aimag(sig(isp,ik)%xdiag)*ryd
           write(6,10) 'Static Correlation', isp, ik
           write(6,'(5g14.6)') real(sig(isp,ik)%scsdiag,dp)*ryd, &
                aimag(sig(isp,ik)%scsdiag)*ryd
           if (sig_in%xc == XC_GW) then
              write(6,10) 'Static Vertex', isp, ik
              write(6,'(5g14.6)') real(sig(isp,ik)%sgsdiag,dp)*ryd, &
                aimag(sig(isp,ik)%sgsdiag)*ryd
           endif
        endif
        if (sig(isp,ik)%noffd_s > 0) then
           write(6,10) 'Exchange, off-diagonal', isp, ik
           write(6,'(5g14.6)') real(sig(isp,ik)%xoffd,dp)*ryd, &
                aimag(sig(isp,ik)%xoffd)*ryd
           write(6,10) 'Static Correlation, off-diagonal', isp, ik
           write(6,'(5g14.6)') real(sig(isp,ik)%scsoffd,dp)*ryd, &
                aimag(sig(isp,ik)%scsoffd)*ryd
           if (sig_in%xc == XC_GW) then
              write(6,10) 'Static Vertex, off-diagonal', isp, ik
              write(6,'(5g14.6)') real(sig(isp,ik)%sgsoffd,dp)*ryd, &
                aimag(sig(isp,ik)%sgsoffd)*ryd
           endif
        endif
        if (sig_in%xc == XC_GW) then
           if (sig(isp,ik)%ndiag > 0) then
              write(6,10) 'Screened Exchange', isp, ik
              write(6,'(3g14.6)') real(sig(isp,ik)%sexdiag,dp)*ryd, &
                aimag(sig(isp,ik)%sexdiag)*ryd
              write(6,10) 'Correlation', isp, ik
              write(6,'(3g14.6)') real(sig(isp,ik)%scdiag,dp)*ryd, &
                aimag(sig(isp,ik)%scdiag)*ryd
              write(6,10) 'Vertex', isp, ik
              write(6,'(3g14.6)') real(sig(isp,ik)%sgdiag,dp)*ryd, &
                aimag(sig(isp,ik)%sgdiag)*ryd
           endif
           if (sig(isp,ik)%noffd > 0) then
              write(6,10) 'Correlation, off-diagonal', isp, ik
              write(6,'(2g14.6)') real(sig(isp,ik)%scoffd,dp)*ryd, &
                aimag(sig(isp,ik)%scoffd)*ryd
              write(6,10) 'Vertex, off-diagonal', isp, ik
              write(6,'(2g14.6)') real(sig(isp,ik)%sgoffd,dp)*ryd, &
                aimag(sig(isp,ik)%sgoffd)*ryd
           endif
        endif
        write(6,*)
#endif
     enddo
  enddo
10 format(/,1x,a,' (isp,ik) = ',i3,1x,i4)

end subroutine save_sigma
!===================================================================
