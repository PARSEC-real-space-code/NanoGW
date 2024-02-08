!===================================================================
!
! If chkpt < 0 at input, exit without doing anything.
! If chkpt >= 0 at input, open file sigma.chkpt.dat and check sanity of
! file. If it passes sanity tests, then the following parameters will
! have value modified at output:
!     chkpt = irp + (nq + iq)*nrep
!     sig(:)%xdiag
!     sig(:)%sexdiag
!     sig(:)%scdiag
!     sig(:)%scsdiag
!     sig(:)%sgdiag
!     sig(:)%sgsdiag
!     sig(:)%xoffd
!     sig(:)%scoffd
!     sig(:)%scsoffd
!     sig(:)%sgoffd
!     sig(:)%sgsoffd
!
! Copyright (C) 2009 Murilo L. Tiago, http://users.ices.utexas.edu/~mtiago
! This file is part of RGWBS. It is distributed under the GPL v1.
!
!-------------------------------------------------------------------
subroutine read_sigma(nkpt,nspin,sig_in,sig,chkpt)

  use typedefs
  implicit none

  ! arguments
  type (siginfo), intent(in) :: sig_in
  type (siginfo), dimension(nspin,nkpt), intent(inout) :: sig
  integer, intent(in) :: nkpt, nspin
  integer, intent(inout) :: chkpt

  ! local variables
  type (siginfo) :: sig_tmp
  integer :: isp, ik, ii, chkpt_in
  logical :: lcheck
  ! number of input unit
  integer, parameter :: iunit = 60

  lcheck = .false.
  if (chkpt < 0) return

  open(iunit,file='sigma.chkpt.dat',form='unformatted',status='old',iostat=ii)
  if (ii == 0) then
     read(iunit) chkpt_in, ii
     call map_check('sigma.chkpt.dat','spin',1,nspin,ii,lcheck)
     if (.not. lcheck) goto 8

     do ik = 1, nkpt
        do isp = 1, nspin
           read(iunit) sig_tmp%nmap, sig_tmp%indxk
           read(iunit) sig_tmp%ndiag_s, sig_tmp%ndiag
           read(iunit) sig_tmp%noffd_s, sig_tmp%noffd
           read(iunit) sig_tmp%nen, sig_tmp%deltae

           call map_check('sigma.chkpt.dat','number of QP wavefunctions', &
                1,sig(isp,ik)%nmap,sig_tmp%nmap,lcheck)
           if (.not. lcheck) goto 8

           call map_check('sigma.chkpt.dat','number of diagonal matrix elements', &
                1,sig(isp,ik)%ndiag_s,sig_tmp%ndiag_s,lcheck)
           if (.not. lcheck) goto 8

           call map_check('sigma.chkpt.dat','number of diagonal matrix elements', &
                1,sig(isp,ik)%ndiag,sig_tmp%ndiag,lcheck)
           if (.not. lcheck) goto 8

           call map_check('sigma.chkpt.dat','number of off-diagonal matrix elements', &
                1,sig(isp,ik)%noffd_s,sig_tmp%noffd_s,lcheck)
           if (.not. lcheck) goto 8

           call map_check('sigma.chkpt.dat','number of off-diagonal matrix elements', &
                1,sig(isp,ik)%noffd,sig_tmp%noffd,lcheck)
           if (.not. lcheck) goto 8

           call map_check('sigma.chkpt.dat','number of energy points', &
                1,sig_in%nen,sig_tmp%nen,lcheck)
           if (.not. lcheck) goto 8

           call vector_check('sigma.chkpt.dat','Energy range (eV)', &
                1,sig_in%deltae,sig_tmp%deltae,lcheck)
           if (.not. lcheck) goto 8

           if (sig(isp,ik)%nmap > 0) then
              allocate(sig_tmp%map(sig_tmp%nmap))
              read(iunit) (sig_tmp%map(ii),ii=1,sig_tmp%nmap)
              call map_check('sigma.chkpt.dat','index of QP wavefunctions', &
                   sig_tmp%nmap,sig(isp,ik)%map,sig_tmp%map,lcheck)
              deallocate(sig_tmp%map)
              if (.not. lcheck) goto 8
           endif

           if (sig(isp,ik)%ndiag_s > 0) then

              allocate(sig_tmp%diag(sig_tmp%ndiag_s))
              read(iunit) (sig_tmp%diag(ii),ii=1,sig_tmp%ndiag_s)
              call map_check('sigma.chkpt.dat','index of diagonal matrix elements', &
                   sig_tmp%ndiag_s,sig(isp,ik)%diag,sig_tmp%diag,lcheck)
              deallocate(sig_tmp%diag)
              if (.not. lcheck) goto 8

              read(iunit)
              read(iunit)
              read(iunit)
              read(iunit) (sig(isp,ik)%xdiag(ii),ii=1,sig(isp,ik)%ndiag_s)
              read(iunit) (sig(isp,ik)%scsdiag(ii),ii=1,sig(isp,ik)%ndiag_s)
              if (sig_in%xc == XC_GW) then
                 read(iunit) (sig(isp,ik)%sgsdiag(ii),ii=1,sig(isp,ik)%ndiag_s)
              endif
           endif

           if (sig_in%xc == XC_GW) then
              if (sig(isp,ik)%ndiag > 0) then
                 read(iunit) (sig(isp,ik)%scdiag(:,ii),ii=1,2*sig(isp,ik)%ndiag)
                 read(iunit) (sig(isp,ik)%sexdiag(:,ii),ii=1,2*sig(isp,ik)%ndiag)
                 read(iunit) (sig(isp,ik)%sgdiag(:,ii),ii=1,2*sig(isp,ik)%ndiag)
              endif
           endif

           if (sig(1,ik)%noffd_s + sig(nspin,ik)%noffd_s > 0) then

              allocate(sig_tmp%off1(sig_tmp%noffd_s))
              read(iunit) (sig_tmp%off1(ii),ii=1,sig_tmp%noffd_s)
              call map_check('sigma.chkpt.dat','index of off-diagonal matrix elements', &
                   sig_tmp%noffd_s,sig(isp,ik)%off1,sig_tmp%off1,lcheck)
              deallocate(sig_tmp%off1)
              if (.not. lcheck) goto 8

              allocate(sig_tmp%off2(sig_tmp%noffd_s))
              read(iunit) (sig_tmp%off2(ii),ii=1,sig_tmp%noffd_s)
              call map_check('sigma.chkpt.dat','index of off-diagonal matrix elements', &
                   sig_tmp%noffd_s,sig(isp,ik)%off2,sig_tmp%off2,lcheck)
              deallocate(sig_tmp%off2)
              if (.not. lcheck) goto 8

              read(iunit)
              read(iunit) (sig(isp,ik)%xoffd(ii),ii=1,sig(isp,ik)%noffd_s)
              read(iunit) (sig(isp,ik)%scsoffd(ii),ii=1,sig(isp,ik)%noffd_s)
              if (sig_in%xc == XC_GW) then
                 read(iunit) (sig(isp,ik)%sgsoffd(ii),ii=1,sig(isp,ik)%noffd_s)
              endif
           endif
           if (sig_in%xc == XC_GW) then
              if (sig(isp,ik)%noffd > 0) then
                 read(iunit) (sig(isp,ik)%scoffd(:,ii),ii=1,2*sig(isp,ik)%noffd)
                 read(iunit) (sig(isp,ik)%sgoffd(:,ii),ii=1,2*sig(isp,ik)%noffd)
              endif
           endif

        enddo
     enddo
  endif

8 continue

  close(iunit)
  if (lcheck) chkpt = chkpt_in
  write(6,'(/,a,i5,/)') ' Checkpoint control = ',chkpt

end subroutine read_sigma
!===================================================================
