!===================================================================
!
! Sets up and solves the Bethe-Salpeter Equation, using input matrix
! elements of the Coulomb potential and TDLDA exchange-correlation
! potential in the basis of one-electron transitions.
! Energy of quasiparticle states, and optionally off-diagonal matrix
! elements of the (Sigma-Vxc) operator are also needed as input, and
! calculated using the "sigma" code.
!
! The interaction kernel has the mixing part (interaction between
! positive energy excitations and negative energy excitations)
! and energy dependence included via a constant energy reference, 
! that should be chosen to be close to the excitation energy of
! interest.
!
! One should figure out if the energy dependence is implemented
! correctly.
!
! Copyright (C) 2009 Murilo Tiago, Univ. of Texas, Austin, TX, USA
! mtiago@ices.utexas.edu
!
! First version written by Murilo Tiago, Univ. of Minnesota, April 2004
!
!    This program is free software; you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, version 1.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program; if not, write to the Free Software
!    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston MA  02110-1301 USA
!
!-------------------------------------------------------------------
program bsesolv

  use typedefs
  use mpi_module
  implicit none
#ifdef MPI
  include 'mpif.h'
#endif
  type (gspace) :: gvec
  type (kptinfo) :: kpt
  type (qptinfo) :: qpt, q_bse
  type (polinfo), dimension(2) :: pol_in, bsepol_in
  type (kernelinfo), dimension(:,:), allocatable :: k_p, k_vv, k_cc
  type (kernelinfo), dimension(:), allocatable :: k_vc, k_x
  type (polinfo), dimension(:,:), allocatable :: pol
  type (polinfo), dimension(:), allocatable :: bsepol
  type (options) :: opt

  character (len=800) :: lastwords
  character (len=40), allocatable :: routnam(:)

  logical :: nolda, tamm_d, writeig, trip_flag, trunc_c, &
       mix_flag, noxchange, snorm, readocc, hf, init_gr
  integer :: nspin, ii, irp, isp, iq, nbuff, lcache, nmap, dft_code
  real(dp) :: tsec(2), mem1, mem2, mem3, xsum, xmax, ecutb2, &
       tdldacut, bsecut, eref

  logical, allocatable :: wmap(:)

  ! W Gao ISDF method
  type (ISDF) :: isdf_in
  logical :: doisdf, verbose
  integer :: n_intp, intp_type, isdf_type, info
  integer, allocatable :: timerlist(:)

  ! WG debug
  integer :: outdbg
  character (len=20) :: dbg_filename

  !-------------------------------------------------------------------
  ! Initialization.
  !
  call header('BSESOLV')
  ! W Gao open dbg files
  write(dbg_filename,"(i7)") peinf%inode
  outdbg = peinf%inode+198812
  dbg_filename = "kernel_dbg"//adjustl(dbg_filename)
  open(outdbg,file=dbg_filename,status='unknown',iostat=info) 

  !-------------------------------------------------------------------
  ! Read input parameters from rgwbs.in.
  !
  call input_g(pol_in,qpt,tdldacut,nbuff,lcache,w_grp%npes, &
          nolda,tamm_d,r_grp%num,dft_code,doisdf,n_intp,intp_type,isdf_type,&
          isdf_in%lessmemory,isdf_in%fastselect,opt%eigsolver,opt%linear_algebra,opt%pdgesv_nbl2d,.false.)
  ! Currently, the ISDF method is not implemented in BSE package
  if (doisdf) then
     if (peinf%master) then
        write(6, *) "ISDF method is not implemented for BSE for now."
        write(6, *) "Proceed without using ISDF"
     endif
     doisdf = .false.
  endif
  call input_b(bsepol_in,q_bse,writeig,trip_flag,trunc_c,mix_flag, &
     noxchange,snorm,readocc,hf,bsecut,eref,ecutb2)
  !-------------------------------------------------------------------
  ! Determine the set of wavefunctions to read: if n-th wavefunction is
  ! needed, then wmap(n) = .true.; otherwise, wmap(n) = .false.
  !
  nmap = 0
  if (max(pol_in(1)%ncond,pol_in(1)%nval, &
       pol_in(2)%ncond,pol_in(2)%nval,bsepol_in(1)%ncond,bsepol_in(1)%nval, &
       bsepol_in(2)%ncond,bsepol_in(2)%nval) > 0) then
     if (pol_in(1)%ncond > 0) nmap = max(nmap,maxval(pol_in(1)%cmap))
     if (pol_in(1)%nval > 0) nmap = max(nmap,maxval(pol_in(1)%vmap))
     if (pol_in(2)%ncond > 0) nmap = max(nmap,maxval(pol_in(2)%cmap))
     if (pol_in(2)%nval > 0) nmap = max(nmap,maxval(pol_in(2)%vmap))
     if (bsepol_in(1)%ncond > 0) nmap = max(nmap,maxval(bsepol_in(1)%cmap))
     if (bsepol_in(1)%nval > 0) nmap = max(nmap,maxval(bsepol_in(1)%vmap))
     if (bsepol_in(2)%ncond > 0) nmap = max(nmap,maxval(bsepol_in(2)%cmap))
     if (bsepol_in(2)%nval > 0) nmap = max(nmap,maxval(bsepol_in(2)%vmap))
     allocate(wmap(nmap))
     wmap = .false.
     do isp = 1, 2
        do ii = 1, pol_in(isp)%ncond
           wmap(pol_in(isp)%cmap(ii)) = .true.
        enddo
        do ii = 1, pol_in(isp)%nval
           wmap(pol_in(isp)%vmap(ii)) = .true.
        enddo
     enddo
     do isp = 1, 2
        do ii = 1, bsepol_in(isp)%ncond
           wmap(bsepol_in(isp)%cmap(ii)) = .true.
        enddo
        do ii = 1, bsepol_in(isp)%nval
           wmap(bsepol_in(isp)%vmap(ii)) = .true.
        enddo
     enddo
  else
     allocate(wmap(1))
  endif
  if (min(pol_in(1)%ncond,pol_in(1)%nval,pol_in(2)%ncond,pol_in(2)%nval) < &
       0) then
     deallocate(wmap)
     allocate(wmap(1))
     nmap = 0
  endif
  if (min(bsepol_in(1)%ncond,bsepol_in(1)%nval,bsepol_in(2)%ncond, &
       bsepol_in(2)%nval) < 0) then
     deallocate(wmap)
     allocate(wmap(1))
     nmap = 0
  endif

  !-------------------------------------------------------------------
  ! Read wave-function file.
  !
  init_gr = .false.
  if ( dft_code == PARATEC ) then
     call paratec_wfn(gvec,kpt,nmap,nspin,wmap,init_gr)
  else
     call parsec_wfn(gvec,kpt,nmap,nspin,wmap,init_gr)
  endif
  deallocate(wmap)
  if (trip_flag .and. nspin == 2) then
     write(lastwords,*) 'ERROR: cannot do BSE triplet kernel with ', &
          'spin polarized wavefunctions. Stop.'
     call die(lastwords)
  endif
  !-------------------------------------------------------------------
  ! Calculate characters of representations.
  !
  if (kpt%lcplx) then
     call zcharac_group(gvec%syms,gvec,kpt,70,nspin,kpt%wfn(1,1)%nstate)
  else
     call dcharac_group(gvec%syms,gvec,kpt,70,nspin,kpt%wfn(1,1)%nstate)
  endif

  !-------------------------------------------------------------------
  ! Initialize pol and k_p structures.
  !
  allocate(pol(gvec%syms%ntrans,qpt%nk))
  pol(:,:)%lv = .false.
  pol(:,:)%ltv = .false.
  pol(:,:)%ntr = 0
  allocate(bsepol(gvec%syms%ntrans))
  bsepol(:)%ntr = 0
  allocate(k_p(gvec%syms%ntrans,qpt%nk))
  k_p(:,:)%ncol = 0
  k_p(:,:)%nbuff = nbuff
  k_p(:,:)%lcache = lcache
  k_p(:,:)%isdf = doisdf
  allocate(k_vv(gvec%syms%ntrans,qpt%nk))
  k_vv(:,:)%ncol = 0
  k_vv(:,:)%nbuff = nbuff
  k_vv(:,:)%lcache = lcache
  k_vv(:,:)%isdf = doisdf
  allocate(k_cc(gvec%syms%ntrans,qpt%nk))
  k_cc(:,:)%ncol = 0
  k_cc(:,:)%nbuff = nbuff
  k_cc(:,:)%lcache = lcache
  k_cc(:,:)%isdf = doisdf
  allocate(k_vc(gvec%syms%ntrans))
  k_vc(:)%ncol = 0
  k_vc(:)%nbuff = nbuff
  k_vc(:)%lcache = lcache
  k_vc(:)%isdf = doisdf
  allocate(k_x(gvec%syms%ntrans))
  k_x(:)%ncol = 0
  k_x(:)%nbuff = nbuff
  k_x(:)%lcache = lcache
  k_x(:)%isdf = doisdf
  call stopwatch(peinf%master,'Calling setup_b')
  call timacc(2,1,tsec)
  if (kpt%lcplx) then
     call zsetup_b(gvec,kpt,qpt,q_bse,pol_in,pol,bsepol_in,bsepol, &
          k_p,k_vv,k_cc,k_vc,k_x,nspin,readocc,tdldacut,bsecut)
  else
     call dsetup_b(gvec,kpt,qpt,q_bse,pol_in,pol,bsepol_in,bsepol, &
          k_p,k_vv,k_cc,k_vc,k_x,nspin,readocc,tdldacut,bsecut)
  endif
  print *, "finish setup_b"
  call timacc(2,2,tsec)

  !-------------------------------------------------------------------
  ! Print out warnings, information etc.
  !
  if (kpt%lcplx) tamm_d = .true.
  if (peinf%master) then
     write(6,'(/,a,/,/,a,/,2a,/)') repeat('-',65), &
          ' BSE input data: ', ' ', repeat('-',15)
     write(6,'(2a)') ' Number of transitions per representation ', &
          'in TDLDA polarizabillity:'
     write(6,'(8i8)') ((pol(ii,iq)%ntr,ii=1,gvec%syms%ntrans),iq=1,qpt%nk)
     write(6,'(a,i10,/)') ' total = ', sum(pol(:,:)%ntr)
     if (tdldacut > zero) then
        write(6,*) 'Energy cutoff applied in TDLDA polarizability = ', &
             tdldacut*ryd, ' eV'
     else
        write(6,*) 'No energy cutoff in TDLDA polarizability'
     endif
     if (nolda) write(6,'(a,/)') &
          ' LDA kernel is not included in polarizability'
     if (tamm_d) then
        write(6,'(2a,/)') ' Calculating TDLDA ', &
             'polarizability within the Tamm-Dancoff approximation.'
     else
        write(6,'(2a,/)') ' Not using the Tamm-Dancoff ', &
             'approximation in TDLDA polarizability.'
     endif
     write(6,'(2a)') ' Number of transitions per representation ', &
          'in BSE polarizabillity: '
     write(6,'(8i8)') (bsepol(ii)%ntr,ii=1,gvec%syms%ntrans)
     write(6,'(a,i10,/)') ' total = ', sum(bsepol(:)%ntr)
     if (bsecut > zero) then
        write(6,'(2a,f20.10,/)') ' Energy cutoff applied in ', &
             'BSE polarizability (eV) = ', bsecut*ryd
     else
        write(6,'(a,/)') ' No energy cutoff in BSE polarizability'
     endif
     if (snorm) then
        write(6,'(a,/)') 'Renormalizing Sum rule'
     else
        write(6,'(a,/)') 'Sum rule not renormalized'
     endif
     if (trip_flag) then
        write(6,'(a,/)') ' BSE for spin triplet excitations'
     else
        write(6,'(a,/)') ' BSE for spin singlet excitations'
     endif
     if (mix_flag) then
        write(6,'(a,/)') ' Including block mixing in BSE'
     else
        write(6,'(a,/)') ' No block mixing in BSE'
     endif
     if (eref > 0) then
        write(6,'(a,f10.4,/)') &
             ' Energy reference for dynamical kernel (eV) = ', eref
     else
        write(6,'(a,a,/)') ' Negative energy reference.', &
             ' Neglecting dynamical screening effects '
     endif
     if (ecutb2 > zero) then
        write(6,'(2a,g10.4,/)') ' Energy resolution in energy poles, ', &
             'BSE (eV) = ', ecutb2 * ryd
     else
        write(6,'(a,/)') ' Using no energy resolution in BSE, energy poles'
     endif
     if (hf) write(6,'(2a,/)') ' Using the ', &
          'Hartree-Fock approximation in self-energy'
     !
     ! WARNING: for periodic boundary conditions, use a truncated Coulomb
     ! potential (see Hanke's article) in diagonalization. That should be
     ! always done.
     !
     write(6,'(a,/)') &
          ' Coulomb potential is being truncated in the BSE equation.'
     !
     ! Estimate memory usage. Take into account only matrices that scale
     ! as N^4 with respect to the number of LDA states.
     ! For diagonalization, we store 4 matrices: hamiltonian/eigenvectors, 
     ! temporary array (in eigensolver), Hartree kernel, and XC kernel.
     ! Parallelized diagonalization also uses a second temporary matrix.
     ! We also need to store all the various kernel and potential matrices
     ! (there are lots fo them!).
     ! BSE diagonalization is not paralellized.
     !
     mem1 = sum(kpt%wfn(:,:)%nmem)*two/real(nspin,dp) * &
          real(nspin*gvec%nr,dp)/two/131072.d0 / real(w_grp%npes,dp)
     if (kpt%lcplx) mem1 = mem1 * two
     write(6,'(/,a,f10.2,a)') ' Memory needed to store wavefunctions : ', &
           mem1,' MB/proc.'
     mem1 = real(5*gvec%nr,dp)/131072.d0
     if (kpt%lcplx) mem1 = mem1 * two
     write(6,'(a,a,f10.2,a)') ' Memory needed to calculate kernel ', &
          'matrix elements : ',mem1,' MB/proc.'
     xmax = 0
     do iq = 1, qpt%nk
        do irp = 1, gvec%syms%ntrans
           xsum = real(pol(irp,iq)%ntr * pol(irp,iq)%ntr,dp)
           if (xmax < xsum) xmax = xsum
        enddo
     enddo
     mem1 = xmax/1024.d0*3.d0/128.d0/r_grp%npes
     if (r_grp%npes > 1) mem1 = mem1*4.d0/3.d0
     xmax = 0
     do iq = 1, qpt%nk
        do irp = 1, gvec%syms%ntrans
           xsum = real(k_x(irp)%nrow * k_x(irp)%ncol,dp) + &
                real(pol(irp,iq)%ntr * pol(irp,iq)%ntr,dp) / &
                real(r_grp%npes,dp)
           if (xsum > xmax) xmax = xsum
           xsum = real(k_vv(irp,iq)%nrow * k_vv(irp,iq)%ncol,dp) * two + &
                real(pol(irp,iq)%ntr * k_vv(irp,iq)%ncol,dp) + &
                real(pol(irp,iq)%ntr * pol(irp,iq)%ntr,dp) / &
                real(r_grp%npes,dp)
           if (xsum > xmax) xmax = xsum
           xsum = real(k_vv(irp,iq)%nrow * k_vv(irp,iq)%ncol,dp) + &
                real(pol(irp,iq)%ntr * k_vv(irp,iq)%ncol,dp) * two + &
                real(pol(irp,iq)%ntr * pol(irp,iq)%ntr,dp) / &
                real(r_grp%npes,dp)
           if (xsum > xmax) xmax = xsum
           xsum = real(k_cc(irp,iq)%nrow * k_cc(irp,iq)%ncol,dp) * two + &
                real(pol(irp,iq)%ntr * k_cc(irp,iq)%ncol,dp) + &
                real(pol(irp,iq)%ntr * k_vv(irp,iq)%ncol,dp) * two + &
                real(pol(irp,iq)%ntr * pol(irp,iq)%ntr,dp) / &
                real(r_grp%npes,dp)
           if (xsum > xmax) xmax = xsum
           xsum = real(k_cc(irp,iq)%nrow * k_cc(irp,iq)%ncol,dp) + &
                real(pol(irp,iq)%ntr * k_cc(irp,iq)%ncol,dp) * two + &
                real(pol(irp,iq)%ntr * k_vv(irp,iq)%ncol,dp) * two + &
                real(pol(irp,iq)%ntr * pol(irp,iq)%ntr,dp) / &
                real(r_grp%npes,dp)
           if (xsum > xmax) xmax = xsum
           xsum = real(k_vc(irp)%nrow * k_vc(irp)%ncol,dp) * two + &
                real(pol(irp,iq)%ntr * k_vc(irp)%ncol,dp) + &
                real(pol(irp,iq)%ntr * pol(irp,iq)%ntr,dp) / &
                real(r_grp%npes,dp)
           if (xsum > xmax) xmax = xsum
           xsum = real(k_vc(irp)%nrow * k_vc(irp)%ncol,dp) + &
                real(pol(irp,iq)%ntr * k_vc(irp)%ncol,dp) * two + &
                real(pol(irp,iq)%ntr * pol(irp,iq)%ntr,dp) / &
                real(r_grp%npes,dp)
           if (xsum > xmax) xmax = xsum
           xmax = xmax + real(bsepol(irp)%ntr * bsepol(irp)%ntr,dp)
        enddo
     enddo
     mem2 = xmax/1024.d0/128.d0/r_grp%npes

     ! Single processor diagonalization.
     xmax = 0
     do irp = 1, gvec%syms%ntrans
        xmax = xmax + real(bsepol(irp)%ntr * bsepol(irp)%ntr,dp)
     enddo
     mem3 = xmax/1024.d0*three/128.d0

     write(6,'(/,a,f10.2,a)') &
          ' Memory needed for polarizability diagonalization : ', &
          mem1, ' MB/proc.'
     write(6,'(a,f10.2,a)') ' Memory needed in potentials calculation : ', &
          mem2, ' MB/proc.'
     write(6,'(a,f10.2,a)') &
          ' Memory needed for single proc. BSE diagonalization : ', &
          mem3, ' MB/proc.'
     write(6,'(/,a,/)') repeat('-',65)
     call flush(6)
  endif

  ! Convert input energy values to rydbergs.
  eref = eref/ryd
  ! We use only the square of energy resolution.
  ecutb2 = ecutb2 * ecutb2

  !-------------------------------------------------------------------
  ! Do calculation.
  !
  if (kpt%lcplx) then
     call zcalculate_bse(gvec,kpt,qpt,q_bse,bsepol_in,k_p,k_vv,k_cc, &
          k_vc,k_x,pol,bsepol,nolda,tamm_d,trip_flag,trunc_c,writeig, &
          mix_flag,noxchange,snorm,hf,nspin,ecutb2,eref,isdf_in,opt)
  else
     call dcalculate_bse(gvec,kpt,qpt,q_bse,bsepol_in,k_p,k_vv,k_cc, &
          k_vc,k_x,pol,bsepol,nolda,tamm_d,trip_flag,trunc_c,writeig, &
          mix_flag,noxchange,snorm,hf,nspin,ecutb2,eref,isdf_in,opt)
  endif

  if (peinf%master) call delete_file(90,'bse_chkpt.dat')
  close(outdbg)
  !-------------------------------------------------------------------
  ! Time accounting.
  !
  irp = 9
  isp = 3
  allocate(routnam(irp+isp))
  allocate(timerlist(irp+isp))
  routnam(1)='SETUP_B:'           ; timerlist(1)=2
  routnam(2)='KERNEL:'            ; timerlist(2)=3
  routnam(3)='DIAG_POL:'          ; timerlist(3)=4
  routnam(4)='POTENTIAL:'         ; timerlist(4)=5
  routnam(5)='EXCHANGE_B:'        ; timerlist(5)=6
  routnam(6)='DIRECT_B:'          ; timerlist(6)=7
  routnam(7)='DIRECT_S:'          ; timerlist(7)=8
  routnam(8)='DIRECT_MIX:'        ; timerlist(8)=9
  routnam(9)='DIAG_BSE:'          ; timerlist(9)=10

  routnam(10)='POISSON_FFT:'      ; timerlist(10)=11
  routnam(11)='EIGENSOLVER:'      ; timerlist(11)=12
  routnam(12)='INTEGRATION:'      ; timerlist(12)=13

  call timacc(1,2,tsec)
  call finalize(peinf%master,peinf%comm,irp,isp,routnam,timerlist)

end program bsesolv
!===================================================================
