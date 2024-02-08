!===================================================================
!
! Read input wavefunctions/parameters from file GWC, created by PARATEC.
!
! Only Master PE reads input, and broadcasts all info to the other PEs.
!
! Besides reading input, this subroutine performs a series of tasks:
!   constructs real-space and reciprocal-space grids;
!   sets up symmetry group and calculates representation
!      (inside subroutines);
!   reads electron density from GWC and constructs the DFT 
!      exchange-correlation potential and kernel (derivative of potential);
!   defines the FFT grid.
!
! Copyright (C) 2009 Murilo L. Tiago, http://users.ices.utexas.edu/~mtiago
! This file is part of RGWBS. It is distributed under the GPL v1.
!
!-------------------------------------------------------------------
subroutine paratec_wfn(gvec,kpt,nmap,nspin,wmap,init_gr)

  use typedefs
  use mpi_module
  use fft_module
  use fd_module
  implicit none

#ifdef MPI
  include 'mpif.h'
#endif

  ! arguments
  ! specification of real-space grid
  type (gspace), intent(inout) :: gvec
  ! k-point structure
  type (kptinfo), intent(inout) :: kpt
  ! number of orbitals to be saved for each spin channel and k-point
  integer, intent(in) :: nmap
  ! number of spin channels (1 or 2)
  integer, intent(out) :: nspin
  ! list of orbitals to be saved for each spin channel and k-point
  logical, intent(in) :: wmap(nmap)
  ! true if spacial derivatives (e.g.: gradient, laplacian) will be performed
  logical, intent(inout) ::  init_gr

  ! Local variables
  integer, parameter :: infile = 25  ! number of input unit
  character (len=800) :: lastwords
  integer :: ii, jj, nn, isp, ikp, icount, info, ndim, nwedge, itmp(3)
  integer :: i1, i2, i3, j1, j2, j3
  real(dp):: norm, rtmp, mtmp(3,3)

  integer, allocatable :: rindex(:), gtable(:,:), invindex(:)
  integer, allocatable :: ifband(:,:), map_found(:,:)
  real(dp), allocatable :: rho_in(:,:)
  complex(dpc), allocatable :: zcg(:), zwfn(:)
#ifdef MPI
  integer :: status(MPI_STATUS_SIZE)
#endif

  !-------------------------------------------------------------------
  ! Define fixed parameters. These parameters should be used whenever
  ! wavefunctions are input in plane-wave basis.
  !
  ! wavefunctions are complex in real space.
  kpt%lcplx = .true.
  ! bulk boundary conditions, periodic along all cartesian directions
  gvec%per = 3
  gvec%rmax = zero
  gvec%shift = zero
  init_gr = .true.

  ! trivial group symmetry
  gvec%syms%ntrans = 1
  allocate(gvec%syms%trans(3,3,gvec%syms%ntrans))
  allocate(gvec%syms%chi(gvec%syms%ntrans,gvec%syms%ntrans))
  gvec%syms%trans = zero
  do ii = 1, 3
     gvec%syms%trans(ii,ii,1) = one
  enddo
  gvec%syms%chi = 1

  fd%norder = -1
  !-------------------------------------------------------------------
  ! Read info about reciprocal space from GWC. Ignore symmetry operations.
  !
  if(peinf%master) then
     write(6,'(/,a,/,a,a,/)') ' Electronic structure data :', &
          ' ', repeat('-',27)
     write(6,*) 'Reading electron wavefunctions from file GWC'
     open(infile,file='GWC',form='unformatted',status='old',iostat=info)
     if (info /= 0) call die('ERROR: GWC file not found')
     read(infile) gvec%bdot
     read(infile) gvec%celvol
     mtmp = gvec%bdot
     call mtrxin(mtmp,rtmp,info)
     if (info /= 0) call die('ERROR in mtrxin. Stop.')
     rtmp = abs( gvec%celvol - eight * pi * pi * pi / sqrt(abs(rtmp)) )
      if ( rtmp > 1.d-6 ) then
        write(lastwords,*) 'ERROR: volume mismatch. ', rtmp
        call die(lastwords)
     endif
     read(infile) nn
     do ii = 1, 2*nn
        read(infile)
     enddo
     read(infile) nspin
     read(infile) kpt%nk
     read(infile)
     read(infile)
     allocate(kpt%weight(kpt%nk))
     read(infile) kpt%weight
     allocate(kpt%fk(3,kpt%nk))
     do ikp = 1, kpt%nk
        read(infile) kpt%fk(:,ikp)
     enddo
     read(infile) icount
     allocate(ifband(kpt%nk,nspin))
     read(infile) ifband
     if (maxval(abs(ifband)) /= 1) then
        write(lastwords,*) 'ERROR: GWC seems to have missing bands. ', ifband
        call die(lastwords)
     endif
     read(infile) ifband
  endif
#ifdef MPI
  call MPI_BCAST(gvec%bdot,9, &
       MPI_DOUBLE_PRECISION,peinf%masterid,peinf%comm,info)
  call MPI_BCAST(gvec%celvol,1, &
       MPI_DOUBLE_PRECISION,peinf%masterid,peinf%comm,info)
  call MPI_BCAST(nspin,1, &
       MPI_INTEGER,peinf%masterid,peinf%comm,info)
  call MPI_BCAST(kpt%nk, 1, &
       MPI_INTEGER,peinf%masterid,MPI_COMM_WORLD,info)
  if (.not. peinf%master) allocate(kpt%fk(3,kpt%nk))
  call MPI_BCAST(kpt%fk, kpt%nk*3, &
       MPI_DOUBLE_PRECISION,peinf%masterid,MPI_COMM_WORLD,info)
  if (.not. peinf%master) allocate(kpt%weight(kpt%nk))
  call MPI_BCAST(kpt%weight, kpt%nk, &
       MPI_DOUBLE_PRECISION,peinf%masterid,MPI_COMM_WORLD,info)
  call MPI_BCAST(icount,1, &
       MPI_INTEGER,peinf%masterid,peinf%comm,info)
#endif

  !-------------------------------------------------------------------
  ! Read DFT eigenvalues and related data.
  !
  allocate(kpt%wfn(nspin,kpt%nk))
  do ikp = 1, kpt%nk
     do isp = 1, nspin
        kpt%wfn(isp,ikp)%nstate = icount
        allocate(kpt%wfn(isp,ikp)%e0(kpt%wfn(isp,ikp)%nstate))
        allocate(kpt%wfn(isp,ikp)%occ0(kpt%wfn(isp,ikp)%nstate))
        if (peinf%master) then
           read(infile) (kpt%wfn(isp,ikp)%e0(ii),ii=1,kpt%wfn(isp,ikp)%nstate)
           do ii = 1, ifband(ikp,isp)
              kpt%wfn(isp,ikp)%occ0(ii) = one
           enddo
           do ii = ifband(ikp,isp) + 1, kpt%wfn(isp,ikp)%nstate
              kpt%wfn(isp,ikp)%occ0(ii) = zero
           enddo
        endif
#ifdef MPI
        call MPI_BCAST(kpt%wfn(isp,ikp)%e0,kpt%wfn(isp,ikp)%nstate, &
             MPI_DOUBLE_PRECISION,peinf%masterid,peinf%comm,info)
        call MPI_BCAST(kpt%wfn(isp,ikp)%occ0,kpt%wfn(isp,ikp)%nstate, &
             MPI_DOUBLE_PRECISION,peinf%masterid,peinf%comm,info)
#endif
        allocate(kpt%wfn(isp,ikp)%e1(kpt%wfn(isp,ikp)%nstate))
        call read_scissors(peinf%master,isp,kpt%wfn(isp,ikp)%nstate, &
             kpt%wfn(isp,ikp)%e0,kpt%wfn(isp,ikp)%e1)
        allocate(kpt%wfn(isp,ikp)%irep(kpt%wfn(isp,ikp)%nstate))
        kpt%wfn(isp,ikp)%irep = 1
     enddo
  enddo
  !-------------------------------------------------------------------
  ! Read reciprocal space grid and define FFT grid. Nfft must be even number!
  ! nfft(4) is used only in real transforms.
  !
  if (peinf%master) then
     read(infile)
     read(infile) ndim
     allocate(gtable(3,ndim))
     do jj = 1, ndim
        read(infile) (gtable(ii,jj),ii=1,3)
     enddo
     itmp = 2 * maxval(gtable,dim=2) + 2
     call adjust_fft(itmp)
     deallocate(gtable)
     ndim = 0

     mtmp = gvec%bdot
     call mtrxin(mtmp,rtmp,info)
     if (info /= 0) call die('ERROR in mtrxin. Stop.')
     do ii = 1, 3
        gvec%step(ii) = two * pi * sqrt(mtmp(ii,ii)) / real(itmp(ii),dp)
     enddo
  endif
#ifdef MPI
  call MPI_BCAST(itmp,3, &
       MPI_INTEGER,peinf%masterid,peinf%comm,info)
  call MPI_BCAST(gvec%step,3, &
       MPI_DOUBLE_PRECISION,peinf%masterid,peinf%comm,info)
#endif
  gvec%hcub = gvec%celvol / real(itmp(1) * itmp(2) * itmp(3))

  do ii = 1, 3
     fft_box%nfft(ii) = itmp(ii)
  enddo
  fft_box%nfft(4) = fft_box%nfft(1) + 2

  !-------------------------------------------------------------------
  ! Define electron density. We compute it from FFTs of the occupied
  ! wavefunctions. Because of grid discretization, there could be an
  ! error in the density!
  ! After it is done, rewind input file.
  !
  nwedge = fft_box%nfft(1)*fft_box%nfft(2)*fft_box%nfft(3)
  allocate(rho_in(nwedge,3))
  rho_in = zero
  icount = 0

  call zinitialize_FFT(peinf%inode,fft_box)
  do ikp = 1, kpt%nk
     if (peinf%master) read(infile) ndim
#ifdef MPI
     call MPI_BCAST(ndim,1,MPI_INTEGER,peinf%masterid,peinf%comm,info)
#endif
     allocate(gtable(3,ndim))
     if (peinf%master) then
        do ii = 1, ndim
           read(infile) (gtable(jj,ii),jj=1,3)
        enddo
     endif
#ifdef MPI
     call MPI_BCAST(gtable,3*ndim,MPI_INTEGER,peinf%masterid,peinf%comm,info)
#endif
     allocate(zcg(ndim*nspin))
     do ii = 1, kpt%wfn(1,ikp)%nstate
        if (peinf%master) then
           read(infile) (zcg(jj),jj=1,ndim*nspin)
        endif
        do isp = 1, nspin
           if (kpt%wfn(isp,ikp)%occ0(ii) > zero) then
#ifdef MPI
              icount = mod(ii-1,peinf%npes)
              if (icount /= peinf%masterid .and. peinf%master) &
                   call MPI_SEND(zcg,ndim*nspin,MPI_DOUBLE_COMPLEX, &
                   icount,icount,peinf%comm,info)
              if (icount == peinf%inode .and. (.not. peinf%master)) &
                   call MPI_RECV(zcg,ndim*nspin,MPI_DOUBLE_COMPLEX, &
                   peinf%masterid,icount,peinf%comm,status,info)
#endif
              if (icount == peinf%inode) then
                 fft_box%zbox = zzero
                 do jj = 1, ndim
                    j3 = gtable(3,jj)
                    i3 = j3 + 1
                    if (j3 < 0) i3 = fft_box%nfft(3) + i3
                    i3 = i3 - fft_box%nfft(3)/2-1
                    j2 = gtable(2,jj)
                    i2 = j2 + 1
                    if (j2 < 0) i2 = fft_box%nfft(2) + i2
                    i2 = i2 - fft_box%nfft(2)/2 - 1
                    j1 = gtable(1,jj)
                    i1 = j1 + 1
                    if (j1 < 0) i1 = fft_box%nfft(1) + i1
                    i1 = i1 - fft_box%nfft(1)/2 - 1
                    fft_box%zbox(i1,i2,i3) = zcg(jj + ndim*(isp-1))
                 enddo
                 call zdo_FFT(1,fft_box)
                 fft_box%zbox = sqrt(fft_box%scale) * fft_box%zbox
                 jj = 0
                 do j3 = -fft_box%nfft(3)/2, fft_box%nfft(3)/2 - 1
                    do j2 = -fft_box%nfft(2)/2, fft_box%nfft(2)/2 - 1
                       do j1 = -fft_box%nfft(1)/2, fft_box%nfft(1)/2 - 1
                          jj = jj + 1
                          rho_in(jj,isp + 1) = rho_in(jj,isp + 1) + &
                                kpt%weight(ikp) * kpt%wfn(isp,ikp)%occ0(ii) * &
                               fft_box%zbox(j1,j2,j3)*conjg(fft_box%zbox(j1,j2,j3))
                       enddo
                    enddo
                 enddo
              endif
           endif
        enddo
#ifdef MPI
        if (mod(ii-1,peinf%npes) == peinf%npes - 1) &
             call MPI_BARRIER(peinf%comm,info)
#endif
     enddo
#ifdef MPI
     call MPI_BARRIER(peinf%comm,info)
#endif
     deallocate(zcg)
     deallocate(gtable)
  enddo

  if (peinf%master) then
     rewind(infile)
     read(infile)
     read(infile)
     read(infile) nn
     do ii = 1, 2*nn
        read(infile)
     enddo
     do ii = 1, 9 + (nspin + 1) * kpt%nk
        read(infile)
     enddo
     read(infile) nn
     do ii = 1, nn
        read(infile)
     enddo
  endif

  rho_in(:,1) = rho_in(:,2) + rho_in(:,3)
  rho_in = rho_in * two / real(nspin,dp) / gvec%hcub
  call dpsum(nwedge*3,peinf%npes,peinf%comm,rho_in)

  allocate(gtable(3,nwedge))
  jj = 0
  do j3 = -fft_box%nfft(3)/2, fft_box%nfft(3)/2 - 1
     do j2 = -fft_box%nfft(2)/2, fft_box%nfft(2)/2 - 1
        do j1 = -fft_box%nfft(1)/2, fft_box%nfft(1)/2 - 1
           jj = jj + 1
           gtable(1,jj) = j1
           gtable(2,jj) = j2
           gtable(3,jj) = j3
        enddo
     enddo
  enddo
  ndim = nwedge

  !-------------------------------------------------------------------
  ! Never truncate the Coulomb potential in a solid.
  !
  fft_box%trunc = -one

  !-------------------------------------------------------------------
  ! Define the long-wavelength part of Coulomb potential.
  !
  gvec%long = vq_average(3,kpt%nk,gvec%celvol,gvec%rmax)

  !---------------------------------------------------------------
  ! Master PE print data out.
  !
  if (peinf%master) then
     if (nspin == 2) then
        write(6,'(/,a,/)') ' Spin polarized calculation.'
     else
        write(6,'(/,a,/)') ' Spin unpolarized calculation.'
     endif
     do ikp = 1, kpt%nk
        do isp = 1, nspin
           write(6,'(a,3f13.6)') ' K-point = ', kpt%fk(:,ikp)
           write(6,*) ' Number of occupied bands, spin ', isp, &
                ' : ', ifband(ikp,isp)
           write(6,*) 'DFT Eigenvalues (eV), spin ', isp, ' :'
           write(6,'(10f7.2)') &
                (kpt%wfn(isp,ikp)%e0(ii)*ryd,ii=1,kpt%wfn(isp,ikp)%nstate)
           write(6,*)
        enddo
     enddo
     write(6,'(2a,e16.8)') ' Maximum value of electron density ', &
          '(a.u.^-3) = ', maxval(rho_in(:,1))
     write(6,'(a,f10.4)') ' Total number of electrons = ', &
          sum(rho_in(:,1)) * gvec%hcub * real(gvec%syms%ntrans,dp)
     deallocate(ifband)
  endif

  !-------------------------------------------------------------------
  ! Reorder grid points in descending value of electron density.
  ! gtable is the full, original grid read from GWC.
  ! gvec%r is the set of grid points in wedge, ordered according to
  ! electron density.
  !
  allocate(rindex(nwedge))
  rindex = 0
  rho_in = -rho_in
  call quicksort(nwedge,rho_in(1,1),rindex)
  rho_in = -rho_in
  gvec%nr = nwedge

  allocate(gvec%r(3,gvec%nr))
  allocate(invindex(nwedge))
  do jj = 1, gvec%nr
     ii = rindex(jj)
     gvec%r(:,jj) = gtable(:,ii)
     invindex(ii) = jj
  enddo
  deallocate(gtable)

  !-------------------------------------------------------------------
  ! Define MPI goups.
  !
  call create_r_grp(gvec%syms%ntrans)
  call create_w_grp()

  !-------------------------------------------------------------------
  ! Define exchange-correlation matrix elements.
  !
  call paratec_xc(nspin,ii,jj)
  call define_xc(nspin,ii,jj,init_gr)

  !-------------------------------------------------------------------
  ! Set distributed grid.
  !
  call paratec_atoms(gvec,peinf%master,info)
  if (info /= 0) call die(' ')
  if (peinf%master) write(6,'(/,a,/,a,/)') ' Grid data :',' -----------'
  allocate(gtable(3,gvec%nr*gvec%syms%ntrans),stat=info)
  call alccheck('gtable','paratec_wfn',3*gvec%nr*gvec%syms%ntrans,info)
  call grid_setup(gvec,gtable)
  deallocate(gtable)

  !-------------------------------------------------------------------
  ! More printouts.
  !
  if (peinf%master) then
        write(6,'(a,/)') ' Gradients and Laplacians calculated with FFT'
        write(6,*) ' BULK BOUNDARY CONDITIONS !'
        write(6,'(/,a)') ' Unit lattice vectors (a.u.) = '
           write(6,'(3f20.15)') gvec%avec
     write(6,'(/,a,3f16.6,/)') ' Grid spacing (a.u.) = ', gvec%step
     write(6,'(a,f16.6)') ' Volume of periodic cell (a.u.^3) = ', gvec%celvol
     write(6,'(/,a,4i4,/)') ' Real-space FFT grid = ', fft_box%nfft
     write(6,'(a,f20.10,/)') ' Long-wavelength Coulomb energy (eV) = ', &
          gvec%long * ryd
     if (fft_box%trunc > zero) then
        write(6,'(a,f20.10)') ' Coulomb truncation length (a.u.) = ', &
             fft_box%trunc
     else
        write(6,*) 'No Coulomb truncation'
     endif
     write(6,'(/,a,i10,a)') ' Full grid has ', ndim, ' points.'
     write(6,'(a,i10,a,/)') ' Irreducible wedge has ', gvec%nr, ' points.'
  endif

  !-------------------------------------------------------------------
  ! Determine which wavefunctions to keep in memory.
  !
  do isp = 1, nspin
     do ikp = 1, kpt%nk
        allocate(kpt%wfn(isp,ikp)%map(kpt%wfn(isp,ikp)%nstate))
        kpt%wfn(isp,ikp)%map = 0
        icount = 0
        if (nmap == 0) then
           icount = kpt%wfn(isp,ikp)%nstate
           do ii = 1, kpt%wfn(isp,ikp)%nstate
              kpt%wfn(isp,ikp)%map(ii) = ii
           enddo
        else
           do ii = 1, nmap
              if (.not. wmap(ii)) cycle
              icount = icount + 1
              kpt%wfn(isp,ikp)%map(ii) = icount
              if (ii > kpt%wfn(isp,ikp)%nstate) then
                 write(lastwords,*) 'ERROR: You requested at least ', ii, &
                      ' states but GWC contains only ', &
                      kpt%wfn(isp,ikp)%nstate, ' Stop.'
                 call die(lastwords)
              endif
           enddo
        endif
        kpt%wfn(isp,ikp)%nmem = icount
        allocate(kpt%wfn(isp,ikp)%imap(kpt%wfn(isp,ikp)%nmem))
        kpt%wfn(isp,ikp)%imap = 0
        do ii = 1, kpt%wfn(isp,ikp)%nstate
           if (kpt%wfn(isp,ikp)%map(ii) /= 0) &
                kpt%wfn(isp,ikp)%imap(kpt%wfn(isp,ikp)%map(ii)) = ii
        enddo
     enddo
  enddo

  !-------------------------------------------------------------------
  ! Save electron density.
  !
  allocate(kpt%rho(gvec%nr,nspin),stat=info)
  call alccheck('kpt%rho','paratec_wfn',gvec%nr*nspin,info)
  kpt%rho = zero
  if (peinf%master) then
     do isp = 1, nspin
        do ii = 1, gvec%nr
           rtmp = rho_in(rindex(ii),isp+1)
           if (rtmp > zero) kpt%rho(ii,isp) = rtmp
        enddo
     enddo
  endif
#ifdef MPI
  call MPI_BCAST(kpt%rho,gvec%nr*nspin, &
       MPI_DOUBLE_PRECISION,peinf%masterid,peinf%comm,info)
#endif
  deallocate(rho_in)

  !-------------------------------------------------------------------
  ! Read in wavefunctions, calculate characters and store the important
  ! ones in memory.
  !
  call stopwatch(peinf%master,'paratec_wfn:  reading wave functions')
  w_grp%nr = gvec%nr
  w_grp%ldn = w_grp%nr / w_grp%npes + 1
  w_grp%mydim = w_grp%nr / w_grp%npes
  nn = mod(w_grp%nr, w_grp%npes)
  w_grp%offset = w_grp%ldn*w_grp%inode
  if (w_grp%inode < nn) w_grp%mydim = w_grp%mydim + 1
  if (w_grp%inode > nn) w_grp%offset = w_grp%offset - w_grp%inode + nn
  do isp = 1, nspin
     do ikp = 1, kpt%nk
        allocate(kpt%wfn(isp,ikp)%zwf(w_grp%ldn,kpt%wfn(isp,ikp)%nmem))
        kpt%wfn(isp,ikp)%zwf = zzero
     enddo
  enddo

  allocate(zwfn(nwedge))
  do ikp = 1, kpt%nk
     if (peinf%master) then
        allocate(map_found(max(kpt%wfn(1,ikp)%nmem,kpt%wfn(nspin,ikp)%nmem),nspin))
        map_found = 0
        read(infile) ndim
        allocate(gtable(3,ndim))
        do ii = 1, ndim
           read(infile) (gtable(jj,ii),jj=1,3)
        enddo
        allocate(zcg(ndim*nspin))
     endif
     itmp = 0
     do ii = 1, kpt%wfn(1,ikp)%nstate
        if (kpt%wfn(1,ikp)%map(ii) == 0 .or. kpt%wfn(nspin,ikp)%map(ii) == 0) then
           if (peinf%master) read(infile)
           cycle
        endif
        if (peinf%master) read(infile) (zcg(jj),jj=1,ndim*nspin)
        do isp = 1, nspin

           if (peinf%master) then
              zwfn = zzero
              fft_box%zbox = zzero
              do jj = 1, ndim
                 j3 = gtable(3,jj)
                 i3 = j3 + 1
                 if (j3 < 0) i3 = fft_box%nfft(3) + i3
                 i3 = i3 - fft_box%nfft(3)/2-1
                 j2 = gtable(2,jj)
                 i2 = j2 + 1
                 if (j2 < 0) i2 = fft_box%nfft(2) + i2
                 i2 = i2 - fft_box%nfft(2)/2 - 1
                 j1 = gtable(1,jj)
                 i1 = j1 + 1
                 if (j1 < 0) i1 = fft_box%nfft(1) + i1
                 i1 = i1 - fft_box%nfft(1)/2 - 1
                 fft_box%zbox(i1,i2,i3) = zcg(jj + ndim*(isp-1))
              enddo
              call zdo_FFT(1,fft_box)
              jj = 0
              do j3 = -fft_box%nfft(3)/2, fft_box%nfft(3)/2 - 1
                 do j2 = -fft_box%nfft(2)/2, fft_box%nfft(2)/2 - 1
                    do j1 = -fft_box%nfft(1)/2, fft_box%nfft(1)/2 - 1
                       jj = jj + 1
                       zwfn(invindex(jj)) = fft_box%zbox(j1,j2,j3) * sqrt(fft_box%scale)
                    enddo
                 enddo
              enddo
           endif

           if (kpt%wfn(isp,ikp)%map(ii) == 0) cycle
           if (peinf%master) then
              map_found(kpt%wfn(isp,ikp)%map(ii),isp) = ii
              norm = real(dot_product(zwfn,zwfn),dp)
#ifdef DEBUG
              write(6,'(a,3i8,g14.5)') ' GWC ', isp, ii, gvec%nr, norm
#endif
              if (abs(norm-one) > 1.d-6) write(6,*) &
                   ' WARNING: The wavefuctions are not normalized ', &
                   ii, gvec%nr, norm, isp
           endif

           itmp(isp) = itmp(isp) + 1
#ifdef MPI
           call MPI_BCAST(zwfn,nwedge, &
                MPI_DOUBLE_COMPLEX,peinf%masterid,peinf%comm,info)
#endif
           call zcopy(w_grp%mydim,zwfn(w_grp%offset+1),1,kpt%wfn(isp,ikp)%zwf(1,itmp(isp)),1)
        enddo
     enddo
     if (peinf%master) then
        do isp = 1, nspin
           write(6,'(/,a,i10,2(a,i4))') ' Read wavefunctions for ', &
                itmp(isp), ' levels, spin ', isp, ' k-point ', ikp
           do ii = 1, kpt%wfn(isp,ikp)%nstate
              if (kpt%wfn(isp,ikp)%map(ii) /= 0) then
                 if (map_found(kpt%wfn(isp,ikp)%map(ii),isp) == 0) then
                    write(lastwords,*) ' ERROR! wavefunction for state ', &
                         ii, ' not found.'
                    call die(lastwords)
                 endif
              endif
           enddo
        enddo
        deallocate(zcg)
        deallocate(gtable)
        deallocate(map_found)
     endif
  enddo
  deallocate(zwfn)

  call zfinalize_FFT(peinf%inode,fft_box)

  deallocate(rindex)
  deallocate(invindex)
  if (peinf%master) close(infile)
  if (peinf%master) write(6,'(/,a,/)') repeat('-',65)

  !-------------------------------------------------------------------
  ! Synchronization.
  !
#ifdef MPI
  call MPI_BARRIER(peinf%comm,info)
#endif

end subroutine paratec_wfn
!===================================================================
