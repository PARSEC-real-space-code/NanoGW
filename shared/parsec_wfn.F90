!===================================================================
!
! Read input wavefunctions/parameters from file parsec.dat
! first lines in input file has creation date.
!
! Only Master PE reads input, and broadcasts all info to the other PEs.
!
! Besides reading input, this subroutine performs a series of tasks:
!   constructs real-space and reciprocal-space grids;
!   sets up symmetry group and calculates representation
!      (inside subroutines);
!   reads electron density from parsec.dat and constructs the DFT 
!      exchange-correlation potential and kernel (derivative of potential);
!   defines the FFT grid.
!
! Copyright (C) 2009 Murilo L. Tiago, http://users.ices.utexas.edu/~mtiago
! This file is part of RGWBS. It is distributed under the GPL v1.
!
!-------------------------------------------------------------------
subroutine parsec_wfn(gvec,kpt,nmap,nspin,wmap,init_gr)

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
  logical :: readwfn
  character (len=26) :: datelabel
  character (len=800) :: lastwords
  integer :: ii, jj, nn, isp, ikp, icount, info, ndim, nwedge, itmp(3)
  integer :: isave
  real(dp):: norm, rtmp, vtmp(3), mtmp(3,3)

  integer, allocatable :: rindex(:), gtable(:,:)
  integer, allocatable :: iord(:), map_found(:)
  real(dp), allocatable :: rho_in(:,:)
  real(dp), allocatable :: dcg(:), dwfn(:)
  complex(dpc), allocatable :: zcg(:), zwfn(:)

  !-------------------------------------------------------------------
  ! Read info for real space grid from parsec.dat.
  !
  if(peinf%master) then
     write(6,'(/,a,/,a,a,/)') ' Electronic structure data :', &
          ' ', repeat('-',27)
     write(6,*) 'Reading electron wavefunctions from file parsec.dat/wfn.dat'
     open(infile,file='parsec.dat',form='unformatted',status='old',iostat=info)
     readwfn = .false.
     if (info /= 0) then
        write(6,*) 'WARNING: parsec.dat file not found'
        write(6,*) '         Looking for file wfn.dat'
        call stopwatch(peinf%master,'parsec_wfn:  reading wfn.dat')
        open(infile,file='wfn.dat',form='unformatted',status='old',iostat=info)
        readwfn = .true.
        if (info /= 0) then
           write(lastwords,*) 'ERROR: wfn.dat file not found'
           call die(lastwords)
        endif
     endif
     read(infile) datelabel
     write(6,'(/,a,a,a)') ' Wavefunction file created on ', datelabel, ' UTC'
     if (readwfn) then
        ! Assume non-periodic system, real wave-functions.
        gvec%per = 0
        kpt%lcplx = .false.
        read(infile) rtmp, vtmp(1), ndim, nwedge, gvec%syms%ntrans, nspin
        gvec%step = rtmp
        gvec%avec = zero
        do ii = 1, 3
           gvec%avec(ii,ii) = one
        enddo
        fd%lap_dir_num = 3
        fd%lap_dir(1:3) = 1
        fd%lap_dir(4:6) = 0
        fd%lap_neig = 0
        fd%lap_neig(1,1) = 1
        fd%lap_neig(2,2) = 1
        fd%lap_neig(3,3) = 1
        fd%lap_dir_step = gvec%step(1)
        fd%b_lap(1:3) = one
        fd%b_lap(4:6) = zero
        kpt%nk = 1
        allocate(kpt%fk(3,kpt%nk))
        kpt%fk = zero
        allocate(kpt%weight(kpt%nk))
        kpt%weight = one
     else
        read(infile) nspin,ii,gvec%per
        if (ii == 0) kpt%lcplx = .false.
        if (ii == 1) kpt%lcplx = .true.
        !
        ! Must modify the value of gvec%per: values at input are:
        !    0  ->  cluster, no periodicity
        !    2  ->  wire, periodicity along x direction only
        !    3  ->  slab, periodicity along x and y directions
        !    1  ->  bulk, periodicity along all cartesian directions
        !
        ii = gvec%per
        select case(ii)
        case(1)
           gvec%per = 3
        case(2)
           gvec%per = 1
        case(3)
           gvec%per = 2
        end select
        if (kpt%lcplx .and. peinf%master) write(6,'(/,a,/,a,/,a,/)') &
             repeat('-*-',23), ' Wavefunctions are complex', repeat('-*-',23)
        if (gvec%per < 0 .or. gvec%per > 3) then
           write(6,*) 'WARNING: boundary conditions not consistent '
           write(6,*) gvec%per
           call die(' STOP. ')
        endif
        if (gvec%per > 0) then
           read(infile) gvec%step
           read(infile) gvec%avec
           read(infile) fd%lap_dir_num
           read(infile) (fd%lap_dir(ii),ii=1,3) !, &
           !     ((fd%lap_neig(ii,jj),ii=1,3),jj=1,3)  ! W.Gao: lap_neig is not written by parsec1.4
           read(infile) (fd%lap_dir_step(ii),ii=1,3) !, fd%b_lap  ! W.Gao: b_lap is not written by parsec1.4
           ii = 3
           do jj = 1, 3
              if (fd%lap_dir(jj) /= 0) then
                 ii = ii + 1
                 fd%lap_neig(ii,:) = fd%lap_neig(jj,:)
                 fd%lap_dir_step(ii) = fd%lap_dir_step(jj)
                 fd%b_lap(ii) = fd%b_lap(jj+3)
              endif
           enddo
           fd%lap_dir_num = fd%lap_dir_num + 3
           if (ii /= fd%lap_dir_num) then
                write(lastwords,*) 'ERROR! Laplacian directions incorrect ', &
                     ii, fd%lap_dir_num, ' Stop.'
                call die(lastwords)
           endif
           fd%lap_dir(fd%lap_dir_num+1:6) = 0
           fd%lap_neig(fd%lap_dir_num+1:6,:) = 0
           fd%lap_dir_step(fd%lap_dir_num+1:6) = zero
           fd%b_lap(fd%lap_dir_num+1:6) = zero
           fd%lap_dir(1:3) = 1
           fd%lap_neig(1:3,:) = 0
           fd%lap_neig(1,1) = 1
           fd%lap_neig(2,2) = 1
           fd%lap_neig(3,3) = 1
           fd%lap_dir_step(1:3) = gvec%step(1:3)
           read(infile) kpt%nk
           read(infile)
           read(infile)
           read(infile)
           allocate(kpt%fk(3,kpt%nk))
           read(infile) kpt%fk
           allocate(kpt%weight(kpt%nk))
           read(infile) kpt%weight
        else
           read(infile) rtmp
           gvec%step = rtmp
           gvec%avec = zero
           do ii = 1, 3
              gvec%avec(ii,ii) = one
           enddo
           fd%lap_dir_num = 3
           fd%lap_dir(1:3) = 1
           fd%lap_dir(4:6) = 0
           fd%lap_neig = 0
           fd%lap_neig(1,1) = 1
           fd%lap_neig(2,2) = 1
           fd%lap_neig(3,3) = 1
           fd%lap_dir_step = gvec%step(1)
           fd%b_lap(1:3) = one
           fd%b_lap(4:6) = zero
           kpt%nk = 1
           allocate(kpt%fk(3,kpt%nk))
           kpt%fk = zero
           allocate(kpt%weight(kpt%nk))
           kpt%weight = one
        endif
        read(infile) ndim
        read(infile) gvec%shift
        read(infile) nwedge, gvec%syms%ntrans
     endif
#ifdef DEBUG
     write(6,*) "Finished reading header."
#endif
  endif
#ifdef MPI
  call MPI_BCAST(readwfn,1, &
       MPI_LOGICAL,peinf%masterid,peinf%comm,info)
  call MPI_BCAST(gvec%per,1, &
       MPI_INTEGER,peinf%masterid,peinf%comm,info)
  call MPI_BCAST(gvec%step,6, &
       MPI_DOUBLE_PRECISION,peinf%masterid,peinf%comm,info)
  call MPI_BCAST(gvec%avec,9, &
       MPI_DOUBLE_PRECISION,peinf%masterid,peinf%comm,info)
  call MPI_BCAST(fd%lap_dir_num,1, &
       MPI_INTEGER,peinf%masterid,peinf%comm,info)
  call MPI_BCAST(fd%lap_dir,6, &
       MPI_INTEGER,peinf%masterid,peinf%comm,info)
  call MPI_BCAST(fd%lap_neig,18, &
       MPI_INTEGER,peinf%masterid,peinf%comm,info)
  call MPI_BCAST(fd%lap_dir_step,6, &
       MPI_DOUBLE_PRECISION,peinf%masterid,peinf%comm,info)
  call MPI_BCAST(fd%b_lap,6, &
       MPI_DOUBLE_PRECISION,peinf%masterid,peinf%comm,info)
  call MPI_BCAST(kpt%nk, 1, &
       MPI_INTEGER,peinf%masterid,MPI_COMM_WORLD,info)
  if (.not. peinf%master) allocate(kpt%fk(3,kpt%nk))
  call MPI_BCAST(kpt%fk, kpt%nk*3, &
       MPI_DOUBLE_PRECISION,peinf%masterid,MPI_COMM_WORLD,info)
  if (.not. peinf%master) allocate(kpt%weight(kpt%nk))
  call MPI_BCAST(kpt%weight, kpt%nk, &
       MPI_DOUBLE_PRECISION,peinf%masterid,MPI_COMM_WORLD,info)
  call MPI_BCAST(gvec%syms%ntrans,1, &
       MPI_INTEGER,peinf%masterid,peinf%comm,info)
  call MPI_BCAST(nspin,1, &
       MPI_INTEGER,peinf%masterid,peinf%comm,info)
  call MPI_BCAST(kpt%lcplx,1, &
       MPI_LOGICAL,peinf%masterid,peinf%comm,info)
#endif
  !
  ! K-vectors are writen in Cartesian components. Must convert them to units
  ! of reciprocal lattice vectors.
  !
  kpt%fk = matmul(transpose(gvec%avec),kpt%fk) / (two *pi)
  !
  ! Read info about symmetry operations.
  !
  allocate(gvec%syms%trans(3,3,gvec%syms%ntrans))
  allocate(gvec%syms%chi(gvec%syms%ntrans,gvec%syms%ntrans))
  if (peinf%master) then
     if (readwfn) read(infile) gvec%shift
     read(infile)
     read(infile) gvec%syms%trans
     read(infile)
     read(infile)
     read(infile) gvec%syms%chi
  endif
#ifdef MPI
  call MPI_BCAST(gvec%shift,3, &
       MPI_DOUBLE_PRECISION,peinf%masterid,peinf%comm,info)
  call MPI_BCAST(gvec%syms%trans,9*gvec%syms%ntrans, &
       MPI_DOUBLE_PRECISION,peinf%masterid,peinf%comm,info)
  call MPI_BCAST(gvec%syms%chi,gvec%syms%ntrans*gvec%syms%ntrans, &
       MPI_INTEGER,peinf%masterid,peinf%comm,info)
#endif

  !-------------------------------------------------------------------
  ! Read real space grid and determine the maximum distance from a grid
  ! point to the origin of Cartesian system. This distance is relevant
  ! only in non-periodic systems. It will be used to define the truncation
  ! length in the Coulomb potential.
  !
  if(peinf%master) then
     allocate(gtable(3,nwedge))
     read(infile) ((gtable(ii,jj),ii=1,3),jj=1,nwedge)
     gvec%rmax = zero
     vtmp = zero
     do nn = 1, nwedge
        do jj = 1, gvec%syms%ntrans
           call unfold(gtable(1,nn),gvec%syms%trans(1,1,jj),gvec%shift(1),itmp(1))
           do ii = gvec%per + 1, 3
              vtmp(ii) = (itmp(ii) + gvec%shift(ii))*gvec%step(ii)
           enddo
           rtmp = dot_product(vtmp,vtmp)
           if ( gvec%rmax < rtmp ) gvec%rmax = rtmp
        enddo
     enddo
     gvec%rmax = sqrt(gvec%rmax)
  endif
#ifdef MPI
  call MPI_BCAST(gvec%rmax,1, &
       MPI_DOUBLE_PRECISION,peinf%masterid,peinf%comm,info)
  call MPI_BCAST(nwedge,1, &
       MPI_INTEGER,peinf%masterid,peinf%comm,info)
  if(.not.peinf%master) allocate(gtable(3,nwedge))
  call MPI_BCAST(gtable,3*nwedge, &
       MPI_INTEGER,peinf%masterid,peinf%comm,info)
#endif

  !-------------------------------------------------------------------
  ! Read DFT eigenvalues and related data.
  ! Must be careful with wfn.dat file because the electron density is
  ! written in different place. After it is done, add charge density
  ! for this component into rho_in.
  !
  allocate(kpt%wfn(nspin,kpt%nk))
  allocate(rho_in(nwedge,3))
  rho_in = zero

  if (peinf%master) then
     isp = 1
     do ikp = 1, kpt%nk
        read(infile) kpt%wfn(isp,ikp)%nstate
        allocate(kpt%wfn(isp,ikp)%irep(kpt%wfn(isp,ikp)%nstate))
        read(infile) (kpt%wfn(isp,ikp)%irep(ii),ii=1,kpt%wfn(isp,ikp)%nstate)
        allocate(kpt%wfn(isp,ikp)%e0(kpt%wfn(isp,ikp)%nstate))
        read(infile) (kpt%wfn(isp,ikp)%e0(ii),ii=1,kpt%wfn(isp,ikp)%nstate)
        allocate(kpt%wfn(isp,ikp)%e1(kpt%wfn(isp,ikp)%nstate))
        call read_scissors(peinf%master,isp,kpt%wfn(isp,ikp)%nstate, &
             kpt%wfn(isp,ikp)%e0,kpt%wfn(isp,ikp)%e1)
        allocate(kpt%wfn(isp,ikp)%occ0(kpt%wfn(isp,ikp)%nstate))
        read(infile) (kpt%wfn(isp,ikp)%occ0(ii),ii=1,kpt%wfn(isp,ikp)%nstate)
     enddo
     read(infile)
     read(infile) (rho_in(ii,2),ii=1,nwedge)

     if (nspin == 2) then
        if (readwfn) then
           do ii = 1, kpt%wfn(isp,1)%nstate
              read(infile)
           enddo
        endif
        isp = 2
        do ikp = 1, kpt%nk
           read(infile) kpt%wfn(isp,ikp)%nstate
           allocate(kpt%wfn(isp,ikp)%irep(kpt%wfn(isp,ikp)%nstate))
           read(infile) (kpt%wfn(isp,ikp)%irep(ii),ii=1,kpt%wfn(isp,ikp)%nstate)
           allocate(kpt%wfn(isp,ikp)%e0(kpt%wfn(isp,ikp)%nstate))
           read(infile) (kpt%wfn(isp,ikp)%e0(ii),ii=1,kpt%wfn(isp,ikp)%nstate)
           allocate(kpt%wfn(isp,ikp)%e1(kpt%wfn(isp,ikp)%nstate))
           call read_scissors(peinf%master,isp,kpt%wfn(isp,ikp)%nstate, &
                kpt%wfn(isp,ikp)%e0,kpt%wfn(isp,ikp)%e1)
           allocate(kpt%wfn(isp,ikp)%occ0(kpt%wfn(isp,ikp)%nstate))
           read(infile) (kpt%wfn(isp,ikp)%occ0(ii),ii=1,kpt%wfn(isp,ikp)%nstate)
        enddo
        read(infile)
        read(infile) (rho_in(ii,3),ii=1,nwedge)
        if (readwfn) then
           rewind(infile)
           do ii = 1, 15
              read(infile)
           enddo
        endif
     endif
     rho_in(:,1) = rho_in(:,2) + rho_in(:,3)
#ifdef DEBUG
     write(6,*) 'Finished reading header and charge-density in parsec.dat'
#endif

  endif ! if (peinf%master)

#ifdef MPI
  do ikp = 1, kpt%nk
     do isp = 1, nspin
        call MPI_BCAST(kpt%wfn(isp,ikp)%nstate,1, &
             MPI_INTEGER,peinf%masterid,peinf%comm,info)
        if (.not. peinf%master) then
           allocate(kpt%wfn(isp,ikp)%irep(kpt%wfn(isp,ikp)%nstate))
           allocate(kpt%wfn(isp,ikp)%e0(kpt%wfn(isp,ikp)%nstate))
           allocate(kpt%wfn(isp,ikp)%e1(kpt%wfn(isp,ikp)%nstate))
           allocate(kpt%wfn(isp,ikp)%occ0(kpt%wfn(isp,ikp)%nstate))
        endif
        call MPI_BCAST(kpt%wfn(isp,ikp)%irep,kpt%wfn(isp,ikp)%nstate, &
             MPI_INTEGER,peinf%masterid,peinf%comm,info)
        call MPI_BCAST(kpt%wfn(isp,ikp)%e0,kpt%wfn(isp,ikp)%nstate, &
             MPI_DOUBLE_PRECISION,peinf%masterid,peinf%comm,info)
        call MPI_BCAST(kpt%wfn(isp,ikp)%e1,kpt%wfn(isp,ikp)%nstate, &
             MPI_DOUBLE_PRECISION,peinf%masterid,peinf%comm,info)
        call MPI_BCAST(kpt%wfn(isp,ikp)%occ0,kpt%wfn(isp,ikp)%nstate, &
             MPI_DOUBLE_PRECISION,peinf%masterid,peinf%comm,info)
     enddo
  enddo
  call MPI_BCAST(rho_in,nwedge*3, &
       MPI_DOUBLE_PRECISION,peinf%masterid,peinf%comm,info)
#endif

  !-------------------------------------------------------------------
  ! Define FFT grid and supercell box. Nfft must be even number!
  ! nfft(4) is used only in real transforms.
  !
  do ii = 1, 3
     if (readwfn) then
        vtmp(ii) = four*gvec%rmax
     else
        if (gvec%per >= ii) then
           vtmp(ii) = sqrt(dot_product(gvec%avec(:,ii),gvec%avec(:,ii)))
        else
           vtmp(ii) = four*gvec%rmax
        endif
     endif
     fft_box%nfft(ii) = nint( vtmp(ii)/gvec%step(ii) )
  enddo
  call adjust_fft(fft_box%nfft)
  do ii = 1, 3
     if (mod(fft_box%nfft(ii),2) /= 0) fft_box%nfft(ii) = fft_box%nfft(ii) + 1
  enddo
  fft_box%nfft(4) = fft_box%nfft(1) + 2

  !---------------------------------------------------------------
  ! Calculate volume, hcub, reciprocal space metric.
  !
  if (gvec%per > 0) then
     mtmp = gvec%avec
     do ii = 1, gvec%per
        mtmp(:,ii) = gvec%step(ii)*mtmp(:,ii)/sqrt(sum(mtmp(:,ii)**2))
     enddo
     call mtrxin(mtmp,gvec%hcub,info)
     if (info /= 0) call die('ERROR in mtrxin. Stop.')
     gvec%hcub = abs(gvec%hcub)
     gvec%bvec = gvec%avec
     call mtrxin(gvec%bvec,gvec%celvol,info)
     if (info /= 0) call die('ERROR in mtrxin. Stop.')
     gvec%celvol = abs(gvec%celvol)
     gvec%bvec = two * pi * gvec%bvec
     gvec%bdot = matmul(gvec%bvec,transpose(gvec%bvec))
     do ii = gvec%per + 1, 3
        gvec%hcub = gvec%hcub * gvec%step(ii)
        gvec%bdot(ii,ii) = one
     enddo
  else
     gvec%hcub = gvec%step(1)**3
     gvec%celvol = one
     gvec%bvec = zero
     gvec%bdot = zero
     do ii = 1, 3
        rtmp = fft_box%nfft(ii)*gvec%step(ii)
        gvec%celvol = gvec%celvol*rtmp
        gvec%bvec(ii,ii) = two * pi / rtmp
        gvec%bdot(ii,ii) = gvec%bvec(ii,ii)**2
     enddo
  endif

  !-------------------------------------------------------------------
  ! Truncate the Coulomb potential if the system is non-periodic or
  ! partially periodic.
  !
  if (gvec%per < 3) then
     fft_box%trunc = two*gvec%rmax
  else
     fft_box%trunc = -one
  endif

  !-------------------------------------------------------------------
  ! Define the long-wavelength part of Coulomb potential.
  !
  gvec%long = vq_average(gvec%per,kpt%nk,gvec%celvol,gvec%rmax)

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
           if (gvec%per > 0) write(6,'(a,3f13.6)') ' K-point = ', kpt%fk(:,ikp)
           write(6,*) ' Number of occupied levels, spin ', isp, &
                ' : ', sum(kpt%wfn(isp,ikp)%occ0)
           write(6,*) 'DFT Eigenvalues (eV), spin ', isp, ' :'
           write(6,'(10f7.2)') &
                (kpt%wfn(isp,ikp)%e0(ii)*ryd,ii=1,kpt%wfn(isp,ikp)%nstate)
           write(6,*) 'Irreducible representations, spin ', isp, ' :'
           write(6,'(10i5)') &
                (kpt%wfn(isp,ikp)%irep(ii),ii=1,kpt%wfn(isp,ikp)%nstate)
           write(6,*)
        enddo
     enddo
     write(6,'(2a,e16.8)') ' Maximum value of electron density ', &
          '(a.u.^-3) = ', maxval(rho_in(:,1))
     write(6,'(a,f10.4)') ' Total number of electrons = ', &
          sum(rho_in(:,1)) * gvec%hcub * real(gvec%syms%ntrans,dp)
  endif

  !-------------------------------------------------------------------
  ! Reorder grid points in descending value of electron density.
  ! gtable is the full, original grid read from parsec.dat/wfn.dat.
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
  do jj = 1, gvec%nr
     ii = rindex(jj)
     gvec%r(:,jj) = gtable(:,ii)
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
  call parsec_xc(nspin,ii,jj)
  call define_xc(nspin,ii,jj,init_gr)
  if (gvec%per > 0) init_gr = .true.

  !-------------------------------------------------------------------
  ! Set distributed grid.
  !
  fd%norder = 0
  if (init_gr) call parsec_atoms(gvec,peinf%master,fd%norder,info)
  if (info /= 0) call die(' ')
  if (init_gr) then
     allocate(gtable(3,gvec%nr*gvec%syms%ntrans),stat=info)
     call alccheck('gtable','parsec_wfn',3*gvec%nr*gvec%syms%ntrans,info)
     call grid_setup(gvec,gtable)
     if (fd%norder > 0) call fd_setup(gvec,gtable)
     deallocate(gtable)
  endif

  !-------------------------------------------------------------------
  ! More printouts.
  !
  if (peinf%master) then
     write(6,'(/,a,/,a,/)') ' Grid data :',' -----------'
     if (fd%norder > 0) then
        write(6,'(2a,i3,/)') ' Gradients and Laplacians calculated ', &
             'with finite differences, order ', fd%norder
        if (abs(sum(fd%b_lap) - three) > 1.d-6) &
             write(6,'(2a,g20.10)') 'WARNING! Trace of laplacian ', &
             'directions is not exactly 3 ! ', sum(fd%b_lap)
        write(6,'(a,i2,/,a)') ' Number of finite difference directions = ', &
             fd%lap_dir_num, ' Finite difference directions:'
        write(6,'(6(3i6,3x,a,f9.6,/))') ((fd%lap_neig(ii,jj),jj=1,3), &
             '  weight = ', fd%b_lap(ii), ii = 1, fd%lap_dir_num)
     else
        write(6,'(a,/)') ' Gradients and Laplacians calculated with FFT'
     endif
     if (gvec%per > 0) then
        select case(gvec%per)
        case (1)
           write(6,'(a,//,a,f16.6)') '  WIRE BOUNDARY CONDITIONS !', &
                ' Length of periodic cell (a.u.) = ', gvec%celvol
        case (2)
           write(6,'(a,//,a,f16.6)') '  SLAB BOUNDARY CONDITIONS !', &
                ' Area of periodic cell (a.u.^2) = ', gvec%celvol
        case (3)
           write(6,'(a,//,a,f16.6)') '  BULK BOUNDARY CONDITIONS !', &
                ' Volume of periodic cell (a.u.^3) = ', gvec%celvol
        end select
        write(6,'(/,a)') ' Unit lattice vectors (a.u.) = '
        do jj = 1, gvec%per
           write(6,'(3f20.15)') (gvec%avec(ii,jj),ii=1,gvec%per)
        enddo
     else
        write(6,'(a,/)') '  ELECTRONIC SYSTEM IS CONFINED !'
     endif
!#ifdef DEBUG
     write(6,'(/,a,3(/,3f20.10))') ' Lattice vectors (a.u.): ', gvec%avec
     write(6,'(a,3(/,3f20.10))') ' Normalized lattice vectors: ', gvec%avec_norm
     write(6,'(a,3(/,3f20.10))') ' Reciprocal lattice vectors (a.u.^-1): ', gvec%bvec
     write(6,'(a,3(/,3f20.10))') ' Reciprocal space metric (a.u.^-2): ', gvec%bdot
     write(6,'(a,g20.8)') ' Numeric volume infinitesimal (a.u.^3): ', gvec%hcub
!#endif
     if (gvec%per < 3) write(6,'(/,a,f16.6)') &
          ' Confining radius (a.u.) = ', gvec%rmax
     write(6,'(/,a,3f16.6)') ' Grid spacing (a.u.) = ', gvec%step
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
        write(6,*) " in parsec_wfn: nmap = ", nmap
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
                      ' states but parsec.dat/wfn.dat contains only ', &
                      kpt%wfn(isp,ikp)%nstate, ' Stop.'
                 call die(lastwords)
              endif
           enddo
        endif
        kpt%wfn(isp,ikp)%nmem = icount
        write(6,*) " in parsec_wfn: nmem = ", kpt%wfn(isp,ikp)%nmem
        allocate(kpt%wfn(isp,ikp)%imap(kpt%wfn(isp,ikp)%nmem))
        kpt%wfn(isp,ikp)%imap = 0
        do ii = 1, kpt%wfn(isp,ikp)%nstate
           if (kpt%wfn(isp,ikp)%map(ii) /= 0) &
                kpt%wfn(isp,ikp)%imap(kpt%wfn(isp,ikp)%map(ii)) = ii
        enddo

        ! W Gao dbg
        if(peinf%master .and. .true.) then
          write(*,'("isp = ",i2," ikp = ",i4)') isp, ikp
          write(*,'("    i   wfn%map(ii) ")')
          do ii = 1, kpt%wfn(isp,ikp)%nstate
            write(*,'(2i8)') ii, kpt%wfn(isp,ikp)%map(ii)
          enddo
          write(*,'("    i   wfn%imap(ii) ")')
          do ii = 1, kpt%wfn(isp,ikp)%nmem
            write(*,'(2i8)') ii, kpt%wfn(isp,ikp)%imap(ii)
          enddo
        endif
        ! dbg
     enddo ! ikp
  enddo ! isp
   
  !-------------------------------------------------------------------
  ! Save electron density.
  !
  allocate(kpt%rho(gvec%nr,nspin),stat=info)
  call alccheck('kpt%rho','parsec_wfn',gvec%nr*nspin,info)
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
  call stopwatch(peinf%master,'parsec_wfn:  reading wave functions')
  w_grp%nr = gvec%nr
  w_grp%ldn = w_grp%nr / w_grp%npes + 1
  w_grp%mydim = w_grp%nr / w_grp%npes
  nn = mod(w_grp%nr, w_grp%npes)
  w_grp%offset = w_grp%ldn*w_grp%inode
  if (w_grp%inode < nn) w_grp%mydim = w_grp%mydim + 1
  if (w_grp%inode > nn) w_grp%offset = w_grp%offset - w_grp%inode + nn
  do isp = 1, nspin
     do ikp = 1, kpt%nk
        if (kpt%lcplx) then
           allocate(kpt%wfn(isp,ikp)%zwf(w_grp%ldn,kpt%wfn(isp,ikp)%nmem))
           kpt%wfn(isp,ikp)%zwf = zzero
        else
           allocate(kpt%wfn(isp,ikp)%dwf(w_grp%ldn,kpt%wfn(isp,ikp)%nmem))
           kpt%wfn(isp,ikp)%dwf = zero
        endif
     enddo
  enddo

  if (kpt%lcplx) then
     allocate(zcg(gvec%nr))
     allocate(zwfn(gvec%nr))
  else
     allocate(dcg(gvec%nr))
     allocate(dwfn(gvec%nr))
  endif
#ifdef DEBUG
  if(peinf%master) then
     write(6,*) " Start reading wavefunctions from parsec.dat. "
  endif
#endif 
  do isp = 1, nspin
     do ikp = 1, kpt%nk
        if (peinf%master) then
#ifdef DEBUG
           write(6,'(A,i2,A,i4)') " Reading wfn with spin ",isp," at kpt ", ikp 
#endif
           allocate(map_found(kpt%wfn(isp,ikp)%nmem))
           map_found = 0
           if (readwfn) then
              if (isp == 2) then
                 do ii = 1, 6
                    read(infile)
                 enddo
              endif
              nn = kpt%wfn(isp,ikp)%nstate
              allocate(iord(nn))
              do ii = 1, nn
                 iord(ii) = ii
              enddo
           else
              do ii = 1, 7
                 read(infile)
              enddo
              read(infile) nn
              if (nn == 0) then
                 write(lastwords,*) &
                      ' ERROR: parsec.dat contains no wavefunctions ! ', jj
                 call die(lastwords)
              endif
              allocate(iord(nn))
              read(infile) iord
           endif
        endif
#ifdef MPI
        call MPI_BCAST(nn,1,MPI_INTEGER,peinf%masterid,peinf%comm,info)
#endif
        icount = 0
        do ii = 1, nn
           if (peinf%master) then
              if (kpt%lcplx) then
                 read(infile) (zcg(jj),jj=1,nwedge)
                 do jj = 1, gvec%nr
                    zwfn(jj) = zcg(rindex(jj))
                 enddo
              else
                 read(infile) (dcg(jj),jj=1,nwedge)
                 do jj = 1, gvec%nr
                    dwfn(jj) = dcg(rindex(jj))
                 enddo
              endif
              isave = kpt%wfn(isp,ikp)%map(iord(ii))
           endif
#ifdef MPI
           call MPI_BCAST(isave,1,MPI_INTEGER,peinf%masterid,peinf%comm,info)
#endif
           ! If we do not need this state, skip it.
           if (isave == 0) cycle
           ! If we do, mark it as found.
           if (peinf%master) then
              map_found(kpt%wfn(isp,ikp)%map(iord(ii))) = iord(ii)
              !  Check normalization.
              if (kpt%lcplx) then
                 norm = real(dot_product(zwfn,zwfn),dp) * real(gvec%syms%ntrans,dp)
              else
                 norm = dot_product(dwfn,dwfn) * real(gvec%syms%ntrans,dp)
              endif
#ifdef DEBUG
              write(6,'(a,3i8,g14.5)') ' wfn ', &
                   isp, iord(ii), gvec%nr, norm
#endif
              if (abs(norm-one) > 1.d-6) write(6,*) &
                   ' WARNING: The wavefuctions are not normalized ', &
                   iord(ii), gvec%nr, norm, isp
           endif
           icount = icount + 1
           if (kpt%lcplx) then
#ifdef MPI
              call MPI_BCAST(zwfn,gvec%nr, &
                   MPI_DOUBLE_COMPLEX,peinf%masterid,peinf%comm,info)
#endif
              call zcopy(w_grp%mydim,zwfn(w_grp%offset+1),1,kpt%wfn(isp,ikp)%zwf(1,icount),1)
           else
#ifdef MPI
              call MPI_BCAST(dwfn,gvec%nr, &
                   MPI_DOUBLE_PRECISION,peinf%masterid,peinf%comm,info)
#endif
              call dcopy(w_grp%mydim,dwfn(w_grp%offset+1),1,kpt%wfn(isp,ikp)%dwf(1,icount),1)
           endif
        enddo                     !ii (loop on states)
        if (peinf%master) then
           if (kpt%lcplx) then
              write(6,'(/,a,i10,2(a,i4))') ' Read wavefunctions for ', &
                   icount, ' levels, spin ', isp, ' k-point ', ikp
           else
              write(6,'(/,a,i10,a,i4)') ' Read wavefunctions for ', &
                   icount, ' levels, spin ', isp
           endif
           do ii = 1, kpt%wfn(isp,ikp)%nstate
              if (kpt%wfn(isp,ikp)%map(ii) /= 0) then
                 if (map_found(kpt%wfn(isp,ikp)%map(ii)) == 0) then
                    write(6,*) ' ERROR! wavefunction for state ', ii, ' not found.'
                    call die(' Stop.')
                 endif
              endif
           enddo
           deallocate(map_found)
           deallocate(iord)
        endif
     enddo
  enddo
  if (kpt%lcplx) then
     deallocate(zcg, zwfn)
  else
     deallocate(dcg, dwfn)
  endif
  deallocate(rindex)
  if (peinf%master) close(infile)
  if (peinf%master) write(6,'(/,a,/)') repeat('-',65)

  !-------------------------------------------------------------------
  ! Synchronization.
  !
#ifdef MPI
  call MPI_BARRIER(peinf%comm,info)
#endif

end subroutine parsec_wfn
!===================================================================
