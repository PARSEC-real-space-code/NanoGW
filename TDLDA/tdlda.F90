!===================================================================
!
! Diagonalizes the TDLDA eigenvalue equations.
!
! Input file is "tdlda.in" (optional). Output
! is writen on screen and on files "eigenvalues_rpa" (RPA-LDA energy
! eigenvalues and oscillator strengths), "eigenvalues_lda" (the same
! for TDLDA) and "pol_eig.dat" (eigenvectors).
!
! This code follows the numerical methodology presented in:
!    M.E. Casida, in "Recent Advances in Density-Functional Methods", ed. D.P. Chong (1995)
!    I. Vasiliev, S. Ogut and J.R. Chelikowsky, Phys. Rev. B 65, 115416 (2002)
!
! Copyright (C) 2009 Murilo Tiago, Univ. of Texas, Austin, TX, USA
! mtiago@ices.utexas.edu
!
! First version written by Murilo Tiago, Univ. of Minnesota, April 2004
!
!    This program is free software; you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation; either version 1, or (at your option)
!    any later version.
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
program tdlda
#ifdef HIPMAGMA
  use magma
#endif
  use typedefs
  use mpi_module
  implicit none
#ifdef MPI
  include 'mpif.h'
#endif
  type(gspace) gvec
  type(kptinfo) :: kpt
  type(qptinfo) :: qpt
  type(polinfo), dimension(2) :: pol_in
  type(polinfo), allocatable :: pol(:, :)
  type(kernelinfo), allocatable :: k_p(:, :)
  type(options) :: opt

  character(len=800) :: lastwords
  character(len=40), allocatable :: routnam(:)
  logical :: nolda, tamm_d, rpaonly, trip_flag, noxchange, trunc_c, init_gr
  integer :: ii, jj, isp, irp, iq, nspin, nmap, nbuff, lcache, dft_code
  real(dp) :: tsec(2), mem1, xsum, xmax, tdldacut
  logical, allocatable :: wmap(:)
#ifdef MPI
  integer :: info
#endif

  ! W Gao dbg
  type(ISDF) :: isdf_in
  logical :: doisdf
  ! the number of val and cond states pairs. Since there could be two spins,
  ! ncv(1) correspond to spin up, and ncv(2) correspond to spin down.
  ! Assumption: for different kpts, the number of states with the same
  !  spin are the same
  integer :: n_intp, n_intp_r, n_dummy, maxicc, maxivv, iv, ic, &
             ivv, icc, maxnc, maxnv, maxncv, igrid, jgrid, &
             ihomo, ikp, intp_type, isdf_type, &
             kflag, vcirp, virp, cirp, ipt
  ! cvt.f90
  integer, allocatable :: intp(:), pairmap(:, :, :, :, :), invpairmap(:, :, :, :, :), icv(:), &
                          ncv(:, :, :), ivlist(:, :, :, :), iclist(:, :, :, :), nv(:, :, :), nc(:, :, :), timerlist(:), &
                          intp_r_tmp(:), intp_r(:), idum(:), sort_idx(:)
  real(dp), allocatable :: rho_at_intp_r(:)
  logical, allocatable  :: not_duplicate(:)

#ifdef DEBUG
  ! WG debug
  integer :: rho_intp_dbg
  character(len=20) :: dbg_filename
#endif
  logical :: verbose

  !-------------------------------------------------------------------
  ! Initialization.
  !
  call header('TDLDA')
  !-------------------------------------------------------------------
  ! Read input info.
  !
  if (peinf%master) write (*, *) " Call input_g"
  if (peinf%master) write (*, *) " dp = ", dp
  call input_g(pol_in, qpt, tdldacut, nbuff, lcache, w_grp%npes, &
               nolda, tamm_d, r_grp%num, dft_code, doisdf, n_intp, intp_type, isdf_type, &
               isdf_in%lessmemory, isdf_in%fastselect, opt%eigsolver, opt%linear_algebra, opt%pdgesv_nbl2d, .false.)
  if (peinf%master) write (*, *) " input_g done"
  call input_t(tamm_d, rpaonly, trip_flag, noxchange, trunc_c)
  if (opt%linear_algebra == 2 .or. opt%eigsolver == 2) then
#ifndef HIPMAGMA
    write (6, *) "Not compiled with '-DHIPMAGMA', use cpu for ", &
      " eigsolver and linear algebra !!"
    opt%linear_algebra = 1
    opt%eigsolver = 1
#endif
  end if

  !-------------------------------------------------------------------
  ! Determine the set of wavefunctions to read: if n-th wavefunction is
  ! needed, then wmap(n) = .true.; otherwise, wmap(n) = .false.
  !
  nmap = 0
  ! if any of pol_in(1)%ncond, .... is larger than 0, then
  if (max(pol_in(1)%ncond, pol_in(1)%nval, &
          pol_in(2)%ncond, pol_in(2)%nval) > 0) then
    ! nmap is equal to the max of pol_in(1)%cmap(:), ...
    if (pol_in(1)%ncond > 0) nmap = max(nmap, maxval(pol_in(1)%cmap))
    if (pol_in(1)%nval > 0) nmap = max(nmap, maxval(pol_in(1)%vmap))
    if (pol_in(2)%ncond > 0) nmap = max(nmap, maxval(pol_in(2)%cmap))
    if (pol_in(2)%nval > 0) nmap = max(nmap, maxval(pol_in(2)%vmap))
    allocate (wmap(nmap))
    wmap = .false. ! Initialization
    do isp = 1, 2
      do ii = 1, pol_in(isp)%ncond
        wmap(pol_in(isp)%cmap(ii)) = .true.
      end do
      do ii = 1, pol_in(isp)%nval
        wmap(pol_in(isp)%vmap(ii)) = .true.
      end do
    end do
    ! W Gao dbg
#ifdef DEBUG
    if (peinf%master) then
      write (*, '(a,i5)') " nmap = ", nmap
      write (*, '("   i    wmap(i)   ")')
      do ii = 1, nmap
        if (wmap(ii)) then
          write (*, '(i8," True")') ii
        else
          write (*, '(i8," False")') ii
        end if
      end do
    end if
#endif
  else
    allocate (wmap(1))
  end if
  if (min(pol_in(1)%ncond, pol_in(1)%nval, pol_in(2)%ncond, pol_in(2)%nval) < 0) then
    deallocate (wmap)
    allocate (wmap(1))
    nmap = 0
  end if

  !-------------------------------------------------------------------
  ! Read wave-function file.
  !
  init_gr = .false.
  if (dft_code == PARATEC) then
    call paratec_wfn(gvec, kpt, nmap, nspin, wmap, init_gr)
  else
    call parsec_wfn(gvec, kpt, nmap, nspin, wmap, init_gr)
  end if
  deallocate (wmap)
  if (trip_flag .and. nspin == 2) then
    write (lastwords, *) 'ERROR: cannot do TDLDA triplet kernel with ', &
      'spin polarized wavefunctions. Stop.'
    call die(lastwords)
  end if

  !-------------------------------------------------------------------
  ! Calculate characters of representations.
  !
  if (kpt%lcplx) then
    call zcharac_group(gvec%syms, gvec, kpt, 70, nspin, kpt%wfn(1, 1)%nstate)
  else
    call dcharac_group(gvec%syms, gvec, kpt, 70, nspin, kpt%wfn(1, 1)%nstate)
  end if

  !-------------------------------------------------------------------
  ! Initialize pol and k_p structures.
  !
  allocate (pol(gvec%syms%ntrans, qpt%nk))
  pol(:, :)%ntr = 0
  allocate (k_p(gvec%syms%ntrans, qpt%nk))
  k_p(:, :)%ncol = 0
  k_p(:, :)%nbuff = nbuff
  k_p(:, :)%lcache = lcache
  k_p(:, :)%isdf = doisdf
  ! Skip construction of kernel if we are doing RPA spectrum only.
  if (rpaonly) k_p(:, :)%ncol = -1

  call stopwatch(peinf%master, 'Calling setup_g')
  call timacc(2, 1, tsec)
  if (kpt%lcplx) then
    call zsetup_g(gvec, kpt, qpt, pol_in, pol, k_p, nspin, tdldacut, .false., isdf_in%lessmemory, .false.)
  else
    call dsetup_g(gvec, kpt, qpt, pol_in, pol, k_p, nspin, tdldacut, .false., isdf_in%lessmemory, .false.)
    call timacc(2, 2, tsec)
  end if

  if (doisdf) then
    ! open DEBUG file
#ifdef DEBUG
    rho_intp_dbg = 100001
    open (rho_intp_dbg, file="rho_intp_dbg.dat", status='unknown', iostat=info)
#endif
    ! --- prepare some inputs for the ISDF method ---
    ! W Gao find the index of highest occupied orbital
    ihomo = 1
    do isp = 1, nspin
      do ikp = 1, kpt%nk
        do ii = 1, kpt%wfn(isp, ikp)%nstate
          if (kpt%wfn(isp, ikp)%occ0(ii) > tol_occ .and. &
              ihomo < ii) then
            ihomo = ii
          end if
        end do ! ii loop
      end do ! ikp loop
    end do ! isp loop
    ! if n_intp can not be found in rgwbs.in or invalid (i.e., less than the
    ! number of occupied states), then set it to a default value
    if (n_intp < ihomo) then
      n_intp = int(2.0*ihomo)
    end if
    allocate (intp(n_intp))
    intp(1:n_intp) = 0
    ! --- find interpolation points for ISDF method ---
    call stopwatch(peinf%master, "before call cvt")
    call timacc(51, 1, tsec)
    if (intp_type == 1) then
      if (peinf%master) then
        write (*, *) " intp_type == 1"
        call cvt(gvec, kpt%rho, nspin, n_intp, intp)
      end if
    elseif (intp_type == 2) then
      if (peinf%master) write (*, *) " intp_type == 2"
      call cvt_wfn(gvec, kpt%wfn, nspin, kpt%nk, n_intp, intp)
    elseif (intp_type == 3) then
      if (peinf%master) write (*, *) " intp_type == 3, call cvt_parallel_sym1"
      call cvt_parallel_sym1(gvec, kpt%rho, nspin, n_intp, intp)
    elseif (intp_type == 4) then
      if (peinf%master) write (*, *) " intp_type == 4, call cvt_parallel_sym2"
      call cvt_parallel_sym2(gvec, kpt%rho, nspin, n_intp, intp)
    elseif (intp_type == 5) then
      if (peinf%master) then
        write (*, *) " intp_type = 4, read interpolation points from intp.dat"
        open (1111, file="intp.dat", form='formatted', status='old')
        do ii = 1, n_intp
          read (1111, *) intp(ii)
        end do
      end if
    else
      write (*, *) 'Type', intp_type, 'method for finding interpolation points is', &
        ' not implememted so far. The default method will be used.'
    end if
    call MPI_BARRIER(peinf%comm, info)
    call timacc(51, 2, tsec)
    call stopwatch(peinf%master, "after call cvt")
    ! broadcast intp to all processors
    call MPI_bcast(intp(1), n_intp, MPI_INTEGER, peinf%masterid, peinf%comm, info)
    ! pick out the interpolation points in reduced zone
    n_intp_r = 0
    do ipt = 1, n_intp
      if (mod(intp(ipt), gvec%syms%ntrans) == 0) n_intp_r = n_intp_r + 1
    end do
    allocate (intp_r_tmp(n_intp_r))
    ii = 0
    do ipt = 1, n_intp
      if (mod(intp(ipt), gvec%syms%ntrans) == 0) then
        ii = ii + 1
        intp_r_tmp(ii) = intp(ipt)
      end if
    end do
    allocate (rho_at_intp_r(n_intp_r))
    allocate (not_duplicate(n_intp_r))
    allocate (sort_idx(n_intp_r))
    not_duplicate = .true.
#ifdef DEBUG
    if (peinf%master) write (rho_intp_dbg, *) &
      "# rho at interpolation points "
#endif
    do ipt = 1, n_intp_r
      igrid = intp_r_tmp(ipt)
      jgrid = (igrid - 1)/gvec%syms%ntrans + 1
      rho_at_intp_r(ipt) = kpt%rho(jgrid, 1)
#ifdef DEBUG
      if (peinf%master) write (rho_intp_dbg, '(i7,i7,f30.23)') ipt, intp_r_tmp(ipt), kpt%rho(jgrid, 1)
#endif
    end do
    call quicksort(n_intp_r, rho_at_intp_r, sort_idx)
    n_dummy = 0
    ! Note: if there are some intp points having the same rho, then
    !       we pick up the last point of these duplicated points.
#ifdef DEBUG
    if (peinf%master) write (rho_intp_dbg, *) &
      "# rho at interpolation points (sorted) "
#endif
    do ipt = 1, n_intp_r - 1
      if (abs(rho_at_intp_r(sort_idx(ipt)) - &
              rho_at_intp_r(sort_idx(ipt + 1))) < 1.0e-14) then
        not_duplicate(sort_idx(ipt)) = .false.
#ifdef DEBUG
        if (peinf%master) write (rho_intp_dbg, '(i7,i7,f30.23,a)') sort_idx(ipt), &
          intp_r_tmp(sort_idx(ipt)), rho_at_intp_r(sort_idx(ipt)), 'F'
#endif
      else
        not_duplicate(sort_idx(ipt)) = .true.
        n_dummy = n_dummy + 1
#ifdef DEBUG
        if (peinf%master) write (rho_intp_dbg, '(i7,i7,f30.23,a)') sort_idx(ipt), &
          intp_r_tmp(sort_idx(ipt)), rho_at_intp_r(sort_idx(ipt)), 'T'
#endif
      end if
    end do
    ! The last point in intp_r_tmp(:)
    not_duplicate(sort_idx(n_intp_r)) = .true.
    n_dummy = n_dummy + 1
#ifdef DEBUG
    if (peinf%master) write (rho_intp_dbg, '(i7,i7,f30.23,a)') &
      sort_idx(n_intp_r), &
      intp_r_tmp(sort_idx(n_intp_r)), &
      rho_at_intp_r(sort_idx(n_intp_r)), 'T'
#endif
    if (peinf%master) then
      write (6, *) "Original n_intp_r =", n_intp_r
      write (6, *) "After removing duplicated points, n_intp_r =", n_dummy
    end if
    allocate (intp_r(n_dummy))
    ii = 0
    do ipt = 1, n_intp_r
      if (not_duplicate(ipt)) then
        ii = ii + 1
        intp_r(ii) = intp_r_tmp(ipt)
      end if
    end do
    n_intp_r = n_dummy
    deallocate (intp_r_tmp)
    deallocate (rho_at_intp_r)
    deallocate (not_duplicate)
    !
    if (peinf%master) write (*, *) " Finding interpolation points successfully. "
    !
    ! ISDF will deal with all the pair products of wave functions as defined in
    ! pol_in(isp)%vmap(:) and pol_in(isp)%cmap(:)
    allocate (ncv(nspin, kpt%nk, gvec%syms%ntrans))
    allocate (nv(nspin, kpt%nk, gvec%syms%ntrans))
    allocate (nc(nspin, kpt%nk, gvec%syms%ntrans))
    allocate (icv(gvec%syms%ntrans))
    allocate (idum(gvec%syms%ntrans))
    ! Initializing with zeros
    ncv = 0
    nc = 0
    nv = 0
    maxivv = 0 ! the absolute index of the highest occupied state
    maxicc = 0 ! the absolute index of the highest unoccupied state
    do isp = 1, nspin
      maxivv = max(maxval(pol_in(isp)%vmap(:)), maxivv)
      maxicc = max(maxval(pol_in(isp)%cmap(:)), maxicc)
    end do
    allocate (pairmap(maxivv, maxicc, nspin, kpt%nk, gvec%syms%ntrans))
    pairmap = 0
    do ikp = 1, kpt%nk
      do isp = 1, nspin
        do iv = 1, pol_in(isp)%nval
          ivv = pol_in(isp)%vmap(iv)
          irp = kpt%wfn(isp, ikp)%irep(ivv)
          nv(isp, ikp, irp) = nv(isp, ikp, irp) + 1
        end do
        do ic = 1, pol_in(isp)%ncond
          icc = pol_in(isp)%cmap(ic)
          irp = kpt%wfn(isp, ikp)%irep(icc)
          nc(isp, ikp, irp) = nc(isp, ikp, irp) + 1
        end do
        icv = 0
        do iv = 1, pol_in(isp)%nval
          do ic = 1, pol_in(isp)%ncond
            ivv = pol_in(isp)%vmap(iv)
            icc = pol_in(isp)%cmap(ic)
            virp = kpt%wfn(isp, ikp)%irep(ivv)
            cirp = kpt%wfn(isp, ikp)%irep(icc)
            vcirp = gvec%syms%prod(virp, cirp)
            icv(vcirp) = icv(vcirp) + 1
            ! pairmap maps the real valence and conduction band index to
            ! the |vc> pair index obtained in isdf.f90
            pairmap(ivv, icc, isp, ikp, vcirp) = icv(vcirp)
          end do
        end do
        ncv(isp, ikp, 1:gvec%syms%ntrans) = icv(1:gvec%syms%ntrans)
        if (peinf%master) write (6, *) "isp", isp, ", ikp", ikp, ", ncv ", &
          (icv(ii), ii=1, gvec%syms%ntrans)
      end do ! isp
    end do ! ikp
    maxnc = maxval(nc)
    maxnv = maxval(nv)
    maxncv = maxval(ncv)
    allocate (invpairmap(2, maxncv, nspin, kpt%nk, gvec%syms%ntrans))
    allocate (ivlist(maxnv, nspin, kpt%nk, gvec%syms%ntrans))
    allocate (iclist(maxnc, nspin, kpt%nk, gvec%syms%ntrans))
    ! Initialization with zeros
    invpairmap = 0
    ivlist = 0
    iclist = 0
    do ikp = 1, kpt%nk
      do isp = 1, nspin
        if (peinf%master) write (6, *) "ikp ", ikp, " isp", isp
        idum = 0
        do iv = 1, pol_in(isp)%nval
          ivv = pol_in(isp)%vmap(iv)
          irp = kpt%wfn(isp, ikp)%irep(ivv)
          idum(irp) = idum(irp) + 1
#ifdef DEBUG
          if (peinf%master) write (6, *) "ivv ", ivv, "irp", irp, &
            "idum(irp)", idum(irp)
#endif
          ivlist(idum(irp), isp, ikp, irp) = ivv
        end do
        idum = 0
        do ic = 1, pol_in(isp)%ncond
          icc = pol_in(isp)%cmap(ic)
          irp = kpt%wfn(isp, ikp)%irep(icc)
          idum(irp) = idum(irp) + 1
#ifdef DEBUG
          if (peinf%master) write (6, *) "icc ", icc, "irp", irp, &
            " idum(irp)", idum(irp)
#endif
          iclist(idum(irp), isp, ikp, irp) = icc
        end do
        icv = 0
        do iv = 1, pol_in(isp)%nval
          do ic = 1, pol_in(isp)%ncond
            ivv = pol_in(isp)%vmap(iv)
            icc = pol_in(isp)%cmap(ic)
            virp = kpt%wfn(isp, ikp)%irep(ivv)
            cirp = kpt%wfn(isp, ikp)%irep(icc)
            vcirp = gvec%syms%prod(virp, cirp)
            icv(vcirp) = icv(vcirp) + 1
            invpairmap(1, icv(vcirp), isp, ikp, vcirp) = ivv
            invpairmap(2, icv(vcirp), isp, ikp, vcirp) = icc
          end do ! ic loop
        end do ! iv loop
      end do ! isp loop
    end do ! ikp loop
#ifdef DEBUG
    if (peinf%inode == 1) then
      do ikp = 1, kpt%nk
        do isp = 1, nspin
          do irp = 1, gvec%syms%ntrans
            write (6, *) "ikp ", ikp, "isp ", isp, "irp ", irp
            write (6, '(10(i6,i3))') (ivlist(ii, isp, ikp, irp), ii=1, nv(isp, ikp, irp))
            write (6, '(10(i6,i3))') (iclist(ii, isp, ikp, irp), ii=1, nc(isp, ikp, irp))
            write (6, '(10("(",2i6,")"))') (invpairmap(1, ii, isp, ikp, irp), &
                                            invpairmap(2, ii, isp, ikp, irp), ii=1, ncv(isp, ikp, irp))
          end do
        end do
      end do
    end if
    if (peinf%master) write (6, *) " Finish getting invpairmap and pairmap. "
#endif
    !

    isdf_in%maxivv = maxivv
    isdf_in%maxicc = maxicc
    isdf_in%maxncv = maxncv
    isdf_in%maxnv = maxnv
    isdf_in%maxnc = maxnc
    !
    allocate (isdf_in%ncv(nspin, kpt%nk, gvec%syms%ntrans))
    isdf_in%ncv = ncv
    allocate (isdf_in%invpairmap(2, maxncv, nspin, kpt%nk, gvec%syms%ntrans))
    isdf_in%invpairmap = invpairmap
    allocate (isdf_in%pairmap(maxivv, maxicc, nspin, kpt%nk, gvec%syms%ntrans))
    isdf_in%pairmap = pairmap
    allocate (isdf_in%nv(nspin, kpt%nk, gvec%syms%ntrans))
    isdf_in%nv = nv
    allocate (isdf_in%ivlist(maxnv, nspin, kpt%nk, gvec%syms%ntrans))
    isdf_in%ivlist = ivlist
    allocate (isdf_in%nc(nspin, kpt%nk, gvec%syms%ntrans))
    isdf_in%nc = nc
    allocate (isdf_in%iclist(maxnc, nspin, kpt%nk, gvec%syms%ntrans))
    isdf_in%iclist = iclist
    !
    !call MPI_BARRIER(peinf%comm, info)
    !stop
    isdf_in%n_slice = 7
    isdf_in%n_intp_r = n_intp_r
    allocate (isdf_in%intp_r(n_intp_r))
    isdf_in%intp_r(1:n_intp_r) = intp_r(1:n_intp_r)
    if (isdf_in%lessmemory) then
      do isp = 1, nspin
        do ikp = 1, kpt%nk
          if (kpt%wfn(isp, ikp)%nmem /= kpt%wfn(1, 1)%nmem) then
            write (6, '(a,i2,a,i5,a)') " kpt%wfn(", isp, ",", ikp, &
              ")%nmem is not equal to kpt%wfn(1,1)%nmem."
            write (6, *) " Can't allocate isdf_in%Psi_intp "
          end if
        end do ! ikp
      end do ! isp
      allocate (isdf_in%Psi_intp(n_intp_r, kpt%wfn(1, 1)%nmem, nspin, kpt%nk))
      isdf_in%Psi_intp = zero
    else ! don't save memory
      ! Cmtrx is intialized to zero in isdf_parallel.f90
      allocate (isdf_in%Cmtrx(n_intp_r, maxncv, nspin, kpt%nk, gvec%syms%ntrans))
      ! Mmtrx is intialized to zero in isdf_parallel.f90
    end if
    allocate (isdf_in%Mmtrx(n_intp_r, n_intp_r, nspin, nspin, kpt%nk, 2, gvec%syms%ntrans))
    !
    deallocate (intp_r)
    deallocate (ncv)
    deallocate (invpairmap)
    deallocate (pairmap)
    deallocate (nv)
    deallocate (ivlist)
    deallocate (nc)
    deallocate (iclist)
    !
    ! --- perform ISDF method to interpolate pair products of wave functions ---
    !
    call stopwatch(peinf%master, "before call ISDF (interpolative separable &
&       density fitting")
    kflag = 1
    if (trip_flag) kflag = 3
    if (noxchange) kflag = 2
    if (nolda) kflag = 0
    if (peinf%inode == 0) then
      write (6, *) "Reduced interpolation points = ", n_intp_r
      write (6, *) "Max number of orbital pairs ", maxncv
      write (6, *) "nspin ", nspin
      write (6, *) "num of kpoints ", kpt%nk
      write (6, *) "num of sym ops ", gvec%syms%ntrans
    end if
    verbose = .false.
    call timacc(52, 1, tsec)
    if (isdf_in%lessmemory) then
      if (peinf%master) write (*, *) " call isdf_parallel_sym_lessmemory"
      call isdf_parallel_sym_lessmemory(gvec, pol_in, kpt, nspin, isdf_in, kflag, opt, verbose)
      if (peinf%master) write (*, *) " done isdf"
    else
      if (peinf%master) write (*, *) " call isdf_parallel_sym"
      call isdf_parallel_sym(gvec, pol_in, kpt, nspin, isdf_in, kflag, verbose)
      if (peinf%master) write (*, *) " done isdf"
    end if
    call timacc(52, 2, tsec)
    call stopwatch(peinf%master, " after call isdf")

#ifdef DEBUG
    close (rho_intp_dbg)
#endif
  end if ! doisdf

  ! The outputs are Cmtrx and Mmtrx, which are used by k_integrate_isdf() for
  ! calculation of K(v,c,v',c') later !!
  ! --- finished ISDF ---

  !-------------------------------------------------------------------
  ! Print out warnings, information etc.
  !
  if (kpt%lcplx) tamm_d = .true.
  if (peinf%master) then
    write (6, '(/,a,/,/,a,/,2a,/)') repeat('-', 65), &
      ' Polarizability input data: ', ' ', repeat('-', 25)
    write (6, '(2a)') ' Number of transitions per representation ', &
      'in TDLDA polarizabillity:'
    write (6, '(8i8)') ((pol(irp, iq)%ntr, irp=1, gvec%syms%ntrans), iq=1, qpt%nk)
    write (6, '(a,i10,/)') ' total = ', sum(pol(:, :)%ntr)
    if (tdldacut > zero) then
      write (6, *) 'Energy cutoff applied in TDLDA polarizability = ', &
        tdldacut*ryd, ' eV'
    else
      write (6, *) 'No energy cutoff in TDLDA polarizability'
    end if
    if (nolda) write (6, '(/,a,/)') &
      ' LDA kernel is not included in polarizability'
    if (tamm_d) then
      write (6, '(2a,/)') ' Calculating TDLDA ', &
        'polarizability within the Tamm-Dancoff approximation.'
    else
      write (6, '(2a,/)') ' Not using the Tamm-Dancoff ', &
        'approximation in TDLDA polarizability.'
    end if
    !
    ! WARNING: for periodic boundary conditions, use a truncated Coulomb
    ! potential (see Hanke's article) in diagonalization. For that, set
    ! is_pbc = .false.
    ! In other applications (e.g., when calculating self-energies), one
    ! should not truncate the Coulomb potential.
    !
    write (6, '(a,/)') &
      ' Coulomb potential is being truncated in the TDLDA equation.'
    if (noxchange) write (6, '(a,/,a,/,a,/)') repeat('*', 65), &
      ' Exchange kernel not included in TDLDA ', repeat('*', 65)
    !
    ! Estimate memory usage.
    ! For diagonalization, we store 4 matrices: hamiltonian/eigenvectors,
    ! temporary array (in eigensolver), Hartree kernel, and XC kernel.
    ! Parallelized diagonalization also uses a second temporary matrix.
    !
    mem1 = sum(kpt%wfn(:, :)%nmem)*two/real(nspin, dp)* &
           real(nspin*gvec%nr, dp)/two/131072.d0/real(w_grp%npes, dp)
    if (kpt%lcplx) mem1 = mem1*two
    write (6, '(/,a,f10.2,a)') ' Memory needed to store wavefunctions : ', &
      mem1, ' MB/proc.'
    mem1 = real(5*gvec%nr, dp)/131072.d0
    if (kpt%lcplx) mem1 = mem1*two
    write (6, '(a,a,f10.2,a)') ' Memory needed to calculate kernel ', &
      'matrix elements : ', mem1, ' MB/proc.'
    xmax = 0
    do iq = 1, qpt%nk
      do irp = 1, gvec%syms%ntrans
        xsum = real(pol(irp, iq)%ntr*pol(irp, iq)%ntr, dp)
        if (xmax < xsum) xmax = xsum
      end do
    end do
    mem1 = xmax/1024.d0*3.d0/128.d0/r_grp%npes
    if (r_grp%npes > 1) mem1 = mem1*4.d0/3.d0
    write (6, '(/,a,f10.2,a/)') ' Memory needed for diagonalization : ', &
      mem1, ' MB/proc.'
    write (6, '(a,/)') repeat('-', 65)
    call flush (6)
  end if

  if (noxchange) trunc_c = .false.

  ii = 0
  xmax = zero
  call stopwatch(peinf%master, "after call calculate_tdlda")
  if (kpt%lcplx) then
    call zcalculate_tdlda(gvec, kpt, qpt, k_p, pol, nspin, ii, &
                          tamm_d, nolda, rpaonly, trip_flag, noxchange, trunc_c, xsum, xmax, &
                          isdf_in, opt)
  else
    call dcalculate_tdlda(gvec, kpt, qpt, k_p, pol, nspin, ii, &
                          tamm_d, nolda, rpaonly, trip_flag, noxchange, trunc_c, xsum, xmax, &
                          isdf_in, opt)
  end if
  call stopwatch(peinf%master, "after call calculate_tdlda")

  ! Deallocate arrays for ISDF method
  if (doisdf) then
    deallocate (isdf_in%ncv)
    deallocate (isdf_in%pairmap)
    deallocate (isdf_in%invpairmap)
    deallocate (isdf_in%ivlist)
    deallocate (isdf_in%iclist)
    deallocate (isdf_in%nv)
    deallocate (isdf_in%nc)
    if (isdf_in%lessmemory) then
      deallocate (isdf_in%Psi_intp)
    else
      deallocate (isdf_in%Cmtrx)
    end if
    deallocate (isdf_in%Mmtrx)
  end if

  !-------------------------------------------------------------------
  ! Time accounting.
  !
  ii = 3
  jj = 11
  allocate (routnam(ii + jj))
  allocate (timerlist(ii + jj))
  routnam(1) = 'SETUP_T:'; timerlist(1) = 2
  routnam(2) = 'KERNEL:'; timerlist(2) = 3
  routnam(3) = 'DIAG_POL:'; timerlist(3) = 4

  routnam(4) = 'POISSON_FFT:'; timerlist(4) = 11
  routnam(5) = 'EIGENSOLVER:'; timerlist(5) = 12
  routnam(6) = 'INTEGRATION:'; timerlist(6) = 13
  routnam(7) = 'Find intp pts:'; timerlist(7) = 51
  routnam(8) = 'ISDF_PARALLEL:'; timerlist(8) = 52
  routnam(9) = 'Calc intp vectors:'; timerlist(9) = 53
  routnam(10) = 'Calc <zeta|K|zeta>:'; timerlist(10) = 54
  routnam(11) = '(in cvt step1):'; timerlist(11) = 61
  routnam(12) = '(in cvt step2):'; timerlist(12) = 62
  routnam(13) = 'k_integrate_isdf:'; timerlist(13) = 63
  routnam(14) = 'kernel k_print:'; timerlist(14) = 64

  call timacc(1, 2, tsec)
!#ifdef HIPMAGMA
!  if (opt%linear_algebra .eq. 2 .or. opt%eigsolver .eq. 2) then
!  call magmaf_finalize()
!  endif
!#endif
  call finalize(peinf%master, peinf%comm, ii, jj, routnam, timerlist)

end program tdlda
!===================================================================
