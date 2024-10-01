
! Weiwei Gao, Feb. 2018
!
! In the density fitting method, the products of wave functions are
! written as linear combinations of interpolation vectors zeta(r).
! See equation 5 in J. Chem. Theory. Comput. 2017, 13, 5420-5431:
!
! \phi_i(r)\psi_j(r) = \sum^{Nu}_{u=1} \zeta_u(r) \phi_i(r_u) \psi_j(r_u)
!
! <==> written in a matrix form: Z = Zeta * C
!  ________________________________________________________________________________
! |  Matrix      |     Shape      |        Contents                                |
! |______________|________________|________________________________________________|
! |   P(r,r_u)   |  Ng * N_u      |    \sum_i \phi_i(r)\psi_i(r_u)                 |
! |              |                |     i=1...Nv, r are full-grid points           |
! |______________|________________|________________________________________________|
! |   Q(r,r_u)   |  Ng * N_u      |    \sum_i \phi_i(r)\psi_i(r_u)                 |
! |              |                |     i=1...Nc, r are full-grid points           |
! |______________|________________|________________________________________________|
! |  Zeta_u(r)   |  Ng *  Nu      |    {..., \zeta_u(r), ...}                      |
! |              |                |     u=1...Nu                                   |
! |______________|________________|________________________________________________|
! |  C(r_u,i,j)  |  Nu * (Nc*Nv)  |    {..., \phi_i(ru)*\psi_j(ru), ...}           |
! |              |                |                                                |
! |______________|________________|________________________________________________|
!  _______________________________________________________________________
! |  Matrix      |                  Parallelization                       |
! |______________|________________________________________________________|
! |   P(r,r_u)   |             Each proc store part of Z                  |
! |              |             P(ngfl, Nc*Nv), ngfl = mydim               |
! |______________|________________________________________________________|
! |   Q(r,r_u)   |             Each proc store part of Z                  |
! |              |             Q(ngfl, Nc*Nv), ngfl = mydim               |
! |______________|________________________________________________________|
! |  Zeta_u(r)   |             Each proc store part of zeta               |
! |              |         zeta(ngfl, Nu), ngfl = mydim*ntrans            |
! |______________|________________________________________________________|
! |  C(r_u,i,j)  |             Every proc store a full copy of C          |
! |              |              Cmtrx ( Nu, Nc*Nv )                       |
! |______________|________________________________________________________|
!
!
! This subroutine calculates the interpolation vectors zeta_u(r)
!
! n_intp   : the number of interpolation vectors or points, or Nu
! n_intp_r : the number of interpolation vectors or points in reduced r-space domain
! intp     : the index of interpolation points in the full grid
! zeta(gvec%nr, n_intp_r) : the interpolation vectors
! kpt%wfn(isp,ikp)%dwf(:,:) : store the wavefunctions \phi_i, \psi_j
!
! For now, this is only partially parallelized
! 1. each processor evaluates part of the matrix Z and part of the matrix C
! 2. collect matrix C from different procs. distribute matrix Z to procs.
! 3. every proc in w_grp solve part of the linear equations to get zeta
!
! This subroutine consists of the following main steps:
!
! Step 1.0: Prepare P(r,r_u)
! Step 1.5: Prepare A=C.C^T=P(r_u,r_v)Q(r_u,r_v), B=Z.C^T=P(r,r_u)Q(r,r_v)
! Step 1.6: Deallocate P and Q
! Step 2.0: solve a linear quation A.zeta=B to get zeta
! Step 3.0: Calculate Mmtrx(1) = <zeta(i)| V_coul(r,r') | zeta(j)> and
!       Mmtrx(2) = <zeta(i)| f^LDA(r) | zeta(j)>
! Step 3.5: Calculate C(r_u,i,j)

subroutine isdf_parallel_sym(gvec, pol_in, kpt, nspin, isdf_in, kflag, &
                             verbose)

  use typedefs
  use mpi_module
  use myconstants
  use fft_module
  use xc_functionals
  implicit none
#ifdef MPI
  include 'mpif.h'
#endif

  type(gspace), intent(in) :: gvec
  type(polinfo), intent(in), dimension(2) :: pol_in
  type(kptinfo), intent(in) :: kpt
  type(ISDF), intent(inout) :: isdf_in
  ! n_intp_r: number of interpolation points in reduced r-space domain
  ! intp_r(n_intp_r): the index of interpolation points in full grids
  ! ncv(2): the number of wave function pairs for two spins
  !
  ! invpairmap maps the index of pair |cv> to c and v
  integer, intent(in) :: nspin, kflag
  ! If verbose is true, then print out additional debug information
  logical, intent(in) :: verbose

  ! ------------------ Local variables ------------------
  ! P(gvec%mydim, n_intp_r, nspin, kpt%nk)
  ! Q(gvec%mydim, n_intp_r, nspin, kpt%nk)
  ! Cmtrx(n_intp_r, Nc*Nv, nspin, kpt%nk)
  ! Issue need to be addressed: if there are N wfn_grp, does it mean these matrix need to be
  ! calculated by all processors at the same time, or can we distribute the workload
  ! later ??
  real(dp), allocatable :: &
    PsiV(:, :, :), PsiV_intp(:, :, :), &   ! PsiV: wfn on reduced domain
    PsiC(:, :, :), PsiC_intp(:, :, :), &
    P(:, :, :), P_intp(:, :, :), &   ! P on reduced domain
    Q(:, :, :), Q_intp(:, :, :), &   ! Q on reduced domain
    zeta(:, :, :, :), tmp_Zmtrx(:), diff_vec(:), &
    fxc(:, :, :), fxc_loc(:, :, :), fzeta(:), acc_B(:, :), &
    ! matrices and vectors used for solving linear equations
    Amtrx(:, :), Bmtrx(:, :), Xmtrx(:, :), tmpmtrx(:, :, :, :), &
    rho_h(:), rho_h_distr(:), &
    tmpCmtrx(:, :)
  real(dp) :: diff, weight, qkt(3), tsec(2), &
              vcharac(gvec%syms%ntrans), ccharac(gvec%syms%ntrans), &
              tmpvec(isdf_in%n_intp_r), matel
  integer, allocatable :: inv_ivlist(:, :, :, :), inv_iclist(:, :, :, :), ipiv(:)
  ! counters and temporary integers
  integer :: ipt, ii, jj, iv, ic, icv, irp, jrp, rsp, csp, i_row, i_col, lrp1, lrp2, &
             IVV, ICC, JVV, JCC, isp, ikp, errinfo, ipe, einfo, &
             ivrp, icrp, n_intp_r, maxncv, &
             maxnv, maxnc, ig, icv1, icv2
  integer :: status, mpi_status(MPI_STATUS_SIZE)
  ! The index of w_grp, r_grp, processor that works on the calculation of
  ! zeta(:,:,nspin,nk)
  integer :: workwgrp, workrgrp, workproc
  ! Each processor in a w_grp store part of the wave function
  ! offset: index of grid point from which a processor start to store the wave function
  ! ncount: number of elements of wave functions that a processor stores
  integer, dimension(0:w_grp%npes - 1) :: offset, ncount
  ! temporary dummy variable
  integer, dimension(0:w_grp%npes - 1) :: idum
  ! the number of grid points in irreducible wedge, ngr = gvec%nr
  integer :: ngr, ngrid
  ! the number of full grid point, ngf = ngr * (# of sym operations)
  integer :: ngfl, iptf, iptr, ioff, ioff1, ioff2, rcond

#ifdef DEBUG
  ! variables for debug and test of accuracy
  character(50) :: dbg_filename = "isdf_dbg.dat"
  integer :: dbgunit = 20171130
  integer :: outdbg
#endif
  ! external functions
  real(dp), external :: ddot
  type(xc_type) :: xc_lda
  !
  workrgrp = 0
  workwgrp = 0
  call MPI_BCAST(workrgrp, 1, MPI_INTEGER, peinf%masterid, peinf%comm, errinfo)
  call MPI_BCAST(workwgrp, 1, MPI_INTEGER, peinf%masterid, peinf%comm, errinfo)

#ifdef DEBUG
  if (peinf%master .and. verbose) then
    !
    write (*, *) "call isdf_parallel(), write debug info to ", dbg_filename
    !
    open (dbgunit, file=dbg_filename, form='formatted', status='replace')
    !
  end if
#endif
  !
  ngrid = w_grp%mydim
  ! the number of real-space grid points in reduced real-space domain
  ngr = gvec%nr
  ! the number of grid points of interpolation vectors zeta(:) stored in each processor
  ngfl = w_grp%mydim*gvec%syms%ntrans
  n_intp_r = isdf_in%n_intp_r
  maxncv = isdf_in%maxncv
  maxnv = isdf_in%maxnv
  maxnc = isdf_in%maxnc
  !
  ! each processor stores part of P, Q and zeta
  !
  ! Allocating intermediate variables
  allocate (P(w_grp%mydim, n_intp_r, gvec%syms%ntrans))
  allocate (Q(w_grp%mydim, n_intp_r, gvec%syms%ntrans))
  allocate (acc_B(w_grp%mydim, n_intp_r))
  allocate (P_intp(n_intp_r, n_intp_r, gvec%syms%ntrans))
  allocate (Q_intp(n_intp_r, n_intp_r, gvec%syms%ntrans))
  allocate (Amtrx(n_intp_r, n_intp_r))
  allocate (Bmtrx(n_intp_r, w_grp%mydim))
  allocate (Xmtrx(n_intp_r, w_grp%mydim))
  allocate (zeta(w_grp%mydim, n_intp_r, nspin, gvec%syms%ntrans)) ! interpolation vectors
  allocate (PsiV(w_grp%mydim, maxnv, gvec%syms%ntrans))
  allocate (PsiV_intp(n_intp_r, maxnv, gvec%syms%ntrans))
  allocate (PsiC(w_grp%mydim, maxnc, gvec%syms%ntrans))
  allocate (PsiC_intp(n_intp_r, maxnc, gvec%syms%ntrans))
  allocate (tmpCmtrx(n_intp_r, maxncv))
  allocate (tmp_Zmtrx(w_grp%mydim))
  allocate (diff_vec(w_grp%mydim))
  !
  ! initialize matrices with zero
  !
  isdf_in%Cmtrx = zero
  zeta = zero
  !
#ifdef DEBUG
  if (verbose) then
    do ipe = 1, peinf%npes
      if (peinf%inode == ipe) then
        write (6, '(a)') " Index of interpolation points in full domain: "
        write (6, *) "n_intp_r", n_intp_r
        write (6, '(5i14)') (isdf_in%intp_r(ii), ii=1, n_intp_r)
      end if
      call MPI_BARRIER(peinf%comm, errinfo)
    end do
  end if
#endif
  !
  idum = 0
  idum(w_grp%inode) = w_grp%offset + 1
  call MPI_ALLREDUCE(idum, offset, w_grp%npes, MPI_INTEGER, MPI_SUM, &
                     w_grp%comm, errinfo)
  !
  idum = 0
  idum(w_grp%inode) = w_grp%mydim
  call MPI_ALLREDUCE(idum, ncount, w_grp%npes, MPI_INTEGER, MPI_SUM, &
                     w_grp%comm, errinfo)
  !
#ifdef DEBUG
  if (verbose .and. peinf%master) then
    !
    write (dbgunit, *) " in isdf() "
    write (dbgunit, *) " w_grp%mygr ", w_grp%mygr, " w_grp%inode = ", w_grp%inode, &
      " workwgrp = ", workwgrp, &
      " workrgrp = ", workrgrp, &
      " offset: ", (offset(ii), ii=0, w_grp%npes - 1), &
      " ncount: ", (ncount(ii), ii=0, w_grp%npes - 1)
    !
  end if
#endif
  !
  ! isdf_in%ivlist(:) maps the index of valence states used in
  !   calculation (i.e., stored in memory) to the real index of valence states
  allocate (inv_ivlist(isdf_in%maxivv, nspin, kpt%nk, gvec%syms%ntrans)) ! For now, only deal with confined system. So we assume kpt%nk=1, and there is no dependence on k here
  allocate (inv_iclist(isdf_in%maxicc, nspin, kpt%nk, gvec%syms%ntrans))
  ! inv_ivlist(:) maps the real index of valence states
  !   to the index of valence states used in the calculation
  inv_ivlist = 0
  inv_iclist = 0
  !
  do ikp = 1, kpt%nk
    do isp = 1, nspin
      !
      do irp = 1, gvec%syms%ntrans
        do iv = 1, isdf_in%nv(isp, ikp, irp)
          IVV = isdf_in%ivlist(iv, isp, ikp, irp)
          inv_ivlist(IVV, isp, ikp, irp) = iv
        end do
        !
        do ic = 1, isdf_in%nc(isp, ikp, irp)
          ICC = isdf_in%iclist(ic, isp, ikp, irp)
          inv_iclist(ICC, isp, ikp, irp) = ic
        end do
        !
      end do
      !
    end do
  end do
  !
  isdf_in%Mmtrx = zero ! Mmtrx dimension: Mmtrx(n_intp_r, n_intp_r, nspin, nspin)
  !
  ! qkt is set to zero. This is only valid for
  !  tests of nonperiodic system, and should be updated later.
  !
  qkt = 0
  allocate (rho_h(ngr))     ! note: gvec%nr is equal to w_grp%nr
  allocate (rho_h_distr(w_grp%ldn*w_grp%npes))
  !
  ! kflag = 0 : calculate kernel K^x  ( Coulomb )
  !         1 :                  K^x + K^f ( Coulomb + 1/2 * F_xc )
  !         2 :                  K^f   ( 1/2 * F_xc )
  !         3 :                  K^t   ( 1/2 * F_xc for spin triplet )
  !
  ! Generate Coulomb potential
  !
  if (kflag < 2) then
    call dinitialize_FFT(peinf%inode, fft_box)
    call dcreate_coul_0D(gvec%bdot, qkt, fft_box)
    if (peinf%master) write (6, *) " finished creating FFT plans and coul_0D"
  end if
  !
  ! Calculate LDA kernel
  !
  if (kflag > 0) then
    !
    allocate (fxc(w_grp%mydim, nspin, nspin), stat=errinfo)
    allocate (fzeta(w_grp%ldn*w_grp%npes), stat=errinfo)
    fxc = zero
    !
    ! Copy the charge density to fxc
    !
    do isp = 1, nspin
      call dcopy(w_grp%mydim, kpt%rho(w_grp%offset + 1, isp), 1, &
                 fxc(1, isp, 1), 1)
    end do
    call xc_init(nspin, XC_LDA_X, XC_LDA_C_PZ, 0, zero, one, .false., xc_lda)
    xc_lda%has_grad = .false.
    call fxc_get(xc_lda, nspin, w_grp%mydim, kflag, fxc)
    call xc_end(xc_lda)
    !
#ifdef DEBUG
    write (6, *) "inode: ", peinf%inode, ": finished initializing xc functional"
#endif
  end if
  !
  do ikp = 1, kpt%nk
    if (peinf%master) write (6, *) " ikp = ", ikp
    !
    do isp = 1, nspin
      if (peinf%master) write (6, *) " isp = ", isp
      !
      ! initialize matrices with zero
      !
      P = zero
      Q = zero
      PsiV = zero
      PsiC = zero
      PsiV_intp = zero
      PsiC_intp = zero
      call timacc(53, 1, tsec)
      !
      ! The following loop calculate P(r,r_u,jrp), Q(r,r_u,jrp) for all representations
      !
#ifdef DEBUG
      if (peinf%master) then
        write (6, *) " kpt%wfn(isp,ikp)%map ", &
          kpt%wfn(isp, ikp)%map(:)
      end if
#endif
      do jrp = 1, gvec%syms%ntrans
        if (peinf%master .and. verbose) write (6, *) " jrp = ", jrp
        !
#ifdef DEBUG
        if (verbose .and. peinf%master) &
          write (dbgunit, *) ' isp = ', isp, ', ikp = ', ikp
#endif
        !
        do iv = 1, isdf_in%nv(isp, ikp, jrp)
          IVV = isdf_in%ivlist(iv, isp, ikp, jrp)
          JVV = kpt%wfn(isp, ikp)%map(IVV)
          ! PsiV_i(r)
          call dcopy(w_grp%mydim, &
                     kpt%wfn(isp, ikp)%dwf(1, JVV), 1, PsiV(1, iv, jrp), 1)
        end do
        !
        do ic = 1, isdf_in%nc(isp, ikp, jrp)
          ICC = isdf_in%iclist(ic, isp, ikp, jrp)
          JCC = kpt%wfn(isp, ikp)%map(ICC)
          ! PsiC_i(r)
          call dcopy(w_grp%mydim, &
                     kpt%wfn(isp, ikp)%dwf(1, JCC), 1, PsiC(1, ic, jrp), 1)
        end do
        !
        ! pick the interpolation points
        !
        do ipt = 1, n_intp_r
          iptf = isdf_in%intp_r(ipt)
          if (gvec%syms%ntrans > 1) then
            iptr = iptf/gvec%syms%ntrans + 1
          else
            iptr = iptf
          end if
          ioff1 = offset(w_grp%inode)
          ioff2 = ioff1 + w_grp%mydim
          if (iptr >= ioff1 &
              .and. &
              iptr < ioff2) then
            jj = iptr - ioff1 + 1
            ! if ( icv .eq. 1 .and. verbose ) write(*, '(4(a,i5))') &
            !    ", wgrp%inode ", w_grp%inode, ", ipt ", ipt, &
            !    ", iptf ", iptf, ", jj ", jj
            PsiV_intp(ipt, 1:isdf_in%nv(isp, ikp, jrp), jrp) = PsiV(jj, 1:isdf_in%nv(isp, ikp, jrp), jrp)
            PsiC_intp(ipt, 1:isdf_in%nc(isp, ikp, jrp), jrp) = PsiC(jj, 1:isdf_in%nc(isp, ikp, jrp), jrp)
          end if
        end do
        !
        call MPI_ALLREDUCE(MPI_IN_PLACE, PsiV_intp(1, 1, jrp), n_intp_r*isdf_in%nv(isp, ikp, jrp), MPI_DOUBLE, MPI_SUM, &
                           w_grp%comm, errinfo)
        call MPI_ALLREDUCE(MPI_IN_PLACE, PsiC_intp(1, 1, jrp), n_intp_r*isdf_in%nc(isp, ikp, jrp), MPI_DOUBLE, MPI_SUM, &
                           w_grp%comm, errinfo)
        !
        ! Prepare P and Q
        !
        ! P(r,r_u,jrp) = \sum_{V~jrp} \PsiV(r) \PsiV(r_u)
        call dgemm('n', 't', w_grp%mydim, n_intp_r, isdf_in%nv(isp, ikp, jrp), one, &
                   PsiV(1, 1, jrp), w_grp%mydim, PsiV_intp(1, 1, jrp), n_intp_r, zero, P(1, 1, jrp), w_grp%mydim)
        ! Q(r,r_u,jrp) = \sum_{C~jrp} \PsiC(r) \PsiC(r_u)
        call dgemm('n', 't', w_grp%mydim, n_intp_r, isdf_in%nc(isp, ikp, jrp), one, &
                   PsiC(1, 1, jrp), w_grp%mydim, PsiC_intp(1, 1, jrp), n_intp_r, zero, Q(1, 1, jrp), w_grp%mydim)
        ! Calculate P_intp(r_u,r_u,jrp), and Q_intp(r_u,r_u,jrp)
        call dgemm('n', 't', n_intp_r, n_intp_r, isdf_in%nv(isp, ikp, jrp), one, &
                   PsiV_intp(1, 1, jrp), n_intp_r, PsiV_intp(1, 1, jrp), n_intp_r, zero, &
                   P_intp(1, 1, jrp), n_intp_r)
        call dgemm('n', 't', n_intp_r, n_intp_r, isdf_in%nc(isp, ikp, jrp), one, &
                   PsiC_intp(1, 1, jrp), n_intp_r, PsiC_intp(1, 1, jrp), n_intp_r, zero, &
                   Q_intp(1, 1, jrp), n_intp_r)
        !
#ifdef DEBUG
        if (verbose .and. peinf%master) then
          write (dbgunit, *) "jrp =", jrp
          write (dbgunit, '("PsiV = ")')
          call printmatrix(PsiV(1, 1, jrp), w_grp%mydim, isdf_in%nv(isp, ikp, jrp), dbgunit)
          write (dbgunit, '("PsiC = ")')
          call printmatrix(PsiC(1, 1, jrp), w_grp%mydim, isdf_in%nc(isp, ikp, jrp), dbgunit)
          write (dbgunit, '("PsiV_intp = ")')
          call printmatrix(PsiV_intp(1, 1, jrp), n_intp_r, isdf_in%nv(isp, ikp, jrp), dbgunit)
          write (dbgunit, '("PsiC_intp = ")')
          call printmatrix(PsiC_intp(1, 1, jrp), n_intp_r, isdf_in%nc(isp, ikp, jrp), dbgunit)
          write (dbgunit, '("P = ")')
          call printmatrix(P(1, 1, jrp), w_grp%mydim, n_intp_r, dbgunit)
          write (dbgunit, '("Q = ")')
          call printmatrix(Q(1, 1, jrp), w_grp%mydim, n_intp_r, dbgunit)
          write (dbgunit, '("P_intp = ")')
          call printmatrix(P_intp(1, 1, jrp), n_intp_r, n_intp_r, dbgunit)
          write (dbgunit, '("Q_intp = ")')
          call printmatrix(Q_intp(1, 1, jrp), n_intp_r, n_intp_r, dbgunit)
        end if
#endif
        !
      end do ! jrp loop
      !
#ifdef DEBUG
      write (6, *) "inode: ", peinf%inode, ", Calculate Cmtrx"
#endif
      !
      do jrp = 1, gvec%syms%ntrans
        !
        ! Prepare Cmtrx
        !
        do icv = 1, isdf_in%ncv(isp, ikp, jrp)
          IVV = isdf_in%invpairmap(1, icv, isp, ikp, jrp)
          ivrp = kpt%wfn(isp, ikp)%irep(IVV)
          ICC = isdf_in%invpairmap(2, icv, isp, ikp, jrp)
          icrp = kpt%wfn(isp, ikp)%irep(ICC)
          isdf_in%Cmtrx(1:n_intp_r, icv, isp, ikp, jrp) = &
            PsiV_intp(1:n_intp_r, inv_ivlist(IVV, isp, ikp, ivrp), ivrp)* &
            PsiC_intp(1:n_intp_r, inv_iclist(ICC, isp, ikp, icrp), icrp) ! element-wise multiplication
        end do
#ifdef DEBUG
        !if (peinf%master) then
        !   write(dbgunit, '(a,i5,a,i5,a,i5)') " jrp ", jrp, " isp ", isp, " ikp ", ikp
        !   write(dbgunit, *) " Cmtrx = "
        !   do icv = 1, isdf_in%ncv(isp,ikp,jrp)
        !      do jj = 1, n_intp_r
        !         write(dbgunit, '(2i7,e43.35)') icv, jj, isdf_in%Cmtrx(jj, icv, isp, ikp, jrp)
        !      enddo
        !   enddo
        !endif
#endif
      end do ! jrp
      !
      ! Calculate zeta(r,n_intp_r,jrp) for all representations
      !
#ifdef DEBUG
      write (6, *) "inode: ", peinf%inode, " Calculate zeta "
#endif
      do jrp = 1, gvec%syms%ntrans/r_grp%num
        irp = r_grp%g_rep(jrp)
        Amtrx = zero
        Bmtrx = zero
        Xmtrx = zero
        ! Calculate A and B
        !
        ! another way to calculate Amtrx
        ! Loop over all the reps.
        do lrp1 = 1, gvec%syms%ntrans
          !
          ! The direct product of lrp1 and lrp2 should be irp: lrp1 * lrp2 = irp
          !
          do lrp2 = 1, gvec%syms%ntrans
            if (gvec%syms%prod(lrp1, lrp2) == irp) exit
          end do
#ifdef DEBUG
          if (peinf%master .and. verbose) then
            write (dbgunit, *) "lrp1 ", lrp1, ", lrp2 ", lrp2, ", irp ", irp
          end if
#endif
          !
          ! ---------------------
          ! For each irp, spin and ikp, calculate zeta
          ! zeta dimension: Ng*Nu
          ! Set A = ( sum_{lrp1} P_intp(r_u,r_u,lrp1).Q_intp(r_u,r_u,lrp2) )^T  dimension: n_intp_r * n_intp_r
          !     B = ( sum_{lrp1} P(r,r_u,lrp1).Q(r,r_u,lrp2) )^T                dimension: n_intp_r * w_grp%mydim
          ! A is a symmetric matrix!
          ! we solve the linear equation A*(zeta^T) = B to get zeta^T for representation irp
          ! zeta^T = A^-1 * B = (C * C^T)^-1 * C * Z^T
          !
          !  Matrix dimensions:
          !  zeta(ngfl, n_intp_r, :, :)
          !  Amtrx(n_intp_r, n_intp_r)       intermediate variable, store C * C^T
          !  Bmtrx(n_intp_r, w_grp%mydim)  intermediate variable, store C * Z^T
          !  Xmtrx(n_intp_r, w_grp%mydim)         intermediate variable, store zeta^T
          ! ---------------------
          !
          ! calculate A = P_intp.Q_intp (Note: This is an element-wise multipliation)
          !
          Amtrx(1:n_intp_r, 1:n_intp_r) = Amtrx(1:n_intp_r, 1:n_intp_r) + &
                                          P_intp(1:n_intp_r, 1:n_intp_r, lrp1)* &
                                          Q_intp(1:n_intp_r, 1:n_intp_r, lrp2)
          !
          ! calculate B = (P.Q)^T (Note: This is an element-wise multiplication)
          !
          !
          ! Comment: this will crash with intel compiler
          !Bmtrx(1:n_intp_r,1:w_grp%mydim) = Bmtrx(1:n_intp_r,1:w_grp%mydim) + &
          !                     transpose( P(1:w_grp%mydim, 1:n_intp_r, lrp1) * &
          !                                Q(1:w_grp%mydim, 1:n_intp_r, lrp2) )
          ! End comment
          acc_B(1:w_grp%mydim, 1:n_intp_r) = P(1:w_grp%mydim, 1:n_intp_r, lrp1)* &
                                             Q(1:w_grp%mydim, 1:n_intp_r, lrp2)
          Bmtrx(1:n_intp_r, 1:w_grp%mydim) = Bmtrx(1:n_intp_r, 1:w_grp%mydim) + &
                                             transpose(acc_B(1:w_grp%mydim, 1:n_intp_r))
          !
#ifdef DEBUG
          if (verbose .and. peinf%master) then
            write (dbgunit, '(" P*Q = ")')
            call printmatrix(P(:, :, lrp1)*Q(:, :, lrp2), w_grp%mydim, n_intp_r, dbgunit)
          end if
#endif
        end do ! lrp1
        !
        ! solve linear equation A * X = B
        !
        !call dlinear_solver( isdf_in%n_intp_r, w_grp%mydim, Amtrx, Bmtrx, Xmtrx, w_grp%inode, verbose, einfo)
        Xmtrx(:, :) = Bmtrx(:, :)
        allocate (ipiv(n_intp_r))
        call dgesv(isdf_in%n_intp_r, w_grp%mydim, Amtrx, n_intp_r, ipiv, Xmtrx, &
                   n_intp_r, einfo)
        deallocate (ipiv)
        call MPI_BARRIER(peinf%comm, errinfo)
        !
        ! Copy Xmtrx to zeta
        !
        ! if ( verbose ) write(*,*) 'zeta', isp,ikp
        do ii = 1, isdf_in%n_intp_r
          do jj = 1, w_grp%mydim
            zeta(jj, ii, isp, irp) = Xmtrx(ii, jj)
          end do ! jj loop
        end do ! ii loop
#ifdef DEBUG
        if (verbose .and. peinf%master) then
          write (dbgunit, '(" irp =", i3, "  Zeta = ")') irp
          call printmatrix(zeta(1:137, 1:isdf_in%n_intp_r, isp, irp), 137, isdf_in%n_intp_r, dbgunit)
          call printmatrix(zeta(138:274, 1:isdf_in%n_intp_r, isp, irp), 137, isdf_in%n_intp_r, dbgunit)
          call printmatrix(zeta(275:411, 1:isdf_in%n_intp_r, isp, irp), 137, isdf_in%n_intp_r, dbgunit)
        end if
#endif
        !
      end do ! jrp loop
      !
      do jrp = 1, gvec%syms%ntrans
        !
        if (verbose) then
          !
#ifdef DEBUG
          if (peinf%master .and. verbose) then
            write (dbgunit, *) "isp", isp, " jrp", jrp, " ncv", &
              isdf_in%ncv(isp, ikp, jrp)
          end if
#endif
          do icv = 1, isdf_in%ncv(isp, ikp, jrp)
            call dgemv('n', w_grp%mydim, isdf_in%n_intp_r, &
                       one, zeta(1, 1, isp, jrp), w_grp%mydim, &
                       isdf_in%Cmtrx(1, icv, isp, ikp, jrp), 1, zero, tmp_Zmtrx, 1)
            !
            IVV = isdf_in%invpairmap(1, icv, isp, ikp, jrp)
            JVV = kpt%wfn(isp, ikp)%map(IVV)
            ICC = isdf_in%invpairmap(2, icv, isp, ikp, jrp)
            JCC = kpt%wfn(isp, ikp)%map(ICC)
            diff_vec(1:w_grp%mydim) = kpt%wfn(isp, ikp)%dwf(1:w_grp%mydim, IVV)* &
                                      kpt%wfn(isp, ikp)%dwf(1:w_grp%mydim, JCC) - tmp_Zmtrx(1:w_grp%mydim)
            diff = ddot(w_grp%mydim, diff_vec, 1, diff_vec, 1)
            weight = ddot(w_grp%mydim, tmp_Zmtrx, 1, tmp_Zmtrx, 1)
            call MPI_ALLREDUCE(MPI_IN_PLACE, diff, 1, MPI_DOUBLE, &
                               MPI_SUM, w_grp%comm, errinfo)
            call MPI_ALLREDUCE(MPI_IN_PLACE, weight, 1, MPI_DOUBLE, &
                               MPI_SUM, w_grp%comm, errinfo)
            diff = sqrt(diff)/sqrt(weight)
#ifdef DEBUG
            if (peinf%master .and. verbose) then
              if (diff > 0.09) write (dbgunit, '(a)', advance='no') " !! "
              write (dbgunit, *) " icv ", icv, " diff ", diff
              write (dbgunit, '(5e12.5)') tmp_Zmtrx(101:105)
            end if
#endif
          end do
        end if ! verbose
      end do ! jrp

    end do ! isp loop
    !
    call timacc(53, 2, tsec)
    !
#ifdef DEBUG
    write (6, *) "inode: ", peinf%inode, " start to calculate Mmtrx"
    write (6, *) "r_grp%num ", r_grp%num
#endif
    call timacc(54, 1, tsec)
    do jrp = 1, gvec%syms%ntrans/r_grp%num
      ! jrp is the index of representation belong to this r_grp
      ! irp is the real index of the representation
      !
      irp = r_grp%g_rep(jrp)
      !
      ! Now calculate <zeta_u(r,ispin)|V(r,r')|zeta_w(r',jspin)>, where u, w = 1, ..., n_intp
      !  and store it in Mmtrx(n_intp_r, n_intp_r)
      !
      do rsp = 1, nspin
        do csp = 1, nspin
#ifdef DEBUG
          write (6, *) "inode: ", peinf%inode, " jrp = ", jrp, " irp = ", irp, " rsp = ", rsp, " csp = ", csp
#endif
          !
          ! Exp: n_intp_r=1000, r_grp%npes=10, w_grp%npes=5
          !      then w_grp%mygr=0 work on ii = 1,2,3,4,5,  11,12,13,14,15, 21,22,23,24,25, ...
          !           w_grp%mygr=1 work on ii = 6,7,8,9,10, 16,17,18,19,20, 26,27,28,29,30, ...
          !
          do i_row = 1, isdf_in%n_intp_r, r_grp%npes
            !
            ! initialize rho_h and rho_h_local
            !
            rho_h = 0.d0
            rho_h_distr = 0.d0
            !
            do ipe = 0, w_grp%npes - 1
              ii = w_grp%mygr*w_grp%npes + i_row + ipe
              ! Note: w_grp%ldn = w_grp%mydim or w_grp%mydim+1
              ioff = w_grp%ldn*ipe + 1
              if (ii > isdf_in%n_intp_r) cycle
              ! Note:  ngrid = w_grp%mydim
              ! rho_h_distr are distributed over w_grp
              call dcopy(ngrid, &
                         zeta(1, ii, rsp, irp), 1, &
                         rho_h_distr(ioff), 1)
              if (kflag > 0) then
                call dcopy(ngrid, zeta(1, ii, rsp, irp), 1, &
                           fzeta(ioff), 1)
                ! calculate fzeta = f_lda * zeta
                call dmultiply_vec(ngrid, fxc(1, rsp, csp), fzeta(ioff))
              end if ! kflag > 0
            end do ! ipe loop
            !
            if (kflag < 2) then
              ii = w_grp%mygr*w_grp%npes + i_row + w_grp%inode
#ifdef DEBUG
              if (peinf%master) write (6, *) "ii ", ii
#endif
              call dgather(1, rho_h_distr, rho_h)
              if (ii <= isdf_in%n_intp_r) then
                !
                ! solve poisson equation to get: rho_h(r) = \int V_c(r,r') zeta_ii(r') dr'
                !
                call dpoisson(gvec, rho_h, irp)
              end if
              call dscatter(1, rho_h_distr, rho_h)
            end if ! kflag < 2
            !
            ! Mmtrx is a symmetric matrix
            !
            do ipe = 0, w_grp%npes - 1
              ! ii is real row index
              ii = w_grp%mygr*w_grp%npes + i_row + ipe
              ioff = w_grp%ldn*ipe + 1
              if (ii > isdf_in%n_intp_r) cycle
#ifdef DEBUG
              if (peinf%master) write (6, *) " calculate Mmtrx"
#endif
              do i_col = 1, isdf_in%n_intp_r
                !
                if (rsp == csp .and. i_col < ii) cycle
                !
                ! calculate: \int zeta_jj(r) V_c(r,r') zeta_ii(r') drdr'
                !          = \int zeta_jj(r) rho_h(r) dr
                !
                ! Note: ngrid = w_grp%mydim
#ifdef DEBUG
                !if ( peinf%master ) write(dbgunit, *) " ii ", ii, &
                !   " i_col", i_col
#endif
                if (kflag > 0) then
#ifdef DEBUG
                  if (peinf%master .and. ii == 1 .and. i_col == 1 .and. verbose) then
                    write (dbgunit, *) " test zeta fzeta irp = ", irp
                    do jj = 1, 10
                      write (dbgunit, *) jj, zeta(jj, 1, 1, irp), fzeta(ioff + jj - 1)
                    end do
                    do jj = 822, 831
                      write (dbgunit, *) jj, zeta(jj, 1, 1, irp), fzeta(ioff + jj - 1)
                    end do
                  end if
#endif
                  isdf_in%Mmtrx(ii, i_col, rsp, csp, ikp, 2, irp) = &
                    ddot(ngrid, zeta(1, i_col, csp, irp), 1, fzeta(ioff), 1)/gvec%hcub* &
                    gvec%syms%ntrans
                  !
                  ! Mmtrx is a symmetric matrix
                  !
                  isdf_in%Mmtrx(i_col, ii, csp, rsp, ikp, 2, irp) = &
                    isdf_in%Mmtrx(ii, i_col, rsp, csp, ikp, 2, irp)
                end if
                if (kflag < 2) then
#ifdef DEBUG
                  if (peinf%master .and. ii == 1 .and. i_col == 1 .and. verbose) then
                    write (dbgunit, *) " test zeta vzeta irp = ", irp
                    do jj = 1, 10
                      write (dbgunit, *) jj, zeta(jj, 1, 1, irp), rho_h_distr(ioff + jj - 1)
                    end do
                    do jj = 822, 831
                      write (dbgunit, *) jj, zeta(jj, 1, 1, irp), rho_h_distr(ioff + jj - 1)
                    end do
                  end if
#endif
                  isdf_in%Mmtrx(ii, i_col, rsp, csp, ikp, 1, irp) = &
                    ddot(ngrid, zeta(1, i_col, csp, irp), 1, rho_h_distr(ioff), 1)/gvec%hcub* &
                    gvec%syms%ntrans
                  !
                  ! Mmtrx is a symmetric matrix
                  !
                  isdf_in%Mmtrx(i_col, ii, csp, rsp, ikp, 1, irp) = &
                    isdf_in%Mmtrx(ii, i_col, rsp, csp, ikp, 1, irp)
                end if
                !
              end do ! i_col loop
              !
            end do ! ipe loop
            !
          end do ! i_row loop
          !
#ifdef DEBUG
          !if ( peinf%master .and. verbose ) then
          !   write(dbgunit,*) " jrp = ", jrp, " irp = ", irp, " rsp = ", rsp, " csp = ", csp
          !   write(dbgunit, '( "Mmtrx = " )' )
          !   call printmatrix ( isdf_in%Mmtrx( 1,1,rsp,csp,ikp,1,jrp ), &
          !     isdf_in%n_intp_r, isdf_in%n_intp_r, dbgunit)
          !endif
          !
#endif
        end do ! csp loop
        !
      end do ! rsp loop
      !
    end do ! jrp loop
    call timacc(54, 2, tsec)
    !
  end do ! ikp loop
  deallocate (PsiV)
  deallocate (PsiV_intp)
  deallocate (PsiC)
  deallocate (PsiC_intp)
  !
  ! clean up all the allocated variables
  !
  if (peinf%master .and. verbose) then
    write (*, *) " DEALLOCATING arrays"
  end if
  deallocate (Amtrx)
  deallocate (Bmtrx)
  deallocate (Xmtrx)
  deallocate (P)
  deallocate (P_intp)
  deallocate (Q)
  deallocate (Q_intp)
  deallocate (acc_B)
  deallocate (inv_ivlist)
  deallocate (inv_iclist)
  !
  if (kflag > 0) then
    deallocate (fxc)
    deallocate (fzeta)
  end if
  deallocate (zeta) ! no longer needed
  deallocate (rho_h)
  deallocate (rho_h_distr)
  !
  ! old code here --
  !if ( kflag > 0 ) then
  !   do jrp = 1, gvec%syms%ntrans
  !      call MPI_ALLREDUCE( MPI_IN_PLACE, isdf_in%Mmtrx(1,1,1,1,1,2,jrp), n_intp_r * n_intp_r * nspin * nspin * kpt%nk, &
  !        MPI_DOUBLE, MPI_SUM, w_grp%comm, errinfo )
  !   enddo
  !endif
  !if ( kflag < 2 ) then
  !   do jrp = 1, gvec%syms%ntrans
  !      call MPI_ALLREDUCE( MPI_IN_PLACE, isdf_in%Mmtrx(1,1,1,1,1,1,jrp), n_intp_r * n_intp_r * nspin * nspin * kpt%nk, &
  !        MPI_DOUBLE, MPI_SUM, w_grp%comm, errinfo )
  !   enddo
  !endif
  ! end old code ---
#ifdef DEBUG
  outdbg = 198812 + peinf%inode
  do jrp = 1, gvec%syms%ntrans
    write (outdbg, '(a,i2,a)') " Large Mmtrx (:, :, rsp=1, csp=1, ikp=1, 1, jrp=", jrp, ") = "
    do ii = 1, n_intp_r
      do jj = 1, n_intp_r
        if (abs(isdf_in%Mmtrx(ii, jj, 1, 1, 1, 1, jrp)) > 1e12) then
          write (outdbg, *) ii, jj, isdf_in%Mmtrx(ii, jj, 1, 1, 1, 1, jrp)
        end if
      end do ! ii
    end do ! jj
  end do
  do jrp = 1, gvec%syms%ntrans
    write (outdbg, '(a,i2,a)') " Large Mmtrx (:, :, rsp=1, csp=1, ikp=1, 2, jrp=", jrp, ") = "
    do ii = 1, n_intp_r
      do jj = 1, n_intp_r
        if (abs(isdf_in%Mmtrx(ii, jj, 1, 1, 1, 2, jrp)) > 1e12) then
          write (outdbg, *) ii, jj, isdf_in%Mmtrx(ii, jj, 1, 1, 1, 2, jrp)
        end if
      end do
    end do
  end do
#endif
  !
  ! new code here ---
  call MPI_ALLREDUCE(MPI_IN_PLACE, isdf_in%Mmtrx(1, 1, 1, 1, 1, 1, 1), n_intp_r*n_intp_r*nspin*nspin*kpt%nk*2*gvec%syms%ntrans, &
                     MPI_DOUBLE, MPI_SUM, w_grp%comm, errinfo)
  ! new code end here
  !
#ifdef DEBUG
  if (peinf%master .and. verbose) then
    do jrp = 1, gvec%syms%ntrans
      write (dbgunit, '(a,i2,a)') " Mmtrx (:, :, rsp=1, csp=1, ikp=1, 1, jrp=", jrp, ") = "
      call printmatrix(isdf_in%Mmtrx(1:isdf_in%n_intp_r, 1:isdf_in%n_intp_r, 1, 1, 1, 1, jrp), isdf_in%n_intp_r, isdf_in%n_intp_r, dbgunit)
    end do
    do jrp = 1, gvec%syms%ntrans
      write (dbgunit, '(a,i2,a)') " Mmtrx (:, :, rsp=1, csp=1, ikp=1, 2, jrp=", jrp, ") = "
      call printmatrix(isdf_in%Mmtrx(1:isdf_in%n_intp_r, 1:isdf_in%n_intp_r, 1, 1, 1, 2, jrp), isdf_in%n_intp_r, isdf_in%n_intp_r, dbgunit)
    end do
  end if
#endif
  !
  ! test calculating < 1 2 | F | 3 4 >
  if (peinf%master .and. verbose) then
    isp = 1
    ikp = 1
    irp = 1
#ifdef DEBUG
    write (dbgunit, *) " test calculate < 1 2 | V_coul | 3 4 > "
#endif
    tmpvec = zero
    do icv1 = 1, 10
      do icv2 = 1, 10
        call dgemv('n', isdf_in%n_intp_r, isdf_in%n_intp_r, one, &
                   isdf_in%Mmtrx(1, 1, isp, isp, ikp, 1, irp), isdf_in%n_intp_r, &
                   isdf_in%Cmtrx(1, icv1, isp, ikp, irp), 1, zero, tmpvec, 1)
        matel = ddot(isdf_in%n_intp_r, tmpvec, 1, &
                     isdf_in%Cmtrx(1, icv2, isp, ikp, irp), 1)
#ifdef DEBUG
        write (dbgunit, *) " icv1 ", icv1, " icv2 ", icv2, matel
#endif
      end do
    end do
#ifdef DEBUG
    write (dbgunit, *) " test calculate < 1 2 | f_xc | 3 4 > "
#endif
    tmpvec = zero
    do icv1 = 1, 10
      do icv2 = 1, 10
        call dgemv('n', isdf_in%n_intp_r, isdf_in%n_intp_r, one, &
                   isdf_in%Mmtrx(1, 1, isp, isp, ikp, 2, irp), isdf_in%n_intp_r, &
                   isdf_in%Cmtrx(1, icv1, isp, ikp, irp), 1, zero, tmpvec, 1)
        matel = ddot(isdf_in%n_intp_r, tmpvec, 1, &
                     isdf_in%Cmtrx(1, icv2, isp, ikp, irp), 1)
#ifdef DEBUG
        write (dbgunit, *) " icv1 ", icv1, " icv2 ", icv2, matel
#endif
      end do
    end do
  end if
#ifdef DEBUG
  if (peinf%master .and. verbose) then
    close (dbgunit)
  end if
#endif
  return
  !
end subroutine isdf_parallel_sym
