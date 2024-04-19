
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
! Step 1.6: DeALLOCATE P and Q
! Step 2.0: solve a linear quation A.zeta=B to get zeta
! Step 3.0: Calculate Mmtrx(1) = <zeta(i)| V_coul(r,r') | zeta(j)> and
!       Mmtrx(2) = <zeta(i)| f^LDA(r) | zeta(j)>
! Step 3.5: Calculate C(r_u,i,j)

! How to reduce memory requirements:
! [ done ]  1. Split zeta into n_slice, and calculate one slice each time, store zeta and
!    vc_zeta in hard drive using HDF5
! [      ]  2. deallocate all the wave functions after using it for
!    constructing PsiV and PsiC
! [      ]  3. store Mmtrx in hard drive using HDF5
! [ working on ]  4. don't calculate Cmtrx, only calculate and store Psi_intp 
subroutine isdf_parallel_sym_UltraLowMem ( gvec, pol_in, kpt, nspin, isdf_in, kflag, &
      opt, verbose )
 
#ifdef HIPMAGMA
  use magma
#endif
  use HDF5
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
  !
  ! n_intp_r: number of interpolation points in reduced r-space domain
  ! intp_r(n_intp_r): the index of interpolation points in full grids
  ! ncv(2): the number of wave function pairs for two spins
  !
  ! invpairmap maps the index of pair |cv> to c and v
  integer, intent(in) :: nspin, kflag
  type(options), intent(in) :: opt
  ! If verbose is true, then print out additional debug information
  logical, intent(in) :: verbose
  !
  ! ------------------ Local variables ------------------
  ! P(gvec%mydim, n_intp_r, nspin, kpt%nk)
  ! Q(gvec%mydim, n_intp_r, nspin, kpt%nk)
  ! Cmtrx(n_intp_r, Nc*Nv, nspin, kpt%nk)
  ! Issue need to be addressed: if there are N wfn_grp, does it mean these matrix need to be
  ! calculated by all processors at the same time, or can we distribute the workload
  ! later ??
  !
  real(dp), allocatable ::               &
     PsiV(:,:), PsiV_intp(:,:,:),        &   ! PsiV: wfn on reduced domain
     PsiC(:,:), PsiC_intp(:,:,:),        &
     P(:,:,:), P_intp(:,:,:),            &   ! P on reduced domain 
     Q(:,:,:), Q_intp(:,:,:),            &   ! Q on reduced domain 
     zeta(:,:,:), vzeta(:,:,:),          &
     fxc(:,:,:), fzeta(:,:),             &
     tmp_array(:,:),                     &
     tmp_Psi_intp_loc(:,:),              &
     tmp_Mmtrx_loc(:,:),                 &
     ! matrices and vectors used for solving linear equations
     Amtrx(:,:,:), Bmtrx(:,:), Xmtrx(:,:), &
     rho_h(:), rho_h_distr(:,:), &
     Amtrx_bl(:), Bmtrx_bl(:)
  real(dp) :: qkt(3), tsec(2), norm_factor
  !real(dp) ::matel, tmpvec(isdf_in%n_intp_r), &
  !   tmpCmtrx1(isdf_in%n_intp_r), &
  !   tmpCmtrx2(isdf_in%n_intp_r)
  !integer, allocatable :: inv_ivlist(:,:,:,:), inv_iclist(:,:,:,:)
  integer, allocatable :: ipiv(:)
  ! counters and temporary integers
  integer :: ipt, ii, jj, iv, ic, icv, irp, &
     jrp, rsp, csp, i_row, i_col, lrp1, lrp2, & 
     IVV, ICC, JVV, JCC, isp, ikp, errinfo, ipe, einfo, &
     ivrp, icrp, n_intp_r, maxncv, & 
     maxnv, maxnc, ig, icv1, icv2, &
     status, mpi_status(MPI_STATUS_SIZE), n_row, &
     ldn_intp_r 
  ! Each processor in a w_grp store part of the wave function
  ! offset: index of grid point from which a processor start to store the wave function
  ! ncount: number of elements of wave functions that a processor stores
  !integer, dimension(0:w_grp%npes-1) :: offset, ncount
  ! temporary dummy variable
  integer, dimension(0:w_grp%npes-1) :: idum
  ! the number of grid points in irreducible wedge, ngr = gvec%nr
  integer :: ngr, ngrid 
  ! the number of full grid point, ngf = ngr * (# of sym operations)
  integer :: ngfl, iptf, iptr, ioff, ioff1, ioff2, rcond, rank, outdbg

  ! variables for debug and test of accuracy 
  character(50) :: dbg_filename = "isdf_dbg.dat"
  integer, parameter :: dbgunit = 20171130, amatunit = 20190220
  ! external functions
  real(dp), external :: ddot
  integer,  external :: numroc
  type(xc_type) :: xc_lda
  !
  ! tmporary hdf5 file "zeta_tmp_iproc.h5", where iproc is the index of a processor
  ! the following variables is used for temporarily storing zeta and Vcoul_zeta
  !
  character(len=40) :: h5filename ! File name
  character(len=7)  :: iproc_string 
  character(len=2)  :: isp_string, irp_string, jrp_string
  character(len=40) :: dset_zeta = "Zeta"
  character(len=40) :: dset_vzeta = "VcZeta"
  integer(HID_T) :: file_id, &  ! file identifier
    dset_zeta_id(nspin, gvec%syms%ntrans), &
    dset_vczeta_id(nspin, gvec%syms%ntrans), & ! dataset identifier
    dspace_zeta, dspace_vczeta, &             ! dataspace identifier
    subdspace
  integer(hsize_t) :: data_dims(2), subdim(2), shift(2), stride(2), block(2)
  integer :: h5err
  ! scalapack and blacs variables
  integer :: context_system, icntxt_1d, descb_1d(9), info, MB_, &
    mycol_1d, myn_intp_r, myrow_1d, ndim_ipiv, npcol_1d, nprow_1d, &
    RSRC_, desca_1d(9), desca_2d(9), descb_2d(9), nprow_2d, npcol_2d, &
    nbl_2d, ldrowA, ldcolA, myrow_2d, mycol_2d, icntxt_2d, ldrowB, ldcolB
  logical, parameter :: use_old_code = .false.
  !
  outdbg=198812+w_grp%inode
  if( w_grp%master ) then
     !
     ! write(*,*) "call isdf_parallel(), write debug info to ", dbg_filename
     !
     open(dbgunit, file=dbg_filename, form='formatted', status='replace')
     !
  endif
  if( w_grp%master ) then
     open(amatunit, file="amat.dat", form='formatted', status='replace')
  endif
  ! the number of real-space grid points stored in current proc
  ngrid= w_grp%mydim
  ! the number of real-space grid points in reduced real-space domain
  ngr = gvec%nr 
  ! the number of grid points of interpolation vectors zeta(:) stored in each processor
  ngfl = w_grp%mydim * gvec%syms%ntrans
  n_intp_r = isdf_in%n_intp_r
  myn_intp_r = w_grp%myn_intp_r
  ldn_intp_r = w_grp%ldn_intp_r
  maxncv   = isdf_in%maxncv
  maxnv    = isdf_in%maxnv
  maxnc    = isdf_in%maxnc
  !
  ! Initialize h5 FORTRAN interface
  !
  call h5open_f(h5err)
  !
  ! create h5 file here
  !
#ifdef DEBUG
  if (w_grp%master .and. verbose) write(6,*) "create h5 file"
#endif
  write(iproc_string, '(I7.7)') w_grp%inode
  h5filename = "zeta_tmp_"//iproc_string//".h5"
  call h5fcreate_f(h5filename, H5F_ACC_TRUNC_F, file_id, h5err)
  !
  ! each processor stores part of P, Q and zeta
  !
#ifdef DEBUG
  if ( verbose ) then
     do ipe = 0, w_grp%npes-1
        if (w_grp%inode .eq. ipe) then
           write(6, '(a)') " Index of interpolation points in full domain: "
           write(6, * ) "n_intp_r", n_intp_r
           write(6, '(5i14)') ( isdf_in%intp_r(ii), ii=1, n_intp_r )
        endif
        call MPI_BARRIER(w_grp%comm, errinfo)
     enddo
  endif
#endif
  !
  ! isdf_in%ivlist(:) maps the index of valence states used in
  !   calculation (i.e., stored in memory) to the real index of valence states
  !ALLOCATE( inv_ivlist(isdf_in%maxivv, nspin, kpt%nk, gvec%syms%ntrans) ) ! For now, only deal with confined system. So we assume kpt%nk=1, and there is no dependence on k here
  !ALLOCATE( inv_iclist(isdf_in%maxicc, nspin, kpt%nk, gvec%syms%ntrans) )
  ALLOCATE( tmp_array(w_grp%mydim, isdf_in%n_intp_r) )
  ! inv_ivlist(:) maps the real index of valence states
  !   to the index of valence states used in the calculation
  !inv_ivlist = 0
  !inv_iclist = 0
  !
  !do ikp = 1, kpt%nk
  !  do isp = 1, nspin
  !    do irp = 1, gvec%syms%ntrans
  !      do iv = 1, isdf_in%nv(isp,ikp,irp)
  !        IVV = isdf_in%ivlist(iv, isp, ikp, irp)
  !        inv_ivlist(IVV, isp, ikp, irp) = iv
  !      enddo
  !      do ic = 1, isdf_in%nc(isp,ikp,irp)
  !        ICC = isdf_in%iclist(ic, isp, ikp, irp)
  !        inv_iclist(ICC, isp, ikp, irp) = ic
  !      enddo
  !     enddo
  !  enddo
  !enddo
  !
  ! qkt is set to zero. This is only valid for
  !  tests of nonperiodic system, and should be updated later.
  !
  qkt = 0
  !
  ! kflag = 0 : calculate kernel K^x  ( Coulomb )
  !         1 :                  K^x + K^f ( Coulomb + 1/2 * F_xc )
  !         2 :                  K^f   ( 1/2 * F_xc )
  !         3 :                  K^t   ( 1/2 * F_xc for spin triplet )
  !
  ! Generate Coulomb potential
  !
  if (kflag < 2 ) then
     call dinitialize_FFT(w_grp%inode, fft_box)
     call dcreate_coul_0D(gvec%bdot, qkt, fft_box)
     if (w_grp%master) write(6,*) " finished creating FFT plans and coul_0D"
  endif
  !
  ! Calculate LDA kernel
  !
  if ( kflag > 0 ) then
     !
     ALLOCATE( fxc( w_grp%mydim, nspin, nspin ), stat = errinfo )
     fxc = zero
     ! 
     ! Copy the charge density to fxc
     !
     do isp = 1, nspin
        call dcopy( w_grp%mydim, kpt%rho(w_grp%offset+1, isp), 1, &
           fxc(1, isp, 1), 1 )
     enddo 
     call xc_init( nspin, XC_LDA_X, XC_LDA_C_PZ, 0, zero, one, .false., xc_lda )
     xc_lda%has_grad = .false.
     call fxc_get( xc_lda, nspin, w_grp%mydim, kflag, fxc )
     call xc_end( xc_lda )
     !
     write(6,*) "inode: ", w_grp%inode, ": finished initializing xc functional"
  endif
  !
  ! create dataspace and dataset here
  !
  rank = 2
  data_dims(1) = w_grp%mydim
  data_dims(2) = n_intp_r
  ! create dataspace
#ifdef DEBUG
  if ( w_grp%master ) write(6,*) "Create dataspace: dspace_zeta"
#endif
  call h5screate_simple_f(rank, data_dims, dspace_zeta, h5err)
  ! create dataspace for vc_zeta
#ifdef DEBUG
  if ( w_grp%master ) write(6,*) "Create dataspace: dspace_vczeta"
#endif
  call h5screate_simple_f(rank, data_dims, dspace_vczeta, h5err)
  ! create hdf5 datasets for zeta and V_c*zeta for each spin and 
  ! group representation
  do isp = 1, nspin
    write(isp_string, '(i2.2)') isp
    do jrp = 1, gvec%syms%ntrans
      write(jrp_string, '(i2.2)') jrp
      dset_zeta = "zeta_isp"//isp_string//"_jrp"//jrp_string
      ! create dataset for each spin and representation
      call h5dcreate_f(file_id, dset_zeta, H5T_NATIVE_DOUBLE, &
        dspace_zeta, dset_zeta_id(isp, jrp), h5err)
    enddo
    do jrp = 1, gvec%syms%ntrans/r_grp%num
       irp = r_grp%g_rep(jrp)
       write(irp_string, '(i2.2)') irp
       dset_vzeta = 'vczeta_isp'//isp_string//"_irp"//irp_string
       !
       call h5dcreate_f(file_id, dset_vzeta, H5T_NATIVE_DOUBLE, &
         dspace_vczeta, dset_vczeta_id(isp, irp), h5err)
    enddo ! jrp loop
  enddo
  !
  isdf_in%Mmtrx_loc = zero ! Dimension: Mmtrx_loc(w_grp%myn_intp_r, n_intp_r, nspin, nspin, kpt%nk, 2, gvec%syms%ntrans)
  !
  norm_factor = 1.d0/gvec%hcub * gvec%syms%ntrans
  ! Here we assume there is only 1 r_grp and 1 w_grp
  call blacs_get(-1, 0, context_system)
  icntxt_1d = context_system
  nprow_1d = 1
  npcol_1d = w_grp%npes
  call blacs_gridinit(icntxt_1d, 'c', nprow_1d, npcol_1d)
  call blacs_gridinfo(icntxt_1d, nprow_1d, npcol_1d, myrow_1d, mycol_1d)
  if (w_grp%inode .ne. mycol_1d) then
    call die("isdf_parallel_sym_UtraLowMem: inode is not equal to mycol_1d")
  endif
  ! SUBROUTINE descinit( DESC, M, N, MB, NB, IRSRC, ICSRC, ICTXT,
  !   $ LLD, INFO )
  ! Initialize the descriptor of the Amtrx in the 1d process grid
  ! Note: scalapack requires desca( mb_ ) .eq. desca( nb_ ) .eq. descab( mb_ )
  ! see source code
  if ( w_grp%npes .ge. 2 ) then
    call descinit(desca_1d, n_intp_r, n_intp_r, &
      w_grp%ldn_intp_r, w_grp%ldn_intp_r, 0, 0, icntxt_1d, n_intp_r, info)
    !print *, "desca_1d ", desca_1d(1:9), " info ", info
    ! Initialize the descriptor of the Bmtrx in the 1d process grid
    call descinit(descb_1d, n_intp_r, gvec%nr, &
      w_grp%ldn_intp_r, 1, 0, 0, icntxt_1d, n_intp_r, info)
    !print *, "descb_1d ", descb_1d(1:9), " info ", info
  else
    call descinit(desca_1d, n_intp_r, n_intp_r, &
      n_intp_r, n_intp_r, 0, 0, icntxt_1d, n_intp_r, info)
    !print *, "desca_1d ", desca_1d(1:9), " info ", info
    ! Initialize the descriptor of the Bmtrx in the 1d process grid
    call descinit(descb_1d, n_intp_r, gvec%nr, &
      n_intp_r, gvec%nr, 0, 0, icntxt_1d, n_intp_r, info)
    !print *, "descb_1d ", descb_1d(1:9), " info ", info
  endif

  RSRC_ = 7
  MB_ = 5

  ! Information for 2D process grid, need for solving linear equation
  nprow_2d = int(sqrt(real(w_grp%npes)+1e-6))
  do ii = nprow_2d, 1, -1
      if (mod(w_grp%npes, ii) == 0) exit
  enddo
  nprow_2d = ii
  npcol_2d = w_grp%npes/nprow_2d
  nbl_2d = min(opt%pdgesv_nbl2d, n_intp_r/(max(nprow_2d, npcol_2d)))
  nbl_2d = max(nbl_2d, 1)
  if (peinf%master) print *, " pdgesv blocksize ", nbl_2d, "x", nbl_2d
  icntxt_2d = context_system
  ! Initialize the SL context for 1d process grid
  call blacs_gridinit(icntxt_2d, 'c', nprow_2d, npcol_2d)
  ! Get the row and col indices for the current process
  call blacs_gridinfo(icntxt_2d, nprow_2d, npcol_2d, myrow_2d, mycol_2d)
  ldrowA = numroc(n_intp_r, nbl_2d, myrow_2d, 0, nprow_2d)
  ldcolA = numroc(n_intp_r, nbl_2d, mycol_2d, 0, npcol_2d)
  ldrowB = numroc(n_intp_r, nbl_2d, myrow_2d, 0, nprow_2d)
  ldcolB = numroc(gvec%nr , nbl_2d, mycol_2d, 0, npcol_2d)
  if (myrow_2d .ne. -1) then
    allocate(Amtrx_bl(ldrowA*ldcolA))
    allocate(Bmtrx_bl(ldrowB*ldcolB))
    call descinit(desca_2d, n_intp_r, n_intp_r, nbl_2d, nbl_2d, 0, 0, &
      icntxt_2d, ldrowA, info)
    call descinit(descb_2d, n_intp_r, gvec%nr,  nbl_2d, nbl_2d, 0, 0, &
      icntxt_2d, ldrowA, info)
  else
    desca_2d(2) = -1
    descb_2d(2) = -1
  endif

  ! assume kpt%wfn(isp,ikp)%nmem is the same for all isp and ikp
  allocate(tmp_Psi_intp_loc(ldn_intp_r, kpt%wfn(1,1)%nmem))
  allocate(tmp_Mmtrx_loc(ldn_intp_r, n_intp_r))
  tmp_Mmtrx_loc = 0.d0
  do ikp = 1, kpt%nk
#ifdef DEBUG
    if (w_grp%master) write(6,*) " ikp = ", ikp
#endif
    !
    ! compute zeta functions (interpolation vectors)
    !
    do isp = 1, nspin
      ! construct isdf%Psi_intp_loc
      call stopwatch(peinf%master, "Start distribute_intp_pts")
      call distribute_intp_pts( kpt%wfn(isp,ikp), isdf_in, &
        gvec%syms%ntrans, isp, ikp )
      call stopwatch(peinf%master, "Finish distribute_intp_pts")
      !do ipe = 0, w_grp%npes-1
      !  if (w_grp%inode .eq. ipe) then
      !    if (ipe .eq. 0) write(6, *) isp, ikp, " Psi_intp_loc "
      !    print *, " ipe = ", ipe
      !    call printmatrix( isdf_in%Psi_intp_loc(1,1,1,1), w_grp%myn_intp_r, &
      !    kpt%wfn(1,1)%nmem, 6 ) 
      !  endif
      !  call MPI_BARRIER(w_grp%comm, info)
      !enddo
      !
      ALLOCATE( P      (n_intp_r, w_grp%mydim,gvec%syms%ntrans ) )
      ALLOCATE( Q      (n_intp_r, w_grp%mydim,gvec%syms%ntrans ) )
      ALLOCATE( P_intp (n_intp_r, myn_intp_r, gvec%syms%ntrans ) )
      ALLOCATE( Q_intp (n_intp_r, myn_intp_r, gvec%syms%ntrans ) )
      ALLOCATE( Amtrx  (n_intp_r, myn_intp_r, gvec%syms%ntrans ) )
      ALLOCATE( Bmtrx  (n_intp_r, w_grp%mydim  ) )
      ALLOCATE( Xmtrx  (w_grp%mydim,  n_intp_r ) )
      !
      P_intp = zero
      Q_intp = zero
      Amtrx = zero
      !
      ! initialize matrices with zero
      !
      P = zero
      Q = zero
      call timacc(53,1,tsec)
      !
      ! The following loop calculate P(r,r_u,jrp), Q(r,r_u,jrp) for all representations 
      !
      !if(w_grp%master) then
      !  write(6,*) " kpt%wfn(isp,ikp)%map ", &
      !    kpt%wfn(isp,ikp)%map(:)
      !endif
      do jrp = 1, gvec%syms%ntrans
#ifdef DEBUG
        if (w_grp%master .and. verbose) write(6,*) " jrp = ", jrp
#endif
        !
        if (verbose .and. w_grp%master) &
           write(dbgunit, *) ' isp = ', isp, ', ikp = ', ikp
        !
        ! ldn_intp_r = max(myn_intp_r) among all the procs
        do ipe = 0, w_grp%npes-1
          if (w_grp%inode .eq. ipe) then
            tmp_Psi_intp_loc(1:myn_intp_r,1:kpt%wfn(isp,ikp)%nmem) = &
              isdf_in%Psi_intp_loc(1:myn_intp_r, 1:kpt%wfn(isp,ikp)%nmem, isp, ikp)
          endif
          !MPI_BCAST(BUFFER, COUNT, DATATYPE, ROOT, COMM, IERROR)
          !        <type> BUFFER(*)  
          !        INTEGER COUNT, DATATYPE, ROOT, COMM, IERROR
          call MPI_BCAST(tmp_Psi_intp_loc(1,1), ldn_intp_r*kpt%wfn(isp,ikp)%nmem, &
            MPI_DOUBLE, ipe, w_grp%comm, info)
          call timacc(78, 1, tsec)
          do iv = 1, isdf_in%nv(isp,ikp,jrp)
            IVV = isdf_in%ivlist(iv,isp,ikp,jrp)
            JVV = kpt%wfn(isp,ikp)%map(IVV)
            ! Compute P(n_intp_r, w_grp%mydim)
            !
            ! old code
            !
            !call dgemm( 'n', 't', w_grp%n_intp_end(ipe)-w_grp%n_intp_start(ipe)+1, &
            !  w_grp%mydim, 1, one, &
            !  tmp_Psi_intp_loc(1,JVV), ldn_intp_r, &
            !  kpt%wfn(isp,ikp)%dwf(1,JVV), w_grp%ldn, one, &
            !  P(w_grp%n_intp_start(ipe),1,jrp), n_intp_r, jrp )
            !call dgemm( 'n', 't', w_grp%n_intp_end(ipe)-w_grp%n_intp_start(ipe)+1, &
            !  myn_intp_r, 1, one, &
            !  tmp_Psi_intp_loc(1,JVV), ldn_intp_r, &
            !  isdf_in%Psi_intp_loc(1,JVV,isp,ikp), myn_intp_r, one, &
            !  P_intp(w_grp%n_intp_start(ipe),1,jrp), n_intp_r, jrp )
            !
            ! new code
            !
            call dgemm_hl( 'n', 't', w_grp%n_intp_end(ipe)-w_grp%n_intp_start(ipe)+1, &
              w_grp%mydim, 1, one, &
              tmp_Psi_intp_loc(1,JVV), ldn_intp_r, &
              kpt%wfn(isp,ikp)%dwf(1,JVV), w_grp%ldn, one, &
              P(w_grp%n_intp_start(ipe),1,jrp), n_intp_r, &
              opt%linear_algebra )
            call dgemm_hl( 'n', 't', w_grp%n_intp_end(ipe)-w_grp%n_intp_start(ipe)+1, &
              myn_intp_r, 1, one, &
              tmp_Psi_intp_loc(1,JVV), ldn_intp_r, &
              isdf_in%Psi_intp_loc(1,JVV,isp,ikp), myn_intp_r, one, &
              P_intp(w_grp%n_intp_start(ipe),1,jrp), n_intp_r, &
              opt%linear_algebra )
            !
          enddo ! iv loop
          !
          do ic = 1, isdf_in%nc(isp,ikp,jrp)
            ICC = isdf_in%iclist(ic,isp,ikp,jrp)
            JCC = kpt%wfn(isp,ikp)%map(ICC)
            ! Compute Q
            !
            ! old code
            !
            !call dgemm( 'n', 't', w_grp%n_intp_end(ipe)-w_grp%n_intp_start(ipe)+1, &
            !  w_grp%mydim, 1, one, &
            !  tmp_Psi_intp_loc(1,JCC), ldn_intp_r, &
            !  kpt%wfn(isp,ikp)%dwf(1,JCC), w_grp%ldn, one, &
            !  Q(w_grp%n_intp_start(ipe),1,jrp), n_intp_r )
            !call dgemm( 'n', 't', w_grp%n_intp_end(ipe)-w_grp%n_intp_start(ipe)+1, &
            !  myn_intp_r, 1, one, &
            !  tmp_Psi_intp_loc(1,JCC), ldn_intp_r, &
            !  isdf_in%Psi_intp_loc(1,JCC,isp,ikp), myn_intp_r, one, &
            !  Q_intp(w_grp%n_intp_start(ipe),1,jrp), n_intp_r )
            !
            ! new code
            !
            call dgemm_hl( 'n', 't', w_grp%n_intp_end(ipe)-w_grp%n_intp_start(ipe)+1, &
              w_grp%mydim, 1, one, &
              tmp_Psi_intp_loc(1,JCC), ldn_intp_r, &
              kpt%wfn(isp,ikp)%dwf(1,JCC), w_grp%ldn, one, &
              Q(w_grp%n_intp_start(ipe),1,jrp), n_intp_r, &
              opt%linear_algebra )
            call dgemm_hl( 'n', 't', w_grp%n_intp_end(ipe)-w_grp%n_intp_start(ipe)+1, &
              myn_intp_r, 1, one, &
              tmp_Psi_intp_loc(1,JCC), ldn_intp_r, &
              isdf_in%Psi_intp_loc(1,JCC,isp,ikp), myn_intp_r, one, &
              Q_intp(w_grp%n_intp_start(ipe),1,jrp), n_intp_r, &
              opt%linear_algebra )
          enddo ! ic loop
          call timacc(78, 2, tsec)
          !
        enddo
        !
        !if ( .true. .and. w_grp%master ) then
        !  write(dbgunit, *) "jrp =", jrp
        !  write(dbgunit, '("P = ")' ) 
        !  call printmatrix ( P(1,1,jrp), n_intp_r, w_grp%mydim, 6 )
        !  write(dbgunit, '("Q = ")' )             
        !  call printmatrix ( Q(1,1,jrp), n_intp_r, w_grp%mydim, 6 )
        !  write(dbgunit, '("P_intp = ")' ) 
        !  call printmatrix ( P_intp(1,1,jrp), n_intp_r, myn_intp_r, 6 )
        !  write(dbgunit, '("Q_intp = ")' ) 
        !  call printmatrix ( Q_intp(1,1,jrp), n_intp_r, myn_intp_r, 6 )
        !endif
        !
      enddo ! jrp loop
      call stopwatch(peinf%master,' Finish constructing P, Q, P_intp, Q_intp. ')
      !
      ! Calculate zeta(r,n_intp_r,jrp) for all representations
      !
      do jrp = 1, gvec%syms%ntrans/r_grp%num
        irp = r_grp%g_rep(jrp)
        Bmtrx = zero
        Xmtrx = zero
        ! Calculate A and B
        !
        ! Loop over all the reps.
        call timacc(80, 1, tsec)
        do lrp1 = 1, gvec%syms%ntrans
          !
          ! The direct product of lrp1 and lrp2 should be irp: lrp1 * lrp2 = irp
          !
          do lrp2 = 1, gvec%syms%ntrans
             if(gvec%syms%prod(lrp1,lrp2) == irp) exit
          enddo
          if (w_grp%master) then
             write(dbgunit,*) "lrp1 ", lrp1, ", lrp2 ", lrp2, ", irp ", irp
          endif
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
          !  Amtrx(n_intp_r, n_intp_r)     intermediate variable, store C * C^T
          !  Bmtrx(n_intp_r, w_grp%mydim)         intermediate variable, store C * Z^T
          !  Xmtrx(n_intp_r, w_grp%mydim)         intermediate variable, store zeta^T
          ! ---------------------
          !
          ! calculate A = P_intp.Q_intp (Note: This is an element-wise multipliation)
          !
          Amtrx(1:n_intp_r, 1:myn_intp_r,irp) = &
            Amtrx (1:n_intp_r, 1:myn_intp_r,irp) + &
              P_intp(1:n_intp_r, 1:myn_intp_r,lrp1) * &
              Q_intp(1:n_intp_r, 1:myn_intp_r,lrp2)
          !
          ! calculate B = (P.Q)^T (Note: This is an element-wise multiplication)
          !
          ! Comment: this will crash with intel compiler
          !Bmtrx(1:n_intp_r,1:w_grp%mydim) = Bmtrx(1:n_intp_r,1:w_grp%mydim) + &
          !                     transpose( P(1:w_grp%mydim, 1:n_intp_r, lrp1) * &
          !                                Q(1:w_grp%mydim, 1:n_intp_r, lrp2) )
          ! End comment
          Bmtrx(1:n_intp_r,1:w_grp%mydim) = &
             Bmtrx(1:n_intp_r,1:w_grp%mydim) + &
              P(1:n_intp_r, 1:w_grp%mydim, lrp1) * &
              Q(1:n_intp_r, 1:w_grp%mydim, lrp2)
          !
        enddo ! lrp1
        call timacc(80, 2, tsec)
        !
        ! solve linear equation A * X = B
        !
        !if (.true. .and. w_grp%master) then
        !  write(dbgunit, '(" irp =", i3, "  Amtrx = ")') irp
        !  call printmatrix ( Amtrx(1:n_intp_r, 1:myn_intp_r, irp), &
        !    n_intp_r,  myn_intp_r, dbgunit )
        !  write(dbgunit, '(" irp =", i3, "  Bmtrx = ")') irp
        !  call printmatrix ( Bmtrx(1:n_intp_r, 1:w_grp%mydim) , &
        !    n_intp_r, w_grp%mydim, dbgunit )
        !endif
        ! Note that dgesv will change Amtrx, so we need a temporary Amtrx1
        !Amtrx1(1:n_intp_r,1:n_intp_r) = Amtrx(1:n_intp_r,1:n_intp_r,irp)
        !call dgesv(isdf_in%n_intp_r, w_grp%mydim, Amtrx1(1,1), n_intp_r, ipiv, Bmtrx, &
        !  n_intp_r, einfo)
        ! Use Scalapack to solve the linear system
        call MPI_BARRIER(w_grp%comm, errinfo)
        if (w_grp%master) print *, " jrp ", jrp
        if (use_old_code) then
          ndim_ipiv = numroc(n_intp_r, w_grp%ldn_intp_r, myrow_1d, desca_1d(RSRC_), nprow_1d) + desca_1d(MB_) 
          print *, "ndim_ipiv ", ndim_ipiv
          ALLOCATE(ipiv(ndim_ipiv))
          !if (w_grp%master) print *, "gvec%nr = ", gvec%nr
          !if (w_grp%master) print *, "desca_1d ", desca_1d(1:9)
          !if (w_grp%master) print *, "LOCc(n_intp_r)", numroc(n_intp_r, 1, 0, 0, 1)
          call stopwatch(w_grp%master,"Call pdgesv " )
          call pdgesv(n_intp_r, gvec%nr, Amtrx(1,1,irp), 1, 1, desca_1d, ipiv, Bmtrx, &
            1, 1, descb_1d, info)
          call stopwatch(w_grp%master,"Finish pdgesv " )
        else
          call stopwatch(w_grp%master,"Call pdgemr2d to redistribute Amtrx to 2D grid" )
          call pdgemr2d(n_intp_r, n_intp_r, Amtrx(1,1,irp), 1, 1, desca_1d, Amtrx_bl, &
            1, 1, desca_2d, icntxt_1d)
          call stopwatch(w_grp%master,"Call pdgemr2d to redistribute Bmtrx to 2D grid" )
          call pdgemr2d(n_intp_r, gvec%nr , Bmtrx, 1, 1, descb_1d, Bmtrx_bl, &
            1, 1, descb_2d, icntxt_1d)
          call stopwatch(w_grp%master,"Call pdgesv " )
          ndim_ipiv = numroc(n_intp_r, nbl_2d, myrow_2d, desca_2d(RSRC_), nprow_2d) + desca_2d(MB_) 
          !print *, "ndim_ipiv ", ndim_ipiv
          if (myrow_2d .ne. -1 .and. mycol_2d .ne. -1) then
          ALLOCATE(ipiv(ndim_ipiv))
          call pdgesv(n_intp_r, gvec%nr, Amtrx_bl, 1, 1, desca_2d, ipiv, Bmtrx_bl, &
            1, 1, descb_2d, info)
          endif
          call stopwatch(w_grp%master,"Call pdgemr2d to redistribute Amtrx to 1D grid" )
          call pdgemr2d(n_intp_r, n_intp_r, Amtrx_bl, 1, 1, desca_2d, Amtrx(1,1,irp), &
            1, 1, desca_1d, icntxt_1d)
          call stopwatch(w_grp%master,"Call pdgemr2d to redistribute Bmtrx to 1D grid" )
          call pdgemr2d(n_intp_r, gvec%nr , Bmtrx_bl, 1, 1, descb_2d, Bmtrx, &
            1, 1, descb_1d, icntxt_1d)
        endif
        !pdgetrs( TRANS, N, NRHS, A, IA, JA, DESCA, IPIV, B,
        !      $                    IB, JB, DESCB, INFO )
        !call pdgetrs('n', n_intp_r, gvec%nr, Amtrx, 1, 1, desca_1d, ipiv, &
        !   Bmtrx, 1, 1, descb_1d, info)
        !if (w_grp%master) print *, info
#ifdef DEBUG
        write(6, *) " inode ", w_grp%inode, " info = ", info
#endif
        if (myrow_2d .ne. -1 .and. mycol_2d .ne. -1) DEALLOCATE(ipiv)
        Xmtrx(1:w_grp%mydim,1:n_intp_r) = transpose(Bmtrx(1:n_intp_r,1:w_grp%mydim)) ! probably we can remove Bmtrx here now
        ! for debug use
        !write(outdbg, *) " Large Xmtrx "
        !do ii = 1, w_grp%mydim
        !  do jj = 1, isdf_in%n_intp_r
        !    if ( abs(Xmtrx(ii, jj)) > 1.0e7 ) then
        !      write(outdbg, *) ii, jj, Xmtrx(ii, jj)
        !    endif
        !  enddo ! jj
        !enddo ! ii 
        if (.true. .and. w_grp%master) then
          write(dbgunit, '(" irp =", i3, "  Zeta = ")') irp
          call printmatrix ( Xmtrx(1:w_grp%mydim, 1:isdf_in%n_intp_r), &
            w_grp%mydim, &
            isdf_in%n_intp_r, dbgunit )
        endif
        !
        ! Copy Xmtrx to zeta
        !
        ! Xmtrx (w_grp%mydim, n_intp_r)
        tmp_array(1:w_grp%mydim, 1:isdf_in%n_intp_r) = &
          Xmtrx(1:w_grp%mydim, 1:isdf_in%n_intp_r)
        !call h5dwrite_f( dset_zeta_id(isp,irp), H5T_NATIVE_DOUBLE, &
        !  Xmtrx(1:w_grp%mydim,1:isdf_in%n_intp_r), &
        !  data_dims, h5err, subdspace, dspace_zeta )
        call h5dwrite_f( dset_zeta_id(isp,irp), H5T_NATIVE_DOUBLE, &
          tmp_array, data_dims, h5err )
        if (peinf%master) write(6,*) "done write hdf5"
        !
      enddo ! jrp loop, jrp = 1, gvec%syms%ntrans/r_grp%num
      !
    enddo ! isp loop
    call MPI_BARRIER(peinf%comm, errinfo)
    if (w_grp%master) write(6,*) "deallocate arrays"
    !
    ! DEALLOCATE arrays
    !
    DEALLOCATE(P)
    DEALLOCATE(Q)
    DEALLOCATE(P_intp)
    DEALLOCATE(Q_intP)
    DEALLOCATE(Amtrx)
    DEALLOCATE(Bmtrx)
    DEALLOCATE(Xmtrx)
    !
    call timacc(53,2,tsec)
    jj = 3
    call timacc(53, jj, tsec)
    if (peinf%master) write(6, '(i3, f15.6, f15.6, a, i12, a)') 53, tsec(1), tsec(2), " sec", jj, " construct P, Q, P_intp, Q_intp"
    jj = 3
    call timacc(78, jj, tsec)
    if (peinf%master) write(6, '(i3, f15.6, f15.6, a, i12, a)') 78, tsec(1), tsec(2), " sec", jj, " dgemm"
    jj = 3
    call timacc(80, jj, tsec)
    if (peinf%master) write(6, '(i3, f15.6, f15.6, a, i12, a)') 80, tsec(1), tsec(2), " sec", jj, " hadamard"
    !
!#ifdef DEBUG
    if (peinf%master) write(6,*) "inode: ", peinf%inode, " start to calculate Mmtrx"
    ! if (peinf%master) write(6,*) "r_grp%num ", r_grp%num
!#endif
    call timacc(54,1,tsec)
    if ( kflag < 2 ) then
      !
      ALLOCATE(rho_h(ngr))     ! note: gvec%nr is equal to w_grp%nr
      ALLOCATE(rho_h_distr(w_grp%ldn, w_grp%npes))
      do jrp = 1, gvec%syms%ntrans/r_grp%num
        ! jrp is the index of representation belong to this r_grp
        ! irp is the real index of the representation
        !
        irp = r_grp%g_rep(jrp)
        !
        ! Now calculate <zeta_u(r,ispin)|V(r,r')|zeta_w(r',jspin)>, where u, w = 1, ..., n_intp
        ! and store it in Mmtrx(n_intp_r, n_intp_r)
        !
        do rsp = 1, nspin
          !
          ! Exp: n_intp_r=1000, r_grp%npes=10, w_grp%npes=5
          !      then w_grp%mygr=0 work on ii = 1,2,3,4,5,  11,12,13,14,15, 21,22,23,24,25, ...
          !           w_grp%mygr=1 work on ii = 6,7,8,9,10, 16,17,18,19,20, 26,27,28,29,30, ...
          !
          ! create a new subdataspace
          rank = 2
          subdim(1) = w_grp%mydim
          subdim(2) = w_grp%npes
          call h5screate_simple_f(rank, subdim, subdspace, h5err)
          ! Warning: currently, only works for: 1 w_grp per each r_grp
          !          this means r_grp%npes == w_grp%npes
          do i_row = 1, isdf_in%n_intp_r, r_grp%npes
            !
            ! locate a slice of dataset
            !
            shift(1)  = 0
            shift(2)  = w_grp%mygr*w_grp%npes+i_row-1
            if ( shift(2)+1 > isdf_in%n_intp_r ) cycle
            if ( shift(2)+subdim(2) > isdf_in%n_intp_r ) then
              subdim(2) = isdf_in%n_intp_r-shift(2)
              call h5sclose_f(subdspace, h5err)
              call h5screate_simple_f(rank, subdim, subdspace, h5err)
            endif
            !
            ! read rho_h_distr from hdf5 file
            !
            stride = (/1,1/)
            block = (/1,1/)
            call h5sselect_hyperslab_f( dspace_zeta, H5S_SELECT_SET_F, &
               shift, subdim, h5err, stride, block )
            call h5dread_f( dset_zeta_id(rsp,irp), H5T_NATIVE_DOUBLE, &
               tmp_array(1:w_grp%mydim, 1:w_grp%npes), data_dims, &
               h5err, subdspace, dspace_zeta )
            rho_h_distr(1:w_grp%mydim, 1:w_grp%npes) = &
               tmp_array(1:w_grp%mydim, 1:w_grp%npes)
            if (w_grp%master) write(dbgunit,*) " i_row ", i_row, ": done h5read_3"
            !
            ! initialize rho_h
            !
            rho_h = 0.d0
            !
            ii = w_grp%mygr*w_grp%npes + i_row + w_grp%inode
            call dgather(1,rho_h_distr,rho_h)
            call timacc(79, 1, tsec)
            if ( ii <= isdf_in%n_intp_r ) then
              !
              ! solve poisson equation to get: rho_h(r) = \int V_c(r,r') zeta_ii(r') dr'
              !
              call dpoisson(gvec, rho_h, irp)
            endif
            call timacc(79, 2, tsec)
            call dscatter(1,rho_h_distr,rho_h)
            !
            ! subdim(1) = w_grp%mydim is the leading dimension of tmp_array
            tmp_array(1:subdim(1), 1:subdim(2)) = rho_h_distr(1:subdim(1), 1:subdim(2))
            call h5sselect_hyperslab_f(dspace_vczeta, H5S_SELECT_SET_F, &
              shift, subdim, h5err, stride, block)
            call h5dwrite_f(dset_vczeta_id(rsp,irp), H5T_NATIVE_DOUBLE, &
              tmp_array(1:subdim(1), 1:subdim(2)), data_dims, &
              h5err, subdspace, dspace_vczeta)
            !if(w_grp%master) write(6,*) "done write hdf5"
            !
          enddo ! i_row loop
          !
        enddo ! rsp loop
        !
      enddo ! jrp loop
      !
      DEALLOCATE(rho_h)
      DEALLOCATE(rho_h_distr)
      !
    endif ! kflag < 2 
    !
    if (peinf%master) write(6,*) "inode: ", w_grp%inode, " allocate zeta" 
    ALLOCATE ( zeta(w_grp%mydim, isdf_in%n_intp_r, nspin) )
    if (kflag>0) ALLOCATE ( fzeta(w_grp%mydim, isdf_in%n_intp_r) )
    if (kflag<2) ALLOCATE ( vzeta(w_grp%mydim, isdf_in%n_intp_r, nspin) )
    !
    do jrp = 1, gvec%syms%ntrans/r_grp%num
      irp = r_grp%g_rep(jrp)
      do isp = 1, nspin
        ! open dataspace
        call h5sclose_f(dspace_zeta, h5err)
        call h5dget_space_f( dset_zeta_id(isp,irp), dspace_zeta, h5err )
        if(kflag<2) then
          call h5sclose_f(dspace_vczeta, h5err)
          call h5dget_space_f( dset_vczeta_id(isp,irp), dspace_vczeta, h5err )
        endif
      enddo

      do isp = 1, nspin
        !call h5dread_f( dset_zeta_id(isp, irp), H5T_NATIVE_DOUBLE, &
        !   zeta( 1:w_grp%mydim,1:isdf_in%n_intp_r, isp ), data_dims, h5err, subdspace, dspace_zeta )
        call h5dread_f( dset_zeta_id(isp, irp), H5T_NATIVE_DOUBLE, &
           tmp_array, data_dims, h5err )
        zeta( 1:w_grp%mydim, 1:isdf_in%n_intp_r, isp ) = tmp_array( 1:w_grp%mydim, 1:isdf_in%n_intp_r )
        !if (w_grp%master) write(dbgunit,*) "done h5read_1"
        !if (w_grp%master) then
        !  write(dbgunit, '("Read zeta() irp",i3," isp ",i3 )') irp, isp
        !  call printmatrix( zeta(1,1,isp), w_grp%mydim, isdf_in%n_intp_r, 6 )
        !endif
        if(kflag<2) then
        call h5dread_f( dset_vczeta_id(isp, irp), H5T_NATIVE_DOUBLE, &
           tmp_array, data_dims, h5err )
        vzeta( 1:w_grp%mydim, 1:isdf_in%n_intp_r, isp ) = tmp_array( 1:w_grp%mydim, 1:isdf_in%n_intp_r )
        endif
        !if (w_grp%master) write(dbgunit,*) "done h5read_2"
      enddo
      ! calculate mtrx now
      do rsp = 1, nspin
        do csp = 1, nspin
          do ipe = 0, w_grp%npes-1
          !print *, "inode", w_grp%inode, kflag
          if ( kflag > 0 ) then
            do ii = 1, isdf_in%n_intp_r
              fzeta(1:w_grp%mydim, ii) = &
                fxc(1:w_grp%mydim, rsp, csp) &
                 * zeta(1:w_grp%mydim, ii, csp)
            enddo ! ii loop
            if (w_grp%master) then
              write(dbgunit,*) " test zeta fzeta irp =", irp
              do jj = 1, 10
                write(dbgunit,*) jj-1, zeta(jj,1,1), fzeta(jj,1)
              enddo
            endif

            n_row = w_grp%n_intp_end(ipe)-w_grp%n_intp_start(ipe)+1
            ! old code
            !call dgemm('T', 'N', n_row, n_intp_r, w_grp%mydim, &
            !  norm_factor, zeta(1,w_grp%n_intp_start(ipe),rsp), w_grp%mydim, &
            !  fzeta(1,1), w_grp%mydim, 1.d0, &
            !  tmp_Mmtrx_loc(1, 1), ldn_intp_r)
            ! new code
            call timacc(77, 1, tsec)
            call dgemm_hl('T', 'N', n_row, n_intp_r, w_grp%mydim, &
               norm_factor, zeta(1,w_grp%n_intp_start(ipe),rsp), w_grp%mydim, &
               fzeta(1,1), w_grp%mydim, 1.d0, &
               tmp_Mmtrx_loc(1, 1), ldn_intp_r, opt%linear_algebra)
            call timacc(77, 2, tsec)
            call MPI_REDUCE(MPI_IN_PLACE, tmp_Mmtrx_loc(1,1), n_row*n_intp_r, &
              MPI_DOUBLE, MPI_SUM, ipe, w_grp%comm, errinfo)
            if (w_grp%inode .eq. ipe) then
              isdf_in%Mmtrx_loc(1:n_row, 1:n_intp_r, rsp, csp, ikp, 2, irp) = & 
               tmp_Mmtrx_loc(1:n_row, 1:n_intp_r)
            endif
          endif ! kflag > 0
          if ( kflag < 2 ) then
            if (w_grp%master) then
              write(dbgunit,*) " test zeta vzeta irp =", irp
              do jj = 1, 10
                write(dbgunit,*) jj-1, zeta(jj,1,1), vzeta(jj,1,1)
              enddo
            endif

            n_row = w_grp%n_intp_end(ipe)-w_grp%n_intp_start(ipe)+1
            !write(6000+w_grp%inode, *) "zeta = "
            !call printmatrix ( zeta(1,1,1), &
            !  w_grp%mydim, n_intp_r, 6000+w_grp%inode )
            !write(6000+w_grp%inode, *) "vzeta = "
            !call printmatrix ( vzeta(1,1,1), &
            !  w_grp%mydim, n_intp_r, 6000+w_grp%inode )

            ! old code
            !call dgemm('T', 'N', n_row, n_intp_r, w_grp%mydim, &
            !  norm_factor, zeta(1,w_grp%n_intp_start(ipe),rsp), w_grp%mydim, &
            !  vzeta(1,1,csp), w_grp%mydim, 0.d0, &
            !  tmp_Mmtrx_loc(1, 1), ldn_intp_r)
            ! new code
            call timacc(77, 1, tsec)
            call dgemm_hl('T', 'N', n_row, n_intp_r, w_grp%mydim, &
              norm_factor, zeta(1,w_grp%n_intp_start(ipe),rsp), w_grp%mydim, &
              vzeta(1,1,csp), w_grp%mydim, 0.d0, &
              tmp_Mmtrx_loc(1, 1), ldn_intp_r, opt%linear_algebra)
            call timacc(77, 2, tsec)
            !write(6000+w_grp%inode, *) "tmp_Mmtrx_loc = "
            !call printmatrix ( tmp_Mmtrx_loc (1,1), &
            !  ldn_intp_r, n_intp_r, 6000+w_grp%inode )
            !call MPI_REDUCE(MPI_IN_PLACE, tmp_Mmtrx_loc(1,1), ldn_intp_r*n_intp_r, &
            !  MPI_DOUBLE, MPI_SUM, ipe, w_grp%comm, errinfo)
            call MPI_REDUCE(tmp_Mmtrx_loc(1,1), &
              isdf_in%Mmtrx_loc(1,1,rsp,csp,ikp,1,irp), n_row*n_intp_r, &
              MPI_DOUBLE, MPI_SUM, ipe, w_grp%comm, errinfo)
            !if (w_grp%inode .eq. ipe) then
            !  isdf_in%Mmtrx_loc(1:n_row, 1:n_intp_r, rsp, csp, ikp, 1, irp) = & 
            !   tmp_Mmtrx_loc(1:n_row, 1:n_intp_r)
            !endif
            !if (ipe .eq. w_grp%inode) then
            !write(6000+w_grp%inode, *) "isdf_in%Mmtrx_loc = "
            !call printmatrix ( isdf_in%Mmtrx_loc (1,1,rsp,csp,ikp,1,irp), &
            !  ldn_intp_r, n_intp_r, 6000+w_grp%inode )
            !endif
            !
          endif ! kflag < 2
          enddo ! ipe loop
        enddo ! csp
      enddo ! rsp 
    enddo ! jrp loop
    call timacc(54,2,tsec)
    jj = 3
    call timacc(54, jj, tsec)
    if (peinf%master) write(6, '(i3, f15.6, f15.6, a, i12, a)') 54, tsec(1), tsec(2), " sec", jj, " building Mmtrx "
    jj = 3
    call timacc(77, jj, tsec)
    if (peinf%master) write(6, '(i3, f15.6, f15.6, a, i12, a)') 77, tsec(1), tsec(2), " sec", jj, " dgemm "
    jj = 3
    call timacc(78, jj, tsec)
    if (peinf%master) write(6, '(i3, f15.6, f15.6, a, i12, a)') 79, tsec(1), tsec(2), " sec", jj, " dpoisson "
    !
    DEALLOCATE( zeta  )
    if (kflag > 0) DEALLOCATE( fzeta )
    if (kflag < 2) DEALLOCATE( vzeta )
  enddo ! ikp loop

  if (myrow_2d .ne. -1) then 
    call blacs_gridexit(icntxt_2d)
    DEALLOCATE(Amtrx_bl)
    DEALLOCATE(Bmtrx_bl)
  endif
  if (myrow_1d .ne. -1) call blacs_gridexit(icntxt_1d)
  ! close dataspaces
  call h5sclose_f( dspace_zeta, h5err )
  call h5sclose_f( dspace_vczeta, h5err )
  call h5sclose_f( subdspace, h5err )
  ! close dataset
  do isp = 1, nspin
    do jrp = 1, gvec%syms%ntrans
      call h5dclose_f(dset_zeta_id(isp, jrp), h5err)
    enddo
    do jrp = 1, gvec%syms%ntrans/r_grp%num
      irp = r_grp%g_rep(jrp)
      call h5dclose_f(dset_vczeta_id(isp, irp), h5err)
    enddo ! jrp loop
  enddo
  ! close file
  call h5fclose_f(file_id, h5err)
  ! close FORTRAN interface
  call h5close_f(h5err)
  !
  ! clean up all the ALLOCATEd variables
  !
  !DEALLOCATE(inv_ivlist)
  !DEALLOCATE(inv_iclist)
  if ( kflag > 0 ) then
    DEALLOCATE(fxc)
  endif
  !
  if ( w_grp%master .and. .True. ) then
    do jrp = 1, gvec%syms%ntrans
      write( dbgunit, '(a,i2,a)' ) " Mmtrx_loc (:, :, rsp=1, csp=1, ikp=1, 1, jrp=", jrp, ") = "
      call printmatrix ( isdf_in%Mmtrx_loc (1,1,1,1,1,1,jrp), &
        w_grp%myn_intp_r, n_intp_r, dbgunit )
    enddo
    do jrp = 1, gvec%syms%ntrans
      write( dbgunit, '(a,i2,a)' ) " Mmtrx_loc (:, :, rsp=1, csp=1, ikp=1, 2, jrp=", jrp, ") = "
      call printmatrix ( isdf_in%Mmtrx_loc (1,1,1,1,1,2,jrp), &
        w_grp%myn_intp_r, n_intp_r, dbgunit )
    enddo
  endif
  DEALLOCATE(tmp_Mmtrx_loc)
  DEALLOCATE(tmp_Psi_intp_loc)
  DEALLOCATE(tmp_array)
  !
  if ( w_grp%master ) then
     close ( dbgunit )
  endif
  if ( w_grp%master ) then
     close ( amatunit )
  endif
  return
  !
end subroutine isdf_parallel_sym_UltraLowMem

! prepare isdf%Psi_intp_loc
subroutine distribute_intp_pts(wfn, isdf_in, ntrans, isp, ikp)

  use mpi_module
  use typedefs
  implicit none 
#ifdef MPI
  include 'mpif.h'
#endif
  type(ISDF), intent(inout) :: isdf_in
  type(wavefunction), intent(in) :: wfn
  integer, intent(in) :: ntrans, isp, ikp
  real(dp) :: dummy(wfn%nmem)
  integer, parameter :: nchecks = 1000
  integer :: inode_dest, inode_source, &
    ipt, iptf, iptr, ii, jj, mpi_err, tag
  INTEGER(KIND=MPI_ADDRESS_KIND) :: max_tag
  logical :: flag
  integer :: ioff1 (0:w_grp%npes-1), ioff2 (0:w_grp%npes-1), &
    mpistatus(MPI_STATUS_SIZE)

  call MPI_Comm_get_attr(w_grp%comm, MPI_TAG_UB, max_tag, &
    flag, mpi_err)
  if ( w_grp%master ) print *, "MPI MAX_TAG = ", MAX_TAG
  ioff1 = 0
  ioff2 = 0
  ioff1(w_grp%inode) = w_grp%offset+1
  ioff2(w_grp%inode) = w_grp%offset+w_grp%mydim
  call MPI_ALLREDUCE(MPI_IN_PLACE, ioff1, w_grp%npes, MPI_INTEGER, &
    MPI_SUM, w_grp%comm, mpi_err)
  call MPI_ALLREDUCE(MPI_IN_PLACE, ioff2, w_grp%npes, MPI_INTEGER, &
    MPI_SUM, w_grp%comm, mpi_err)
  !if (w_grp%master) then
  !  print *, " ioff1 ", ioff1(0:w_grp%npes-1)
  !  print *, " ioff2 ", ioff2(0:w_grp%npes-1)
  !endif
  inode_dest = 0 ! Start from the first proc
  !print *, "w_grp%ldn_intp_r", w_grp%ldn_intp_r
  !print *
  do ipt = 1, isdf_in%n_intp_r
    !if (ipt.eq.1)  print *, " w_grp%n_intp_end(inode_dest) ", &
    !  w_grp%n_intp_end(inode_dest)
    if (ipt .gt. w_grp%n_intp_end(inode_dest)) then
      inode_dest = inode_dest + 1 ! update the destination proc
    endif ! 
    ! determine whether this interpolation point
    ! is located in the current process
    iptf = isdf_in%intp_r(ipt)
    if ( ntrans .gt. 1 ) then
      iptr = iptf / ntrans + 1
    else
      iptr = iptf
    endif
    inode_source = iptr / (w_grp%ldn-1)
    ! make sure the starting guess of inode_source is not out-of-bound
    do ii = 1, nchecks
      if (inode_source .ge. w_grp%npes-1) then
        inode_source = inode_source - 1
      else 
        exit 
      endif
    enddo
    ! look for the inode_source
    do ii = 1, nchecks
      if ( iptr .gt. ioff2(inode_source) ) then
          inode_source = inode_source + 1
          continue 
      endif
      if ( iptr .lt. ioff1(inode_source) ) then
          inode_source = inode_source - 1
          continue
      endif
      if ( inode_source .gt. w_grp%npes - 1 .or. &
           inode_source .lt. 0 ) then
        if ( w_grp%master ) print *, " ipt ", ipt, " iptf ", iptf, " iptr ", iptr, &
          " ioff2(w_grp%npes-1) ", ioff2(w_grp%npes-1)
        call die("iptr is out of bound, see sigma.out or tdlda.out.")
      endif
      if ( iptr .ge. ioff1(inode_source) .and. &
           iptr .le. ioff2(inode_source) ) then
          exit
      endif
    enddo
    if ( w_grp%inode .eq. inode_source ) then 
      ! This interpolation point is located in the
      ! current process
      jj = iptr - ioff1(inode_source) + 1
      if (inode_source .eq. inode_dest) then
         ! w_grp%inode .eq. inode_dest
         isdf_in%Psi_intp_loc(ipt-w_grp%n_intp_start(inode_dest)+1, &
           1:wfn%nmem, isp, ikp) &
           = wfn%dwf(jj, 1:wfn%nmem)
      else ! source and dest process are different  
         dummy(1:wfn%nmem) = wfn%dwf(jj,1:wfn%nmem)
         tag = ipt
         call MPI_SEND( dummy(1), &
           wfn%nmem, MPI_DOUBLE, inode_dest, tag, &
           w_grp%comm, mpi_err )
         ! print *, mpi_err
      endif 
    endif
    if ( w_grp%inode .eq. inode_dest .and. &
         inode_dest .ne. inode_source ) then
      ! if the current process is the destination process
      ! and the source process is not the destination process
      tag = ipt
      call MPI_RECV( dummy(1), &
        wfn%nmem, MPI_DOUBLE, inode_source, tag, &
        w_grp%comm, mpistatus, mpi_err )
      ! print *, mpi_err
      isdf_in%Psi_intp_loc(ipt-w_grp%n_intp_start(inode_dest)+1, &
        1:wfn%nmem, isp, ikp) &
        = dummy(1:wfn%nmem)
    endif
  enddo ! ipt

end subroutine 

