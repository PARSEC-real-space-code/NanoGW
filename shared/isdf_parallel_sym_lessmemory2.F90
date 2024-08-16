
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
subroutine isdf_parallel_sym_lessmemory2 ( gvec, pol_in, kpt, nspin, isdf_in, kflag, &
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
     acc_B(:,:), tmp_array(:,:),         &
     ! matrices and vectors used for solving linear equations
     Amtrx(:,:,:), Bmtrx(:,:), Xmtrx(:,:), &
     rho_h(:), rho_h_distr(:,:), Amtrx1(:,:)
  real(dp) :: qkt(3), tsec(2), &
     tmpvec(isdf_in%n_intp_r), matel, norm_factor, &
     tmpCmtrx1(isdf_in%n_intp_r), &
     tmpCmtrx2(isdf_in%n_intp_r)
  integer, allocatable :: inv_ivlist(:,:,:,:), inv_iclist(:,:,:,:), ipiv(:)
  ! counters and temporary integers
  integer :: ipt, ii, jj, iv, ic, icv, irp, jrp, rsp, csp, i_row, i_col, lrp1, lrp2, & 
             IVV, ICC, JVV, JCC, isp, ikp, errinfo, ipe, einfo, &
             ivrp, icrp, n_intp_r, maxncv, & 
             maxnv, maxnc, ig, icv1, icv2, &
             status, mpi_status(MPI_STATUS_SIZE)
  ! Each processor in a w_grp store part of the wave function
  ! offset: index of grid point from which a processor start to store the wave function
  ! ncount: number of elements of wave functions that a processor stores
  integer, dimension(0:w_grp%npes-1) :: offset, ncount
  ! temporary dummy variable
  integer, dimension(0:w_grp%npes-1) :: idum
  ! the number of grid points in irreducible wedge, ngr = gvec%nr
  integer :: ngr, ngrid 
  ! the number of full grid point, ngf = ngr * (# of sym operations)
  integer :: ngfl, iptf, iptr, ioff, ioff1, ioff2, rcond, rank

#ifdef DEBUG
  ! variables for debug and test of accuracy 
  character(50) :: dbg_filename = "isdf_dbg.dat"
  integer, parameter :: dbgunit = 20171130
  integer :: outdbg
#endif

  ! external functions
  real(dp), external :: ddot
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

#ifdef DEBUG
  outdbg=198812+peinf%inode
  if( peinf%master ) then
     ! write(*,*) "call isdf_parallel(), write debug info to ", dbg_filename
     open(dbgunit, file=dbg_filename, form='formatted', status='replace')
  endif
#endif

  ! the number of real-space grid points stored in current proc
  ngrid= w_grp%mydim
  ! the number of real-space grid points in reduced real-space domain
  ngr = gvec%nr 
  ! the number of grid points of interpolation vectors zeta(:) stored in each processor
  ngfl = w_grp%mydim * gvec%syms%ntrans
  n_intp_r = isdf_in%n_intp_r
  maxncv   = isdf_in%maxncv
  maxnv    = isdf_in%maxnv
  maxnc    = isdf_in%maxnc
  
  ! Initialize h5 FORTRAN interface
  !
  call h5open_f(h5err)
  !
  ! create h5 file here
  !
#ifdef DEBUG
  if (peinf%master .and. verbose) write(6,*) "create h5 file"
#endif
  write(iproc_string, '(I7.7)') peinf%inode
  h5filename = "zeta_tmp_"//iproc_string//".h5"
  call h5fcreate_f(h5filename, H5F_ACC_TRUNC_F, file_id, h5err)
  !
  ! each processor stores part of P, Q and zeta
  !
#ifdef DEBUG
  if ( verbose ) then
     do ipe = 1, peinf%npes
        if (peinf%inode .eq. ipe) then
           write(6, '(a)') " Index of interpolation points in full domain: "
           write(6, * ) "n_intp_r", n_intp_r
           write(6, '(5i14)') ( isdf_in%intp_r(ii), ii=1, n_intp_r )
        endif
        call MPI_BARRIER(peinf%comm, errinfo)
     enddo
  endif
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
  if ( verbose .and. peinf%master ) then
    write(dbgunit, *) " In isdf() "
    write(dbgunit, *) " w_grp%mygr ", w_grp%mygr, " w_grp%inode = ", w_grp%inode
    write(dbgunit, *) "    ii      offset(ii)       ncount(ii) "
    do ii = 0, w_grp%npes-1
      write(dbgunit, '(i7,2i9)') ii,offset(ii),ncount(ii)
    enddo
  endif
#endif
  !
  ! isdf_in%ivlist(:) maps the index of valence states used in
  !   calculation (i.e., stored in memory) to the real index of valence states
  ALLOCATE( inv_ivlist(isdf_in%maxivv, nspin, kpt%nk, gvec%syms%ntrans) ) ! For now, only deal with confined system. So we assume kpt%nk=1, and there is no dependence on k here
  ALLOCATE( inv_iclist(isdf_in%maxicc, nspin, kpt%nk, gvec%syms%ntrans) )
  ! inv_ivlist(:) maps the real index of valence states
  !   to the index of valence states used in the calculation
  inv_ivlist = 0
  inv_iclist = 0
  !
  do ikp = 1, kpt%nk
    do isp = 1, nspin
      !
      do irp = 1, gvec%syms%ntrans
        do iv = 1, isdf_in%nv(isp,ikp,irp)
          IVV = isdf_in%ivlist(iv, isp, ikp, irp)
          inv_ivlist(IVV, isp, ikp, irp) = iv
        enddo
        !
        do ic = 1, isdf_in%nc(isp,ikp,irp)
          ICC = isdf_in%iclist(ic, isp, ikp, irp)
          inv_iclist(ICC, isp, ikp, irp) = ic
        enddo
        !
      enddo
      !
    enddo
  enddo
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
     call dinitialize_FFT(peinf%inode, fft_box)
     call dcreate_coul_0D(gvec%bdot, qkt, fft_box)
     if (peinf%master) write(6,*) " finished creating FFT plans and coul_0D"
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
     write(6,*) "inode: ", peinf%inode, ": finished initializing xc functional"
  endif
  !
  ! create dataspace and dataset here
  !
  rank = 2
  data_dims(1) = w_grp%mydim
  data_dims(2) = n_intp_r
  ! create dataspace
#ifdef DEBUG
  if ( peinf%master ) write(6,*) "Create dataspace: dspace_zeta"
#endif
  call h5screate_simple_f(rank, data_dims, dspace_zeta, h5err)
  ! create dataspace for vc_zeta
#ifdef DEBUG
  if ( peinf%master ) write(6,*) "Create dataspace: dspace_vczeta"
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
  isdf_in%Mmtrx = zero ! Mmtrx dimension: Mmtrx(n_intp_r, n_intp_r, nspin, nspin)
  !
  norm_factor = 1.d0/gvec%hcub * gvec%syms%ntrans
  do ikp = 1, kpt%nk
#ifdef DEBUG
    if (peinf%master) write(6,*) " ikp = ", ikp
#endif
    !
    ALLOCATE( P      (w_grp%mydim, n_intp_r, gvec%syms%ntrans ) )
    ALLOCATE( Q      (w_grp%mydim, n_intp_r, gvec%syms%ntrans ) )
    ALLOCATE( acc_B  (w_grp%mydim, n_intp_r ) )
    ALLOCATE( P_intp (n_intp_r,    n_intp_r, gvec%syms%ntrans ) )
    ALLOCATE( Q_intp (n_intp_r,    n_intp_r, gvec%syms%ntrans ) )
    ALLOCATE( PsiV   (w_grp%mydim, maxnv ) )
    ALLOCATE( PsiV_intp(n_intp_r,  maxnv, gvec%syms%ntrans ) )
    ALLOCATE( PsiC   (w_grp%mydim, maxnc ) )
    ALLOCATE( PsiC_intp(n_intp_r,  maxnc, gvec%syms%ntrans ) )
    ALLOCATE( Amtrx  (n_intp_r,    n_intp_r, gvec%syms%ntrans ) )
    ALLOCATE( Amtrx1 (n_intp_r,    n_intp_r ) )
    ALLOCATE( Bmtrx  (n_intp_r,    w_grp%mydim) )
    ALLOCATE( Xmtrx  (w_grp%mydim, n_intp_r ) )
    ! compute zeta functions (interpolation vectors)
    do isp = 1, nspin
      ! construct isdf%Psi_intp
      do ii = 1, kpt%wfn(isp,ikp)%nmem
        !
        ! pick the interpolation points
        !
        do ipt = 1, n_intp_r
          iptf = isdf_in%intp_r(ipt)
          if ( gvec%syms%ntrans > 1 ) then 
            iptr = iptf / gvec%syms%ntrans + 1
          else 
            iptr = iptf
          endif
          ioff1 = offset(w_grp%inode)
          ioff2 = ioff1+w_grp%mydim
          if ( iptr .ge. ioff1 &
               .and. &
               iptr .lt. ioff2 ) then
            jj = iptr - ioff1 + 1
            ! if ( icv .eq. 1 .and. verbose ) write(*, '(4(a,i5))') &
            !    ", wgrp%inode ", w_grp%inode, ", ipt ", ipt, &
            !    ", iptf ", iptf, ", jj ", jj
            isdf_in%Psi_intp(ipt, ii, isp, ikp) = kpt%wfn(isp,ikp)%dwf(jj, ii)
          endif
        enddo ! ipt = 1, n_intp_r
      enddo ! ii = 1, kpt%wfn(isp,ikp)%nmem
      call MPI_ALLREDUCE(MPI_IN_PLACE, isdf_in%Psi_intp(1,1,isp,ikp), n_intp_r*kpt%wfn(isp,ikp)%nmem, MPI_DOUBLE, MPI_SUM, &
        w_grp%comm, errinfo)
      if ( w_grp%inode .eq. 0 ) then
        write(6, *) isp, ikp, " Psi_intp = " 
        call printmatrix(isdf_in%Psi_intp(1,1,1,1), n_intp_r, kpt%wfn(isp,ikp)%nmem, 6) 
      endif 
      !stop 
      !
      P_intp = zero
      Q_intp = zero
      PsiV_intp = zero
      PsiC_intp = zero
      Amtrx = zero
        !
        ! initialize matrices with zero
        !
        P = zero
        Q = zero
        PsiV  = zero
        PsiC  = zero
        call timacc(53,1,tsec)
        !
        ! The following loop calculate P(r,r_u,jrp), Q(r,r_u,jrp) for all representations 
        !
        !if(peinf%master) then
        !  write(6,*) " kpt%wfn(isp,ikp)%map ", &
        !    kpt%wfn(isp,ikp)%map(:)
        !endif
        do jrp = 1, gvec%syms%ntrans
#ifdef DEBUG
          if (peinf%master) write(6,*) " jrp = ", jrp
          if (verbose .and. peinf%master) &
             write(dbgunit, *) ' isp = ', isp, ', ikp = ', ikp
#endif
          do iv = 1, isdf_in%nv(isp,ikp,jrp)
            IVV = isdf_in%ivlist(iv,isp,ikp,jrp)
            JVV = kpt%wfn(isp,ikp)%map(IVV)
            ! PsiV_i(r)
            call dcopy(w_grp%mydim, & 
              kpt%wfn(isp,ikp)%dwf(1, JVV),1,PsiV(1,iv),1)
              !
              ! pick the interpolation points
              !
              do ipt = 1, n_intp_r
                iptf = isdf_in%intp_r(ipt)
                if ( gvec%syms%ntrans > 1 ) then 
                  iptr = iptf / gvec%syms%ntrans + 1
                else 
                  iptr = iptf
                endif
                ioff1 = offset(w_grp%inode)
                ioff2 = ioff1+w_grp%mydim
                if ( iptr .ge. ioff1 &
                     .and. &
                     iptr .lt. ioff2 ) then
                  jj = iptr - ioff1 + 1
                  ! if ( icv .eq. 1 .and. verbose ) write(*, '(4(a,i5))') &
                  !    ", wgrp%inode ", w_grp%inode, ", ipt ", ipt, &
                  !    ", iptf ", iptf, ", jj ", jj
                  PsiV_intp(ipt, iv, jrp) = kpt%wfn(isp,ikp)%dwf(jj, JVV)
                endif
              enddo
          enddo ! iv loop
          !
          do ic = 1, isdf_in%nc(isp,ikp,jrp)
            ICC = isdf_in%iclist(ic,isp,ikp,jrp)
            JCC = kpt%wfn(isp,ikp)%map(ICC)
            ! PsiC_i(r)
            call dcopy(w_grp%mydim, & 
              kpt%wfn(isp,ikp)%dwf(1, JCC),1,PsiC(1,ic),1)
            !
            ! pick the interpolation points
            !
            do ipt = 1, n_intp_r
              iptf = isdf_in%intp_r(ipt)
              if ( gvec%syms%ntrans > 1 ) then 
                iptr = iptf / gvec%syms%ntrans + 1
              else 
                iptr = iptf
              endif
              ioff1 = offset(w_grp%inode)
              ioff2 = ioff1+w_grp%mydim
              if ( iptr .ge. ioff1 &
                   .and. &
                   iptr .lt. ioff2 ) then
                jj = iptr - ioff1 + 1
                ! if ( icv .eq. 1 .and. verbose ) write(*, '(4(a,i5))') &
                !    ", wgrp%inode ", w_grp%inode, ", ipt ", ipt, &
                !    ", iptf ", iptf, ", jj ", jj
                PsiC_intp(ipt, ic, jrp) = kpt%wfn(isp,ikp)%dwf(jj, JCC)
              endif 
            enddo ! ipt 
          enddo ! ic loop
          !
          call MPI_ALLREDUCE(MPI_IN_PLACE, PsiV_intp(1,1,jrp), n_intp_r*isdf_in%nv(isp,ikp,jrp), MPI_DOUBLE, MPI_SUM, &
            w_grp%comm, errinfo)
          call MPI_ALLREDUCE(MPI_IN_PLACE, PsiC_intp(1,1,jrp), n_intp_r*isdf_in%nc(isp,ikp,jrp), MPI_DOUBLE, MPI_SUM, &
            w_grp%comm, errinfo)
          ! 
          ! Prepare P_intp and Q_intp
          !
          ! Calculate P(r_u,r_u,jrp), and Q(r_u,r_u,jrp)

          ! -- obsolete code --
          !call dgemm('n','t',n_intp_r,n_intp_r, isdf_in%nv(isp,ikp,jrp), one, &
          !  PsiV_intp(1,1,jrp), n_intp_r, PsiV_intp(1,1,jrp), n_intp_r, zero, &
          !  P_intp(1,1,jrp), n_intp_r)
          !call dgemm('n','t',n_intp_r,n_intp_r, isdf_in%nc(isp,ikp,jrp), one, &
          !  PsiC_intp(1,1,jrp), n_intp_r, PsiC_intp(1,1,jrp), n_intp_r, zero, &
          !  Q_intp(1,1,jrp), n_intp_r)
           
          ! -- new code --
          call dgemm_hl('n','t',n_intp_r,n_intp_r, isdf_in%nv(isp,ikp,jrp), one, &
            PsiV_intp(1,1,jrp), n_intp_r, PsiV_intp(1,1,jrp), n_intp_r, zero, &
            P_intp(1,1,jrp), n_intp_r, opt%linear_algebra)
          !if (peinf%inode .eq. 0) write(6, *) "after gemm_hl ", &
          !    PsiV_intp(1,1,jrp), PsiV_intp(1,1,jrp), P_intp(1,1,jrp)
          call dgemm_hl('n','t',n_intp_r,n_intp_r, isdf_in%nc(isp,ikp,jrp), one, &
            PsiC_intp(1,1,jrp), n_intp_r, PsiC_intp(1,1,jrp), n_intp_r, zero, &
            Q_intp(1,1,jrp), n_intp_r, opt%linear_algebra)
          !if (peinf%inode .eq. 0) write(6, *) "after gemm_hl ", &
          !    PsiC_intp(1,1,jrp), PsiC_intp(1,1,jrp), Q_intp(1,1,jrp)
          !
          ! Prepare P and Q
          !
          ! P(r,r_u,jrp) = \sum_{V~jrp} \PsiV(r) \PsiV(r_u)
          ! -- old code --
          !call dgemm('n','t',w_grp%mydim, n_intp_r, isdf_in%nv(isp,ikp,jrp), one, &
          !  PsiV(1,1), w_grp%mydim, PsiV_intp(1,1,jrp), n_intp_r, zero, P(1,1,jrp), w_grp%mydim)

          ! -- new code --
          call dgemm_hl('n','t',w_grp%mydim, n_intp_r, isdf_in%nv(isp,ikp,jrp), one, &
             PsiV(1,1), w_grp%mydim, PsiV_intp(1,1,jrp), n_intp_r, zero, P(1,1,jrp), w_grp%mydim, &
             opt%linear_algebra)
          !if (peinf%inode .eq. 0) write(6, *) "after gemm_hl ", &
          !    PsiV(1,1), PsiV_intp(1,1,jrp), P(1,1,jrp)

          ! Q(r,r_u,jrp) = \sum_{C~jrp} \PsiC(r) \PsiC(r_u)
          ! -- old code -- 
          !call dgemm('n','t',w_grp%mydim, n_intp_r, isdf_in%nc(isp,ikp,jrp), one, &
          !  PsiC(1,1), w_grp%mydim, PsiC_intp(1,1,jrp), n_intp_r, zero, Q(1,1,jrp), w_grp%mydim)

          ! -- new code --
          call dgemm_hl('n','t',w_grp%mydim, n_intp_r, isdf_in%nc(isp,ikp,jrp), one, &
             PsiC(1,1), w_grp%mydim, PsiC_intp(1,1,jrp), n_intp_r, zero, Q(1,1,jrp), w_grp%mydim, &
             opt%linear_algebra)
          !if (peinf%inode .eq. 0) write(6, *) "after gemm_hl ", &
          !    PsiC(1,1), PsiC_intp(1,1,jrp), Q(1,1,jrp)

#ifdef DEBUG
          if ( .true. .and. peinf%master ) then
            write(dbgunit, *) "jrp =", jrp
            !write(dbgunit, '("PsiV = ")' ) 
            !call printmatrix ( PsiV(1,1), w_grp%mydim, isdf_in%nv(isp,ikp,jrp), dbgunit )
            !write(dbgunit, '("PsiC = ")' ) 
            !call printmatrix ( PsiC(1,1), w_grp%mydim, isdf_in%nc(isp,ikp,jrp), dbgunit )
            !write(dbgunit, '("PsiV_intp = ")' ) 
            !call printmatrix ( PsiV_intp(1,1,jrp), n_intp_r, isdf_in%nv(isp,ikp,jrp), dbgunit )
            !write(dbgunit, '("PsiC_intp = ")' ) 
            !call printmatrix ( PsiC_intp(1,1,jrp), n_intp_r, isdf_in%nc(isp,ikp,jrp), dbgunit )
            write(dbgunit, '("P = ")' ) 
            call printmatrix ( P(1,1,jrp), w_grp%mydim, n_intp_r, dbgunit )
            write(dbgunit, '("Q = ")' ) 
            call printmatrix ( Q(1,1,jrp), w_grp%mydim, n_intp_r, dbgunit )
            write(dbgunit, '("P_intp = ")' ) 
            call printmatrix ( P_intp(1,1,jrp), n_intp_r, n_intp_r, dbgunit )
            write(dbgunit, '("Q_intp = ")' ) 
            call printmatrix ( Q_intp(1,1,jrp), n_intp_r, n_intp_r, dbgunit )
          endif
#endif
          !
        enddo ! jrp loop
        !stop
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
          do lrp1 = 1, gvec%syms%ntrans
            !
            ! The direct product of lrp1 and lrp2 should be irp: lrp1 * lrp2 = irp
            !
            do lrp2 = 1, gvec%syms%ntrans
               if(gvec%syms%prod(lrp1,lrp2) == irp) exit
            enddo
#ifdef DEBUG
            if (peinf%master) then
               write(dbgunit,*) "lrp1 ", lrp1, ", lrp2 ", lrp2, ", irp ", irp
            endif
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
            !  Amtrx(n_intp_r, n_intp_r)     intermediate variable, store C * C^T
            !  Bmtrx(n_intp_r, w_grp%mydim)         intermediate variable, store C * Z^T
            !  Xmtrx(n_intp_r, w_grp%mydim)         intermediate variable, store zeta^T
            ! ---------------------
              !
              ! calculate A = P_intp.Q_intp (Note: This is an element-wise multipliation)
              !
              Amtrx(1:n_intp_r,1:n_intp_r,irp) = Amtrx (1:n_intp_r,1:n_intp_r,irp) + &
                                             P_intp(1:n_intp_r,1:n_intp_r,lrp1) * &
                                             Q_intp(1:n_intp_r,1:n_intp_r,lrp2)
            !
            ! calculate B = (P.Q)^T (Note: This is an element-wise multiplication)
            !
            ! Comment: this will crash with intel compiler
            !Bmtrx(1:n_intp_r,1:w_grp%mydim) = Bmtrx(1:n_intp_r,1:w_grp%mydim) + &
            !                     transpose( P(1:w_grp%mydim, 1:n_intp_r, lrp1) * &
            !                                Q(1:w_grp%mydim, 1:n_intp_r, lrp2) )
            ! End comment
            acc_B(1:w_grp%mydim, 1:n_intp_r) = P(1:w_grp%mydim, 1:n_intp_r, lrp1) * &
                                            Q(1:w_grp%mydim, 1:n_intp_r, lrp2)
            Bmtrx(1:n_intp_r,1:w_grp%mydim) = Bmtrx(1:n_intp_r,1:w_grp%mydim) + &
              transpose( acc_B(1:w_grp%mydim, 1:n_intp_r) )
            !
#ifdef DEBUG
            if (verbose .and. peinf%master) then
               write(dbgunit, '(" P*Q = ")')
               call printmatrix ( acc_B(1:w_grp%mydim, 1:n_intp_r), w_grp%mydim, n_intp_r, dbgunit)
            endif
#endif
          enddo ! lrp1
          !
          ! solve linear equation A * X = B
          !
#ifdef DEBUG
          if (.true. .and. peinf%master) then
            write(dbgunit, '(" irp =", i3, "  Amtrx = ")') irp
            call printmatrix ( Amtrx(1:isdf_in%n_intp_r, 1:isdf_in%n_intp_r, irp), &
              isdf_in%n_intp_r,  isdf_in%n_intp_r, dbgunit )
            write(dbgunit, '(" irp =", i3, "  Bmtrx = ")') irp
            call printmatrix ( Bmtrx(1:isdf_in%n_intp_r, 1:w_grp%mydim) , &
              isdf_in%n_intp_r, w_grp%mydim, dbgunit )
          endif
#endif
          ALLOCATE(ipiv(n_intp_r))
          ! Note that dgesv will change Amtrx, so we need a temporary Amtrx1
          Amtrx1(1:n_intp_r,1:n_intp_r) = Amtrx(1:n_intp_r,1:n_intp_r,irp)
          call dgesv(isdf_in%n_intp_r, w_grp%mydim, Amtrx1(1,1), n_intp_r, ipiv, Bmtrx, &
            n_intp_r, einfo)
          !call magmaf_dgesv( )
#ifdef DEBUG
          write(6, *) " inode ", peinf%inode, " einfo = ", einfo
#endif
          DEALLOCATE(ipiv)
          Xmtrx(1:w_grp%mydim,1:n_intp_r) = transpose(Bmtrx(1:n_intp_r,1:w_grp%mydim)) ! probably we can remove Bmtrx here now
#ifdef DEBUG
          ! for debug use
          write(outdbg, *) " Large Xmtrx "
          do ii = 1, w_grp%mydim
            do jj = 1, isdf_in%n_intp_r
              if ( abs(Xmtrx(ii, jj)) > 1.0e7 ) then
                write(outdbg, *) ii, jj, Xmtrx(ii, jj)
              endif
            enddo ! jj
          enddo ! ii
#endif
          call MPI_BARRIER(peinf%comm, errinfo)
#ifdef DEBUG
          if (.true. .and. peinf%master) then
            write(dbgunit, '(" irp =", i3, "  Zeta = ")') irp
            call printmatrix ( Xmtrx(1:w_grp%mydim, 1:isdf_in%n_intp_r), &
              w_grp%mydim, isdf_in%n_intp_r, dbgunit )
          endif
#endif
          !
          ! Copy Xmtrx to zeta
          !
          ! Xmtrx (w_grp%mydim, n_intp_r)
          ALLOCATE(tmp_array(w_grp%mydim, isdf_in%n_intp_r))
          tmp_array(1:w_grp%mydim, 1:isdf_in%n_intp_r) = &
            Xmtrx(1:w_grp%mydim, 1:isdf_in%n_intp_r)
          !call h5dwrite_f( dset_zeta_id(isp,irp), H5T_NATIVE_DOUBLE, &
          !  Xmtrx(1:w_grp%mydim,1:isdf_in%n_intp_r), &
          !  data_dims, h5err, subdspace, dspace_zeta )
          call h5dwrite_f( dset_zeta_id(isp,irp), H5T_NATIVE_DOUBLE, &
            tmp_array, data_dims, h5err )
          if(peinf%master) write(6,*) "done write hdf5"
          DEALLOCATE(tmp_array)
          !
        enddo ! jrp loop, jrp = 1, gvec%syms%ntrans/r_grp%num
        !
    enddo ! isp loop
    !
#ifdef DEBUG
    if (peinf%master) write(6,*) "deallocate arrays"
#endif
    !
    ! DEALLOCATE arrays
    !
    DEALLOCATE(P)
    DEALLOCATE(Q)
    DEALLOCATE(acc_B)
    DEALLOCATE(P_intp)
    DEALLOCATE(Q_intP)
    DEALLOCATE(PsiV)
    DEALLOCATE(PsiC)
    DEALLOCATE(PsiV_intp)
    DEALLOCATE(PsiC_intp)
    DEALLOCATE(Amtrx)
    DEALLOCATE(Amtrx1)
    DEALLOCATE(Bmtrx)
    DEALLOCATE(Xmtrx)
    !
    call timacc(53,2,tsec)
    !
#ifdef DEBUG
    write(6,*) "inode: ", peinf%inode, " start to calculate Mmtrx"
    write(6,*) "r_grp%num ", r_grp%num
#endif
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
               rho_h_distr(1:w_grp%mydim, 1:w_grp%npes), data_dims, &
               h5err, subdspace, dspace_zeta )
#ifdef DEBUG
            if (peinf%master) write(dbgunit,*) " i_row ", i_row, ": done h5read_3"
#endif
            !
            ! initialize rho_h
            !
            rho_h = 0.d0
            !
            ii = w_grp%mygr*w_grp%npes + i_row + w_grp%inode
            call dgather(1,rho_h_distr,rho_h)
            if ( ii <= isdf_in%n_intp_r ) then
              !
              ! solve poisson equation to get: rho_h(r) = \int V_c(r,r') zeta_ii(r') dr'
              !
              call dpoisson(gvec, rho_h, irp)
            endif
            call dscatter(1,rho_h_distr,rho_h)
            !
            call h5sselect_hyperslab_f(dspace_vczeta, H5S_SELECT_SET_F, &
              shift, subdim, h5err, stride, block)
            call h5dwrite_f(dset_vczeta_id(rsp,irp), H5T_NATIVE_DOUBLE, &
              rho_h_distr(1:subdim(1), 1:subdim(2)), data_dims, &
              h5err, subdspace, dspace_vczeta)
            if(peinf%master) write(6,*) "done write hdf5"
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
        ALLOCATE(tmp_array(1:w_grp%mydim,1:isdf_in%n_intp_r))
        do isp = 1, nspin
          !call h5dread_f( dset_zeta_id(isp, irp), H5T_NATIVE_DOUBLE, &
          !   zeta( 1:w_grp%mydim,1:isdf_in%n_intp_r, isp ), data_dims, h5err, subdspace, dspace_zeta )
          call h5dread_f( dset_zeta_id(isp, irp), H5T_NATIVE_DOUBLE, &
             tmp_array, data_dims, h5err )
          zeta( 1:w_grp%mydim, 1:isdf_in%n_intp_r, isp ) = tmp_array( 1:w_grp%mydim, 1:isdf_in%n_intp_r )
#ifdef DEBUG
          if (peinf%master) write(dbgunit,*) "done h5read_1"
          if (peinf%master) then
            write(dbgunit, '("Read zeta() irp",i3," isp ",i3 )') irp, isp
            call printmatrix( zeta(1,1,isp), w_grp%mydim, isdf_in%n_intp_r, &
              dbgunit)
          endif
#endif
          if(kflag<2) then
          call h5dread_f( dset_vczeta_id(isp, irp), H5T_NATIVE_DOUBLE, &
             tmp_array, data_dims, h5err )
          vzeta( 1:w_grp%mydim, 1:isdf_in%n_intp_r, isp ) = tmp_array( 1:w_grp%mydim, 1:isdf_in%n_intp_r )
          endif
#ifdef DEBUG
          if (peinf%master) write(dbgunit,*) "done h5read_2"
#endif
        enddo
        DEALLOCATE(tmp_array)
        ! calculate mtrx now
        do rsp = 1, nspin
          do csp = 1, nspin
            if ( kflag > 0 ) then
              do ii = 1, isdf_in%n_intp_r
                fzeta(1:w_grp%mydim, ii) = &
                  fxc(1:w_grp%mydim, rsp, csp) &
                   * zeta(1:w_grp%mydim, ii, csp)
              enddo ! ii loop
#ifdef DEBUG
              if (peinf%master) then
                write(dbgunit,*) " test zeta fzeta irp =", irp
                do jj = 1, 10
                  write(dbgunit,*) jj-1, zeta(jj,1,1), fzeta(jj,1)
                enddo
              endif
#endif
              ! -- old code -- 
              !call dgemm('T', 'N', isdf_in%n_intp_r, isdf_in%n_intp_r, w_grp%mydim, &
              !  norm_factor, zeta(1,1,rsp), w_grp%mydim, fzeta(1,1), &
              !  w_grp%mydim, 1.d0, &
              !  isdf_in%Mmtrx(1, 1, rsp, csp, ikp, 2, irp), &
              !  isdf_in%n_intp_r)

              ! -- new code --
              call dgemm_hl('t', 'n', isdf_in%n_intp_r, isdf_in%n_intp_r, w_grp%mydim, &
                norm_factor, zeta(1,1,rsp), w_grp%mydim, fzeta(1,1), w_grp%mydim, 1.d0, &
                isdf_in%Mmtrx(1, 1, rsp, csp, ikp, 2, irp), &
                isdf_in%n_intp_r, opt%linear_algebra)

              !
            endif ! kflag > 0
            if ( kflag < 2 ) then
#ifdef DEBUG
              if (peinf%master) then
                write(dbgunit,*) " test zeta vzeta irp =", irp
                do jj = 1, 10
                  write(dbgunit,*) jj-1, zeta(jj,1,1), vzeta(jj,1,1)
                enddo
              endif
#endif
              ! -- old code --
              !call dgemm('T', 'N', isdf_in%n_intp_r, isdf_in%n_intp_r, w_grp%mydim, &
              !  norm_factor, zeta(1,1,rsp), w_grp%mydim, vzeta(1,1,csp), w_grp%mydim, 1.d0, &
              !  isdf_in%Mmtrx(1, 1, rsp, csp, ikp, 1, irp), &
              !  isdf_in%n_intp_r)

              ! -- new code --
              call dgemm_hl('t', 'n', isdf_in%n_intp_r, isdf_in%n_intp_r, w_grp%mydim, &
                norm_factor, zeta(1,1,rsp), w_grp%mydim, vzeta(1,1,csp), w_grp%mydim, 1.d0, &
                isdf_in%Mmtrx(1, 1, rsp, csp, ikp, 1, irp), &
                isdf_in%n_intp_r, opt%linear_algebra)

            endif
          enddo ! csp
        enddo ! rsp 
    enddo ! jrp loop
    call timacc(54,2,tsec)
    !
    DEALLOCATE( zeta  )
    if (kflag > 0) DEALLOCATE( fzeta )
    if (kflag < 2) DEALLOCATE( vzeta )
  enddo ! ikp loop
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
  DEALLOCATE(inv_ivlist)
  DEALLOCATE(inv_iclist)
  if ( kflag > 0 ) then
    DEALLOCATE(fxc)
  endif
  !
  call MPI_ALLREDUCE( MPI_IN_PLACE, isdf_in%Mmtrx(1,1,1,1,1,1,1), n_intp_r * n_intp_r * nspin * nspin * kpt%nk * 2 * gvec%syms%ntrans, &
    MPI_DOUBLE, MPI_SUM, w_grp%comm, errinfo )
  !
#ifdef DEBUG
  if ( peinf%master .and. .True. ) then
    do jrp = 1, gvec%syms%ntrans
      write( dbgunit, '(a,i2,a)' ) " Mmtrx (:, :, rsp=1, csp=1, ikp=1, 1, jrp=", jrp, ") = "
      !call printmatrix ( isdf_in%Mmtrx (1,1,1,1,1,2,jrp), n_intp_r, n_intp_r, dbgunit ) 
      call printmatrix ( isdf_in%Mmtrx (1:isdf_in%n_intp_r,1:isdf_in%n_intp_r,1,1,1,1,jrp), isdf_in%n_intp_r, isdf_in%n_intp_r, dbgunit )
    enddo
    do jrp = 1, gvec%syms%ntrans
      write( dbgunit, '(a,i2,a)' ) " Mmtrx (:, :, rsp=1, csp=1, ikp=1, 2, jrp=", jrp, ") = "
      !call printmatrix ( isdf_in%Mmtrx (1,1,1,1,1,2,jrp), n_intp_r, n_intp_r, dbgunit ) 
      call printmatrix ( isdf_in%Mmtrx (1:isdf_in%n_intp_r,1:isdf_in%n_intp_r,1,1,1,2,jrp), isdf_in%n_intp_r, isdf_in%n_intp_r, dbgunit )
    enddo
  endif
#endif
  !
  !stop
  !
  ! test calculating < 1 2 | F | 3 4 > 
  if ( peinf%master ) then
    isp = 1
    ikp = 1
    irp = 1
#ifdef DEBUG
    write (dbgunit, *) " test calculate < 1 2 | V_coul | 3 4 > "
#endif
    tmpvec = zero
    do icv1 = 1, 10
      ivv = isdf_in%invpairmap(1, icv1, isp, ikp, irp)
      icc = isdf_in%invpairmap(2, icv1, isp, ikp, irp)
      ii = kpt%wfn(isp,ikp)%map(ivv)
      jj = kpt%wfn(isp,ikp)%map(icc)
      tmpCmtrx1(1:isdf_in%n_intp_r) = &
       isdf_in%Psi_intp(1:isdf_in%n_intp_r, ii, isp, ikp) * &
       isdf_in%Psi_intp(1:isdf_in%n_intp_r, jj, isp, ikp)
      do icv2 = 1, 10
        ivv = isdf_in%invpairmap(1, icv2, isp, ikp, irp)
        icc = isdf_in%invpairmap(2, icv2, isp, ikp, irp)
        ii = kpt%wfn(isp,ikp)%map(ivv)
        jj = kpt%wfn(isp,ikp)%map(icc)
        tmpCmtrx2(1:isdf_in%n_intp_r) = &
         isdf_in%Psi_intp(1:isdf_in%n_intp_r, ii, isp, ikp) * &
         isdf_in%Psi_intp(1:isdf_in%n_intp_r, jj, isp, ikp)
        call dgemv('n', isdf_in%n_intp_r, isdf_in%n_intp_r, one, &
          isdf_in%Mmtrx (1, 1, isp, isp, ikp, 1, irp), isdf_in%n_intp_r, &
          tmpCmtrx1(1), 1, zero, tmpvec, 1)
        matel = ddot(isdf_in%n_intp_r, tmpvec, 1, &
          tmpCmtrx2(1), 1)
#ifdef DEBUG
        write (dbgunit, *) " icv1 ", icv1, " icv2 ", icv2, matel
#endif
      enddo ! icv2
    enddo ! icv1
#ifdef DEBUG
    write (dbgunit, *) " test calculate < 1 2 | f_xc | 3 4 > "
#endif
    tmpvec = zero
    do icv1 = 1, 10
      ivv = isdf_in%invpairmap(1, icv1, isp, ikp, irp)
      icc = isdf_in%invpairmap(2, icv1, isp, ikp, irp)
      ii = kpt%wfn(isp,ikp)%map(ivv)
      jj = kpt%wfn(isp,ikp)%map(icc)
      tmpCmtrx1(1:isdf_in%n_intp_r) = &
       isdf_in%Psi_intp(1:isdf_in%n_intp_r, ii, isp, ikp) * &
       isdf_in%Psi_intp(1:isdf_in%n_intp_r, jj, isp, ikp)
      do icv2 = 1, 10
        ivv = isdf_in%invpairmap(1, icv2, isp, ikp, irp)
        icc = isdf_in%invpairmap(2, icv2, isp, ikp, irp)
        ii = kpt%wfn(isp,ikp)%map(ivv)
        jj = kpt%wfn(isp,ikp)%map(icc)
        tmpCmtrx2(1:isdf_in%n_intp_r) = &
         isdf_in%Psi_intp(1:isdf_in%n_intp_r, ii, isp, ikp) * &
         isdf_in%Psi_intp(1:isdf_in%n_intp_r, jj, isp, ikp)
        call dgemv('n', isdf_in%n_intp_r, isdf_in%n_intp_r, one, &
          isdf_in%Mmtrx (1, 1, isp, isp, ikp, 2, irp), isdf_in%n_intp_r, &
          tmpCmtrx1(1), 1, zero, tmpvec, 1)
        matel = ddot(isdf_in%n_intp_r, tmpvec, 1, &
          tmpCmtrx2(1), 1)
#ifdef DEBUG
        write (dbgunit, *) " icv1 ", icv1, " icv2 ", icv2, matel
#endif
      enddo ! icv2
    enddo ! icv1
  endif
#ifdef DEBUG
  if ( peinf%master ) then
     close ( dbgunit )
  endif
#endif
  return
  !
end subroutine isdf_parallel_sym_lessmemory2
