subroutine cvt_wfn(gvec, wfn, nspin, nk, n_intp, intp)

 use typedefs
 use mpi_module
 implicit none
#ifdef MPI
  include 'mpif.h'
#endif

 ! number of spin
 integer, intent(in) :: nspin, nk ! probably we don't need nk here !!!
 ! gvec stores info about the real-space grids 
 type (gspace), intent(in) :: gvec 
 ! charge density in real-space grids, stored in irreducible wedge
 type (wavefunction), intent(in) :: wfn(nspin,nk)
 ! input parameter, the number of interpolation points needed
 integer, intent(in) :: n_intp
 ! outputs, the index of intp pts in the full grid points
 integer, intent(out) :: intp(n_intp)

 ! counters for different loops
 integer :: isp, istate, errinfo
 ! full grid points and charge density on full grid points
 real(dp), allocatable :: rho(:,:), tmprho(:,:), tmpvec(:)

 ! get the charge density on the irreducible grid
 allocate(rho(gvec%nr,nspin))
 allocate(tmprho(gvec%nr,nspin))
 allocate(tmpvec(w_grp%mydim))
  
 ! This only works for confined system for now, so I assume there is only one kpoint.
 rho = 0.d0
 tmprho = 0.d0
 do isp = 1, nspin
    do istate = 1, wfn(isp,1)%nmem
       tmpvec(1:w_grp%mydim) = wfn(isp,1)%dwf(1:w_grp%mydim,istate)
       ! calculate: tmpvec = wfn(isp,1)%dwf ** 2
       call dmultiply_vec(w_grp%mydim, wfn(isp,1)%dwf(1,istate), tmpvec(1))
       tmprho(w_grp%offset+1:w_grp%offset+w_grp%mydim,isp) = &
         tmprho(w_grp%offset+1:w_grp%offset+w_grp%mydim,isp) + &
         tmpvec(1:w_grp%mydim)
    enddo ! istate loop
 enddo ! isp loop
 call MPI_AllREDUCE(tmprho, rho, gvec%nr*nspin, MPI_DOUBLE, MPI_SUM, &
    w_grp%comm, errinfo)
 deallocate(tmpvec)
 deallocate(tmprho)
 if(peinf%master) then
    call cvt(gvec, rho, nspin, n_intp, intp)
 endif
 deallocate(rho)

end subroutine cvt_wfn
