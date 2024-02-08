! This doesn't work
! 
! INPUT:   charge density rho(:,nspin) for two different spins
! OUTPUT:  index of interpolation points at full grids
!
! Issue need to be resolved: 
! - How do you make sure there are no duplicated interpolation points??

subroutine cvt_parallel_sym(gvec, rho, nspin, n_intp, intp)

 use typedefs
 use mpi_module
 USE IFPORT   ! This is required to call subroutine rand() if intel compiler is used
 implicit none
#ifdef MPI
  include 'mpif.h'
#endif

 ! -- PARAMETERS --
 ! number of maximum iterations to find interpolation points
 integer, parameter :: max_iter = 100000
 ! random seed for initializing random interpolation points 
 ! integer, parameter :: rseed = 12348
 integer :: rseed
 ! convergence threshold 
 real(dp), parameter :: tol_conv = 1e-3

 ! number of spin
 integer, intent(in) :: nspin
 ! gvec stores info about the real-space grids 
 type (gspace), intent(in) :: gvec 
 ! charge density in real-space grids, stored in irreducible wedge
 real(dp), intent(in) :: rho(gvec%nr, nspin)
 ! input parameter, the number of interpolation points needed
 integer, intent(in) :: n_intp
 ! outputs, the index of intp pts in the full grid points
 integer, intent(out) :: intp(n_intp)
 ! the real-space coordinates (in a.u.) of intp pts
 real(dp) :: newpts(3,n_intp), oldpts(3,n_intp)
 ! counters for different loops
 integer :: iter, ipt, itran, igrid, jgrid, ii, &
   ig, info, nn, mm, gnn, gmm, istart,  &
   iend, igstart, igend, ipes
   
 integer :: outdbg = 12345 ! unit of file for debugging
 ! full grid points and charge density on full grid points
 real(dp), allocatable :: fullgrid(:,:)
 integer,  allocatable :: icenter(:)
 real(dp) :: bounds(2,3), dist, weightedpos(3), weights, & 
  diff, vtmp(3), dist_tmp, mindist
 integer :: pt_tmp(3), values(8)

 open(outdbg, file="cvt_debug.dat", form='formatted', status='unknown')
 ! generate all points in the unfolded real space grid
 allocate(fullgrid(3,gvec%nr * gvec%syms%ntrans))
 ! get the charge density on the full grid
 allocate(icenter(gvec%nr * gvec%syms%ntrans))
 ! 
 nn = mod( n_intp, peinf%npes )
 mm = n_intp / peinf%npes
 istart = 0
 iend   = 0
 ! peinf%inode is from 0 to peinf%npes
 do ipes = 0, peinf%inode
   istart = iend + 1
   iend = istart + mm - 1
   if ( ipes < nn) iend = iend + 1
 enddo
 if (.True.) then
   do ipes = 0, peinf%inode 
     if (peinf%inode == ipes) then
       write(6, '(A,i6,2(A,i6))') "inode: ", peinf%inode, &
        ", istart ", istart, ", iend ", iend
     endif
   enddo
 endif
 !
 gnn = mod( gvec%nr, peinf%npes )
 gmm = gvec%nr / peinf%npes 
 igstart = 0
 igend   = 0
 ! peinf%inode is from 0 to peinf%npes
 do ipes = 0, peinf%inode 
   igstart = igend + 1
   igend   = igstart + gmm - 1
   if (ipes < gnn) igend = igend + 1
 enddo 
 if (.True.) then
   do ipes = 0, peinf%inode 
     if (peinf%inode == ipes) then
       write(6, '(A,i6,2(A,i6))') "inode: ", peinf%inode, &
        ", igstart ", igstart, ", igend ", igend
     endif
   enddo
 endif
 !
 fullgrid = 0
 igrid = (igstart-1)*gvec%syms%ntrans ! counter for full-grid points
 do ig = igstart, igend 
    do itran = 1, gvec%syms%ntrans
       igrid = igrid + 1
       call unfold(gvec%r(1,ig), &
         gvec%syms%trans(1,1,itran),gvec%shift(1),pt_tmp(1))
       do ii = 1,3
         fullgrid(ii,igrid) = (pt_tmp(ii) + gvec%shift(ii)) * gvec%step(ii)
       enddo
    enddo
 enddo
 call MPI_ALLREDUCE( MPI_IN_PLACE, fullgrid(1,1), 3*gvec%nr*gvec%syms%ntrans, &
    MPI_DOUBLE, MPI_SUM, peinf%comm, info)
 !
 ! find the lower and higher bounds of full-grid points
 do ii = 1,3
    bounds(1,ii) = minval(fullgrid(ii,:))/2.0
    bounds(2,ii) = maxval(fullgrid(ii,:))/2.0
    call MPI_ALLREDUCE( MPI_IN_PLACE, bounds(1,ii), 1, MPI_DOUBLE, MPI_MIN, &
       peinf%comm, info )
    call MPI_ALLREDUCE( MPI_IN_PLACE, bounds(2,ii), 1, MPI_DOUBLE, MPI_MAX, &
       peinf%comm, info )
 enddo
 !
 call date_and_time(VALUES=values)
 rseed = values(6)*11+values(7)*1000+values(8)
 !rseed = 1518543090
 if (peinf%master) write(6,*) "# rseed for random number generator is: ", rseed
 call srand(rseed)
 ! generate initial guess of interpolation points
 ! write(outdbg,'(a)') "# Initial guess of interpolation points "
 do ipt = 1, n_intp
    do ii = 1,3
        ! generate some random points in the full grids
        ! multiply by 0.9 to make sure these random points are inside the
        ! boundary
        newpts(ii,ipt) = bounds(1,ii) + rand(0)*0.9*(bounds(2,ii)-bounds(1,ii)) 
    enddo
    ! print out the intial random interpolation points
    if (peinf%master) write(outdbg,'("# ",3f10.4)') newpts(1:3,ipt)
    if (peinf%inode == 2) write(outdbg,'("# ",3f10.4)') newpts(1:3,ipt)
 enddo
 !
 ! perform centroidal voronoi tesselation (CVT) algorithm
 ! For more details, see arXiv: 1711.01531v1 by Kun Dong, Wei Hu and Lin Lin
 ! (2017)
 oldpts = newpts
 do iter = 1, max_iter
    ! for each point in the full grid, find which interpolation points is 
    ! the closest it. Then put the index of intp point to icenter(:)
    icenter = 0
    igrid = (igstart-1)*gvec%syms%ntrans
    do ig = igstart, igend
       do itran = 1, gvec%syms%ntrans
          igrid = igrid + 1
          ! set a very large initial value for dist
          dist = 100.0 * ( maxval(bounds(2,:)) - minval(bounds(1,:)) )**2  
          ! find which intp point is closest to the current grid point
          do ipt = 1, n_intp
             vtmp = newpts(:,ipt) - fullgrid(:,igrid)
             dist_tmp = dot_product(vtmp, vtmp)
             if (dist_tmp < dist) then
                dist = dist_tmp
                icenter(igrid) = ipt
             endif
          enddo ! ipt
       enddo ! itran
    enddo ! ipt
    call MPI_ALLREDUCE( MPI_IN_PLACE, icenter, gvec%nr*gvec%syms%ntrans, MPI_INTEGER, MPI_SUM, &
       peinf%comm, info ) 
    
    ! Now update the interpolation pts
    diff = 0.0
    newpts = 0
    do ipt = istart, iend
       weightedpos(1:3) = 0.0
       weights = 0.0
       do igrid = 1, gvec%nr * gvec%syms%ntrans
          if (icenter(igrid) .ne. ipt) cycle
          jgrid = igrid / gvec%syms%ntrans + 1
          weightedpos = weightedpos + fullgrid(:,igrid) * & 
            ( rho(jgrid,1) + rho(jgrid,nspin) )
          weights = weights + ( rho(jgrid,1) + rho(jgrid,nspin) )
       enddo
       ! update the new intp points with the centroid
       ! The following is a simple method to avoid the case of weights == 0. It
       ! may be changed later with a better method.
       ! Simple minded fix, just generate a new randome points
       if(weights .lt. 1.0e-10) then
          newpts(:,ipt) = oldpts(:,ipt)
       else
          newpts(:,ipt) = weightedpos(:) / weights
       endif
       vtmp = newpts(1:3,ipt) - oldpts(1:3,ipt)
       diff = diff + sqrt(dot_product(vtmp,vtmp))
    enddo ! loop ipt = istart, iend
    call MPI_ALLREDUCE( MPI_IN_PLACE, diff, 1, MPI_DOUBLE, MPI_SUM, &
       peinf%comm, info )
    call MPI_ALLREDUCE( MPI_IN_PLACE, newpts, 3*n_intp, MPI_DOUBLE, MPI_SUM, &
       peinf%comm, info )
    if (peinf%inode == 0) then 
        write(*,'(15(f8.3))') ((newpts(ii,ipt),ii=1,3),ipt=1,n_intp)
        write(*,'(i8,a,f18.12)') iter, " diff (a.u.) ", diff/n_intp
    endif
    if (diff/n_intp < tol_conv) then ! conv threshold is satisfied, break the loop??
       exit 
    endif
    oldpts = newpts
 enddo ! iter

 ! Loop over all the grid points to find a grid point that is closest to newpts(:,ipt)
 intp = 0
 do ipt = istart, iend
    mindist = 1.e9 ! initialize a very large number
    do igrid = 1, gvec%nr * gvec%syms%ntrans
       vtmp = newpts(1:3,ipt) - fullgrid(1:3,igrid)
       dist = sqrt(dot_product(vtmp,vtmp)) 
       if (dist < mindist) then
          mindist = dist
          intp(ipt) = igrid
       endif
    enddo ! igrid
 enddo ! ipt
 call MPI_ALLREDUCE( MPI_IN_PLACE, intp, n_intp, MPI_INTEGER, MPI_SUM, &
     peinf%comm, info )
 !
 if (peinf%master) then
 do ipt = 1, n_intp
    write(outdbg,'(i10,i15,3f8.3," H ",3f8.3)') (intp(ipt)-1)/gvec%syms%ntrans,intp(ipt),fullgrid(1:3,intp(ipt)), &
    fullgrid(1:3,intp(ipt))*0.529177+ 4.1011
 enddo
 endif
 !
 deallocate(fullgrid)
 deallocate(icenter)

 close(outdbg) ! close the dbg file
end subroutine cvt_parallel_sym
