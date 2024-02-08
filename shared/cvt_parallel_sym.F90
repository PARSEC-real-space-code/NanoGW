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
 use esdf
 use psp_module
#ifdef INTEL
 USE IFPORT   ! This is required to call subroutine rand() if intel compiler is used
#endif
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
 real(dp), parameter :: tol_conv = 1e-5

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
 real(dp) :: initialpts(3,n_intp), newpts(3,n_intp), oldpts(3,n_intp), &
   in_mindist(2), out_mindist(2)
 ! counters for different loops
 integer :: iter, ipt, jpt, itran, igrid, jgrid, ii, &
   ig, info, nn, mm, gnn, gmm, istart,  &
   iend, igstart, igend, ipes, ity, iat, jj, &
   flag, minig
 integer :: outdbg = 12345 ! unit of file for debugging
 integer :: outdbg2 = 123456
 ! full grid points and charge density on full grid points
 real(dp), allocatable :: fullgrid(:,:)
 integer,  allocatable :: icenter(:)
 real(dp) :: bounds(2,3), dist, weightedpos(3), weights, & 
  diff, vtmp(3), dist_tmp, mindist, xmin, ymin, zmin, &
  xmax, ymax, zmax, xmov, ymov, zmov

 integer :: pt_tmp(3), values(8), select_grid(n_intp)

 ! generate all points in the unfolded real space grid
 allocate(fullgrid(3,gvec%nr * gvec%syms%ntrans))
 ! get the charge density on the full grid
 allocate(icenter(gvec%nr * gvec%syms%ntrans))
 ! 
 nn = mod( n_intp, peinf%npes )
 mm = n_intp / peinf%npes
 istart = 0
 iend   = 0
 ! peinf%inode [0, peinf%npes]
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
 !
 do ii = 1,3
    bounds(1,ii) = minval(fullgrid(ii,:))
    bounds(2,ii) = maxval(fullgrid(ii,:))
    call MPI_ALLREDUCE( MPI_IN_PLACE, bounds(1,ii), 1, MPI_DOUBLE, MPI_MIN, &
       peinf%comm, info )
    call MPI_ALLREDUCE( MPI_IN_PLACE, bounds(2,ii), 1, MPI_DOUBLE, MPI_MAX, &
       peinf%comm, info )
 enddo
 !
 call date_and_time(VALUES=values)
 rseed = values(6)*11+values(7)*1000+values(8)
 ! rseed = 1518543090
 if (peinf%master) write(6,*) "# rseed for random number generator is: ", rseed
 call srand(rseed)
 ! generate initial guess of interpolation points
 ! write(outdbg,'(a)') "# Initial guess of interpolation points "
 if (peinf%master) open(outdbg, file="initial_pts.vasp", form='formatted', status='unknown')
 do ipt = 1, n_intp
    do ii = 1,3
       ! generate some random points in the full grids
       ! multiply by 0.9 to make sure these random points are inside the
       ! boundary
       newpts(ii,ipt) = bounds(1,ii) + &
                 0.25*(bounds(2,ii)-bounds(1,ii)) + &
          rand(0)*0.5*(bounds(2,ii)-bounds(1,ii)) 
    enddo
    ! save initial pts for printing out later
    if (peinf%master) then
       initialpts(1:3,ipt) = newpts(1:3,ipt)  
    endif
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
          jgrid = (igrid-1) / gvec%syms%ntrans + 1
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
 if (.True.) then
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
 else !
    select_grid = 0
    do ipt = 1, n_intp
       mindist = 1.e9 ! initialize a very large number
       ! Loop over all the grid points to find a grid point that is closest to newpts(:,ipt)
       do igrid = igstart, igend ! igrid is the index in the reduced domain
          flag = 0
          ! check if this grid point has already exist in select_grid()
          ! if it does, flag becomes 1
          if(ipt .ge. 2) then
             do jpt = 1, ipt-1
                if(igrid .eq. select_grid(jpt)) then
                   flag = 1
                   exit
                endif
             enddo
          endif
          ! if this grid point (in reduced grid) already exist, then skip it
          if (flag .eq. 1 ) cycle
          ! if this grid point is not in select_grid(), then proceed
          do itran = 1, gvec%syms%ntrans
             ! ig is the index in the full domain
             ig = (igrid-1)*gvec%syms%ntrans + itran
             vtmp = newpts(1:3,ipt) - fullgrid(1:3,ig)
             dist = sqrt(dot_product(vtmp,vtmp))
             if (dist < mindist) then
                mindist = dist
                minig = ig
             endif
          enddo ! itrans
       enddo ! igrid 
       in_mindist(1) = mindist ! local minimum dist found in igstart to igend
       in_mindist(2) = minig   ! index of the minimum grid stored in procs
       do ipes = 0, peinf%npes-1
         if(peinf%inode == ipes) &
            write(*,*) peinf%inode, in_mindist(1), int(in_mindist(2))
         call MPI_BARRIER(peinf%comm,info)
       enddo
       call MPI_ALLREDUCE( in_mindist, out_mindist, 1, MPI_2DOUBLE_PRECISION, &
         MPI_MINLOC, peinf%comm, info )
       do ipes = 0, peinf%npes-1
         if(peinf%inode == ipes) &
            write(*,*) peinf%inode, out_mindist(1), int(out_mindist(2))
         call MPI_BARRIER(peinf%comm,info)
       enddo   
       intp(ipt) = int( out_mindist(2) )
       select_grid(ipt) = intp(ipt) / gvec%syms%ntrans + 1
    enddo ! ipt
 endif
 !
 if (peinf%master) then

    xmin =  99999.0
    xmax = -99999.0
    ymin =  99999.0
    ymax = -99999.0
    zmin =  99999.0
    zmax = -99999.0

    do ity = 1, size(psp)
      do iat = 1, psp(ity)%natmi
        if (xmin > psp(ity)%ratm(1,iat)) xmin = psp(ity)%ratm(1,iat)
        if (xmax < psp(ity)%ratm(1,iat)) xmax = psp(ity)%ratm(1,iat)
        if (ymin > psp(ity)%ratm(2,iat)) ymin = psp(ity)%ratm(2,iat)
        if (ymax < psp(ity)%ratm(2,iat)) ymax = psp(ity)%ratm(2,iat)
        if (zmin > psp(ity)%ratm(3,iat)) zmin = psp(ity)%ratm(3,iat)
        if (zmax < psp(ity)%ratm(3,iat)) zmax = psp(ity)%ratm(3,iat)
      enddo
    enddo
     
    xmov = 0.5*(xmin+xmax)
    ymov = 0.5*(ymin+ymax)
    zmov = 0.5*(zmin+zmax)
     
    open(outdbg2, file="final_pts.vasp", form='formatted', status='unknown')
    write(outdbg2, '( "Fr atoms are interpolation points " )')
    write(outdbg2, '( "1.0")')
    write(outdbg2, '( 3f12.6 )') (bounds(2,1)-bounds(1,1))*0.52917721, 0.0, 0.0 ! lattice vectors
    write(outdbg2, '( 3f12.6 )') 0.0, (bounds(2,2)-bounds(1,2))*0.52917721, 0.0 ! lattice vectors
    write(outdbg2, '( 3f12.6 )') 0.0, 0.0, (bounds(2,3)-bounds(1,3))*0.52917721 ! lattice vectors
    do ity = 1, size(psp)
      write(outdbg2, '("   ", A, "  ")', advance='no') psp(ity)%name
    enddo
    write(outdbg2, '("  Fr ")')
    do ity = 1, size(psp)
      write(outdbg2, '(i5)', advance='no') psp(ity)%natmi
    enddo
    write(outdbg2, '(i5)') n_intp
    write(outdbg2, '(" Cartesian ")')
    do ity = 1, size(psp)
      do iat=1, psp(ity)%natmi
        if (sqrt(xmov*xmov+ymov*ymov+zmov*zmov) > 0.1) then
          write(outdbg2,'(3f12.6)') \
             (psp(ity)%ratm(1,iat)-xmov+0.5*(bounds(2,1)-bounds(1,1)))*0.52917721,\
             (psp(ity)%ratm(2,iat)-ymov+0.5*(bounds(2,2)-bounds(1,2)))*0.52917721,\
             (psp(ity)%ratm(3,iat)-zmov+0.5*(bounds(2,3)-bounds(1,3)))*0.52917721
        else
          write(outdbg2,'(3f12.6)') ( (psp(ity)%ratm(jj,iat) + \
          0.5*(bounds(2,jj)-bounds(1,jj)))*0.52917721, jj=1,3)
        endif
      enddo
    enddo
    do ipt = 1, n_intp
      write(outdbg2,'(3f12.6)') (( fullgrid(jj,intp(ipt)) + \
         0.5*(bounds(2,jj)-bounds(1,jj)) )*0.52917721, jj=1,3)
      !write(outdbg2,'(3f12.6)') (( newpts(1:3, ipt) + \
      !  0.5*(bounds(2,jj)-bounds(1,jj)) )*0.52917721, jj=1,3)
    enddo ! ipt
    close(outdbg2)

    open(outdbg2, file="final_pts.vasp", form='formatted', status='unknown')
    write(outdbg2, '( "Fr atoms are interpolation points " )')
    write(outdbg2, '( "1.0")')
    write(outdbg2, '( 3f12.6 )') (bounds(2,1)-bounds(1,1))*0.52917721, 0.0, 0.0 ! lattice vectors
    write(outdbg2, '( 3f12.6 )') 0.0, (bounds(2,2)-bounds(1,2))*0.52917721, 0.0 ! lattice vectors
    write(outdbg2, '( 3f12.6 )') 0.0, 0.0, (bounds(2,3)-bounds(1,3))*0.52917721 ! lattice vectors
    do ity = 1, size(psp)
      write(outdbg2, '("   ", A, "  ")', advance='no') psp(ity)%name
    enddo
    write(outdbg2, '("  Fr ")')
    do ity = 1, size(psp)
      write(outdbg2, '(i5)', advance='no') psp(ity)%natmi
    enddo
    write(outdbg2, '(i5)') n_intp
    write(outdbg2, '(" Cartesian ")')
    do ity = 1, size(psp)
      do iat=1, psp(ity)%natmi
        if (sqrt(xmov*xmov+ymov*ymov+zmov*zmov) > 0.1) then
          write(outdbg2,'(3f12.6)') \
             (psp(ity)%ratm(1,iat)-xmov+0.5*(bounds(2,1)-bounds(1,1)))*0.52917721,\
             (psp(ity)%ratm(2,iat)-ymov+0.5*(bounds(2,2)-bounds(1,2)))*0.52917721,\
             (psp(ity)%ratm(3,iat)-zmov+0.5*(bounds(2,3)-bounds(1,3)))*0.52917721
        else
          write(outdbg2,'(3f12.6)') ( (psp(ity)%ratm(jj,iat) + \
          0.5*(bounds(2,jj)-bounds(1,jj)))*0.52917721, jj=1,3)
        endif
      enddo
    enddo
    do ipt = 1, n_intp
      write(outdbg2,'(3f12.6)') (( fullgrid(jj,intp(ipt)) + \
         0.5*(bounds(2,jj)-bounds(1,jj)) )*0.52917721, jj=1,3)
      !write(outdbg2,'(3f12.6)') (( newpts(1:3, ipt) + \
      !  0.5*(bounds(2,jj)-bounds(1,jj)) )*0.52917721, jj=1,3)
    enddo ! ipt
    close(outdbg2)

 endif
 !
 deallocate(fullgrid)
 deallocate(icenter)

end subroutine cvt_parallel_sym
