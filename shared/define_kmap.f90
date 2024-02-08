!===================================================================
!
! Determine the shape of kernel matrices: which quadruplets of
! orbitals belong to the current representation (irp). Input
! arrays mapcol_in and maprow_in contain the list of allowed quadruplets
! on the columns and rows of the kernel, respectively. This
! subroutine goes through those lists and selects the quadruplets
! that belong to representation irp. The select quadruplets are
! stored in kernel structure.
!
! OUTPUT:
!    kernel%nn
!    kernel%ncol
!    kernel%ncol_up
!    kernel%col
!    kernel%nrow
!    kernel%nrow_up
!    kernel%row
!
! Copyright (C) 2009 Murilo L. Tiago, http://users.ices.utexas.edu/~mtiago
! This file is part of RGWBS. It is distributed under the GPL v1.
!
!-------------------------------------------------------------------
subroutine define_kmap(syms,kpt,kernel,irp,ncol_in,mapcol_in, &
  nrow_in,maprow_in,ncol_up,nrow_up)

  use typedefs
  implicit none

  ! arguments
  type (symmetries), intent(in) :: syms
  type (kptinfo), intent(in) :: kpt
  type(kernelinfo), intent(inout) :: kernel
  integer, intent(in) :: irp, ncol_in, nrow_in, ncol_up, nrow_up
  integer, intent(in) :: mapcol_in(4,ncol_in), maprow_in(4,nrow_in)

  ! local variables
  integer :: ii, jj, iv, ic, ik, jk, rprod, isp

  !-------------------------------------------------------------------
  ! Start with columns of the kernel.
  !
  kernel%ncol_up = 0
  jj = 0
  do ii = 1, ncol_in
     iv = mapcol_in(1,ii)
     ic = mapcol_in(2,ii)
     jk = mapcol_in(3,ii)
     ik = mapcol_in(4,ii)
     isp = 1
     if (ii > ncol_up) isp = 2
     rprod = syms%prod(kpt%wfn(isp,jk)%irep(iv),kpt%wfn(isp,ik)%irep(ic))
     if (rprod == irp) jj = jj + 1
     if ( ii <= ncol_up ) kernel%ncol_up = jj
  enddo
  kernel%ncol = jj
  if ( kernel%ncol /= 0 ) then
     allocate(kernel%col(4,kernel%ncol))
     jj = 0
     do ii = 1, ncol_in
        iv = mapcol_in(1,ii)
        ic = mapcol_in(2,ii)
        jk = mapcol_in(3,ii)
        ik = mapcol_in(4,ii)
        isp = 1
        if (ii > ncol_up) isp = 2
        rprod = syms%prod(kpt%wfn(isp,jk)%irep(iv),kpt%wfn(isp,ik)%irep(ic))
        if (rprod == irp) then
           jj = jj + 1
           kernel%col(:,jj) = mapcol_in(:,ii)
        endif
     enddo
  endif
  !-------------------------------------------------------------------
  ! Now, work with rows.
  !
  kernel%nrow_up = 0
  jj = 0
  do ii = 1, nrow_in
     iv = maprow_in(1,ii)
     ic = maprow_in(2,ii)
     jk = maprow_in(3,ii)
     ik = maprow_in(4,ii)
     isp = 1
     if (ii > nrow_up) isp = 2
     rprod = syms%prod(kpt%wfn(isp,jk)%irep(iv),kpt%wfn(isp,ik)%irep(ic))
     if (rprod == irp) jj = jj + 1
     if ( ii <= nrow_up ) kernel%nrow_up = jj
  enddo
  kernel%nrow = jj
  if ( kernel%nrow /= 0 ) then
     allocate(kernel%row(4,kernel%nrow))
     jj = 0
     do ii = 1, nrow_in
        iv = maprow_in(1,ii)
        ic = maprow_in(2,ii)
        jk = maprow_in(3,ii)
        ik = maprow_in(4,ii)
        isp = 1
        if (ii > nrow_up) isp = 2
        rprod = syms%prod(kpt%wfn(isp,jk)%irep(iv),kpt%wfn(isp,ik)%irep(ic))
        if (rprod == irp) then
           jj = jj + 1
           kernel%row(:,jj) = maprow_in(:,ii)
        endif
     enddo
  endif
  kernel%nn = 0

end subroutine define_kmap
!===================================================================
