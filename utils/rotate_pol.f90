!===================================================================
!
! Read polarizability eigenvectors and perform the DFT -> QP rotation.
!
! Input file:   pol_in.dat   : input eigenvectors
!                              (same as pol_diag.dat or bse_diag.dat)
! Output file : pol_out.dat  : output eigenvectors
!
! Copyright (C) 2009 Murilo Tiago, Univ. of Texas, Austin, TX, USA
! mtiago@ices.utexas.edu
!
! First version written by Murilo Tiago, Oak Ridge National Laboratory, February 2008.
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
program rotate_pol

  implicit none
  
  ! constants
  integer, parameter :: dp = kind(1.0d0)
  integer, parameter :: dpc = kind((1.0d0,1.0d0))
  real(dp), parameter :: zero = 0.0_dp
  real(dp), parameter :: ryd = 13.60569173_dp
  
  ! LDA wavefunctions
  type wavefunction
     ! number of LDA states (same as syms%nstate), read from wfn.dat
     integer :: nstate
     ! sorting index
     integer, pointer :: indx(:), inv_indx(:)
     ! rotation matrix
     real(dp), pointer :: u_rot(:,:)
  end type wavefunction

  type (wavefunction), pointer :: wfn(:)

  character (len=26) :: datelabel  !  date/time tag
  integer :: ii, itr, jtr, isp, irp, iv, ic, jv, jc
  integer :: nstt(2), nstate_in, neig, ntr, n_up, nspin, nrep
  real(dp) :: xsum, xstep
  integer, allocatable :: isort(:), tr_old(:,:), tr_new(:,:)
  real(dp), dimension(:), allocatable :: eqp, eig, v_old, v_new
  real(dp), dimension(:,:), allocatable :: rtmp
  integer, parameter :: infile = 66, outfile = 77, qpfile = 25

  !---------------------------------------------------------------
  ! Read header of pol_in.dat and initialize file pol_out.dat.
  !
  open(infile,file='pol_in.dat',form='unformatted',status='old',iostat=ii)
  write(6,*) ' Performing orbital rotation. Output file is pol_out.dat.'
  read(infile) datelabel
  write(6,*) ' File pol_in.dat written on ',datelabel
  read(infile) nstate_in, nrep

  if (ii /= 0) then
     write(6,*) 'File pol_in.dat not found.'
     stop
  endif

  !---------------------------------------------------------------
  ! Open file qp_mtxel.dat with quasi-particle rotation matrices and
  ! read matrices.
  !
  open(qpfile,file='qp_mtxel.dat',form='formatted')
  read(qpfile,*) nspin,(nstt(ii),ii=1,nspin)
  allocate(wfn(nspin))

  do isp = 1, nspin
     if (nstt(isp) /= nstate_in) then
        write(6,*) ' ERROR! level mismatch in pol_in.dat: ', &
             nstt(isp), nstate_in, isp,' STOP.'
        stop
     endif
     !
     ! Read rotation matrix for this spin component.
     ! It could be that QP levels do not have the same ordering as DFT
     ! levels. We must sort them and update the eigenvalues in wfn structure.
     !
     wfn(isp)%nstate = nstt(isp)
     allocate(wfn(isp)%u_rot(nstt(isp),nstt(isp)))
     allocate(eqp(nstt(isp)))
     do ii = 1, nstt(isp)
        read(qpfile,*) eqp(ii)
     enddo
     do ii = 1, nstt(isp)
        do itr = 1, nstt(isp)
           read(qpfile,*) wfn(isp)%u_rot(itr,ii)
        enddo
     enddo
     allocate(wfn(isp)%indx(nstt(isp)))
     call quicksort(nstt(isp),eqp,wfn(isp)%indx)
     allocate(wfn(isp)%inv_indx(nstt(isp)))
     do ii = 1, nstt(isp)
        wfn(isp)%inv_indx(wfn(isp)%indx(ii)) = ii
     enddo

     deallocate(eqp)

     allocate(rtmp(nstt(isp),nstt(isp)))
     rtmp = wfn(isp)%u_rot
     do ii = 1, nstt(isp)
        wfn(isp)%u_rot(:,ii) = rtmp(:,wfn(isp)%indx(ii))
     enddo
     deallocate(rtmp)

  enddo
  close(qpfile)

  !---------------------------------------------------------------
  ! Initialize file pol_out.dat.
  !
  open(outfile,file='pol_out.dat',form='unformatted',status='unknown')
  call get_date(datelabel)
  write(6,*) ' File pol_out.dat written on ',datelabel
  write(outfile) datelabel
  write(outfile) nstate_in, nrep

  !-------------------------------------------------------------------
  !  Read data in pol_in.dat and update polarizability eigenvectors.
  !
  do irp = 1, nrep
     read(infile) ii, neig, ntr, n_up
     write(outfile) ii, neig, ntr, n_up
     write(6,*) ' Read representation ',irp, ii, neig, ntr, n_up

     allocate(tr_old(2,ntr))
     allocate(tr_new(2,ntr))
     read(infile) ((tr_old(ii,itr),ii=1,2),itr=1,ntr)
     do itr = 1, n_up
        tr_new(1,itr) = wfn(1)%inv_indx(tr_old(1,itr))
        tr_new(2,itr) = wfn(1)%inv_indx(tr_old(2,itr))
     enddo
     do itr = n_up + 1, ntr
        tr_new(1,itr) = wfn(2)%inv_indx(tr_old(1,itr))
        tr_new(2,itr) = wfn(2)%inv_indx(tr_old(2,itr))
     enddo
     ! Must reorder transitions according to ascending (v,c) order.
     allocate(v_old(ntr))
     allocate(isort(ntr))
     xstep = maxval(tr_new(1,:))
     xstep = 1.d0/xstep
     do itr = 1, ntr
        v_old(itr) = tr_new(2,itr) + tr_new(1,itr)*xstep
     enddo
     call quicksort(ntr,v_old,isort)
     deallocate(v_old)
     write(outfile) ((tr_new(ii,isort(itr)),ii=1,2),itr=1,ntr)

     allocate(eig(neig))
     read(infile) (eig(ii),ii=1,neig)
     write(outfile) (eig(ii),ii=1,neig)
     deallocate(eig)

     allocate(v_old(ntr))
     allocate(v_new(ntr))
     do ii = 1, neig
        read(infile) (v_old(itr),itr=1,ntr)
        do itr = 1, ntr
           xsum = zero
           iv = tr_new(1,itr)
           ic = tr_new(2,itr)
           do jtr = 1, n_up
              jv = tr_old(1,jtr)
              jc = tr_old(2,jtr)
              xsum = xsum + &
                   wfn(1)%u_rot(jv,iv) * wfn(1)%u_rot(jc,ic) * v_old(jtr)
           enddo
           do jtr = n_up + 1, ntr
              jv = tr_old(1,jtr)
              jc = tr_old(2,jtr)
              xsum = xsum + &
                   wfn(2)%u_rot(jv,iv) * wfn(2)%u_rot(jc,ic) * v_old(jtr)
           enddo
           v_new(isort(itr)) = xsum
        enddo
        write(outfile) (v_new(itr),itr=1,ntr)
     enddo
     deallocate(v_old)
     deallocate(v_new)

     deallocate(tr_old)
     deallocate(tr_new)
     deallocate(isort)
  enddo

  close(infile)
  close(outfile)

end program rotate_pol
!===================================================================
