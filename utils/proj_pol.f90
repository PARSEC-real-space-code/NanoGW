!===================================================================
!
! Projects selected excitations onto selected cv space. The user
! must select a number of TDLDA (or BSE) excitations and a number of
! occupied and unoccupied orbitals. This code prints out the
! projection of the selected excitations on the subspace spanned by
! the selected occupied and unoccupied orbitals.
!
! INPUT FILES:
!     pol_diag.dat (bse_diag.dat can also be used once it is renamed pol_diag.dat)
!
! SCREEN INPUT:
!     ns : number of excitations to be projected (integer)
!     eigv : list of excitation energies (scalar, from 1 to ns)
!     nrp : list of representations of excitations (scalar, from 1 to ns)
!     nv : number of occupied orbitals in the cv space (integer)
!     nc : number of unoccupied orbitals in the cv space (integer)
!     iv : list of occupied orbitals (integer, from 1 to nv)
!     ic : list of unoccupied orbitals (integer, from 1 to nc)
!
! Copyright (C) 2009 Murilo Tiago, Univ. of Texas, Austin, TX, USA
! mtiago@ices.utexas.edu
!
! First version written by Murilo Tiago, Oak Ridge National Laboratory, September 2007.
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
program proj_pol

  implicit none

  ! constants
  integer, parameter :: dp = kind(1.0d0)
  integer, parameter :: dpc = kind((1.0d0,1.0d0))
  real(dp), parameter :: ryd = 13.60569173_dp

  character (len=26) :: datelabel
  logical :: vflag, cflag, sflag
  integer :: n1, nrep, irp, ii, jj, ll, jv, jc, icount, nv, nc, ns, &
       icplx, iq, nqpt
  real(dp) :: psum, norm
  integer, allocatable :: iv(:), ic(:), nrp(:)
  real(dp), allocatable :: eigv(:)
  real(dp), parameter :: e_tol = 1.d-8

  ! polarizability data
  integer :: ntr
  integer, allocatable :: tr(:,:)
  real(dp), allocatable :: eig(:), dv(:)
  complex(dpc), allocatable :: zv(:)

  !-------------------------------------------------------------------
  ! Initialize eigenvector file.
  !
  open(10,file='pol_diag.dat',form='unformatted',status='old')
  read(10) datelabel
  write(6,*) ' File written on ', datelabel
  read(10) n1, nrep, nqpt, icplx
  write(6,*) ' This file has ', n1, ' Kohn-Sham orbitals and ', nrep, &
       ' representations.'
  if (icplx == 2) write(6,*) 'File has real data.'
  if (icplx == 1) write(6,*) 'File has complex data, ', nqpt, ' q-vectors.'

  !-------------------------------------------------------------------
  ! Read specifications for the projection.
  !
  write(6,*) 'How many excitations to project ?'
  read(5,*) ns
  allocate(nrp(ns))
  allocate(eigv(ns))
  write(6,*) 'List of excitation energies (in eV) and representations ?'
  do ii = 1, ns
     read(5,*) eigv(ii), nrp(ii)
  enddo
  eigv = eigv / ryd
  write(6,*) 'Number of occupied and unoccupied levels in the cv space ?'
  read(5,*) nv, nc
  allocate(iv(nv))
  allocate(ic(nc))
  write(6,*) 'List of occupied levels in the cv space ?'
  read(5,*) (iv(ii),ii=1,nv)
  write(6,*) 'List of unoccupied levels in the cv space ?'
  read(5,*) (ic(ii),ii=1,nc)

  !-------------------------------------------------------------------
  ! Search for selected excitations in eigenvector file.
  !
  do ll = 1, nqpt * nrep
     read(10) irp, iq, ntr
     ! tr stores the list of transitions: n-th single-electron transition in
     ! representation irp is a transition from occupied K-S orbital
     ! tr(n,1) to unoccupied K-S orbital tr(n,2).
     ! n runs from 1 to ntr (=number of transitions for this representation.
     ! Warning: the shape of tr has changed recently from tr(1:2,1:ntr) to
     ! tr(1:ntr,1:2).
     allocate(tr(ntr,2))
     ! eig and v store the list of TDLDA eigenvalues and corresponding
     ! eigenvectors.
     ! Notice that eigenvalues are listed from top to bottom.
     ! Each eigenvector is written on a separate line.
     ! s-th line in the file contains this vector:
     !
     ! zv(i)     = X^s_{vc} ( {e_c - e_v} / omega_s )^{1/2}
     ! tr(i,1)   = v
     ! tr(i,2)   = c
     ! eig(s)    = omega_s, in Ry
     !
     ! Coefficients X^s_{vc} are defined in Tiago & Chelikowsky,
     ! PRB 73, 205334, (2006), Eq. 3.
     !
     allocate(eig(ntr))
     allocate(zv(ntr))
     read(10) (tr(ii,:),ii=1,ntr)
     read(10) (eig(ii),ii=1,ntr)
     do ii = 1, ntr
        sflag = .false.
        do jj = 1, ns
           if (abs(eigv(jj) - eig(ii)) < e_tol .and. &
                nrp(jj) == irp) sflag = .true.
        enddo
        if (sflag) then
           if (icplx == 1) then
              allocate(dv(ntr))
              read(10) (dv(jj),jj=1,ntr)
              zv = dv              
              deallocate(dv)
           else
              read(10) (zv(jj),jj=1,ntr)
           endif
           write(6,'(a,i2,a,f20.10,a)') ' Found excitation representation ', &
                irp, ' energy = ', eig(ii) * ryd, ' eV.'
           psum = 0.d0
           icount = 0
           do jj = 1, ntr
              vflag = .false.
              cflag = .false.
              do jv = 1, nv
                 if (iv(jv) == tr(jj,1)) vflag = .true.
              enddo
              do jc = 1, nc
                 if (ic(jc) == tr(jj,2)) cflag = .true.
              enddo
              if (vflag .and. cflag) then
                 icount = icount + 1
                 psum = psum + abs(zv(jj))**2
              endif
           enddo
           norm = dot_product(zv(1:ntr),zv(1:ntr))
           write(6,*) ' projection = ' , psum, ' out of ', norm
           psum = psum / norm
           write(6,*) ' check sum ', icount, ' = ', nv*nc
           write(6,'(a,i5,2f16.8)') ' Projection ', irp, eig(ii) * ryd, psum
        else
           read(10)
        endif
     enddo
     deallocate(tr, eig, zv)
  enddo
  close(10)

  deallocate(iv, ic, nrp, eigv)

end program proj_pol
!===================================================================
