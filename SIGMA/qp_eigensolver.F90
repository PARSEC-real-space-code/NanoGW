!===================================================================
!
! Diagonalize QP Hamiltonians (one Hamiltonian for each spin and
! representation). At the end, master processor broadcasts data.
! This communication step ensures that the output of this subroutine
! remains global.
!
! Quantities modified in this subroutine:
!    q_p%vqp : eigenvectors of the symmetrized Hamiltonian
!    q_p%eqp : eigenvalues of the symmetrized Hamiltonian
!    q_p%hqp : non-Hermitian projection
!
! The non-Hermitian projection is defined as:
!
!    q_p%vqp * q_p%eqp * adjoint(q_p%vqp) - q_p%hqp
!
! If q_p%hqp is input Hermitian, then the non-Hermitian projection is null.
!
! Copyright (C) 2009 Murilo L. Tiago, http://users.ices.utexas.edu/~mtiago
! This file is part of RGWBS. It is distributed under the GPL v1.
!
!-------------------------------------------------------------------
subroutine qp_eigensolver(q_p,verbose,hqp_sym,ierr)

  use typedefs
  use mpi_module
  implicit none
#ifdef MPI
  include 'mpif.h'
#endif

  ! arguments
  ! quasi-particle data for the current k-point, representation and spin
  type(qpinfo), intent(inout) :: q_p
  logical, intent(in) :: &
       verbose, &           ! output flag
       hqp_sym              ! true if QP Hamiltonian is symmetrized
  integer, intent(out) :: ierr   ! error flag

  ! local variables
  integer :: ii, jj, nn
  complex(dpc) :: v_max
  complex(dpc), allocatable :: h_diff(:,:), rtmp(:,:)

  nn = q_p%neig

  ii = 0
  do jj = 1, nn
     if (q_p%hqp(jj,jj) == zzero) ii = ii + 1
  enddo

  if (ii > 0 .and. verbose) then
     write(6,'(/,a,/,a,/)') repeat('*',65), repeat('*',65)
     write(6,*) ' WARNING !!!! There are missing quasi-particle levels.'
     write(6,*) ' Self-energy was not calculated for all electronic levels.'
     write(6,*) ' Quasi-particle eigenvectors may be incorrect.'
     write(6,*) ' Number of missing electronic levels = ', ii
     write(6,'(/,a,/,a,/)') repeat('*',65), repeat('*',65)
  endif

  do ii = 1, nn
     do jj = 1, ii - 1
        if (hqp_sym) then
           q_p%vqp(ii,jj) = half*conjg(q_p%hqp(jj,ii)) + half*q_p%hqp(ii,jj)
        else
           q_p%vqp(ii,jj) = q_p%hqp(ii,jj)
        endif
        q_p%vqp(jj,ii) = conjg(q_p%vqp(ii,jj))
     enddo
     q_p%vqp(ii,ii) = q_p%hqp(ii,ii)
  enddo

  call zeigensolver(.false.,0,1,0,nn,nn,q_p%vqp,q_p%eqp,ierr)
#ifdef MPI
  call MPI_BCAST(q_p%eqp,nn,MPI_DOUBLE_PRECISION,peinf%masterid,peinf%comm,ii)
  call MPI_BCAST(q_p%hqp,nn*nn,MPI_DOUBLE_COMPLEX,peinf%masterid,peinf%comm,ii)
  call MPI_BCAST(q_p%vqp,nn*nn,MPI_DOUBLE_COMPLEX,peinf%masterid,peinf%comm,ii)
  call MPI_BCAST(ierr,1,MPI_INTEGER,peinf%masterid,peinf%comm,ii)
#endif
  if (ierr /= 0) return
  !
  ! Upright eigenvectors. For each eigenvector, set the phase of its highest
  ! component to zero.
  !
  do ii = 1, nn
     v_max = Zzero
     do jj = 1, nn
        if (abs(v_max) < abs(q_p%vqp(jj,ii))) v_max = q_p%vqp(jj,ii)
     enddo
     v_max = v_max/abs(v_max)
     call Zscal(nn,v_max,q_p%vqp(1,ii),1)
  enddo
  !
  ! Calculate non-Hermitian projection.
  !
  allocate(h_diff(nn,nn),stat=ii)
  call alccheck('h_diff','qp_eigensolver',nn*nn,ii)
  h_diff = zzero
  do ii = 1, nn
     h_diff(ii,ii) = -zone * q_p%eqp(ii)
  enddo
  allocate(rtmp(nn,nn),stat=ii)
  call alccheck('rtmp','qp_eigensolver',nn*nn,ii)
  call zgemm('N','C',nn,nn,nn,zone,h_diff,nn,q_p%vqp,nn,zzero,rtmp,nn)
  call zgemm('N','N',nn,nn,nn,zone,q_p%vqp,nn,rtmp,nn,zone,q_p%hqp,nn)

  deallocate(h_diff, rtmp)

end subroutine qp_eigensolver
!===================================================================
