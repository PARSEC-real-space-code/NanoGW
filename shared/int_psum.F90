#include "../shared/mycomplex.h"
!===================================================================
!
! Global sum routine using the MPI_ALLREDUCE primitive.
!
! Copyright (C) 2009 Murilo L. Tiago, http://users.ices.utexas.edu/~mtiago
! This file is part of RGWBS. It is distributed under the GPL v1.
!
!-------------------------------------------------------------------

subroutine int_psum(ndim, npes, comm, buffer)
  use myconstants
  implicit none
#ifdef MPI
  include 'mpif.h'
#endif
  ! dimension of array, number of PEs, MPI communicator
  integer, intent(in) :: ndim, npes, comm
  ! input: array to be summed, local
  ! output: array summed across all PEs, global
  integer, intent(inout) :: buffer(ndim)

#ifdef MPI
  integer :: work(MAXSIZE_MPI)
  integer :: info, ioff, nblock, mdim, iblock
#endif

  !-------------------------------------------------------------------
  if (npes == 1) return
#ifdef MPI
  if (ndim < 1) return
  ioff = 1
  mdim = MAXSIZE_MPI
  if (ndim > MAXSIZE_MPI) then
    nblock = ndim/MAXSIZE_MPI
    do iblock = 1, nblock
      call MPI_ALLREDUCE(buffer(ioff), work, mdim, &
                         MPI_INTEGER, MPI_SUM, comm, info)
      call Zcopy(mdim, work, 1, buffer(ioff), 1)
      ioff = ioff + mdim
    end do
  end if
  mdim = ndim - ioff + 1
  call MPI_ALLREDUCE(buffer(ioff), work, mdim, &
                     MPI_INTEGER, MPI_SUM, comm, info)
  call Zcopy(mdim, work, 1, buffer(ioff), 1)
#endif

end subroutine int_psum
