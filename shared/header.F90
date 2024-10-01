!===================================================================
!
! Initialize processors, open files and printout header information.
!
! Copyright (C) 2009 Murilo L. Tiago, http://users.ices.utexas.edu/~mtiago
! This file is part of RGWBS. It is distributed under the GPL v1.
!
!-------------------------------------------------------------------
subroutine header(namelabel)

#ifdef _CUDA
! Need to include a couple of modules:
!   cublas: required to use generic BLAS interface
!   cudafor: required to use CUDA runtime API routines (e.g. cudaDeviceSynchronize)
!            not explicitly required if file has *.cuf suffix
  use cublas
  use cudafor
#endif
  use myconstants
  use mpi_module
  implicit none
  !
#ifdef MPI
  include 'mpif.h'
#endif
  ! arguments
  ! name of program
  character(len=*), intent(in) :: namelabel

  ! local variables
  character(len=26) :: datelabel
  real(dp) :: tsec(2)
#ifdef MPI
  integer :: info
#endif

  !-------------------------------------------------------------------
  ! Initialize world structure for MPI involving all processors.
  !
  peinf%num = 1
  peinf%mygr = 0
#ifdef MPI
  call MPI_INIT(info)
  if (info /= MPI_SUCCESS) then
    write (6, *) 'MPI initialization failed!'
    stop
  end if
  peinf%comm = MPI_COMM_WORLD
  call MPI_COMM_RANK(peinf%comm, peinf%inode, info)
  call MPI_COMM_SIZE(peinf%comm, peinf%npes, info)
  call MPI_COMM_GROUP(peinf%comm, peinf%handle, info)
#else
  peinf%npes = 1
  peinf%inode = 0
#endif

  ! Define master PE (peinf%master is local).
  peinf%masterid = 0
  if (peinf%inode == peinf%masterid) then
    peinf%master = .true.
  else
    peinf%master = .false.
  end if

  ! Initialize clocks.
  call timacc(1, 1, tsec)
  call stopwatch(.false., ' ')

  ! Write header.
  call get_date(datelabel)
  if (peinf%master) then
    write (6, '(/,a)') repeat('=', 65)
    write (6, '(/,1x,a,3x,a,a)') namelabel, datelabel, ' UTC'
#ifdef MPI
    write (6, '(/,a,i4,a,/,a,/)') ' RUNNING ON ', &
      peinf%npes, ' PROCESSORS', ' Running MPI version (parallel)'
#else
    write (6, '(/,a,/)') ' Running serial version (no MPI)'
#endif
    write (6, '(a)') 'Version date : VERSIONDATE cut-offs'
    write (6, '(a)') 'architecture = MACH'
    write (6, '(a)') 'pre-processing options = CPPOPT'
    write (6, '(a)') 'compilation options = OPTS'
    write (6, '(/,a,/)') repeat('=', 65)
  end if

  return
end subroutine header
!===================================================================
