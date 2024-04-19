program test_lanczos
 use typedefs
 implicit none
 !
 ! variable types
 !
 integer, parameter :: dp = kind(1.0d0)
 integer, parameter :: io = 14444, ndim = 104, nz = 50
 integer :: dummy, niter, ii, jj, iz

 type(dmat) :: H
 real(8) :: v0(ndim), eps
 real(8) :: zz(nz)
 complex(8) :: zzi(nz), spectra(nz)
 
 eps = 1.0d-1
 niter = 5
 ! print *, io
 open(unit=io, file='H_aux_v0.txt', form='formatted', status='old')

 H%n = ndim
 allocate (H%h(ndim, ndim))
 ! Read H
 read(io, *) ! skip
 do ii = 1, ndim
   do jj = 1, ndim
     read(io, *) dummy, dummy, H%h(ii, jj)
   enddo
 enddo
 read(io, *) ! skip
 do ii = 1, ndim
   read(io, *) dummy, v0(ii)
 enddo
 read(io, *) ! skip
 do iz = 1, nz
   read(io, *) zz(iz)
 enddo
 close(io)

 zzi = zz + (0.d0, 1.d0)*eps 
 ! print *, "v0 = ", v0
 call lanczos_spectra_isdf(v0, ndim, zzi, nz, niter, H, spectra)
 do iz = 1, nz
   print *, zz(iz), real(spectra(iz)), aimag(spectra(iz))
 enddo

end program test_lanczos
