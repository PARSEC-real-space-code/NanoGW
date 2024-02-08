!===================================================================
!
! Read checkpoint file sigma.chkpt.dat and convert it to ASCII format.
! Output file has lots of comments. All energy energy quantities are
! given in eV.
!
! INPUT: sigma.chkpt.dat
! OUTPUT: sigma.chkpt.dat.asc
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
program chkpt_bin_asc

  implicit none

  integer, parameter :: dp = kind(1.0d0)
  integer, parameter :: dpc = kind((1.0d0,1.0d0))
  integer, parameter :: XC_GW = 1               ! GW exchange-correlation
  real(dp), parameter :: ryd = 13.60569173_dp
  real(dp), parameter :: invpi = 0.318309886183791_dp

  integer, parameter :: iunit = 60, ounit = 70
  integer :: isp, ii, jj, ien, ik, chkpt, nspin, nkpt, xcflag
  integer :: ndiag_s, ndiag, noffd_s, noffd, nen, nmap, static_type
  real(dp) :: deltae
  real(dp), allocatable :: sigma(:,:,:), sigma_vxc(:,:,:), e0(:), eqp(:)
  integer, dimension(:) , allocatable :: map, diag, off1, off2

  complex(dpc) :: z_spectr
  complex(dpc), allocatable :: zsigma(:,:,:)
  complex(dpc), allocatable :: static_remainder(:)
  character(len=26) :: spectral_file
  character(len=26) :: sigma_file
  character(len=1)  :: arg
  integer :: narg, arglength, argstatus

  open(iunit,file='sigma.chkpt.dat',form='unformatted')
  open(ounit,file='sigma.chkpt.dat.asc',form='formatted')
  read(iunit) chkpt, nspin, nkpt, xcflag
  write(ounit,*) chkpt, nspin

  !write(6,*) " What type of static remainder is this? "
  !read(*,'(i8)') static_type

  narg = command_argument_count()
  !write(6, *) "narg = ", narg
  if (narg == 0) then
    static_type = 2 ! The default value
  else if (narg .ge. 1) then
    call get_command_argument(1, arg, arglength, argstatus) 
    if (argstatus > 0) then
      ! problem
      write(6,*) "Warning: the 1st argment has problem, use the default static_type = 2."
      static_type = 2
    else
      write(static_type,'(i1)') arg
      if ( static_type .ne. 1 .or. static_type .ne. 2 ) then
        write(6,*) " static_type can only be 1 or 2, you provide static_type = ", static_type
        write(6,*) " We use the default static_type = 2."
        static_type = 2
      endif
    endif
  endif
  write(6, *) "static_type = ", static_type

  do ik = 1, nkpt
     do isp = 1, nspin
        read(iunit) nmap
        read(iunit) ndiag_s, ndiag
        read(iunit) noffd_s, noffd
        read(iunit) nen, deltae
        deltae = deltae * ryd

        write(ounit,*) nen, deltae, nmap
        write(ounit,*) ndiag_s, ndiag, noffd_s, noffd

        if (nmap > 0) then
           allocate(map(nmap))
           read(iunit) (map(ii),ii=1,nmap)

           if (ndiag_s > 0) then
              allocate(diag(ndiag_s))
              allocate(static_remainder(ndiag_s))
              read(iunit) (diag(ii),ii=1,ndiag_s)

              allocate(sigma(ndiag_s,1,1))
              sigma = 0.d0
              allocate(zsigma(ndiag_s,1,1))
              zsigma = cmplx(0.d0,0.d0)

              allocate(e0(ndiag_s))
              read(iunit) (e0(ii),ii=1,ndiag_s)
              read(iunit)
              e0 = e0 * ryd
              write(ounit,14) 'Initial (DFT) energy',isp,'Energy'
              write(ounit,24) (map(diag(ii)),e0(ii),ii=1,ndiag_s)

              read(iunit) (zsigma(ii,1,1),ii=1,ndiag_s)
              sigma = real(zsigma,dp) * ryd
              write(ounit,14) 'Initial self-energy (VXC) diagonal',isp, &
                   '<n|V_xc|n>'
              write(ounit,24) (map(diag(ii)),sigma(ii,1,1),ii=1,ndiag_s)
              if (ndiag > 0) then
                 allocate(sigma_vxc(nen,ndiag,2))
                 sigma_vxc = 0.d0
                 do ii = 1, ndiag
                    sigma_vxc(:,ii,1) = e0(ii) - sigma(ii,1,1)
                 enddo
              endif

              read(iunit) (zsigma(ii,1,1),ii=1,ndiag_s)
              sigma = real(zsigma,dp) * ryd
              write(ounit,14) 'Bare exchange diagonal',isp,'<n|Sigma_x|n>'
              write(ounit,24) (map(diag(ii)),sigma(ii,1,1),ii=1,ndiag_s)
              if (ndiag > 0) then
                 do ii = 1, ndiag
                    sigma_vxc(:,ii,1) = sigma_vxc(:,ii,1) + sigma(ii,1,1)
                 enddo
              endif

              read(iunit) (zsigma(ii,1,1),ii=1,ndiag_s)
              sigma = real(zsigma,dp) * ryd
              static_remainder = zsigma(:,1,1) * ryd
              write(ounit,14) 'Static correlation diagonal',isp, &
                   '<n|Sigma_c|n>'
              write(ounit,24) (map(diag(ii)),sigma(ii,1,1),ii=1,ndiag_s)
              if (ndiag > 0) then
                 do ii = 1, ndiag
                    sigma_vxc(:,ii,1) = sigma_vxc(:,ii,1) + &
                    sigma(ii,1,1)/real(static_type,dp)
                 enddo
              endif

              if (xcflag == XC_GW) then
                 read(iunit) (zsigma(ii,1,1),ii=1,ndiag_s)
                 sigma = real(zsigma,dp) * ryd
                 write(ounit,14) 'Static vertex diagonal',isp, &
                      '<n|Sigma_g|n>'
                 write(ounit,24) (map(diag(ii)),sigma(ii,1,1),ii=1,ndiag_s)
                 if (ndiag > 0) then
                    do ii = 1, ndiag
                       sigma_vxc(:,ii,1) = sigma_vxc(:,ii,1) + &
                       sigma(ii,1,1)/real(static_type,dp)
                    enddo
                 endif
              endif

              deallocate(sigma)
              deallocate(zsigma)

              if (ndiag > 0 .and. xcflag == XC_GW) then
                 allocate(sigma(nen,ndiag,2))
                 sigma = 0.d0
                 allocate(zsigma(nen,ndiag,2))
                 zsigma = cmplx(0.d0,0.d0)

                 allocate(eqp(nen))

                 read(iunit) (((zsigma(ien,ii,jj),ien=1,nen),ii=1,ndiag),jj=1,2)
                 sigma = real(zsigma,dp) * ryd
                 write(ounit,16) 'Correlation',isp,'<n|Sigma_c(E)|n>, complex'
                 do ii = 1, ndiag
                    do ien = 1, nen
                       eqp(ien) = e0(ii) - deltae*0.5d0 + &
                            deltae / real(nen - 1,dp) * real(ien - 1,dp)
                    enddo
                    write(ounit,28) (map(diag(ii)), &
                         eqp(ien),sigma(ien,ii,:),ien=1,nen)
                 enddo

                 sigma_vxc = sigma_vxc + sigma

                 read(iunit) (((zsigma(ien,ii,jj),ien=1,nen),ii=1,ndiag),jj=1,2)
                 sigma = real(zsigma,dp) * ryd
                 write(ounit,16) 'Screened exchange',isp, &
                      '<n|Sigma_sx(E)|n>, complex'
                 do ii = 1, ndiag
                    do ien = 1, nen
                       eqp(ien) = e0(ii) - deltae*0.5d0 + &
                            deltae / real(nen - 1,dp) * real(ien - 1,dp)
                    enddo
                    write(ounit,28) (map(diag(ii)), &
                         eqp(ien),sigma(ien,ii,:),ien=1,nen)
                 enddo

                 read(iunit) (((zsigma(ien,ii,jj),ien=1,nen),ii=1,ndiag),jj=1,2)
                 sigma = real(zsigma,dp) * ryd
                 write(ounit,16) 'Vertex',isp,'<n|Sigma_g(E)|n>, complex'
                 do ii = 1, ndiag
                    do ien = 1, nen
                       eqp(ien) = e0(ii) - deltae*0.5d0 + &
                            deltae / real(nen - 1,dp) * real(ien - 1,dp)
                    enddo
                    write(ounit,28) (map(diag(ii)), &
                         eqp(ien),sigma(ien,ii,:),ien=1,nen)
                 enddo

                 sigma_vxc = sigma_vxc + sigma

                 write(ounit,16) 'Sigma - V_xc + E0',isp, &
                      '<n|Sigma(E) - V_xc + E0|n>, complex'
                 do ii = 1, ndiag
                    write(sigma_file,*) map(diag(ii))
                    if(isp.eq.1) then
                        sigma_file = "Eqp"//trim(adjustl(sigma_file))//"-spin1.dat"
                    else
                        sigma_file = "Eqp"//trim(adjustl(sigma_file))//"-spin2.dat"
                    endif
                    open(44,file=trim(sigma_file),form="formatted",status="unknown")
                    do ien = 1, nen
                       eqp(ien) = e0(ii) - deltae*0.5d0 + &
                            deltae / real(nen - 1,dp) * real(ien - 1,dp)
                    enddo
                    write(ounit,28) (map(diag(ii)), &
                         eqp(ien),sigma_vxc(ien,ii,:),ien=1,nen)
                    !write(44, '(A)') "band  E  E0+Re[Sigma(E)]-Vxc Im[Sigma(E)]  No staic remainder:E0+Re(Sigma(E))-Vxc Im(Sigma(E))"
                    write(44,'(i5,5g20.10)') (map(diag(ii)), &
                         eqp(ien),sigma_vxc(ien,ii,1), &
                         sigma_vxc(ien,ii,2),&
                         sigma_vxc(ien,ii,1)-dble(static_remainder(ii)),&
                         sigma_vxc(ien,ii,2)-aimag(static_remainder(ii)),&
                         ien=1,nen)
                    close(44)
                 enddo

                 write(ounit,*) ' Spectral Function, (1/Pi) * Im{1/(E - ',  &
                      'Sigma + V_xc - E0 )}, diagonal, eV^-1, spin ',isp
                 write(ounit,*) '   n',repeat(' ',9),'E',repeat(' ',14), &
                      '|<n|Sigma(E) - V_xc + E0|n>|'
                 do ii = 1, ndiag
                    write(spectral_file,*) map(diag(ii))
                    if(isp.eq.1) then
                       spectral_file = "Spectral"//trim(adjustl(spectral_file))//"-spin1.dat"
                    else
                       spectral_file = "Spectral"//trim(adjustl(spectral_file))//"-spin2.dat"
                    endif
                    open(45,file=trim(spectral_file),form="formatted",status="unknown")
                    do ien = 1, nen
                       eqp(ien) = e0(ii) - deltae*0.5d0 + &
                            deltae / real(nen - 1,dp) * real(ien - 1,dp)
                       z_spectr = cmplx(sigma_vxc(ien,ii,1),sigma_vxc(ien,ii,2))
                       z_spectr = eqp(ien) - z_spectr
                       z_spectr = invpi / z_spectr
                       sigma_vxc(ien,ii,1) = abs(aimag(z_spectr))
                    enddo
                    write(ounit,'(i5,2g20.10)') (map(diag(ii)), &
                         eqp(ien),sigma_vxc(ien,ii,1),ien=1,nen)
                    write(45,'(i5,2g20.10)') (map(diag(ii)), &
                         eqp(ien),sigma_vxc(ien,ii,1),ien=1,nen)
                    close(45)
                 enddo

                 deallocate(eqp)
                 deallocate(sigma)
                 deallocate(zsigma)
                 deallocate(sigma_vxc)
              endif
              deallocate(e0)
              deallocate(diag)
              deallocate(static_remainder)
           endif

           if (noffd_s > 0) then
              allocate(off1(noffd_s))
              allocate(off2(noffd_s))
              read(iunit) (off1(ii),ii=1,noffd_s)
              read(iunit) (off2(ii),ii=1,noffd_s)
              allocate(sigma(noffd_s,1,1))
              sigma = 0.d0
              allocate(zsigma(noffd_s,1,1))
              zsigma = cmplx(0.d0,0.d0)

              read(iunit) (zsigma(ii,1,1),ii=1,noffd_s)
              sigma = real(zsigma,dp) * ryd
              write(ounit,18) 'Initial self-energy (VXC)',isp,'<n1|Sigma_x|n2>'
              write(ounit,34) (map(off1(ii)),map(off2(ii)), &
                   sigma(ii,1,1),ii=1,noffd_s)

              read(iunit) (zsigma(ii,1,1),ii=1,noffd_s)
              sigma = real(zsigma,dp) * ryd
              write(ounit,18) 'Bare exchange',isp,'<n1|Sigma_x|n2>'
              write(ounit,34) (map(off1(ii)),map(off2(ii)), &
                   sigma(ii,1,1),ii=1,noffd_s)

              read(iunit) (zsigma(ii,1,1),ii=1,noffd_s)
              sigma = real(zsigma,dp) * ryd
              write(ounit,18) 'Static correlation',isp,'<n1|Sigma_c|n2>'
              write(ounit,34) (map(off1(ii)),map(off2(ii)), &
                   sigma(ii,1,1),ii=1,noffd_s)

              read(iunit) (zsigma(ii,1,1),ii=1,noffd_s)
              sigma = real(zsigma,dp) * ryd
              write(ounit,18) 'Static vertex',isp,'<n1|Sigma_g|n2>'
              write(ounit,34) (map(off1(ii)),map(off2(ii)), &
                   sigma(ii,1,1),ii=1,noffd_s)

              deallocate(sigma)
              deallocate(zsigma)
              if (noffd > 0 .and. xcflag == XC_GW) then
                 allocate(sigma(2,noffd,2))
                 sigma = 0.d0
                 allocate(zsigma(2,noffd,2))
                 zsigma = cmplx(0.d0,0.d0)

                 read(iunit) ((zsigma(:,ii,jj),ii=1,noffd),jj=1,2)
                 sigma = real(zsigma,dp) * ryd
                 write(ounit,18) 'Correlation',isp,'<n1|Sigma_g(E_1)|n2>, complex'
                 write(ounit,36) (map(off1(ii)),map(off2(ii)), &
                      (sigma(1,ii,jj),jj=1,2),ii=1,noffd)
                 write(ounit,36) (map(off2(ii)),map(off1(ii)), &
                      (sigma(2,ii,jj),jj=1,2),ii=1,noffd)

                 read(iunit) ((zsigma(:,ii,jj),ii=1,noffd),jj=1,2)
                 sigma = real(zsigma,dp) * ryd
                 write(ounit,18) 'Vertex',isp,'<n1|Sigma_g(E_1)|n2>, complex'
                 write(ounit,36) (map(off1(ii)),map(off2(ii)), &
                      (sigma(1,ii,jj),jj=1,2),ii=1,noffd)
                 write(ounit,36) (map(off2(ii)),map(off1(ii)), &
                      (sigma(2,ii,jj),jj=1,2),ii=1,noffd)

                 deallocate(sigma)
                 deallocate(zsigma)
              endif
              deallocate(off2)
              deallocate(off1)
           endif
           deallocate(map)
        endif

     enddo  ! isp = 1, nspin

  enddo ! ik = 1, nkpt

  close(iunit)
  close(ounit)

14 format(2x,a,', eV, spin ',i2,/,4x,'n',6x,a)
16 format(2x,a,' diagonal, eV, spin ',i2,/,4x,'n',9x,'E',18x,a)
18 format(2x,a,' off-diagonal, eV, spin ',i2,/,3x,'n1',4x,'n2',4x,a)
24 format(i5,g20.10)
28 format(i5,3g20.10)
34 format(2i5,g20.10)
36 format(2i5,2g20.10)

end program chkpt_bin_asc
!===================================================================
