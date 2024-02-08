!===================================================================
!
! Print out the numerical sum rule and static polarizability from
! a TDLDA calculation. On output, one gets the ratio between numerical and
! exact sum rules (assuming sum(wfn%occ0) to be the number of occupied
! LDA states in the system). Energy eigenvalues are also writen on disk:
!
! 'eigenvalues_rpa' : RPA-LDA eigenvalues + oscill. strengths
! 'eigenvalues_lda' : TDLDA (polarizability) eigenvalues + oscill. strengths
!
! Copyright (C) 2009 Murilo L. Tiago, http://users.ices.utexas.edu/~mtiago
! This file is part of RGWBS. It is distributed under the GPL v1.
!
!-------------------------------------------------------------------
subroutine sum_rule(rpaonly,nrep,nspin,per,kpt,pol,celvol,xsum,epsinv)

  use typedefs
  implicit none

  ! arguments
  ! true if only RPA eigenvalues are available; false otherwise
  logical, intent(in) :: rpaonly
  ! number of representations, number of spin components, number of periodic
  ! directions (from 0 to 3)
  integer, intent(in) :: nrep, nspin, per
  ! k-point specification
  type (kptinfo), intent(in) :: kpt
  ! polarizability eigenvalues, oscillator strengths
  type (polinfo), dimension(nrep), intent(in) :: pol
  ! volume of unit cell, relevant in bulk systems only
  real(dp), intent(in) :: celvol
  ! averaged sum rule, RPA inverse dielectric function
  real(dp), intent(out) :: xsum, epsinv

  ! local variables
  integer :: ii, i1, i2, k1, k2, icol, ip, isp
  integer :: irp, neig, ioffset
  real(dp) :: rsum(3), zvalence, occdif, spol(3,2), polrpa
  integer, allocatable :: repr(:), indx(:)
  real(dp), allocatable :: ene0(:), ostr0(:,:), ene(:), ostr(:,:)
  complex(dpc), allocatable :: ostrz(:,:)

  !-------------------------------------------------------------------
  ! Initialize arrays. All arrays are collected over representations.
  !
  neig = sum(pol(:)%ntr)
  allocate(ene0(neig))
  allocate(ostr0(neig,3))
  allocate(ostrz(neig,3))
  allocate(ene(neig))
  allocate(ostr(neig,3))
  allocate(repr(neig))
  ene0 = zero
  ene = zero
  ostr0 = zero
  ostrz = zzero
  ostr = zero
  ioffset = 0

  !-------------------------------------------------------------------
  ! Collect transitions over all representations.
  !
  do irp = 1, nrep
     !
     ! RPA sum rule: oscillator strengths are just dipole matrix elements
     !  and energy factors are differences between LDA eigenvalues, with
     !  occupancy factors included.
     !
     if ( pol(irp)%ntr == 0 ) cycle
     do ii = 1, pol(irp)%ntr
        repr(ii + ioffset) = irp
        i1 = pol(irp)%tr(1,ii)
        i2 = pol(irp)%tr(2,ii)
        k1 = pol(irp)%tr(3,ii)
        k2 = pol(irp)%tr(4,ii)
        isp = 1
        if (ii > pol(irp)%n_up) isp = 2
        occdif = kpt%wfn(isp,1)%occ0(i1) - kpt%wfn(isp,1)%occ0(i2)
        ene0(ii+ioffset) = kpt%wfn(isp,k2)%e1(i2) - kpt%wfn(isp,k1)%e1(i1)
        if (kpt%lcplx) then
           ostrz(ii+ioffset,:) = pol(irp)%zdipole(ii,:) * sqrt(occdif)
        else
           ostrz(ii+ioffset,:) = Zone * pol(irp)%ddipole(ii,:) * sqrt(occdif)
        endif
        do ip = 1, 3
           ostr0(ii+ioffset,ip) = two * ene0(ii+ioffset) * &
                abs(ostrz(ii+ioffset,ip))**2
        enddo
     enddo
     !
     ! LDA sum rule: calculate oscillator strengths from polarizability
     ! eigenvectors. Energy factors are polarizability eigenvalues.
     !
     do icol = 1, pol(irp)%ntr
        do ip = 1, 3
           ostr(icol+ioffset,ip) = pol(irp)%ostr(icol,ip)
        enddo
        ene(icol+ioffset) = pol(irp)%eig(icol)
     enddo
     ioffset = ioffset + pol(irp)%ntr
  enddo

  !-------------------------------------------------------------------
  ! Reorder transitions and calculate polarizabilities,
  ! oscillator strengths.
  !
  rsum = zero
  spol = zero
  do ii = 1, neig
     do ip = 1, 3
        rsum(ip) = rsum(ip) + ostr0(ii,ip)
        spol(ip,1) = spol(ip,1) + four*ostr0(ii,ip)/ene0(ii)/ene0(ii)
     enddo
  enddo
  allocate(indx(neig))
  call quicksort(neig,ene0,indx)

  open(68,file='eigenvalues_rpa',form='formatted')
  write(68,*) neig
  do ii = 1, neig
     write(68,'(f12.8,2x,g16.8,2x,g16.8,2x,g16.8,i4)') &
          ene0(indx(ii))*ryd, (ostr0(indx(ii),ip),ip=1,3), repr(indx(ii))
  enddo
  close(68)

  zvalence = zero
  do ii = 1, kpt%nk
     zvalence = zvalence + sum(kpt%wfn(1,ii)%occ0) + sum(kpt%wfn(nspin,ii)%occ0)
  enddo

  write(6,'(/,a)') repeat('-',65)
  write(6,'(/,a)') ' Checking oscillator strength sum rule'
  write(6,'(a)') ' Ratio between numerical and exact values :'
  write(6,'(25x,a)') '   Field polarization'
  write(6,'(25x,a)') '  x  ----  y  ----  z  ---- average'
  write(6,'(a,3f9.4,3x,f9.4)') '  RPA-LDA Sum Rule  =', &
       rsum / zvalence, sum(rsum) / zvalence / three

  rsum = zero
  do ii = 1, neig
     do ip = 1, 3
        rsum(ip) = rsum(ip) + ostr(ii,ip)
        spol(ip,2) = spol(ip,2) + four*ostr(ii,ip)/ene(ii)/ene(ii)
     enddo
  enddo
  spol = spol / real(kpt%nk,dp)
  if (.not. rpaonly) write(6,'(5x,a,3f9.4,f12.4,/)') 'ALDA Sum Rule  =', &
       rsum / zvalence, sum(rsum) / zvalence / three
  write(6,'(a,/)') repeat('-',65)
  xsum = sum(rsum) / zvalence / three
  call quicksort(neig,ene,indx)
  !-------------------------------------------------------------------
  ! Save RPA polarizability.
  !
  polrpa = sum(spol(:,1)) / three
  !-------------------------------------------------------------------
  ! Print out static polarizabilities.
  ! Units are (bohr)^3 = 0.1482 angstroms^3.
  ! 
  write(6,'(a)') ' Static polarizability (a.u.) '
  write(6,'(30x,a)') '   Field polarization '
  write(6,'(30x,a)') '  x  ----  y  ----  z  ---- average'
  write(6,'(a,3(1x,f8.2),1x,f12.3)') ' RPA-LDA Polarizability =', &
       spol(:,1), polrpa

  if (rpaonly) then
     write(6,'(a,/)') repeat('-',65)
  else
     open(68,file='eigenvalues_lda',form='formatted')
     write(68,*) neig
     do ii = 1, neig
        write(68,'(f12.8,2x,g16.8,2x,g16.8,2x,g16.8,i4)') &
             ene(indx(ii))*ryd, (ostr(indx(ii),ip),ip=1,3), repr(indx(ii))
     enddo
     close(68)

     write(6,'(a,3(1x,f8.2),1x,f12.3,/)') &
          '    ALDA Polarizability =', spol(:,2), sum(spol(:,2)) / three
     write(6,*) ' Average ALDA sum rule = ', xsum
     write(6,*) ' Average static RPA-LDA polarizability = ', polrpa

     write(6,'(a,/)') repeat('-',65)

     write(6,*) ' Lowest energy eigenvalues in polarizability', &
          ', all representations:'
     write(6,*) 'Order      Energy (eV)  Representation'
     write(6,'(i5,f16.6,5x,i5)') (ii, ene(indx(ii))*ryd, repr(indx(ii)), &
          ii = 1, min(10,neig))
  endif
  !-------------------------------------------------------------------
  ! In bulk systems, calculate and print out the RPA dielectric
  ! constant.
  !
  epsinv = one
  if (per == 3) then
     epsinv = one / ( one + &
          polrpa * four * pi / celvol * two / real(nspin,dp) )
     write(6,'(/,a)') repeat('-',65)
     write(6,*) 'Dielectric constant (RPA, no local fields) = ', one/epsinv
     write(6,'(a,/)') repeat('-',65)
  endif
  if (xsum > zero) xsum = sqrt(xsum)

  deallocate(indx)
  deallocate(ene0)
  deallocate(ostr0)
  deallocate(ostrz)
  deallocate(ene)
  deallocate(ostr)
  deallocate(repr)

end subroutine sum_rule
!===================================================================
