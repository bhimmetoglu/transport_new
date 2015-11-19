      program tauef
      !
      ! Author: Burak Himmetoglu
      ! This program computes the effective scattering rate at the Fermi level 
      ! The idea is to extract a single descriptor for electron LO-phonon 
      ! scattering strength, which is defined as the density of states 
      ! times inverse tau at E_F. 
      ! The code uses standard a2Fsave inputs.  

!$    USE omp_lib        

      USE smearing_mod

      implicit none

      integer :: i,j,k,nu,ik,ikk,nk1fit,nk2fit,nk3fit,nkfit,            &
     &           nbnd, nksfit, npk, nsym,                               &
     &           s(3,3,48),ns,nrot,ibnd,io,phband_i, phband_f,          &
     &           nphband, n, nn, jbnd, ibnd_ph, ind_k, cbm_i,nmod
      !
      double precision :: wk, at(3,3), bg(3,3), efermi, alat, shift,    &
     &                    T, wo(50), al(3), xk(3), invtau,aa,cut,deg
      ! 
      logical :: lsoc, lscissors
      !
      double precision, allocatable :: xkfit(:,:),etfit(:,:),wkfit(:),  &
     &                                 tauk_ef(:),dfk(:,:,:),vk(:,:,:)
      !
      integer, allocatable :: eqkfit(:), sfit(:), iflag(:,:), nkeff(:)
      !
      ! OMP variables
      double precision :: t0, ts, tf
      !
      integer :: nthreads
      !
      character*20 :: fil_a2F, fil_info
      !
      double precision, PARAMETER :: Rytocm1 = 109737.57990d0,          &
     &                               RytoGHz = 3.289828D6,              &
     &                               RytoTHz = RytoGHz/1000.d0,         &
     &                               RytoeV = 13.60569253d0,            &
     &                               tpi = 6.283185307d0,               &
     &                               convsig = 2.89d5,                  &
     &                               KtoRy = 1.d0/38.681648/300/RytoeV
                                     ! convsig: conversion factor for 
                                     ! conductivity. From au to (Ohm cm)-1 
      !
      namelist /input/ fil_info, fil_a2F, cut, T, efermi, alat,         &
     &                 aa, phband_i, phband_f, nthreads, lscissors,     &
     &                 cbm_i, shift

      read(5,input)

      !Set number of threads
      !$ call omp_set_num_threads(nthreads)
      call cpu_time(ts)

      npk = 40000
      !
      ! Convert to Ryd
      T = T * KtoRy
      efermi = efermi / RytoeV
      ! 
      ! Read the POP frequencies in cm-1 (just Gamma should be)
      !
      open(11,file='wo.in',status='unknown')
          read(11,*) (wo(i),i=1,3)
      close(11)
      !
      ! Convert to Ryd
      wo = wo / Rytocm1
      !
      ! LO-phonon couplings (Ryd*Ryd) 
      open(11,file='Cnu.txt',status='unknown')
      read(11,*) nmod
      do nu=1,nmod
         read(11,*) al(nu)
      end do
      close(11)
      !
      ! Total number of bands of interest (usually the number of relevant conduction bands)
      nphband = phband_f - phband_i + 1
      !
      ! Read a2Fsave
      open(11,file=fil_a2F,status='unknown')
      !
      read(11,*) nbnd, nksfit
      !
      allocate(etfit(nbnd,nksfit), xkfit(3,nksfit), wkfit(nksfit))
      !
      ! Read band structure and k-point data
      read(11,*) etfit
      read(11,*) ((xkfit(i,ik), i=1,3), ik=1,nksfit)
      read(11,*) wkfit
      read(11,*) nk1fit, nk2fit, nk3fit
      read(11,* ) nsym
      do ns=1,nsym
         read(11,*)  ((s(i,j,ns),j=1,3),i=1,3)
      enddo
      !
      close(11)
      !
      ! Regular grid size
      nkfit=nk1fit*nk2fit*nk3fit
      !
      ! Read info file on k-points (and lattice)
      open(11,file=fil_info,status='unknown')
      !
      read(11,*)
      read(11,*) ((at(i,j),i=1,3),j=1,3)
      !
      read(11,*)
      read(11,*)
      !
      read(11,*) ((bg(i,j),i=1,3),j=1,3)
      ! 
      close(11)
      !
      allocate (eqkfit(nkfit),sfit(nkfit))
      allocate (dfk(nphband,nkfit,3),vk(nphband,nkfit,3))
      allocate (nkeff(nphband), iflag(nphband,nkfit))
      allocate (tauk_ef(nphband))
      !
      ! Find the map between IBZ and the full-grid
      call lint ( nsym, s, .true., at, bg, npk, 0,0,0,                  &
     &  nk1fit,nk2fit,nk3fit, nksfit, xkfit, 1, nkfit, eqkfit, sfit)
      !
      ! Deallocate unnecessary variables
      deallocate(wkfit,sfit)
      !
      ! Weights for regular grid
      wk = 2.0/nkfit
      !
      ! If lscissors is true, then shift band energies (move conduction bands higher in energy)
      if (lscissors) then
        etfit(cbm_i:nbnd,:) = etfit(cbm_i:nbnd,:) + shift/RytoeV
        !cbm = cbm + shift
      end if
      !
      ! Call band velocities and forward derivatives
      call vband_ibz( nk1fit,nk2fit,nk3fit,nphband,nksfit,etfit(phband_i:phband_f,:),eqkfit,at, vk, dfk)
      !
      ! Include the 2pi/a factor
      vk = vk / tpi * alat
      !
      ! Determine the reduced grid for each band
      call reducegrid(nkfit,nksfit,nphband,etfit(phband_i:phband_f,:),eqkfit,efermi, &
     &                cut,T, nkeff,iflag)
      ! 
      tauk_ef = 0.d0
      ! Main do-loop, parallelized
      !$ t0 = omp_get_wtime()
      !
      !$omp parallel default(shared) &
      !$omp private(ik,ikk,ind_k,xk,invtau,deg,ibnd,ibnd_ph) 
      do ibnd=phband_i,phband_f 
         !
         ibnd_ph = ibnd - phband_i + 1
         ! 
         !$omp do reduction(+: tauk_ef)
         do ik=1,nkeff(ibnd_ph)
            !
            ikk = iflag(ibnd_ph,ik)  ! ikk is in full-grid (just reduced)
            ind_k = eqkfit(ikk)      ! ind_k is in IBZ
            xk(:) = xkfit(:,ind_k)   ! xk is in IBZ
            !
            !Compute inverse tau at (ibnd_ph,ind_k) 
            call invtau_nk ( nk1fit,nk2fit,nk3fit,nphband,nksfit,nmod,  &
     &           etfit(phband_i:phband_f,:),                            &
     &           xkfit,xk,vk,dfk,ind_k,ibnd_ph,eqkfit,al,wo,efermi,T,   &
     &           at, bg, aa, invtau )
            !
            ! Call smearing parameter at ibnd_ph and ikk 
            deg = sig0(nk1fit,nk2fit,nk3fit,dfk(ibnd_ph,ikk,:),aa)
            !
            !Compute inverse tau weighted by density of states 
            tauk_ef(ibnd_ph) = tauk_ef(ibnd_ph) + invtau * wk *         &
     &                    w0gauss((efermi-etfit(ibnd,ind_k))/deg)/deg
            !
         end do ! ik
         !$omp end do
         !
      end do ! ibnd
      !$omp end parallel
      !
      !Total integration time
      !$ t0 = omp_get_wtime() - t0
      call cpu_time(tf)
      !
      ! Write into file
      open(11,file='invtauef.out',status='unknown')
      write(11,"(4f14.6)") efermi * RytoeV, (tauk_ef(i) * RytoeV,i=1,nphband)
      close(11)
      !
      ! Free memory
      deallocate(xkfit,etfit,eqkfit,dfk,vk,nkeff,iflag,tauk_ef)

      end program tauef
