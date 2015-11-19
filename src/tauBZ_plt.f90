      program tauBZ_plt
      !
      ! Author: Burak Himmetoglu
      ! This code is for visualizing the inverse tau in the BZ
      ! The code uses standard a2Fsave inputs.  

!$    USE omp_lib        

      USE smearing_mod

      implicit none

      integer :: i,j,k,nu,ik,ikk,nk1fit,nk2fit,nk3fit,nkfit,            &
     &           nbnd, nksfit, npk, nsym,                               &
     &           s(3,3,48),ns,nrot,ibnd,io,phband_i, phband_f,          &
     &           nphband, n, nn, jbnd, ibnd_ph, ind_k, counter,nmod
      !
      double precision :: wk, at(3,3), bg(3,3), efermi, alat,           &
     &                    T, wo(3), al(3), xk(3),aa,cut,deg
      ! 
      logical :: lsoc
      !
      double precision, allocatable :: xkfit(:,:),etfit(:,:),wkfit(:),  &
     &                                 dfk(:,:,:),vk(:,:,:),            &
     &                                 xkg(:,:), invtau(:)
      !
      integer, allocatable :: eqkfit(:), sfit(:)
      !
      ! OMP variables
      double precision :: t0
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
     &                 aa, phband_i, phband_f

      read(5,input)

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
      allocate (invtau(nphband))
      !
      ! Find the map between IBZ and the full-grid
      call lint ( nsym, s, .true., at, bg, npk, 0,0,0,                  &
     &  nk1fit,nk2fit,nk3fit, nksfit, xkfit, 1, nkfit, eqkfit, sfit)
      !
      ! Deallocate unnecessary variables, allocate new ones
      deallocate(wkfit,sfit)
      allocate(xkg(nkfit,3))
      !
      ! Weights for regular grid
      wk = 2.0/nkfit
      !
      ! Call band velocities and forward derivatives
      call vband_ibz( nk1fit,nk2fit,nk3fit,nphband,nksfit,etfit(phband_i:phband_f,:),eqkfit,at, vk, dfk)
      !
      ! Include the 2pi/a factor
      vk = vk / tpi * alat
      !
      ! Generate the uniform k-grid
      do i=1,nk1fit
         do j=1,nk2fit
            do k=1,nk3fit
               n = (k-1) + (j-1)*nk3fit + (i-1)*nk2fit*nk3fit + 1
               ! 
               ! crystal coordinates
               ! 
               xkg(n,1) = dble(i-1)/nk1fit
               xkg(n,2) = dble(j-1)/nk2fit
               xkg(n,3) = dble(k-1)/nk3fit
            end do
         end do
      end do
      ! 
      call cryst_to_cart(nkfit,xkg,bg,1) ! pass to cart coords
      ! 
      ! Serial, due to IO operations
      open(11,file='invtau_BZ.out',status='unknown')
      !
      counter = 0
      do ik=1,nkfit
         !
         do ibnd=phband_i,phband_f
            !
            ibnd_ph = ibnd - phband_i + 1
            !
            ind_k = eqkfit(ik)       ! ind_k is in IBZ
            xk(:) = xkfit(:,ind_k)   ! xk is in IBZ
            !
            !Compute inverse tau at (ibnd_ph,ind_k) 
            call invtau_nk ( nk1fit,nk2fit,nk3fit,nphband,nksfit,nmod,  &
     &           etfit(phband_i:phband_f,:),                            &
     &           xkfit,xk,vk,dfk,ind_k,ibnd_ph,eqkfit,al,wo,efermi,T,   &
     &           at, bg, aa, invtau(ibnd_ph) )
            !
         end do ! ibnd 
         !
         ! Z = 0 plane of BZ
         if ( xkg(ik,3) .lt. 1e-8 ) then
            !
            write(11,'(5f14.6)') xkg(ik,1), xkg(ik,2), (invtau(i)*RytoeV,i=1,3)
            counter = counter + 1
            if ( counter .eq. nk1fit ) then
               write(11,*)
               counter = 0
            end if
            !
         end if
         !
      end do ! ik 
      !  
      close(11)
      !
      ! Free memory
      deallocate(xkfit,etfit,eqkfit,dfk,vk,invtau,xkg)

      end program tauBZ_plt
