      subroutine vband_kq ( nk1,nk2,nk3,nbnd,nksfit,xk,vk,dfk, vkq,dfkq )
      !
      ! Written by Burak Himmetoglu (burakhmmtgl@gmail.com)
      ! Uses some parts of the PW distribution
      !
      ! Computes band velocities and forward derivatives using finite
      ! differences in the k+q grid.
      !

!$    use omp_lib

      implicit none

      integer, intent(in) :: nk1,nk2,nk3, nbnd, nksfit

      double precision, intent(in) :: vk(nbnd,nk1*nk2*nk3,3),dfk(nbnd,nk1*nk2*nk3,3), &
     &                                xk(3)

      double precision, intent(out) :: vkq(nbnd,nk1*nk2*nk3,3), dfkq(nbnd,nk1*nk2*nk3,3)
                                       ! vkq is band velocity
                                       ! dfkq is forward derivative (used in adaptive smearing)

      integer :: i,j,k,kx,ky,kz,iq,n,ik,jk,kk,ibnd,nktot,nm,np

      ! Create the k-points, then compute velocities   

      nktot = nk1*nk2*nk3
    
      ! k vector to be added on q 
      kx = nint(nk1*xk(1))
      if (kx < 0) kx = kx + nk1
      if (kx > nk1) kx = kx - nk1
      !
      ky = nint(nk2*xk(2))
      if (ky < 0) ky = ky + nk2
      if (ky > nk2) ky = ky - nk2
      !
      kz = nint(nk3*xk(3))
      if (kz < 0) kz = kz + nk3
      if (kz > nk3) kz = kz - nk3
      !
      !
      do ibnd=1,nbnd
         do i=1,nk1
            do j=1,nk2
               do k=1,nk3
                  !
                  ! Index in the regular grid 
                  iq = k-1 + (j-1)*nk3 + (i-1)*nk2*nk3 + 1
                  !
                  ! Shift
                  ik = i+kx
                  if (ik > nk1) ik = ik - nk1
                  jk = j+ky
                  if (jk > nk2) jk = jk - nk2
                  kk = k+kz
                  if (kk > nk3) kk = kk - nk3
                  n = (kk-1)+(jk-1)*nk3+(ik-1)*nk2*nk3 + 1
                  !   
                  vkq(ibnd,iq,:) = vk(ibnd,n,:) 
                  dfkq(ibnd,iq,:) = dfk(ibnd,n,:)
                  !
               end do  ! k
            end do ! j
         end do ! i
      end do ! ibnd
      !
      !
      end subroutine vband_kq
