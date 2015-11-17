      subroutine invtau_nk (nk1,nk2,nk3,nbnd,nksfit,nmod,etfit,         &
                  xqfit,xk,vk,dfk,ind_k,ind_ph,eqk,al,wo,ef,T,at,bg,    &
                  aa, tauk )
      !
      ! Author: Burak Himmetoglu
      ! This subroutine computes the el-ph scattering rate for LO-phonons 
      ! The subroutine takes the following variables as input:
      ! nk1,nk2,nk3 : Full grid size
      ! nbnd,nksfit : # bands (of interest) and # qpoints in IBZ
      ! nmod        : # LO modes of interest
      ! etfit       : Band energies in IBZ 
      ! xqfit       : Grid in IBZ (in crystal coordinates)
      ! ind_k       : index (in the full grid) of the k-point
      ! ind_ph      : band index
      ! xk          : k-point coordinates where invtau is computed (in crystal coordinates) at ind_k
      ! dfk         : Forward derivative of energy at the k-grid
      ! vk          : Band velocity at the k-grid
      ! eqk         : Map of IBZ to full grid
      ! al          : el-ph coupling parameters
      ! wo          : LO phonon frequencies
      ! ef          : Fermi energy
      ! T           : Temperature
      ! at          : real space basis
      ! bg          : reciprocal space basis
      ! aa          : Parameter for auto-smearing calculation ( 0.8 < aa < 1.4)  
      !
      ! Output is tauk


!$    USE omp_lib

      USE smearing_mod
      
      implicit none

      ! INTENT(IN)
      integer, intent(in) :: nk1,nk2,nk3,nbnd,nksfit,nmod,eqk(nk1*nk2*nk3),&
     &                       ind_k,ind_ph

      double precision, intent(in) :: etfit(nbnd,nksfit), xqfit(3,nksfit),&
     &                                xk(3), al(nmod),wo(nmod), ef, T,  &
     &                                aa, at(3,3), bg(3,3),             &
     &                                vk(nbnd,nk1*nk2*nk3,3),           &
     &                                dfk(nbnd,nk1*nk2*nk3,3)

      !INTENT(OUT)
      double precision, intent(out) :: tauk

      ! LOCAL VARIABLES
      integer :: i,j,k,n,iq,jbnd,nu,nkfit,eqq(nk1*nk2*nk3)

      double precision :: wq, q2(nk1*nk2*nk3), dfkq(nbnd,nk1*nk2*nk3,3),&
     &                    vkq(nbnd,nk1*nk2*nk3,3),                      &
     &                    deg,temp1,temp2,w0g1,w0g2,be,ek,ekq,norm_vk,  &
     &                    norm_vkq,vfac,vkk(3),dfkk(3)

      double precision, parameter :: tpi = 6.283185307d0
      !
      !
      nkfit=nk1*nk2*nk3 
      !
      ! Set the q-pt weight
      wq=1.0/nkfit
      !
      ! Compute q**2 in the Gamma centered grid (Change to invq2)
      call cryst_to_cart(nksfit,xqfit,bg,1) ! to Cartesian coords
      q2 = 0.d0
      do iq=1,nkfit 
         do i=1,3
            q2(iq) = q2(iq) + xqfit(i,eqk(iq))**2.0
         end do
      end do
      !
      call cryst_to_cart(nksfit,xqfit,at,-1) ! back to crystal coords
      !
      ! Mapping between IBZ and full-grid for k+q (for a given xk)
      ! Results in eqq which is the map of the shifted grid
      call kplusq(nk1,nk2,nk3,nksfit,eqk,xk, eqq)
      !
      ! Compute forward derivatives and velocities for k+q grid: dfkq
      call vband_kq( nk1,nk2,nk3,nbnd,nksfit,xk,vk,dfk, vkq,dfkq )
      !
      ! Compute inverse tau
      tauk = 0.d0
      !
      ! band energy e_{i,k} and |vk|^2
      ek = etfit(ind_ph,ind_k) 
      vkk(:) = vk(ind_ph,ind_k,:)
      dfkk(:) = dfk(ind_ph,ind_k,:)
      norm_vk = abs(vkk(1)**2+vkk(2)**2+vkk(3)**2) 
      !
      do iq=2,nkfit ! Avoid q=0 divergence
         !
         do jbnd=1,nbnd
            !
            ! band energy e_{j,k+q}
            ekq = etfit(jbnd,eqq(iq)) 
            ! 
            ! Sum over LO modes
            do nu=1,nmod
               !
               ! BE distribution
               be = 1.d0 / ( exp(wo(nu)/T) - 1.d0 )   
               !
               ! n_{q\nu} + f_{j,k+q}
               temp1 = be + w1gauss((ef-ekq)/T,1) 
               !
               ! 1 - f_{j,k+q} + n_{q\nu}
               temp2 = 1.d0 - w1gauss((ef-ekq)/T,1) + be
               !
               ! Energy conserving delta-functions 
               ! First, compute the nk-dependent smearing parameter
               deg = sig_nk(nk1,nk2,nk3,dfkk,dfkq(jbnd,iq,:),aa)
               !deg = sig_nk(nk1,nk2,nk3,vkk,vkq(jbnd,iq,:),aa)
               !
               ! If deg=0, this means that ek-ekq=0, so the delta function must vanish... 
               if ( deg .le. 1d-10 ) then     
                    w0g1 = 0.d0
                    w0g2 = 0.d0
               else
                  w0g1 = w0gauss((ekq-ek-wo(nu))/deg)/deg
                  w0g2 = w0gauss((ekq-ek+wo(nu))/deg)/deg
               end if
               !
               ! Velocity factor
               vfac = vkk(1)*vkq(jbnd,iq,1)+vkk(2)*vkq(jbnd,iq,2)+vkk(3)*vkq(jbnd,iq,3)
               !
               norm_vkq = abs(vkq(jbnd,iq,1)**2+vkq(jbnd,iq,2)**2+vkq(jbnd,iq,3)**2) 
               if ( norm_vk .le. 1d-10 .or. norm_vkq .le. 1d-10) then
                  vfac = 1.d0
               else
                  vfac = 1.d0 - vfac / sqrt(norm_vk * norm_vkq)
               end if   
               !
               !write(6,"(2I4,7f14.6)") ind_ph, jbnd, vk(:), vkq(jbnd,iq,:), vfac 
               ! Compute tauk 
               !vfac = 1.0 ! DEBUG TEST
               tauk = tauk + tpi * wq * al(nu) / q2(iq) * vfac * ( temp1 * w0g1 + temp2 * w0g2 )
               !
            end do ! nu 
            ! 
         end do ! jbnd
         !
      end do ! iq

      end subroutine invtau_nk
