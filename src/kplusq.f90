      subroutine kplusq (nk1,nk2,nk3,nksfit,eqk,xk, eqq)

      implicit none

      integer, intent(in) :: nk1,nk2,nk3,nksfit,eqk(nk1*nk2*nk3)

      double precision, intent(in) :: xk(3)

      integer, intent(out) :: eqq(nk1*nk2*nk3)

      ! Locals

      integer :: i,j,k,qx,qy,qz,iq,jq,kq,ik,n,nn

      double precision :: xkk(3)

      xkk = xk ! Local copy

      !call cryst_to_cart (1, xkk, at, -1) ! xk alredy in crystal coords, no need for transformation

      qx = nint(nk1*xkk(1))
      if (qx < 0) qx = qx + nk1
      if (qx > nk1) qx = qx - nk1
      !
      qy = nint(nk2*xkk(2))
      if (qy < 0) qy = qy + nk2
      if (qy > nk2) qy = qy - nk2
      !
      qz = nint(nk3*xkk(3))
      if (qz < 0) qz = qz + nk3
      if (qz > nk3) qz = qz - nk3
      !
      ! Find the map
      eqq = 0
      do i=1,nk1
         do j=1,nk2
            do k=1,nk3
               ik = k-1 + (j-1)*nk3 + (i-1)*nk2*nk3 + 1
               iq = i+qx
               if (iq > nk1) iq = iq - nk1
               jq = j+qy
               if (jq > nk2) jq = jq - nk2
               kq = k+qz
               if (kq > nk3) kq = kq - nk3
               nn = (kq-1)+(jq-1)*nk3+(iq-1)*nk2*nk3 + 1
               eqq(ik) = eqk(nn)
            enddo
         enddo
      enddo
      !
      !
      end subroutine kplusq
