c
c     ===========================
      subroutine upwelling(u,v,w)
c     ===========================
c
c     calculate vertical mass advection: pbt*dz/dt
c     upward is positive
c
      implicit none
      include 'param.h'
      include 'pconst.h'
      include 'grdvar.h'
      include 'cvbc.h'
      include 'scalar.h'
c
      real    a,b,dp,t1
      real    u(imt,jmt,km),v(imt,jmt,km),w(imt,jmt,kmp1)
      real    c(imt,jmt,km)
c
      common /works/ a(imt,jmt),b(imt,jmt),dp(imt,jmt)
c
c
      do j=2,jmm
      do i=2,imm
      dp(i,j) = c0
      enddo
      enddo
c
      do k=1,km
c
      do j=2,jmm
      do i=1,imm
      a(i,j) = u(i,j,k)+u(i,j-1,k)
      enddo
      enddo
      do j=1,jmm
      do i=2,imm
      b(i,j) = v(i,j,k)+v(i-1,j,k)
      enddo
      enddo
c
      do j=2,jmm
      do i=2,imm
      t1 = p5*(rdxt(j)*(a(i,j)-a(i-1,j))+rdyt(j)*(b(i,j)-b(i,j-1)))
      dp(i,j)  = dp(i,j) - t1*dz(k)
      c(i,j,k) = t1
      enddo
      enddo
c
      enddo
c
c
c     calculate pbt*dz/dt on the surface of T cell
c
#ifdef snbc
       do j=2,jmm
       do i=2,imm
       w(i,j,1) = emp(i,j) ! note that w is positive upward
       dp(i,j)  = dp(i,j) - emp(i,j)
       enddo
       enddo
#endif
c
      do k=2,km
      do j=2,jmm
      do i=2,imm
      w(i,j,k)=(w(i,j,k-1)+dz(k-1)*(dp(i,j)/pn(i,j)+c(i,j,k-1)))
     &         *tmask(i,j,k)
      enddo
      enddo
      enddo
c
      do k=1,km
      do j=2,jmm
      w(1  ,j,k) = w(imm,j,k)
      w(imt,j,k) = w(2  ,j,k)
      enddo
      enddo
c
      return
      end
