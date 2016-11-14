!
!     ===========================
!BOP
!
! !MODULE: upwelling.f90
! !DESCRIPTION: \input{sections/code-upwelling}
!
! !INTERFACE:
!
      subroutine upwelling(u,v,w,rdxt,rdyt,tmask,dz,pn,imt,jmt,km,imm,jmm,kmp1,  &
                           west,east,north,south,snbc,emp)
!EOP
!-------------------------------------------------------------------------------
!     ===========================
!
!     calculate vertical mass advection: pbt*dz/dt
!     upward is positive
!
      implicit none
      include 'pconst.h'
      include 'mpif.h'
!
      integer imt,jmt,km,imm,jmm,kmp1,i,j,k,snbc
      real tmask(imt,jmt,km)
      real rdxt(jmt),dz(km),pn(imt,jmt),rdyt(jmt),emp(imt,jmt)
      real    t1
      real    u(imt,jmt,km),v(imt,jmt,km),w(imt,jmt,kmp1)
      real    c(imt,jmt,km)
      real    a(imt,jmt),b(imt,jmt),dp(imt,jmt)
      
      integer west,east,north,south
!
!
      do j=2,jmm
      do i=2,imm
      dp(i,j) = c0
      enddo
      enddo
!
      do k=1,km
!
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
!
      do j=2,jmm
      do i=2,imm
      t1 = p5*(rdxt(j)*(a(i,j)-a(i-1,j))+rdyt(j)*(b(i,j)-b(i,j-1)))
      dp(i,j)  = dp(i,j) - t1*dz(k)
      c(i,j,k) = t1
      enddo
      enddo
!
      enddo
!
!
!     calculate pbt*dz/dt on the surface of T cell
!
      if (snbc==1) then
      do j=2,jmm
      do i=2,imm
      w(i,j,1) = emp(i,j) ! note that w is positive upward
      dp(i,j)  = dp(i,j) - emp(i,j)
      enddo
      enddo
      end if
!
      do k=2,km
      do j=2,jmm
      do i=2,imm
      w(i,j,k)=(w(i,j,k-1)+dz(k-1)*(dp(i,j)/pn(i,j)+c(i,j,k-1)))  &
               *tmask(i,j,k)
      enddo
      enddo
      enddo
!
!      do k=1,km
!      do j=2,jmm
!      w(1  ,j,k) = w(imm,j,k)
!      w(imt,j,k) = w(2  ,j,k)
!      enddo
!      enddo
      call swap_array_real3d(w,imt,jmt,kmp1,west,east,north,south)
!
      return
      end
