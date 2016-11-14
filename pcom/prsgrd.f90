!
!     =================
!BOP
!
! !MODULE: prsgrd.f90
! !DESCRIPTION: \input{sections/code-prsgrd}
!
! !INTERFACE:
!
      subroutine prsgrd(bcp,rdxt,rdy,pax,pay,itn,ivn,rho,rhodp,z,dz,rzu,pbxn,pbxs,  &
                        pbye,pbyw,pcxn,pcxs,pcye,pcyw,pdxn,pdxs,pdye,pdyw,imt,jmt,  &
                        km,imm,jmm,kmp1,west,east,north,south)
!EOP
!-------------------------------------------------------------------------------
!     =================
!
!     pressure gradient terms
!
!     rhodp = int(rho), from
!     PCOM: z - pn at T grids, rho is reciprocal of density
!     BCOM: 0 - z  at T grids, rho is density
!
!     PCOM x-component
!     1) pax  = atmospheric pressure gradient
!     2) pbxn = int(rhodp_xbar+z*rho_xbar) at ivn(i,j)
!     3) pbxs = int(rhodp_xbar+z*rho_xbar) at ivn(i,j-1)
!     4) pcxn = int((rhodp)x)              at ivn(i,j)
!     5) pcxs = int((rhodp)x)              at ivn(i,j-1)
!
!     PCOM y-component
!     1) pay  = atmospheric pressure gradient
!     2) pbye = int(rhodp_ybar+z*rho_ybar) at ivn(i,j)
!     3) pbyw = int(rhodp_ybar+z*rho_ybar) at ivn(i-1,j)
!     4) pcye = int((rhodp)y)              at ivn(i,j)
!     5) pcyw = int((rhodp)y)              at ivn(i-1,j)
!
      implicit none
      include 'pconst.h'
      include 'mpif.h'
!
      integer imt,jmt,km,imm,jmm,kmp1,i,j,k
      real wka(imt,jmt),wkb(imt,jmt),wkc(imt,jmt),wk3(imt,jmt,km)
!
      real wkz(kmp1)
      real rdy,rdxt(jmt)
      real bcp(imt,jmt)
      real pax(imt,jmt),pay(imt,jmt)
      real pbxn(imt,jmt),pbxs(imt,jmt)
      real pcxn(imt,jmt),pcxs(imt,jmt)
      real pdxn(imt,jmt),pdxs(imt,jmt)
      real pbye(imt,jmt),pbyw(imt,jmt)
      real pcye(imt,jmt),pcyw(imt,jmt)
      real pdye(imt,jmt),pdyw(imt,jmt)
      real rhodp(imt,jmt,km),rho(imt,jmt,km)
      integer itn(imt,jmt),ivn(imt,jmt)
      real dz(km),rzu(imt,jmt),z(km)
      
      integer west,east,north,south
!
!
!-----------------------------------------------------------------------
!     atmospheric pressure gradients
!-----------------------------------------------------------------------
!     NOTE: set density=const for both PCOM and BCOM
!
      do j=2,jmt
      do i=2,imm
      wka(i,j) = rrho_0*(bcp(i+1,j)-bcp(i,j))*rdxt(j)
      enddo
      enddo
      do j=2,jmm
      do i=2,imt
      wkb(i,j) = rrho_0*(bcp(i,j+1)-bcp(i,j))*rdy
      enddo
      enddo 
      do j=2,jmm
      do i=2,imm
      pax(i,j) = (wka(i,j)+wka(i,j+1))*p5
      pay(i,j) = (wkb(i,j)+wkb(i+1,j))*p5
      enddo
      enddo
!
!
!-----------------------------------------------------------------------
!     PCOM: rhodp(i,j,k) = int(1/density) from z to pn at T grids
!     BCOM: rhodp(i,j,k) = int(density)   from 0 to z  at T grids
!-----------------------------------------------------------------------
      do j=2,jmm
      do i=2,imm
        wkz(1)   = c0
        do k=1,itn(i,j)
        wkz(k+1) = wkz(k) + rho(i,j,k)*dz(k)
        enddo
        do k=1,itn(i,j)
        rhodp(i,j,k) = wkz(itn(i,j)+1) - (wkz(k+1)+wkz(k))*p5
        enddo
      enddo
      enddo
!
!      do k=1,km
!      do j=2,jmm
!      rhodp(1  ,j,k) = rhodp(imm,j,k)
!      rhodp(imt,j,k) = rhodp(2  ,j,k)
!      enddo
!      enddo
      call swap_array_real3d(rhodp,imt,jmt,km,west,east,north,south)
!
!
!-----------------------------------------------------------------------
!     (rhodp)xbar & (rho)xbar*z at (i+1/2,j,k)
!-----------------------------------------------------------------------
      do k=1,km
      do j=2,jmm
      do i=2,imm
      wk3(i,j,k) = (rhodp(i+1,j,k)+rhodp(i,j,k))*p5 &
                 + (rho(i,j,k)+rho(i+1,j,k))*p5*z(k)
      enddo
      enddo
      enddo
      call vinteg_ns(wk3,pbxn,pbxs,ivn,dz,rzu,imt,jmt,km,imm,jmm)
!
!-----------------------------------------------------------------------
!     (rhodp)x at (i+1/2,j,k)
!-----------------------------------------------------------------------
      do k=1,km
      do j=2,jmm
      do i=2,imm
      wk3(i,j,k) = (rhodp(i+1,j,k)-rhodp(i,j,k))*rdxt(j)
      enddo
      enddo
      enddo
      call vinteg_ns(wk3,pcxn,pcxs,ivn,dz,rzu,imt,jmt,km,imm,jmm)
!
!-----------------------------------------------------------------------
!     (rhodp)ybar & (rho)ybar at (i,j+1/2,k)
!-----------------------------------------------------------------------
      do k=1,km
      do j=2,jmm
      do i=2,imm
      wk3(i,j,k) = (rhodp(i,j+1,k)+rhodp(i,j,k))*p5  &
                 + (rho(i,j,k)+rho(i,j+1,k))*p5*z(k)
      enddo
      enddo
      enddo
      call vinteg_ew(wk3,pbye,pbyw,ivn,dz,rzu,imt,jmt,km,imm,jmm)
!
!-----------------------------------------------------------------------
!     (rhodp)y at (i,j+1/2,k)
!-----------------------------------------------------------------------
      do k=1,km
      do j=2,jmm
      do i=2,imm
      wk3(i,j,k) = (rhodp(i,j+1,k)-rhodp(i,j,k))*rdy
      enddo
      enddo
      enddo
      call vinteg_ew(wk3,pcye,pcyw,ivn,dz,rzu,imt,jmt,km,imm,jmm)
!
!
!
      
      do j=2,jmm
      do i=2,imm
      pbxn(i,j) = pbxn(i,j)*rdxt(j)
      pbxs(i,j) = pbxs(i,j)*rdxt(j)
      pbye(i,j) = pbye(i,j)*rdy
      pbyw(i,j) = pbyw(i,j)*rdy
!
      pdxn(i,j) = pdxn(i,j)*rdxt(j)
      pdxs(i,j) = pdxs(i,j)*rdxt(j)
      pdye(i,j) = pdye(i,j)*rdy
      pdyw(i,j) = pdyw(i,j)*rdy
!
      pcxn(i,j) = pcxn(i,j)*p5
      pcxs(i,j) = pcxs(i,j)*p5
      pcye(i,j) = pcye(i,j)*p5
      pcyw(i,j) = pcyw(i,j)*p5
      enddo
      enddo
!
!
!      do j=1,jmt
!
!      pax(1  ,j) = pax(imm,j)
!      pax(imt,j) = pax(2  ,j)
!      pay(1  ,j) = pay(imm,j)
!      pay(imt,j) = pay(2  ,j)
!
!      pbxn(1  ,j) = pbxn(imm,j)
!      pbxn(imt,j) = pbxn(2  ,j)
!      pbxs(1  ,j) = pbxs(imm,j)
!      pbxs(imt,j) = pbxs(2  ,j)
!      pcxn(1  ,j) = pcxn(imm,j)
!      pcxn(imt,j) = pcxn(2  ,j)
!      pcxs(1  ,j) = pcxs(imm,j)
!      pcxs(imt,j) = pcxs(2  ,j)
!      pdxn(1  ,j) = pdxn(imm,j)
!      pdxn(imt,j) = pdxn(2  ,j)
!      pdxs(1  ,j) = pdxs(imm,j)
!      pdxs(imt,j) = pdxs(2  ,j)
!
!      pbye(1  ,j) = pbye(imm,j)
!      pbye(imt,j) = pbye(2  ,j)
!      pbyw(1  ,j) = pbyw(imm,j)
!      pbyw(imt,j) = pbyw(2  ,j)
!      pcye(1  ,j) = pcye(imm,j)
!      pcye(imt,j) = pcye(2  ,j)
!      pcyw(1  ,j) = pcyw(imm,j)
!      pcyw(imt,j) = pcyw(2  ,j)
!      pdye(1  ,j) = pdye(imm,j)
!      pdye(imt,j) = pdye(2  ,j)
!      pdyw(1  ,j) = pdyw(imm,j)
!      pdyw(imt,j) = pdyw(2  ,j)
!
!      enddo
       call swap_array_real2d(pax,imt,jmt,west,east,north,south)
       call swap_array_real2d(pay,imt,jmt,west,east,north,south)
       call swap_array_real2d(pbxn,imt,jmt,west,east,north,south)
       call swap_array_real2d(pbxs,imt,jmt,west,east,north,south)
       call swap_array_real2d(pcxn,imt,jmt,west,east,north,south)
       call swap_array_real2d(pcxs,imt,jmt,west,east,north,south)
       call swap_array_real2d(pdxn,imt,jmt,west,east,north,south)
       call swap_array_real2d(pdxs,imt,jmt,west,east,north,south)
       call swap_array_real2d(pbye,imt,jmt,west,east,north,south)
       call swap_array_real2d(pbyw,imt,jmt,west,east,north,south)
       call swap_array_real2d(pcye,imt,jmt,west,east,north,south)
       call swap_array_real2d(pcyw,imt,jmt,west,east,north,south)
       call swap_array_real2d(pdye,imt,jmt,west,east,north,south)
       call swap_array_real2d(pdyw,imt,jmt,west,east,north,south)
!
      return
      end
