c
c     =================
      subroutine prsgrd
c     =================
c
c     pressure gradient terms
c
c     rhodp = int(rho), from
c     PCOM: z - pn at T grids, rho is reciprocal of density
c     BCOM: 0 - z  at T grids, rho is density
c
c     PCOM x-component
c     1) pax  = atmospheric pressure gradient
c     2) pbxn = int(rhodp_xbar+z*rho_xbar) at ivn(i,j)
c     3) pbxs = int(rhodp_xbar+z*rho_xbar) at ivn(i,j-1)
c     4) pcxn = int((rhodp)x)              at ivn(i,j)
c     5) pcxs = int((rhodp)x)              at ivn(i,j-1)
c
c     PCOM y-component
c     1) pay  = atmospheric pressure gradient
c     2) pbye = int(rhodp_ybar+z*rho_ybar) at ivn(i,j)
c     3) pbyw = int(rhodp_ybar+z*rho_ybar) at ivn(i-1,j)
c     4) pcye = int((rhodp)y)              at ivn(i,j)
c     5) pcyw = int((rhodp)y)              at ivn(i-1,j)
c
      implicit none
      include 'param.h'
      include 'pconst.h'
      include 'scalar.h'
      include 'grdvar.h'
      include 'prog.h'
      include 'cvbc.h'
c
      real wka,wkb,wkc,wk3
c
      common /works/ wka(imt,jmt),wkb(imt,jmt),wkc(imt,jmt)
      common /works/ wk3(imt,jmt,km)
      real           wkz(kmp1)
c
c
c-----------------------------------------------------------------------
c     atmospheric pressure gradients
c-----------------------------------------------------------------------
c     NOTE: set density=const for both PCOM and BCOM
c
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
c
c
c-----------------------------------------------------------------------
c     PCOM: rhodp(i,j,k) = int(1/density) from z to pn at T grids
c     BCOM: rhodp(i,j,k) = int(density)   from 0 to z  at T grids
c-----------------------------------------------------------------------
      do j=2,jmm
      do i=2,imm
        wkz(1)   = c0
        do k=1,itn(i,j)
        wkz(k+1) = wkz(k) + rho(i,j,k)*dz(k)
        enddo
        do k=1,itn(i,j)
#ifdef boussinesq
        rhodp(i,j,k) = (wkz(k+1)+wkz(k))*p5
#else
        rhodp(i,j,k) = wkz(itn(i,j)+1) - (wkz(k+1)+wkz(k))*p5
#endif
        enddo
      enddo
      enddo
c
      do k=1,km
      do j=2,jmm
      rhodp(1  ,j,k) = rhodp(imm,j,k)
      rhodp(imt,j,k) = rhodp(2  ,j,k)
      enddo
      enddo
c
c
c-----------------------------------------------------------------------
c     (rhodp)xbar & (rho)xbar*z at (i+1/2,j,k)
c-----------------------------------------------------------------------
      do k=1,km
      do j=2,jmm
      do i=2,imm
      wk3(i,j,k) = (rhodp(i+1,j,k)+rhodp(i,j,k))*p5
#ifdef boussinesq
     &           - (rho(i,j,k)+rho(i+1,j,k))*p5*z(k)
#else
     &           + (rho(i,j,k)+rho(i+1,j,k))*p5*z(k)
#endif
      enddo
      enddo
      enddo
      call vinteg_ns(wk3,pbxn,pbxs)
c
c-----------------------------------------------------------------------
c     (rhodp)x at (i+1/2,j,k)
c-----------------------------------------------------------------------
      do k=1,km
      do j=2,jmm
      do i=2,imm
      wk3(i,j,k) = (rhodp(i+1,j,k)-rhodp(i,j,k))*rdxt(j)
      enddo
      enddo
      enddo
      call vinteg_ns(wk3,pcxn,pcxs)
c
c-----------------------------------------------------------------------
c     (rhodp)ybar & (rho)ybar at (i,j+1/2,k)
c-----------------------------------------------------------------------
      do k=1,km
      do j=2,jmm
      do i=2,imm
      wk3(i,j,k) = (rhodp(i,j+1,k)+rhodp(i,j,k))*p5
#ifdef boussinesq
     &           - (rho(i,j,k)+rho(i,j+1,k))*p5*z(k)
#else
     &           + (rho(i,j,k)+rho(i,j+1,k))*p5*z(k)
#endif
      enddo
      enddo
      enddo
      call vinteg_ew(wk3,pbye,pbyw)
c
c-----------------------------------------------------------------------
c     (rhodp)y at (i,j+1/2,k)
c-----------------------------------------------------------------------
      do k=1,km
      do j=2,jmm
      do i=2,imm
      wk3(i,j,k) = (rhodp(i,j+1,k)-rhodp(i,j,k))*rdy
      enddo
      enddo
      enddo
      call vinteg_ew(wk3,pcye,pcyw)
c
c
#ifdef boussinesq
      do k=1,km
      do j=2,jmm
      do i=2,imm
      wk3(i,j,k) = (rho(i,j,k)+rho(i+1,j,k))*p5
      enddo
      enddo
      enddo
      call vinteg_ns(wk3,pdxn,pdxs)
c
      do k=1,km
      do j=2,jmm
      do i=2,imm
      wk3(i,j,k) = (rho(i,j,k)+rho(i,j+1,k))*p5
      enddo
      enddo
      enddo
      call vinteg_ew(wk3,pdye,pdyw)
#endif
c
      do j=2,jmm
      do i=2,imm
      pbxn(i,j) = pbxn(i,j)*rdxt(j)
      pbxs(i,j) = pbxs(i,j)*rdxt(j)
      pbye(i,j) = pbye(i,j)*rdy
      pbyw(i,j) = pbyw(i,j)*rdy
c
      pdxn(i,j) = pdxn(i,j)*rdxt(j)
      pdxs(i,j) = pdxs(i,j)*rdxt(j)
      pdye(i,j) = pdye(i,j)*rdy
      pdyw(i,j) = pdyw(i,j)*rdy
c
      pcxn(i,j) = pcxn(i,j)*p5
      pcxs(i,j) = pcxs(i,j)*p5
      pcye(i,j) = pcye(i,j)*p5
      pcyw(i,j) = pcyw(i,j)*p5
      enddo
      enddo
c
c
      do j=1,jmt
c
      pax(1  ,j) = pax(imm,j)
      pax(imt,j) = pax(2  ,j)
      pay(1  ,j) = pay(imm,j)
      pay(imt,j) = pay(2  ,j)
c
      pbxn(1  ,j) = pbxn(imm,j)
      pbxn(imt,j) = pbxn(2  ,j)
      pbxs(1  ,j) = pbxs(imm,j)
      pbxs(imt,j) = pbxs(2  ,j)
      pcxn(1  ,j) = pcxn(imm,j)
      pcxn(imt,j) = pcxn(2  ,j)
      pcxs(1  ,j) = pcxs(imm,j)
      pcxs(imt,j) = pcxs(2  ,j)
      pdxn(1  ,j) = pdxn(imm,j)
      pdxn(imt,j) = pdxn(2  ,j)
      pdxs(1  ,j) = pdxs(imm,j)
      pdxs(imt,j) = pdxs(2  ,j)
c
      pbye(1  ,j) = pbye(imm,j)
      pbye(imt,j) = pbye(2  ,j)
      pbyw(1  ,j) = pbyw(imm,j)
      pbyw(imt,j) = pbyw(2  ,j)
      pcye(1  ,j) = pcye(imm,j)
      pcye(imt,j) = pcye(2  ,j)
      pcyw(1  ,j) = pcyw(imm,j)
      pcyw(imt,j) = pcyw(2  ,j)
      pdye(1  ,j) = pdye(imm,j)
      pdye(imt,j) = pdye(2  ,j)
      pdyw(1  ,j) = pdyw(imm,j)
      pdyw(imt,j) = pdyw(2  ,j)
c
      enddo
c
      return
      end
