c
c     =================
      subroutine readyc
c     =================
c     momentum advections, viscosities & atmpospheric pressure terms
c
      implicit none
      include 'param.h'
      include 'pconst.h'
      include 'scalar.h'
      include 'grdvar.h'
      include 'prog.h'
      include 'cvbc.h'
c
      real ubar,vbar,abc,uvmag,t1,t2,t3
      real a,b,c,u,v,ux,vx,wua,wub,wva,wvb
c
      common /works/ a(imt,jmt),b(imt,jmt),c(imt,jmt)
      common /works/ u(imt,jmt),v(imt,jmt)
      common /works/ ux(imt,jmt),vx(imt,jmt)
      common /works/ wua(imt,jmt),wub(imt,jmt)
      common /works/ wva(imt,jmt),wvb(imt,jmt)
c
c
c-----------------------------------------------------------------------
c     initialize time-averaged pbt (used in baroclinic Eq)
c-----------------------------------------------------------------------
c
      do j=1,jmt
      do i=1,imt
      pmum(i,j) = pmup(i,j)
      pmup(i,j) = pbt (i,j,tau)
      enddo
      enddo
c
c
c-----------------------------------------------------------------------
c     calculate vertical mass advection (dz/dt * pbt)
c-----------------------------------------------------------------------
c
      do j=2,jmm
      do i=2,imm
      if(ivn(i,j).gt.0)then
      spbt(i,j) = p5*sqrt(pbt(i,j  ,tau) + pbt(i+1,j  ,tau) +
     &                    pbt(i,j+1,tau) + pbt(i+1,j+1,tau))
      endif
      enddo
      spbt(1  ,j) = spbt(imm,j)
      spbt(imt,j) = spbt(2  ,j)
      enddo
c
c     du,dv are used temporarily for mass advection
c
      do k=1,km
      do j=1,jmt
      do i=1,imt
      if(umask(i,j,k).gt.c0)then
      du(i,j,k) = spbt(i,j)*up(i,j,k,tau)
      dv(i,j,k) = spbt(i,j)*vp(i,j,k,tau)*cosu(j)
      else
      du(i,j,k) = c0
      dv(i,j,k) = c0
      endif
      enddo
      enddo
      enddo
c
      call upwelling(du,dv,w)
c
c     calculate vertical velocity on the surface of U cell
c
      do k=1,km
      do j=2,jmt
      do i=2,imt
      a(i,j)   = w(i,j,k)/pbt(i,j,tau)
      enddo
      enddo
      do j=2,jmm
      do i=2,imm
      w(i,j,k) = p25*(a(i,j)+a(i+1,j)+a(i,j+1)+a(i+1,j+1))
      enddo
      enddo
      enddo
c
c
c-----------------------------------------------------------------------
c     calculate advections & horizontal viscosities
c-----------------------------------------------------------------------
c
      do 100 k=1,km
c
c     calculate horizontal advection velocities
c
      do j=1,jmt
      do i=1,imt
      if(umask(i,j,k).gt.c0)then
      u(i,j) = up(i,j,k,tau)/spbt(i,j)
      v(i,j) = vp(i,j,k,tau)/spbt(i,j)
      else
      u(i,j) = c0
      v(i,j) = c0
      endif
      enddo
      enddo
c
c     U-advection
c
      do j=2,jmm
      do i=2,imt
      ubar   = p5*(u(i,j)+u(i-1,j))
      a(i,j) = p5*(up(i,j,k,tau)+up(i-1,j,k,tau))*ubar
      b(i,j) = p5*(vp(i,j,k,tau)+vp(i-1,j,k,tau))*ubar
      c(i,j) = ubar
      enddo
      enddo
c
      do j=2,jmm
      do i=2,imm
      if(umask(i,j,k).gt.c0)then
      abc     = c(i+1,j)-c(i,j)
      ux(i,j) = rdxu(j)*((a(i+1,j)-a(i,j)) - p5*up(i,j,k,tau)*abc)
      vx(i,j) = rdxu(j)*((b(i+1,j)-b(i,j)) - p5*vp(i,j,k,tau)*abc)
      endif
      enddo
      enddo
c
c
c     V-advection
c
      do j=2,jmt
      do i=2,imm
      vbar   = p5*(v(i,j)*cosu(j)+v(i,j-1)*cosu(j-1))
      a(i,j) = p5*(up(i,j,k,tau)+up(i,j-1,k,tau))*vbar
      b(i,j) = p5*(vp(i,j,k,tau)+vp(i,j-1,k,tau))*vbar
      c(i,j) = vbar
      enddo
      enddo
c
      do j=2,jmm
      do i=2,imm
      if(umask(i,j,k).gt.c0)then
      abc     = c(i,j+1)-c(i,j)
      ux(i,j) = ux(i,j)+rdyu(j)*((a(i,j+1)-a(i,j))-p5*up(i,j,k,tau)*abc)
      vx(i,j) = vx(i,j)+rdyu(j)*((b(i,j+1)-b(i,j))-p5*vp(i,j,k,tau)*abc)
      endif
      enddo
      enddo
c
c
c     W-advection
c
      if(k.eq.1)then
        do j=2,jmm
        do i=2,imm
        wua(i,j) = c0
        wva(i,j) = c0
        enddo
        enddo
      else
        do j=2,jmm
        do i=2,imm
        wua(i,j) = wub(i,j)
        wva(i,j) = wvb(i,j)
        enddo
        enddo
      endif
c
      do j=2,jmm
      do i=2,imm
      if(k.ge.ivn(i,j))then
        wub(i,j) = c0
        wvb(i,j) = c0
      else 
        wub(i,j) = w(i,j,k+1)*(up(i,j,k,tau)+up(i,j,k+1,tau))
        wvb(i,j) = w(i,j,k+1)*(vp(i,j,k,tau)+vp(i,j,k+1,tau))
      endif
      enddo
      enddo
c
      do j=2,jmm
      do i=2,imm
      if(umask(i,j,k).gt.c0)then
      abc     = w(i,j,k) - w(i,j,k+1)
      ux(i,j) = ux(i,j)+p5*rdz(k)*(wua(i,j)-wub(i,j)-up(i,j,k,tau)*abc)
      vx(i,j) = vx(i,j)+p5*rdz(k)*(wva(i,j)-wvb(i,j)-vp(i,j,k,tau)*abc)
      endif
      enddo
      enddo
c
c
c-----------------------------------------------------------------------
c     + advection + Pa + viscosities(t=taum)
c-----------------------------------------------------------------------
c
      if(leapfrog_c)then
        do j=2,jmm
        do i=2,imm
        du(i,j,k) = - ux(i,j) - pax(i,j)*spbt(i,j) + diffu(i,j,k)
        dv(i,j,k) = - vx(i,j) - pay(i,j)*spbt(i,j) + diffv(i,j,k)
        enddo
        enddo
      else
        do j=2,jmm
        do i=2,imm
        du(i,j,k) = - ux(i,j) - pax(i,j)*spbt(i,j)
        dv(i,j,k) = - vx(i,j) - pay(i,j)*spbt(i,j)
        enddo
        enddo
      endif
c
c
c---------------------------------------------------------------------
c     calculate horizontal viscosities in tau time level
c---------------------------------------------------------------------
c
      do j=2,jmm
      do i=2,imt
      abc     = pbt(i,j,tau)+pbt(i,j+1,tau)
      ux(i,j) = abc*(u(i,j)-u(i-1,j))
      vx(i,j) = abc*(v(i,j)-v(i-1,j))
      enddo
      enddo
      do j=2,jmt
      do i=2,imm
      abc     = pbt(i,j,tau)+pbt(i+1,j,tau)
      a(i,j)  = abc*(u(i,j)-u(i,j-1))
      b(i,j)  = abc*(v(i,j)-v(i,j-1))
      enddo
      enddo
c
      do j=2,jmm
      do i=2,imm
      if(umask(i,j,k).gt.c0)then
      diffu(i,j,k) = am*( (sdxu(j)*(ux(i+1,j)-ux(i,j))+r1c(j)*a(i,j+1)
     &                     -r1d(j)*a(i,j))/spbt(i,j)
     &             + cv1(j)*up(i,j,k,tau)
     &             - cv2(j)*(v(i+1,j)-v(i-1,j))*spbt(i,j) )
c
      diffv(i,j,k) = am*( (sdxu(j)*(vx(i+1,j)-vx(i,j))+r1c(j)*b(i,j+1)
     &                     -r1d(j)*b(i,j))/spbt(i,j)
     &             + cv1(j)*vp(i,j,k,tau)
     &             + cv2(j)*(u(i+1,j)-u(i-1,j))*spbt(i,j) )
      endif
      enddo
      enddo
c
100   continue
c
c
c---------------------------------------------------------------------
c     NOTE: u(i,j) & v(i,j) = u&v at bottom (k=ivn(i,j))
c---------------------------------------------------------------------
c
c
c---------------------------------------------------------------------
c     calculate vertical viscosities in tau time level
c---------------------------------------------------------------------
c
c     2-D working array
      do j=2,jmm
      do i=2,imm
      if(ivn(i,j).gt.0) a(i,j)=kappa_m*gravr/spbt(i,j)**2
      enddo
      enddo
c
c
      do 200 k=1,km
c
      if(k.eq.1)then  
c       windstress bcu&bcv
        do j=2,jmm
        do i=2,imm
        wua(i,j) = bcu(i,j)*gravr/spbt(i,j)*rrho_0
        wva(i,j) = bcv(i,j)*gravr/spbt(i,j)*rrho_0
        enddo
        enddo
      else
        do j=2,jmm
        do i=2,imm
        wua(i,j) = wub(i,j)
        wva(i,j) = wvb(i,j)
        enddo
        enddo
      endif
c
      do j=2,jmm
      do i=2,imm
      if(k.eq.ivn(i,j))then
c       bottom drag (assume rho=c1 for both PCOM and BCOM)
        uvmag    = sqrt(u(i,j)**2+v(i,j)**2)
        abc      = grav/spbt(i,j)*cdbot*uvmag
        wub(i,j) = u(i,j)*abc
        wvb(i,j) = v(i,j)*abc
      else if(k.lt.ivn(i,j))then
        wub(i,j) = a(i,j)*(up(i,j,k,tau)-up(i,j,k+1,tau))*rdzw(k)
        wvb(i,j) = a(i,j)*(vp(i,j,k,tau)-vp(i,j,k+1,tau))*rdzw(k)
      else
        wub(i,j) = c0
        wvb(i,j) = c0
      endif
      enddo
      enddo
c
c
      do j=2,jmm
      do i=2,imm
      diffu(i,j,k) = diffu(i,j,k) + rdz(k)*(wua(i,j)-wub(i,j))
      diffv(i,j,k) = diffv(i,j,k) + rdz(k)*(wva(i,j)-wvb(i,j))
      enddo
      enddo
c
c
c-----------------------------------------------------------------------
c     + diffusion (at tau time level) on euler forward time step
c-----------------------------------------------------------------------
c
      if(.not.leapfrog_c)then
      do j=2,jmm
      do i=2,imm
      du(i,j,k) = du(i,j,k) + diffu(i,j,k)
      dv(i,j,k) = dv(i,j,k) + diffv(i,j,k)
      enddo
      enddo
      endif
c
200   continue
c
c
c-----------------------------------------------------------------------
c     vertical integration of du&dv
c-----------------------------------------------------------------------
c
      call vinteg(du,dub)
      call vinteg(dv,dvb)
c
c
      return
      end
