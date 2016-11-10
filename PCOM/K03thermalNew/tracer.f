c
c     =================
      subroutine tracer
c     =================
c     compute potential temperature and salinity at tau+1 time level
c
      implicit none
      include 'param.h'
      include 'pconst.h'
      include 'scalar.h'
      include 'grdvar.h'
      include 'prog.h'
      include 'cvbc.h'
      real a,b,c,wa,wb,flxa,flxb,dt,abc,aidif,dtts2
c
#ifdef gm90
      include 'isopyc.h'
#endif
c
      real wka(imt,jmt,km),wkb(imt,jmt,km),pbar(imt,jmt),stf(imt,jmt)
c
      common /works/ a(imt,jmt),b(imt,jmt),c(imt,jmt)
      common /works/ wa(imt,jmt),wb(imt,jmt)
      common /works/ flxa(imt,jmt),flxb(imt,jmt),dt(imt,jmt)
c
      real fx,fy
      common /wisop/ fx(imt,jmt),fy(imt,jmt)
c
c
#if defined gm90 || defined implicitvmix
c 2006-1-10 #ifdef gm90
      aidif = p5
#else
      aidif = c1
#endif
c
c-----------------------------------------------------------------------
c     calculate time-averaged pbt and mass advections
c-----------------------------------------------------------------------
      do j=1,jmt
      do i=1,imt
      pmtp(i,j) = pmtp(i,j)*onbc
      enddo
      enddo
      do k=1,km
      do j=1,jmt
      do i=1,imt
      ump(i,j,k) = ump(i,j,k)*oncc
      vmp(i,j,k) = vmp(i,j,k)*oncc
      enddo
      enddo
      enddo
c
      if(leapfrog_t)then
       do j=1,jmt
       do i=1,imt
        pbar(i,j) = (pmtp(i,j) + pmtm(i,j))*p5
       enddo
       enddo
       do k=1,km
       do j=1,jmt
       do i=1,imt
        wka(i,j,k) = (ump(i,j,k)+umm(i,j,k))*p5
        wkb(i,j,k) = (vmp(i,j,k)+vmm(i,j,k))*p5
       enddo
       enddo
       enddo
      else
       do j=1,jmt
       do i=1,imt
        pbar(i,j) = pmtp(i,j)
       enddo
       enddo
       do k=1,km
       do j=1,jmt
       do i=1,imt
        wka(i,j,k) = ump(i,j,k)
        wkb(i,j,k) = vmp(i,j,k)
       enddo
       enddo
       enddo
      endif
c
c
c---------------------------------------------------------------------
c     calculate vertical mass advection
c---------------------------------------------------------------------
      call upwelling(wka,wkb,w)
c
c
c---------------------------------------------------------------------
c     calculate 1/2 mass advection across W & S face of T cells
c---------------------------------------------------------------------
c     du & dv are used temporary as working arrays
      do k=1,km
      do j=2,jmt
      do i=2,imt
      du(i,j,k) = p25*(wka(i-1,j,k)+wka(i-1,j-1,k))
      dv(i,j,k) = p25*(wkb(i,j-1,k)+wkb(i-1,j-1,k))
      enddo
      enddo
      enddo
c
c
c---------------------------------------------------------------------
c     calculate  2*xbar(pbar) & 2*ybar(pbar)
c---------------------------------------------------------------------
      do j=1,jmt
      do i=2,imt
      fx(i,j) = (pbar(i,j)+pbar(i-1,j))*p5
      enddo
      enddo
      do j=2,jmt
      do i=1,imt
      fy(i,j) = (pbar(i,j)+pbar(i,j-1))*p5
      enddo
      enddo
c
c
#ifdef gm90
c---------------------------------------------------------------------
c     diffusion coefficients
c---------------------------------------------------------------------
      call isopyc
c
      do k=1,km
      do j=1,jmt
      do i=1,imt
      wkb(i,j,k) = kappa_h(i,j,k) + ahisop*k3(i,k,j,3)
      enddo
      enddo
      enddo
#endif
c
      if(leapfrog_t) then
        dtts2 = c2dtts
      else
        dtts2 = dtts
      endif
c
c---------------------------------------------------------------------
c     solve for one tracer at a time
c---------------------------------------------------------------------
c     n = 1 => temperature
c     n = 2 => salinity
c
c
      do 1000 n=1,nt
c
c
c---------------------------------------------------------------------
c     advections
c---------------------------------------------------------------------
      do 100 k=1,km
c
      do j=2,jmm
      do i=2,imt
      a(i,j) = du(i,j,k)*(t(i,j,k,n,tau)-t(i-1,j,k,n,tau))
      enddo
      enddo
c
      do j=2,jmt
      do i=2,imm
      b(i,j) = dv(i,j,k)*(t(i,j,k,n,tau)-t(i,j-1,k,n,tau))
      enddo
      enddo
c
      if(k.eq.1)then
        do j=2,jmm
        do i=2,imm
        wa(i,j) = c0
        enddo
        enddo
      else
        do j=2,jmm
        do i=2,imm
        wa(i,j) = wb(i,j)
        enddo
        enddo
      endif
c
      do j=2,jmm
      do i=2,imm
      if(k.ge.itn(i,j))then
        wb(i,j) = c0
      else 
        wb(i,j) = w(i,j,k+1)*p5*(t(i,j,k,n,tau)-t(i,j,k+1,n,tau))
      endif
      enddo
      enddo
c
      do j=2,jmm
      do i=2,imm
      wka(i,j,k) = - rdxt(j)*(a(i+1,j)+a(i,j))
     &             - rdyt(j)*(b(i,j+1)+b(i,j))
     &             - rdz(k) *(wa(i,j)+wb(i,j))
      enddo
      enddo
c
100   continue
c
c
c-----------------------------------------------------------------------
c     diffusions
c-----------------------------------------------------------------------
c
#ifdef gm90
c
c     compute the isopycnal/dipycnal mixing
c     xz and yz isopycnal diffusive flux are solved explicitly;
c     while zz component will be solved implicitly.
c
      call isoflux (wka,n)
c
#else
c
c     horizontal diffusion
c
      do 200 k=1,km
      do j=2,jmm
      do i=2,imt
      a(i,j) = fx(i,j)*(t(i,j,k,n,taum)-t(i-1,j,k,n,taum))*
     &         tmask(i,j,k)*tmask(i-1,j,k)
      enddo
      enddo
c
      do j=2,jmt
      do i=2,imm
      b(i,j) = fy(i,j)*(t(i,j,k,n,taum)-t(i,j-1,k,n,taum))*
     &         tmask(i,j,k)*tmask(i,j-1,k)
      enddo
      enddo
c
      do j=2,jmm
      do i=2,imm
      wka(i,j,k) = wka(i,j,k) + ah*( sdxt(j)*(a(i+1,j)-a(i,j))+
     &                               r1a(j)*b(i,j+1)-r1b(j)*b(i,j) )
      enddo
      enddo
200   continue
c
#endif
c
c-----------------------------------------------------------------------
c     vertical diffusion
c-----------------------------------------------------------------------
c
      do 300 k=1,km
      if(k.eq.1)then
       do j=2,jmm
       do i=2,imm
       flxa(i,j) = c0
       enddo
       enddo
      else
       do j=2,jmm
       do i=2,imm
       flxa(i,j) = flxb(i,j)
       enddo
       enddo
      endif
c
      do j=2,jmm
      do i=2,imm
      if(k.eq.itn(i,j))then
        flxb(i,j) = c0
      else if(k.lt.itn(i,j))then
        flxb(i,j) = (t(i,j,k,n,taum)-t(i,j,k+1,n,taum))*
#ifdef gm90
     &              gravr*wkb(i,j,k)*rdzw(k)
#else
     &              gravr*kappa_h(i,j,k)   *rdzw(k)
#endif
      else
        flxb(i,j) = c0
      endif
      enddo
      enddo
c
      do j=2,jmm
      do i=2,imm
      wka(i,j,k) = wka(i,j,k) + rdz(k)*(flxa(i,j)-flxb(i,j))*aidif
      enddo
      enddo
300   continue
c
c
c-----------------------------------------------------------------------
c     + surface forcing
c-----------------------------------------------------------------------
c
c
      if(n.eq.1)then
c-----------------------------------------------------------------------
c     set sea surface heat flux b.c.
c
c     dT/dt = D(T*-T)/Cp/g/rho/dz
c
c-----------------------------------------------------------------------
      do j=2,jmm
      do i=2,imm
ccc    stf(i,j)   = ddd(i,j)*(bct(i,j)-t(i,j,1,n,tau))
       stf(i,j)   = gamma_t *(bct(i,j)-t(i,j,1,n,tau))
       wka(i,j,1) = wka(i,j,1) + stf(i,j)*rdz(1)*aidif
      enddo
      enddo
c
      else
c-----------------------------------------------------------------------
c     set natural or restoring b.c. for salinity
c-----------------------------------------------------------------------
      do j=2,jmm
      do i=2,imm
#ifdef snbc
       stf(i,j)   = t(i,j,1,n,tau)*emp(i,j)
#else
       stf(i,j)   = gamma_s*(bcs(i,j)-t(i,j,1,n,tau))*pbar(i,j)*dz(1)
#endif
       wka(i,j,1) = wka(i,j,1) + stf(i,j)*rdz(1)*aidif
      enddo
      enddo
c
      endif
c
c
c-----------------------------------------------------------------------
c     compute T&S at tau+1 time level. Doesn't include implicit diffusion
c-----------------------------------------------------------------------
      do k=1,km
      do j=2,jmm
      do i=2,imm
      wka(i,j,k) = t(i,j,k,n,taum) +
     &             wka(i,j,k)/pbar(i,j)*tmask(i,j,k)*dtts2
      enddo
      enddo
      enddo
c
c
c-----------------------------------------------------------------------
c     add the component due to implicit diffusion
c-----------------------------------------------------------------------
c
#ifdef gm90
      call invtri (wka,stf,wkb,aidif,dtts2,pbar)
#endif
#if defined implicitvmix
      call invtri (wka,stf,kappa_h,aidif,dtts2,pbar)
#endif
c
c
c-----------------------------------------------------------------------
c     set T > -1.5 temporarily since no seaice model
c-----------------------------------------------------------------------
      if(n.eq.1)then
      do k=1,km
      do j=2,jmm
      do i=2,imm
      wka(i,j,k) = dmax1(tbice, wka(i,j,k))
      enddo
      enddo
      enddo
      endif
c
c
c-----------------------------------------------------------------------
c     T&S=wka  & Filter & periodic b.c
c-----------------------------------------------------------------------
c
      if(leapfrog_t) then
        do k=1,km
        do j=2,jmm
        do i=2,imm
        t(i,j,k,n,taum) = t(i,j,k,n,tau)
#ifdef asselin_t
     &                  *aft2 +  aft1*(t(i,j,k,n,taum)+wka(i,j,k))
#endif
        enddo
        enddo
        enddo
      endif
c
      do k=1,km
      do j=2,jmm
      do i=2,imm
      t(i,j,k,n,tau) = wka(i,j,k)
      enddo
      enddo
c
      do j=2,jmm
      t(1  ,j,k,n,tau) = t(imm,j,k,n,tau)
      t(imt,j,k,n,tau) = t(2  ,j,k,n,tau)
      t(1  ,j,k,n,taum)= t(imm,j,k,n,taum)
      t(imt,j,k,n,taum)= t(2  ,j,k,n,taum)
      enddo
      enddo
c
1000  continue
c
      if(.not.leapfrog_t) leapfrog_t = .true.
c
      return
      end
