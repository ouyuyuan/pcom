c     =================
      subroutine inirun
c     =================
c     initialization
c
      implicit none
      include 'param.h'
      include 'pconst.h'
      include 'scalar.h'
      include 'grdvar.h'
      include 'prog.h'
      include 'cvbc.h'
      include 'calendar.h'
      include 'diag.h'
c
      real   a,b,c,d,t0
      common /works/ a(imt,jmt),b(imt,jmt),c(imt,jmt)
      common /works/ d(imt,jmt,km)
      real fx,fy
      common /wisop/ fx(imt,jmt),fy(imt,jmt)
c
c-----------------------------------------------------------------------
c     asselin temporal filter parameter
c-----------------------------------------------------------------------
      afb2 = c1-c2*afb1
      afc2 = c1-c2*afc1
      aft2 = c1-c2*aft1
c
c
c-----------------------------------------------------------------------
c     set contral parameter of integration
c-----------------------------------------------------------------------
      nss = int(24.0/dtts)
      ncc = int(dtts/dtuv)
      nbb = int(dtuv*c60/dtsf)
c
#ifdef accelerate
      ncc = 1
#endif
c
      onbb = c1/float(nbb+1)
      oncc = c1/float(ncc+1)
      onbc = c1/float(nbb*ncc+1)
c
      write(6,801)nss,ncc,nbb
801   format(/1x,'nss=',i3,2x,'ncc=',i3,2x,'nbb=',i3/)
c
c
c-----------------------------------------------------------------------
c     set time step to the unit of second
c-----------------------------------------------------------------------
      dtsf   = dtsf * c60
      dtuv   = dtuv * c3600
      dtts   = dtts * c3600
c
      c2dtsf = dtsf * c2
      c2dtuv = dtuv * c2
      c2dtts = dtts * c2
c
c
c-----------------------------------------------------------------------
c     parameters used for semi-implicitly handle of Coriolis term
c-----------------------------------------------------------------------
      do j=1,jmt
      t0      = p5*ff(j)*dtuv
      epea(j) = c1/(c1+t0*t0)
      epeb(j) = t0/(c1+t0*t0)
      enddo
c
      do j=1,jmt
      t0      = p5*ff(j)*c2dtuv
      epla(j) = c1/(c1+t0*t0)
      eplb(j) = t0/(c1+t0*t0)
      enddo
c
      do j=1,jmt
      t0      = p5*ff(j)*dtsf
      ebea(j) = c1/(c1+t0*t0)
      ebeb(j) = t0/(c1+t0*t0)
      enddo
c
      do j=1,jmt
      t0      = p5*ff(j)*c2dtsf
      ebla(j) = c1/(c1+t0*t0)
      eblb(j) = t0/(c1+t0*t0)
      enddo
c
c-----------------------------------------------------------------------
c     surface forcing fields
c-----------------------------------------------------------------------
c     bcf(i,j,12,1) :   sea surface zonal windstres       (dynes/cm**2)
c     bcf(i,j,12,2) :   sea surface meridional windstres  (dynes/cm**2)
c     bcf(i,j,12,3) :   sea surface air temperature       (celsius)
c     bcf(i,j,12,4) :   sea surface air presure           (dynes/cm**2)
c     bcf(i,j,12,5) :   sea surface salinity              (psu)
c     bcf(i,j,12,6) :   rate of evaporation minus precipitation (cm/s)
c     bcf(i,j,12,7) :   coefficient for calculation of HF (w/m2/c)
c
      open(80,file='sbcf.data',form='unformatted',status='old')
      read(80) bcf
      close(80)
c
c
c-----------------------------------------------------------------------
c     read Levitus annual mean temperature and salinity
c-----------------------------------------------------------------------
      open(81,file='tsobs.data',status='old',form='unformatted')
      read(81)pt,ps
      close(81)
c
      do k=1,km
      do j=1,jmt
      pt(1  ,j,k) = pt(imm,j,k)
      pt(imt,j,k) = pt(2  ,j,k)
      ps(1  ,j,k) = ps(imm,j,k)
      ps(imt,j,k) = ps(2  ,j,k)
      enddo
      enddo
c
      do k=1,km
      do j=1,jmt
      do i=1,imt
      t(i,j,k,1,tau ) = pt(i,j,k)
      t(i,j,k,2,tau ) = ps(i,j,k)
      t(i,j,k,1,taum) = t(i,j,k,1,tau)
      t(i,j,k,2,taum) = t(i,j,k,2,tau)
      enddo
      enddo
      enddo
c
c
c-----------------------------------------------------------------------
c     initialize pbt & rho & fixp (for BCOM)
c-----------------------------------------------------------------------
      call setpbt
c
      do j=1,jmt
      do i=1,imt
      pbt(i,j,taum) = pbt(i,j,tau)
      enddo
      enddo
c
c
c-----------------------------------------------------------------------
c     initialize  velocities
c-----------------------------------------------------------------------
      do n=1,2
      do k=1,km
      do j=1,jmt
      do i=1,imt
      up(i,j,k,n) = c0
      vp(i,j,k,n) = c0
      enddo
      enddo
      enddo
      do j=1,jmt
      do i=1,imt
      upb(i,j,n)  = c0
      vpb(i,j,n)  = c0
      enddo
      enddo
      enddo
c
c---------------------------------------------------------------------
c     initialize diagnostic variables in prog.h
c---------------------------------------------------------------------
      do j=1,jmt
      do i=1,imt
      spbt(i,j) = c1
      enddo
      enddo
c
      do k=1,kmp1
      do j=1,jmt
      do i=1,imt
      w(i,j,k) = c0
      enddo
      enddo
      enddo
c
c---------------------------------------------------------------------
c     initialize wroking arrays in prog.h
c---------------------------------------------------------------------
      do j=1,jmt
      do i=1,imt
      dub(i,j) = c0
      dvb(i,j) = c0
      enddo
      enddo
c
      do k=1,km
      do j=1,jmt
      do i=1,imt
      du(i,j,k)    = c0
      dv(i,j,k)    = c0
      diffu(i,j,k) = c0
      diffv(i,j,k) = c0
      pt(i,j,k)    = c0
      ps(i,j,k)    = c0
      enddo
      enddo
      enddo
c
      do j=1,jmt
      do i=1,imt
      pmup(i,j) = pbt(i,j,tau)
      pmtp(i,j) = pbt(i,j,tau)
      pmum(i,j) = pbt(i,j,tau)
      pmtm(i,j) = pbt(i,j,tau)
      enddo
      enddo
c
      do k=1,km
      do j=1,jmt
      do i=1,imt
      ump(i,j,k) = c0
      vmp(i,j,k) = c0
      umm(i,j,k) = c0
      vmm(i,j,k) = c0
      enddo
      enddo
      enddo
c
c-----------------------------------------------------------------------
c     initialize  common /works/
c-----------------------------------------------------------------------
      do k=1,km
      do j=1,jmt
      do i=1,imt
      a(i,j)   = c0
      b(i,j)   = c0
      c(i,j)   = c0
      d(i,j,k) = c0
      enddo
      enddo
      enddo
      do j=1,jmt
      do i=1,imt
      fx(i,j)   = c1
      fy(i,j)   = c1
      enddo
      enddo
c
c-----------------------------------------------------------------------
c     initialize  pressure terms
c-----------------------------------------------------------------------
      do j=1,jmt
      do i=1,imt
      pax (i,j) = c0
      pay (i,j) = c0
      pbxn(i,j) = c0
      pcxn(i,j) = c0
      pdxn(i,j) = c0
      pbxs(i,j) = c0
      pcxs(i,j) = c0
      pdxs(i,j) = c0
      pbye(i,j) = c0
      pcye(i,j) = c0
      pdye(i,j) = c0
      pbyw(i,j) = c0
      pcyw(i,j) = c0
      pdyw(i,j) = c0
      do k=1,km
      rhodp(i,j,k) = c0
      enddo
      enddo
      enddo
c
      do j=1,jmt
      do i=1,imt
      phibx(i,j) = c0
      phiby(i,j) = c0
      enddo
      enddo
      do j=2,jmm
      do i=2,imm
      phibx(i,j) = (phib(i+1,j  )-phib(i,j  ))*rdxt(j  ) +
     &             (phib(i+1,j+1)-phib(i,j+1))*rdxt(j+1)
      phiby(i,j) = (phib(i,j+1)-phib(i,j)+phib(i+1,j+1)-phib(i+1,j))*rdy
      enddo
      enddo
c
      do j=2,jmm
      phibx(1  ,j) = phibx(imm,j)
      phibx(imt,j) = phibx(2  ,j)
      phiby(1  ,j) = phiby(imm,j)
      phiby(imt,j) = phiby(2  ,j)
      enddo
c
c
c-----------------------------------------------------------------------
c     initialize the cumulative monthes of the integration
c-----------------------------------------------------------------------
      month = 1
c
c-----------------------------------------------------------------------
c     read in data for restarting integration
c-----------------------------------------------------------------------
      if(restrt)then
      open(22,file='restr.data',form='unformatted',status='old')
      read(22)month,t,pbt,up,vp,upb,vpb,pmtp,ump,vmp
      close(22)
cxin      month = month + 1
      month =  1    !this will print out, strating from month 1
      endif
c
      return
      end
