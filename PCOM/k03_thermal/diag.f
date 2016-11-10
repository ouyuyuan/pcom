c
c     ===============
      subroutine diag
c     ===============
c
      implicit none
      include 'param.h'
      include 'pconst.h'
      include 'scalar.h'
      include 'prog.h'
      include 'grdvar.h'
      include 'diag.h'
      include 'calendar.h'
c
      character fname*10,ftail*5
      real psi(imt,jmt)
      integer bous
      real t1,h,hx,ek,ep,ea,eb,et,es,msf,maxmsf,pij1,pij2
     *     ,psouth,pnorth,ev,psisou,psinor,dpsi,minmsf
c
      common /works/ msf(jmt,kmp1)
c
c
c     for assemble the file name for output
      encode(5,'(i5)',ftail) year+10000
c
c
      if(mod(year,io_tsuvp).ne.0) goto 100
c
c-----------------------------------------------------------------------
c     initialize common /cdiag/
c-----------------------------------------------------------------------
#ifdef iomth
      if(day.eq.1)then
#else
      if(day.eq.1.and.mth.eq.1) then
#endif
        do k=1,km
        do j=1,jmt
        do i=1,imt
        tmn(i,j,k) = c0
        smn(i,j,k) = c0
        umn(i,j,k) = c0
        vmn(i,j,k) = c0
        enddo
        enddo
        enddo
        do j=1,jmt
        do i=1,imt
        pmn(i,j) = c0
        enddo
        enddo
      endif
c
c
c-----------------------------------------------------------------------
c     accumulate for monthly/annual mean output
c-----------------------------------------------------------------------
      do k=1,km
      do j=1,jmt
      do i=1,imt
      tmn(i,j,k) = tmn(i,j,k) + t(i,j,k,1,tau)
      smn(i,j,k) = smn(i,j,k) + t(i,j,k,2,tau)
      umn(i,j,k) = umn(i,j,k) + up(i,j,k,tau)
      vmn(i,j,k) = vmn(i,j,k) + vp(i,j,k,tau)
      enddo
      enddo
      enddo
c
      do j=1,jmt
      do i=1,imt
      pmn(i,j) = pmn(i,j) + pbt(i,j,tau)
      enddo
      enddo
c
c
      if(day.ne.daypm(mth)) go to 100
c
c-----------------------------------------------------------------------
c     monthly/annual mean output
c-----------------------------------------------------------------------
#ifdef iomth
      fname='M'//ftail(2:5)//monname(mth)
      open(87,file=fname,form='unformatted',status='unknown')
      write(87) month,tmn,smn,umn,vmn,pmn
      close(87)
#else
      if (mth.eq.12) then
      fname='Y'//ftail(2:5)
      open(87,file=fname,form='unformatted',status='unknown')
      write(87) year,tmn,smn,umn,vmn,pmn
      close(87)
      endif
#endif
c
c
100   continue
c
c
      if(day.eq.daypm(mth)) then
c
c-----------------------------------------------------------------------
c     save datafile for restarting integration
c-----------------------------------------------------------------------
      if(mth.eq.12.and.mod(year,io_restr).eq.0)then
      fname='S'//ftail(2:5)
      open(87,file=fname,form='unformatted',status='unknown')
      write(87)month,t,pbt,up,vp,upb,vpb,pmtp,ump,vmp
      close(87)
      endif
c
#ifdef  daily
      endif
#endif
c
c
c-----------------------------------------------------------------------
c     daily/monthly output
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c     total K.E.
c-----------------------------------------------------------------------
      ek = c0
      do k=1,km
      do j=2,jmm
      do i=2,imm
      ek = ek + dz(k)*dxdyu(j)*(up(i,j,k,tau)**2+vp(i,j,k,tau)**2)
      enddo
      enddo
      enddo
      ek = ek*p5/grav
c
c
c-----------------------------------------------------------------------
c     global mean sst & sss (psu)
c-----------------------------------------------------------------------
      ea = c0
      eb = c0
      do j=2,jmm
      do i=2,imm
      ea = ea + t(i,j,1,1,tau)*dxdyt(j)*tmask(i,j,1)
      eb = eb + t(i,j,1,2,tau)*dxdyt(j)*tmask(i,j,1)
      enddo
      enddo
      ea = ea/area
      eb = eb/area
c
c
c-----------------------------------------------------------------------
c     global mean temperature & salinity (psu)
c-----------------------------------------------------------------------
      et = c0
      es = c0
      ep = c0
      do k=1,km
      do j=2,jmm
      do i=2,imm
      t1 = dz(k)*dxdyt(j)*pbt(i,j,tau)*tmask(i,j,k)
      et = et + t(i,j,k,1,tau)*t1
      es = es + t(i,j,k,2,tau)*t1
      ep = ep + t1
      enddo
      enddo
      enddo
      et = et/ep
      es = es/ep
c
c
c-----------------------------------------------------------------------
c     meridional mass transport   (rho cm**3/s)
c-----------------------------------------------------------------------
      do j=2,jmm
      msf(j,1)    = c0
      msf(j,kmp1) = c0
      enddo
c
      do k=2,km
      do j=2,jmm
        hx = dz(k)/grav/rdxt(j)*c1em12
        t1 = c0
        do i=2,imm
        t1 = t1 + vp(i,j,k-1,tau)*spbt(i,j)*hx
        enddo
        msf(j,k)=msf(j,k-1) + t1
      enddo
      enddo
c
      maxmsf = c0
      minmsf = c0
      do k=1,km
      do j=2,jmm
      if(msf(j,k).gt.maxmsf) maxmsf = msf(j,k)
      if(msf(j,k).lt.minmsf) minmsf = msf(j,k)
      enddo
      enddo
c
c     surface elevation
c---------------------------------------------------------------------
      bous=1
      do j=1,jmt
      do i=1,imt
      if(tmask(i,j,1).gt.c0)then
      if(bous.eq.1)then
        psi(i,j) = phib(i,j)/grav
        do k=1,itn(i,j)
        psi(i,j) = psi(i,j)+dz(k)*pbt(i,j,tau)*rho(i,j,k)/grav
        enddo
      endif
      else
       psi(i,j) = c0
      endif
      enddo
      enddo
c
      ev = c0
      do j=2,jmm
      do i=2,imm
      ev = ev + psi(i,j)*dxdyt(j)*tmask(i,j,1)/area
      enddo
      enddo
      do j=2,jmm
      do i=2,imm
      psi(i,j) = psi(i,j) - ev
      enddo
      enddo
c
      psisou=0.0
      psinor=0.0
      do i=2,31
      psisou=psisou+psi(i,2)/30.0d0
      psinor=psinor+psi(i,31)/30.0d0
      enddo
      dpsi=psisou-psinor
c   -----------------------
c
      psouth=0.0
      pnorth=0.0
      do i=2,31
      psouth=psouth+pbt(i,2,tau)*pn(i,2)/30.0d0/grav/1.03
      pnorth=pnorth+pbt(i,31,tau)*pn(i,31)/30.0d0/grav/1.03
      enddo
         
      pij1=pbt(2,31,tau)-pbt(2,2,tau)
      pij2=pbt(31,31,tau)-pbt(31,2,tau)
c
      write(6,299) month,day,ek,ea,et,maxmsf,minmsf,(pnorth-psouth)
     *     ,dpsi,psisou,psinor,ev
      write(16,299) month,day,ek,ea,et,maxmsf,(pnorth-psouth)
     *     ,dpsi,psisou,psinor,ev
c 299  format (i6,4f20.16)
 299  format(i5,i3,e13.6,2f8.4,f9.4,8f10.4)
c      write(6,200) month,day,ek,ea,eb,et,es,maxmsf
c200   format(i5,i3,e13.6,4f10.6,f6.2)
c
#if !defined  daily
      endif
#endif
c
      return
      end
