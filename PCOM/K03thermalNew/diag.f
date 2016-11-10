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
     *     ,psouth,pnorth,ev,psisou,psinor,dpsi,zsum,phfmax,phf(jmt)
     *     ,hpsi(imt,jmt),u_sum(imt,jmt),v_sum(imt,jmt),tusum,tvsum,huv
     *     ,zsf(imt,kmp1),zaxmsf
     *     ,amsf,azsf,ahpsi,au_sum,av_sum
c
      common /works/ msf(jmt,kmp1)
      common /sums/ amsf(jmt,kmp1),azsf(imt,kmp1),ahpsi(imt,jmt)
     *              ,au_sum(imt,jmt),av_sum(imt,jmt)
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
      do k=1,km
      do j=2,jmm
      if(msf(j,k).gt.maxmsf) maxmsf = msf(j,k)
      enddo
      enddo
c
c-----------------------------------------------------------------------
c     zonal mass transport   (rho cm**3/s)
c-----------------------------------------------------------------------
      do i=2,imm
      zsf(i,1)    = c0
      zsf(i,kmp1) = c0
      enddo
c
      do k=2,km
      do i=2,imm
        hx = dz(k)/grav/rdy*c1em12
        t1 = c0
        do j=2,jmm
        t1 = t1 + vp(i,j,k-1,tau)*spbt(i,j)*hx
        enddo
        zsf(i,k)=zsf(i,k-1) + t1
      enddo
      enddo
c
      zaxmsf = c0
      do k=1,km
      do i=2,imm
      if(zsf(i,k).lt.zaxmsf) zaxmsf = zsf(i,k)
      enddo
      enddo
c-----------------------------------------------------------------------
c     Horizontal mass transport   (rho cm**3/s)
c-----------------------------------------------------------------------
      do j=1,jmt
      hpsi(imt,j)    = c0
      enddo
c
      do j=1,jmt-2
      do i=imt-1,2,-1
        t1 = c0
ccc        do k=1,km
        do k=1,km                            
        hx = dz(k)/grav/rdxt(j)*c1em12        ! in unit of Sv
        t1 = t1 + vp(i,j,k,tau)*spbt(i,j)*hx
        enddo
        hpsi(i,j)=hpsi(i+1,j) + t1
      enddo
      enddo
      do j=2,jmt-1
      do i=2,imt-1
        tusum=c0
        tvsum=c0
ccc        do k=1,km
ccc        do k=1,27                              !The upper 980meter
        do k=28,60                              !The upper 980meter
        huv = dz(k)/grav*1.0d-4                  ! in unit of m^2/s
        tusum = tusum + up(i,j,k,tau)*spbt(i,j)*huv
        tvsum = tvsum + vp(i,j,k,tau)*spbt(i,j)*huv
        enddo
        u_sum(i,j)=tusum
        v_sum(i,j)=tvsum
      enddo
      enddo
c------------------------------------------------
c                         The poleward heat flux, in PW
      do j=2,jmm
      do k=1,km
        hx = dz(k)/grav/rdxt(j)*c1em12*4.180d-3
        t1 = c0
        do i=2,imm
        t1 = t1 + t(i,j,k,1,tau)*vp(i,j,k-1,tau)*spbt(i,j)*hx
        enddo
        phf(j)=phf(j) + t1
      enddo
      enddo

      phfmax=c0
      do j=1,jmt
         if (phf(j) .gt. phfmax) phfmax=phf(j)
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
      psisou=psisou+t(i,2,1,2,tau)/30.0d0
      psinor=psinor+t(i,31,1,2,tau)/30.0d0
cc      psisou=psisou+psi(i,2)/30.0d0
cc      psinor=psinor+psi(i,31)/30.0d0
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
      write(6,299) month,day,ek,ea,eb,et,es,maxmsf,zaxmsf,dpsi
c                 ,psisou,psinor,ev
      write(16,299)month,day,ek,ea,eb,et,es,maxmsf,dpsi,psisou,psinor,ev
 299  format(i6,i3,e13.6,2f8.4,f9.4,8f10.4)
      
      if (day .eq. 30) then
cccc -------------------------------this for sum up only
      do i=1,imm,2
      do j=1,jmm,2
cc      write (51,99) u_sum(i,j)
cc      write (52,99) v_sum(i,j)
      au_sum(i,j)=au_sum(i,j)+u_sum(i,j)
      av_sum(i,j)=av_sum(i,j)+v_sum(i,j)
      enddo
      enddo

      do i=1,imm,2
      do j=1,jmm,2
      write (81,99) u_sum(i,j)
      write (82,99) v_sum(i,j)
      enddo
      enddo

      do i=2,imm
      do j=2,jmm
ccc      write (53,99) hpsi(i,j)
      ahpsi(i,j)=ahpsi(i,j)+hpsi(i,j)
      enddo
      enddo

      do i=2,imm
      do j=2,jmm
      write (60,99) psi(i,j)
      write (61,99) t(i,j,1,1,tau)
      write (62,99) t(i,j,5,1,tau)
      write (63,99) t(i,j,10,1,tau)
      write (64,99) t(i,j,15,1,tau)
      write (74,99) (pbt(i,j,tau)-1.0)*pn(i,j)
      enddo
      enddo

      do j=2,jmm
         do k=1,km
            write (70,99) t(2,j,k,1,tau)
            write (71,99) t(10,j,k,1,tau)
            write (73,99) t(20,j,k,1,tau)
      enddo
      enddo

      do i=1,imt
         do k=1,km+1
cccc            write (68,99) zsf(i,k)
            azsf(i,k)=azsf(i,k)+zsf(i,k)
      enddo
      enddo
      do j=1,jmt
         do k=1,km+1
cccc            write (69,99) msf(j,k)
            amsf(j,k)=amsf(j,k)+msf(j,k)
      enddo
      enddo
      do j=1,jmm
         write (75,99) phf(j)
      enddo

      do j=1,jmt
         do k=1,km
         write (76,99) vp(2,j,k,tau)*spbt(2,j)/grav
      enddo
      enddo
c      zsum=0.0
c       do k=1,km
c            write (67,99) z0(k)*0.01
c       enddo
      end if

#endif
 98   format (2f12.6)
 99   format (f12.6)
 999  format (10f12.6)
c
#if !defined  daily
      endif
#endif
c
      return
      end
