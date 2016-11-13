!
!     ===============
      subroutine diag(tmn,smn,pmn,umn,vmn,wmn,t,pbt,spbt,up,vp,upb,vpb,pmtp,ump,vmp,   &
                      phib,rho,year,month,mth,day,daypm,io_tsuvp,io_restr,pn,dz,rdxt,  &
                      dxdyt,dxdyu,tmask,umask,itn,area,imt,jmt,km,nt,imm,jmm,kmp1,nss,  &
                      reczhy,myid,ncpux,ncpuy,west,east,north,south,mat_myid,simt,sjmt)
!     ===============
!
      implicit none
      include 'pconst.h'
      include 'mpif.h'
!

      integer imt,jmt,km,nt,imm,jmm,kmp1,i,j,k,i2,j2,nss
      character fname*10,ftail*5
      integer myid,ncpux,ncpuy,west,east,north,south,simt,sjmt
      integer   daypm(12)
      integer   month,year,mth,day
      integer bous,reczhy
      real*4 uzhy(simt-2,sjmt,km),vzhy(simt-2,sjmt,km),wzhy(simt-2,sjmt,km)
      real*4 tzhy(simt-2,sjmt,km),szhy(simt-2,sjmt,km),pzhy(simt-2,sjmt),undefzhy
      real t1,h,hx,ek,ep,ea,eb,et,es,msf(jmt,kmp1),maxmsf,pij1,pij2
      real psouth,pnorth,ev,psisou,psinor,dpsi,minmsf,mnfactor
           
      real tmn(imt,jmt,km),smn(imt,jmt,km),umn(imt,jmt,km),vmn(imt,jmt,km)
      real pmn(imt,jmt),t(imt,jmt,km,nt,2),pbt(imt,jmt,2),spbt(imt,jmt)
      real upb(imt,jmt,2),vpb(imt,jmt,2),up(imt,jmt,km,2),vp(imt,jmt,km,2)
      real pmtp(imt,jmt),ump(imt,jmt,km),vmp(imt,jmt,km)
      real phib(imt,jmt),rho(imt,jmt,km),tmask(imt,jmt,km),umask(imt,jmt,km)
      real dz(km),rdxt(jmt),dxdyt(jmt),dxdyu(jmt),area
      integer itn(imt,jmt)
      integer io_tsuvp,io_restr
      real wmn(imt,jmt,km),pn(imt,jmt)
      
      integer mat_myid(ncpux+2,ncpuy)
      real stmn(simt,sjmt,km),ssmn(simt,sjmt,km),sumn(simt,sjmt,km)
      real svmn(simt,sjmt,km),spmn(simt,sjmt)
      real st(simt,sjmt,km,nt,2),sup(simt,sjmt,km,2),svp(simt,sjmt,km,2)
      real pbts(simt,sjmt,2),supb(simt,sjmt,2),svpb(simt,sjmt,2)
      real spmtp(simt,sjmt),sump(simt,sjmt,km),svmp(simt,sjmt,km)
      real sw(simt,sjmt,km),smsf(sjmt,kmp1)
      real sdxdyt(sjmt),spsi(simt,sjmt)
      real stmask(simt,sjmt,km),sumask(simt,sjmt,km)
!
!     for assemble the file name for output
      encode(5,'(i5)',ftail) year+10000
!

      call gath_array_real3d(pbts,pbt,mat_myid,ncpux,ncpuy,simt,sjmt,  &
                             2,imt,jmt,myid)
      call gath_array_real3d(stmask,tmask,mat_myid,ncpux,ncpuy,simt,sjmt,  &
                             km,imt,jmt,myid)
      call gath_array_real3d(sumask,umask,mat_myid,ncpux,ncpuy,simt,sjmt,  &
                             km,imt,jmt,myid)
      call gath_array_real1d(sdxdyt,dxdyt,mat_myid,ncpux,ncpuy,sjmt,jmt,myid)

!-----------------------------------------------------------------------
!     accumulate for monthly/annual mean output
!-----------------------------------------------------------------------

      
      if (day.eq.daypm(mth)) then
      mnfactor=c1/(daypm(mth)*nss)
      
      do k=1,km
      do j=1,jmt
      do i=1,imt
      tmn(i,j,k) = tmn(i,j,k)*mnfactor
      smn(i,j,k) = smn(i,j,k)*mnfactor
      umn(i,j,k) = umn(i,j,k)*mnfactor
      vmn(i,j,k) = vmn(i,j,k)*mnfactor
      wmn(i,j,k) = wmn(i,j,k)*mnfactor
      enddo
      enddo
      enddo
      
      do j=1,jmt
      do i=1,imt
      pmn(i,j) = pmn(i,j)*mnfactor
      enddo
      enddo
      
      call gath_array_real2d(spsi,pmn,mat_myid,ncpux,ncpuy,simt,sjmt, &
                             imt,jmt,myid)
                             
      if (myid==0) then
      ev = c0
      do j=2,sjmt-1
      do i=2,simt-1
      ev = ev + spsi(i,j)*sdxdyt(j)*stmask(i,j,1)/area
      enddo
      enddo
      do j=2,sjmt-1
      do i=2,simt-1
      spsi(i,j) = spsi(i,j) - ev
      enddo
      enddo
      end if
!
!-----------------------------------------------------------------------
!     monthly mean output
!-----------------------------------------------------------------------
      call gath_array_real3d(sw,wmn,mat_myid,ncpux,ncpuy,simt,sjmt,  &
                             km,imt,jmt,myid)
      call gath_array_real3d(stmn,tmn,mat_myid,ncpux,ncpuy,simt,sjmt,  &
                             km,imt,jmt,myid)
      call gath_array_real3d(ssmn,smn,mat_myid,ncpux,ncpuy,simt,sjmt,  &
                             km,imt,jmt,myid)
      call gath_array_real3d(sumn,umn,mat_myid,ncpux,ncpuy,simt,sjmt,  &
                             km,imt,jmt,myid)
      call gath_array_real3d(svmn,vmn,mat_myid,ncpux,ncpuy,simt,sjmt,  &
                             km,imt,jmt,myid)
                             
      if (myid==0) then
      undefzhy=999999.9
      uzhy=undefzhy
      vzhy=undefzhy
      wzhy=undefzhy
      tzhy=undefzhy
      szhy=undefzhy
      pzhy=undefzhy
      do i=2,simt-1
        do j=1,sjmt
          do k=1,km
          if (sumask(i,j,k).gt.0) then
            uzhy(i-1,j,k)=sumn(i,j,k)/100.0
            vzhy(i-1,j,k)=svmn(i,j,k)/100.0
            wzhy(i-1,j,k)=sw(i,j,k)/100.0
          end if
          if (stmask(i,j,k).gt.0) then
            tzhy(i-1,j,k)=stmn(i,j,k)
            szhy(i-1,j,k)=ssmn(i,j,k)
          end if
         end do
         if (stmask(i,j,1).gt.0) then
            pzhy(i-1,j)=spsi(i,j)
         end if
        end do
      end do
      write(17,rec=reczhy) uzhy,vzhy,wzhy,tzhy,szhy,pzhy
      reczhy=reczhy+1
      end if
!
!-----------------------------------------------------------------------
!     save datafile for restarting integration
!-----------------------------------------------------------------------
      if(mth.eq.12.and.mod(year,io_restr).eq.0.and.day.eq.daypm(mth))then
      fname='S'//ftail(2:5)
      call gath_array_real3d(supb,upb,mat_myid,ncpux,ncpuy,simt,sjmt,  &
                             2,imt,jmt,myid)
      call gath_array_real3d(svpb,vpb,mat_myid,ncpux,ncpuy,simt,sjmt,  &
                             2,imt,jmt,myid)
      call gath_array_real3d(sump,ump,mat_myid,ncpux,ncpuy,simt,sjmt,  &
                             km,imt,jmt,myid)
      call gath_array_real3d(svmp,vmp,mat_myid,ncpux,ncpuy,simt,sjmt,  &
                             km,imt,jmt,myid)
      call gath_array_real2d(spmtp,pmtp,mat_myid,ncpux,ncpuy,simt,sjmt, &
                             imt,jmt,myid)
      call gath_array_real4d(sup,up,mat_myid,ncpux,ncpuy,simt,sjmt,  &
                             km,2,imt,jmt,myid)
      call gath_array_real4d(svp,vp,mat_myid,ncpux,ncpuy,simt,sjmt,  &
                             km,2,imt,jmt,myid)
      call gath_array_real5d(st,t,mat_myid,ncpux,ncpuy,simt,sjmt,  &
                             km,nt,2,imt,jmt,myid)
      if (myid==0) then
      open(87,file=fname,form='unformatted',status='unknown')
      write(87)month,st,pbts,sup,svp,supb,svpb,spmtp,sump,svmp
      close(87)
      end if
      end if
!
      end if
!
!
!-----------------------------------------------------------------------
!     daily/monthly output
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!     total K.E.
!-----------------------------------------------------------------------
      ek = c0
      do k=1,km
      do j=2,jmm
      do i=2,imm
      ek = ek + dz(k)*dxdyu(j)*(up(i,j,k,tau)**2+vp(i,j,k,tau)**2)
      enddo
      enddo
      enddo
      call gath_var_real(ek,mat_myid,ncpux,ncpuy,myid)
      ek = ek*p5/grav
!-----------------------------------------------------------------------
!     global mean sst & sss (psu)
!-----------------------------------------------------------------------
      ea = c0
      eb = c0
      do j=2,jmm
      do i=2,imm
      ea = ea + t(i,j,1,1,tau)*dxdyt(j)*tmask(i,j,1)
      eb = eb + t(i,j,1,2,tau)*dxdyt(j)*tmask(i,j,1)
      enddo
      enddo
      call gath_var_real(ea,mat_myid,ncpux,ncpuy,myid)
      call gath_var_real(eb,mat_myid,ncpux,ncpuy,myid)
      if (area.gt.0) then
      ea = ea/area
      eb = eb/area
      end if
!
!-----------------------------------------------------------------------
!     global mean temperature & salinity (psu)
!-----------------------------------------------------------------------
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
      call gath_var_real(et,mat_myid,ncpux,ncpuy,myid)
      call gath_var_real(es,mat_myid,ncpux,ncpuy,myid)
      call gath_var_real(ep,mat_myid,ncpux,ncpuy,myid)
      if (ep.gt.0) then
      et = et/ep
      es = es/ep
      end if
!
!-----------------------------------------------------------------------
!     meridional mass transport   (rho cm**3/s)
!-----------------------------------------------------------------------
      do j=2,jmm
      msf(j,1)    = c0
      msf(j,kmp1) = c0
      enddo
!
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
!     
      
      call gath_merid_real2d(smsf,msf,mat_myid,ncpux,ncpuy,sjmt, &
                             km,jmt,myid)
      if (myid==0) then
      maxmsf = c0
      minmsf = c0
      do k=1,km
      do j=2,sjmt-1
      if(smsf(j,k).gt.maxmsf) maxmsf = smsf(j,k)
      if(smsf(j,k).lt.minmsf) minmsf = smsf(j,k)
      enddo
      enddo

      psisou=0.0
      psinor=0.0
      do i=2,31
      psisou=psisou+spsi(i,2)/30.0d0
      psinor=psinor+spsi(i,31)/30.0d0
      enddo
      dpsi=psisou-psinor
!   -----------------------
!
      psouth=0.0
      pnorth=0.0
      do i=2,31
      psouth=psouth+pbt(i,2,tau)*pn(i,2)/30.0d0/grav/1.03
      pnorth=pnorth+pbt(i,31,tau)*pn(i,31)/30.0d0/grav/1.03
      enddo
         
      pij1=pbt(2,31,tau)-pbt(2,2,tau)
      pij2=pbt(31,31,tau)-pbt(31,2,tau)
!
      write(66,"(i5,i3,e13.6,2f9.4,2f12.4,5f15.4)") month,day,ek,  &
           ea,et,maxmsf,minmsf,(pnorth-psouth),dpsi,psisou,psinor,ev
      end if
!
      return
      end
