!     =================
      subroutine inirun(afb1,afc1,aft1,afb2,afc2,aft2,dtts,dtuv,dtsf,nss,ncc,nbb,  &
                        onbb,oncc,onbc,c2dtsf,c2dtuv,c2dtts,epea,epeb,epla,eplb,   &
                        ebea,ebeb,ebla,eblb,bcf,pt,ps,t,pbt,up,vp,upb,vpb,spbt,w,  &
                        dub,dvb,du,dv,diffu,diffv,pmup,pmtp,pmum,pmtm,ump,vmp,umm, &
                        vmm,pax,pay,pbxn,pcxn,pdxn,pbxs,pcxs,pdxs,pbye,pcye,pdye,  &
                        pbyw,pcyw,pdyw,rhodp,phibx,phiby,phib,rdxt,rdy,ff,month,   &
                        restrt,rho,itn,imt,jmt,km,nt,phis,imm,jmm,kmp1,dz,decibar, &
                        myid,ncpux,ncpuy,west,east,north,south,mat_myid,simt,sjmt, &
                        unesco,boussinesq,fixp)
!     =================
!     initialization
!
      implicit none
      include 'pconst.h'
      include 'mpif.h'
!
      integer imt,jmt,km,nt,imm,jmm,kmp1,i,j,k,n,i2,j2,unesco,boussinesq
      real phis
      real a(imt,jmt),b(imt,jmt),c(imt,jmt),d(imt,jmt,km),t0
      real fx(imt,jmt),fy(imt,jmt),fixp(imt,jmt,km)
      integer   month

      real    afb1,afc1,aft1,afb2,afc2,aft2
      real    dtts,dtuv,dtsf,c2dtts,c2dtuv,c2dtsf
      real    onbb,oncc,onbc
      integer nss,ncc,nbb
      real ebea(jmt),ebeb(jmt),ebla(jmt),eblb(jmt)
      real epea(jmt),epeb(jmt),epla(jmt),eplb(jmt)
      real bcf(imt,jmt,12,7)
      real t(imt,jmt,km,nt,2),up(imt,jmt,km,2),vp(imt,jmt,km,2) 
      real pbt(imt,jmt,2),spbt(imt,jmt)
      real upb(imt,jmt,2),vpb(imt,jmt,2),w(imt,jmt,kmp1)
      real du(imt,jmt,km),dv(imt,jmt,km),dub(imt,jmt),dvb(imt,jmt)
      real diffu(imt,jmt,km),diffv(imt,jmt,km),pt(imt,jmt,km),ps(imt,jmt,km)
      real pmup(imt,jmt),pmum(imt,jmt),pmtp(imt,jmt),pmtm(imt,jmt)
      real ump(imt,jmt,km),umm(imt,jmt,km),vmp(imt,jmt,km),vmm(imt,jmt,km)
      real rhodp(imt,jmt,km),rho(imt,jmt,km)
      real pax(imt,jmt),pay(imt,jmt)
      real pbxn(imt,jmt),pbxs(imt,jmt)
      real pcxn(imt,jmt),pcxs(imt,jmt)
      real pdxn(imt,jmt),pdxs(imt,jmt)
      real pbye(imt,jmt),pbyw(imt,jmt)
      real pcye(imt,jmt),pcyw(imt,jmt)
      real pdye(imt,jmt),pdyw(imt,jmt)
      real phibx(imt,jmt),phiby(imt,jmt)
      real phib(imt,jmt)
      real rdxt(jmt),ff(jmt),rdy,dz(km),decibar
      logical restrt
      integer itn(imt,jmt)
      
      integer myid,ncpux,ncpuy,west,east,north,south,simt,sjmt
      integer mat_myid(ncpux+2,ncpuy)
      real sbcf(simt,sjmt,12,7)
      real spt(simt,sjmt,km),sps(simt,sjmt,km)
      real st(simt,sjmt,km,nt,2),sup(simt,sjmt,km,2),svp(simt,sjmt,km,2)
      real supb(simt,sjmt,2),svpb(simt,sjmt,2),pbts(simt,sjmt,2)
      real sump(simt,sjmt,km),svmp(simt,sjmt,km),spmtp(simt,sjmt)

!
!-----------------------------------------------------------------------
!     asselin temporal filter parameter
!-----------------------------------------------------------------------
      afb2 = c1-c2*afb1
      afc2 = c1-c2*afc1
      aft2 = c1-c2*aft1
!
!
!-----------------------------------------------------------------------
!     set contral parameter of integration
!-----------------------------------------------------------------------
      nss = int(24.0/dtts)
      ncc = int(dtts/dtuv)
      nbb = int(dtuv*c60/dtsf)

      onbb = c1/real(nbb+1)
      oncc = c1/real(ncc+1)
      onbc = c1/real(nbb*ncc+1)
!
      if (myid==0) then
      write(66,801) nss,ncc,nbb
801   format(/1x,'nss=',i3,2x,'ncc=',i3,2x,'nbb=',i3/)
      end if
!
!
!-----------------------------------------------------------------------
!     set time step to the unit of second
!-----------------------------------------------------------------------
      dtsf   = dtsf * c60
      dtuv   = dtuv * c3600
      dtts   = dtts * c3600
!
      c2dtsf = dtsf * c2
      c2dtuv = dtuv * c2
      c2dtts = dtts * c2
!
!
!-----------------------------------------------------------------------
!     parameters used for semi-implicitly handle of Coriolis term
!-----------------------------------------------------------------------
      do j=1,jmt
      t0      = p5*ff(j)*dtuv
      epea(j) = c1/(c1+t0*t0)
      epeb(j) = t0/(c1+t0*t0)
      enddo
!
      do j=1,jmt
      t0      = p5*ff(j)*c2dtuv
      epla(j) = c1/(c1+t0*t0)
      eplb(j) = t0/(c1+t0*t0)
      enddo
!
      do j=1,jmt
      t0      = p5*ff(j)*dtsf
      ebea(j) = c1/(c1+t0*t0)
      ebeb(j) = t0/(c1+t0*t0)
      enddo
!
      do j=1,jmt
      t0      = p5*ff(j)*c2dtsf
      ebla(j) = c1/(c1+t0*t0)
      eblb(j) = t0/(c1+t0*t0)
      enddo
      call swap_array_real1d(epea,jmt,north,south)
      call swap_array_real1d(epeb,jmt,north,south)
      call swap_array_real1d(epla,jmt,north,south)
      call swap_array_real1d(eplb,jmt,north,south)
      call swap_array_real1d(ebea,jmt,north,south)
      call swap_array_real1d(ebeb,jmt,north,south)
      call swap_array_real1d(ebla,jmt,north,south)
      call swap_array_real1d(eblb,jmt,north,south)
!
!-----------------------------------------------------------------------
!     surface forcing fields
!-----------------------------------------------------------------------
!     bcf(i,j,12,1) :   sea surface zonal windstres       (dynes/cm**2)
!     bcf(i,j,12,2) :   sea surface meridional windstres  (dynes/cm**2)
!     bcf(i,j,12,3) :   sea surface air temperature       (celsius)
!     bcf(i,j,12,4) :   sea surface air presure           (dynes/cm**2)
!     bcf(i,j,12,5) :   sea surface salinity              (psu)
!     bcf(i,j,12,6) :   rate of evaporation minus precipitation (cm/s)
!     bcf(i,j,12,7) :   coefficient for calculation of HF (w/m2/c)
!
      if (myid==0) then
!
!-----------------------------------------------------------------------
!     read Levitus annual mean temperature and salinity
!-----------------------------------------------------------------------
      open(81,file='tsobs.data',form='unformatted',access='direct', &
           recl=simt*sjmt*km*8,status='old')
      read(81,rec=1) spt
      read(81,rec=2) sps
      close(81)
      end if
      
      call div_array_real3d(spt,pt,mat_myid,ncpux,ncpuy,simt,sjmt,  &
                            km,imt,jmt,myid)
      call div_array_real3d(sps,ps,mat_myid,ncpux,ncpuy,simt,sjmt,  &
                            km,imt,jmt,myid)
!
!      do k=1,km
!      do j=1,jmt
!      pt(1  ,j,k) = pt(imm,j,k)
!      pt(imt,j,k) = pt(2  ,j,k)
!      ps(1  ,j,k) = ps(imm,j,k)
!      ps(imt,j,k) = ps(2  ,j,k)
!      enddo
!      enddo
!
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
!
!
!-----------------------------------------------------------------------
!     initialize pbt & rho & fixp (for BCOM)
!-----------------------------------------------------------------------
      call setpbt(rho,pbt,itn,pt,ps,dz,decibar,imt,jmt,km,kmp1,nt,unesco,boussinesq,fixp)
!
      do j=1,jmt
      do i=1,imt
      pbt(i,j,taum) = pbt(i,j,tau)
      enddo
      enddo
!
!
!-----------------------------------------------------------------------
!     initialize  velocities
!-----------------------------------------------------------------------
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
!
!---------------------------------------------------------------------
!     initialize diagnostic variables in prog.h
!---------------------------------------------------------------------
      do j=1,jmt
      do i=1,imt
      spbt(i,j) = c1
      enddo
      enddo
!
      do k=1,kmp1
      do j=1,jmt
      do i=1,imt
      w(i,j,k) = c0
      enddo
      enddo
      enddo
!
!---------------------------------------------------------------------
!     initialize wroking arrays in prog.h
!---------------------------------------------------------------------
      do j=1,jmt
      do i=1,imt
      dub(i,j) = c0
      dvb(i,j) = c0
      enddo
      enddo
!
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
!
      do j=1,jmt
      do i=1,imt
      pmup(i,j) = pbt(i,j,tau)
      pmtp(i,j) = pbt(i,j,tau)
      pmum(i,j) = pbt(i,j,tau)
      pmtm(i,j) = pbt(i,j,tau)
      enddo
      enddo
!
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
!
!-----------------------------------------------------------------------
!     initialize  common /works/
!-----------------------------------------------------------------------
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
!
!-----------------------------------------------------------------------
!     initialize  pressure terms
!-----------------------------------------------------------------------
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
!
      do j=1,jmt
      do i=1,imt
      phibx(i,j) = c0
      phiby(i,j) = c0
      enddo
      enddo
      do j=2,jmm
      do i=2,imm
      phibx(i,j) = (phib(i+1,j  )-phib(i,j  ))*rdxt(j  ) +  &
                   (phib(i+1,j+1)-phib(i,j+1))*rdxt(j+1)
      phiby(i,j) = (phib(i,j+1)-phib(i,j)+phib(i+1,j+1)-phib(i+1,j))*rdy
      enddo
      enddo
      call swap_array_real2d(phibx,imt,jmt,west,east,north,south)
      call swap_array_real2d(phiby,imt,jmt,west,east,north,south)
!
!      do j=2,jmm
!      phibx(1  ,j) = phibx(imm,j)
!      phibx(imt,j) = phibx(2  ,j)
!      phiby(1  ,j) = phiby(imm,j)
!      phiby(imt,j) = phiby(2  ,j)
!      enddo
!
!
!-----------------------------------------------------------------------
!     initialize the cumulative monthes of the integration
!-----------------------------------------------------------------------
      month = 1
!
!-----------------------------------------------------------------------
!     read in data for restarting integration
!-----------------------------------------------------------------------
      if(restrt) then
      if (myid==0) then
      open(22,file='restr.data',form='unformatted',status='old')
      read(22) month,st,pbts,sup,svp,supb,svpb,spmtp,sump,svmp
      close(22)
      end if
!xin      month = month + 1
      month =  1  !this will print out, strating from month 1
      call div_array_real5d(st,t,mat_myid,ncpux,ncpuy,simt,sjmt,  &
                            km,nt,2,imt,jmt,myid)
      call div_array_real4d(sup,up,mat_myid,ncpux,ncpuy,simt,sjmt,  &
                            km,2,imt,jmt,myid)
      call div_array_real4d(svp,vp,mat_myid,ncpux,ncpuy,simt,sjmt,  &
                            km,2,imt,jmt,myid)
      call div_array_real3d(sump,ump,mat_myid,ncpux,ncpuy,simt,sjmt,  &
                            km,imt,jmt,myid)
      call div_array_real3d(svmp,vmp,mat_myid,ncpux,ncpuy,simt,sjmt,  &
                            km,imt,jmt,myid)
      call div_array_real3d(supb,upb,mat_myid,ncpux,ncpuy,simt,sjmt,  &
                            2,imt,jmt,myid)
      call div_array_real3d(svpb,vpb,mat_myid,ncpux,ncpuy,simt,sjmt,  &
                            2,imt,jmt,myid)
      call div_array_real3d(pbts,pbt,mat_myid,ncpux,ncpuy,simt,sjmt,  &
                            2,imt,jmt,myid)
      call div_array_real2d(spmtp,pmtp,mat_myid,ncpux,ncpuy,simt,sjmt,  &
                            imt,jmt,myid)
      end if
      
!
      return
      end
