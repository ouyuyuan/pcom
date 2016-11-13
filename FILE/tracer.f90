!
!     =================
      subroutine tracer(t,w,pmtp,pmtm,umm,vmm,ump,vmp,onbc,oncc,sdxt,rdxt,rdyt,    &
                        rdz,rdzw,r1a,r1b,tmask,itn,dz,pn,bct,bcs,gamma_t,gamma_s,  &
                        du,dv,ah,kappa_h,gravr,leapfrog_t,c2dtts,dtts,imt,jmt,km,  &
                        nt,imm,jmm,kmp1,west,east,north,south,myid,asselin_t,gm90,  &
                        implicitvmix,snbc,aft1,aft2,emp,dz0,z0,rdy,cosu,K1,K2,K3,   &
                        adv_vetiso,adv_vntiso,adv_vbtiso,rhoi,e,slmxr,ahisop,    &
                        athkdf,fzisop,kref,rdz0,unesco,umn,vmn,wmn,pmn,rho,phib)
!     =================
!     compute potential temperature and salinity at tau+1 time level
!
      implicit none
      include 'pconst.h'
!      include 'isopyc.h'
      include 'mpif.h'
      integer imt,jmt,km,nt,imm,jmm,kmp1,i,j,k,n
      integer asselin_t,implicitvmix,gm90,snbc,unesco,bous
      real abc,aidif,dtts2,aft1,aft2
      real a(imt,jmt),b(imt,jmt),c(imt,jmt),wa(imt,jmt),wb(imt,jmt)
      real flxa(imt,jmt),flxb(imt,jmt),dt(imt,jmt),fx(imt,jmt),fy(imt,jmt)
      real wka(imt,jmt,km),wkb(imt,jmt,km),pbar(imt,jmt),stf(imt,jmt)
      
      real c2dtts,dtts,gravr,onbc,oncc,gamma_t,gamma_s,rdy
      real ah(imt,jmt,km),kappa_h(imt,jmt,km)
      real t(imt,jmt,km,nt,2),w(imt,jmt,kmp1)
      real pmtp(imt,jmt),pmtm(imt,jmt),emp(imt,jmt)
      real ump(imt,jmt,km),umm(imt,jmt,km),vmp(imt,jmt,km),vmm(imt,jmt,km)
      real rdxt(jmt),rdyt(jmt),r1a(jmt),r1b(jmt)
      real tmask(imt,jmt,km),rdz(km),rdzw(km),dz(km),pn(imt,jmt),sdxt(jmt)
      integer itn(imt,jmt)
      real bct(imt,jmt),bcs(imt,jmt)
      real du(imt,jmt,km),dv(imt,jmt,km)
      logical leapfrog_t
      
      real slmxr,ahisop,athkdf,fzisop(km),kref(km),rdz0(km),dz0(km),z0(km),cosu(jmt)
      real K1(imt,km,jmt,3:3),K2(imt,km,jmt,3:3),K3(imt,km,jmt,1:3)
      real adv_vetiso(imt,km,jmt),adv_vntiso(imt,km,jmt),adv_vbtiso(imt,0:km,jmt)
      real rhoi(imt,km,jmt,1:3),e(imt,kmp1,jmt,3)
      
      real umn(imt,jmt,km),vmn(imt,jmt,km),wmn(imt,jmt,km),pmn(imt,jmt),rho(imt,jmt,km)
      real psi(imt,jmt),phib(imt,jmt)
      integer west,east,north,south,myid
!
      if ((gm90==1).or.(implicitvmix==1)) then
      aidif = p5
      else
      aidif = c1
      end if
!
!-----------------------------------------------------------------------
!     calculate time-averaged pbt and mass advections
!-----------------------------------------------------------------------
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
!
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
!
!
!---------------------------------------------------------------------
!     calculate vertical mass advection
!---------------------------------------------------------------------
      call upwelling(wka,wkb,w,rdxt,rdyt,tmask,dz,pn,imt,jmt,km,imm,jmm,kmp1,  &
                     west,east,north,south,snbc,emp)
!
      do k=1,km
      do j=1,jmt
      do i=1,imt
      umn(i,j,k)=umn(i,j,k)+wka(i,j,k)/pbar(i,j)
      vmn(i,j,k)=vmn(i,j,k)+wkb(i,j,k)/(pbar(i,j)*cosu(j))
      wmn(i,j,k)=wmn(i,j,k)+w(i,j,k)*rho(i,j,k)/(pbar(i,j)*grav)
      end do
      end do
      end do
      
      bous=1
      do j=1,jmt
      do i=1,imt
      if(itn(i,j).gt.0)then
      if(bous.eq.1)then
        psi(i,j) = phib(i,j)/grav
        do k=1,itn(i,j)
        psi(i,j) = psi(i,j)+dz(k)*pbar(i,j)*rho(i,j,k)/grav
        enddo
      endif
      else
        psi(i,j) = c0
      endif
      enddo
      enddo
      
      do j=1,jmt
      do i=1,imt
      pmn(i,j)=pmn(i,j)+psi(i,j)
      end do
      end do
!
!---------------------------------------------------------------------
!     calculate 1/2 mass advection across W & S face of T cells
!---------------------------------------------------------------------
!     du & dv are used temporary as working arrays
      do k=1,km
      do j=2,jmt
      do i=2,imt
      du(i,j,k) = p25*(wka(i-1,j,k)+wka(i-1,j-1,k))
      dv(i,j,k) = p25*(wkb(i,j-1,k)+wkb(i-1,j-1,k))
      enddo
      enddo
      enddo
!
!
!---------------------------------------------------------------------
!     calculate  2*xbar(pbar) & 2*ybar(pbar)
!---------------------------------------------------------------------
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
!
!---------------------------------------------------------------------
!     diffusion coefficients
!---------------------------------------------------------------------
      if (gm90==1) then
      call isopyc(imt,jmt,km,kmp1,imm,jmm,nt,itn,tmask,kref,fzisop,rdz0,dz0,z0,  &
                  rdxt,rdyt,rdy,rdzw,dz,pn,cosu,t,slmxr,athkdf,K1,K2,K3,  &
                  adv_vetiso,adv_vntiso,adv_vbtiso,rhoi,e,fx,fy,west,east,north, &
                  south,unesco)
!
      do k=1,km
      do j=1,jmt
      do i=1,imt
      wkb(i,j,k) = kappa_h(i,j,k) + ahisop*k3(i,k,j,3)
      enddo
      enddo
      enddo
      end if
!
      if(leapfrog_t) then
        dtts2 = c2dtts
      else
        dtts2 = dtts
      endif
!
!---------------------------------------------------------------------
!     solve for one tracer at a time
!---------------------------------------------------------------------
!     n = 1 => temperature
!     n = 2 => salinity
!
!
      do 1000 n=1,nt
!
!
!---------------------------------------------------------------------
!     advections
!---------------------------------------------------------------------
      do 100 k=1,km
!
      do j=2,jmm
      do i=2,imt
      a(i,j) = du(i,j,k)*(t(i,j,k,n,tau)-t(i-1,j,k,n,tau))
      enddo
      enddo
!
      do j=2,jmt
      do i=2,imm
      b(i,j) = dv(i,j,k)*(t(i,j,k,n,tau)-t(i,j-1,k,n,tau))
      enddo
      enddo
!
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
!
      do j=2,jmm
      do i=2,imm
      if(k.ge.itn(i,j))then
        wb(i,j) = c0
      else 
        wb(i,j) = w(i,j,k+1)*p5*(t(i,j,k,n,tau)-t(i,j,k+1,n,tau))
      endif
      enddo
      enddo
!
      do j=2,jmm
      do i=2,imm
      wka(i,j,k) = - rdxt(j)*(a(i+1,j)+a(i,j)) &
                   - rdyt(j)*(b(i,j+1)+b(i,j)) &
                   - rdz(k) *(wa(i,j)+wb(i,j))
      enddo
      enddo
!
100   continue
!
!
!-----------------------------------------------------------------------
!     diffusions
!-----------------------------------------------------------------------
!
!     compute the isopycnal/dipycnal mixing
!     xz and yz isopycnal diffusive flux are solved explicitly;
!     while zz component will be solved implicitly.
!
      if (gm90==1) then
      call isoflux(wka,n,imt,jmt,imm,jmm,nt,km,kmp1,itn,gravr,rdz0,dz0,rdz,     &
                   rdzw,rdxt,rdyt,rdy,cosu,tmask,t,ah,K1,K2,K3,adv_vetiso,  &
                   adv_vntiso,adv_vbtiso,fx,fy)
      else
!
!
!     horizontal diffusion
!
      do 200 k=1,km
      do j=2,jmm
      do i=2,imt
      a(i,j) = fx(i,j)*(t(i,j,k,n,taum)-t(i-1,j,k,n,taum))* &
               tmask(i,j,k)*tmask(i-1,j,k)
      enddo
      enddo
!
      do j=2,jmt
      do i=2,imm
      b(i,j) = fy(i,j)*(t(i,j,k,n,taum)-t(i,j-1,k,n,taum))* &
               tmask(i,j,k)*tmask(i,j-1,k)
      enddo
      enddo
!
      do j=2,jmm
      do i=2,imm
      wka(i,j,k) = wka(i,j,k) + ah(i,j,k)*(sdxt(j)*(a(i+1,j)-a(i,j))+  &
                                r1a(j)*b(i,j+1)-r1b(j)*b(i,j))
      enddo
      enddo
200   continue

      end if
!
!
!-----------------------------------------------------------------------
!     vertical diffusion
!-----------------------------------------------------------------------
!
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
!
      if (gm90==1) then
      do j=2,jmm
      do i=2,imm
      if(k.eq.itn(i,j))then
        flxb(i,j) = c0
      else if(k.lt.itn(i,j))then
        flxb(i,j) = (t(i,j,k,n,taum)-t(i,j,k+1,n,taum))* &
                    gravr*wkb(i,j,k)*rdzw(k)
      else
        flxb(i,j) = c0
      endif
      enddo
      enddo
      else
      do j=2,jmm
      do i=2,imm
      if(k.eq.itn(i,j))then
        flxb(i,j) = c0
      else if(k.lt.itn(i,j))then
        flxb(i,j) = (t(i,j,k,n,taum)-t(i,j,k+1,n,taum))* &
                    gravr*kappa_h(i,j,k)*rdzw(k)
      else
        flxb(i,j) = c0
      endif
      enddo
      enddo
      end if
!
      do j=2,jmm
      do i=2,imm
      wka(i,j,k) = wka(i,j,k) + rdz(k)*(flxa(i,j)-flxb(i,j))*aidif
      enddo
      enddo
300   continue
!
!
!-----------------------------------------------------------------------
!     + surface forcing
!-----------------------------------------------------------------------
!
!
      if(n.eq.1)then
!-----------------------------------------------------------------------
!     set sea surface heat flux b.c.
!
!     dT/dt = D(T*-T)/Cp/g/rho/dz
!
!-----------------------------------------------------------------------
      do j=2,jmm
      do i=2,imm
!cc    stf(i,j)   = ddd(i,j)*(bct(i,j)-t(i,j,1,n,tau))
       stf(i,j)   = gamma_t *(bct(i,j)-t(i,j,1,n,tau))
       wka(i,j,1) = wka(i,j,1) + stf(i,j)*rdz(1)*aidif
      enddo
      enddo
!
      else
!-----------------------------------------------------------------------
!     set natural or restoring b.c. for salinity
!-----------------------------------------------------------------------
      if (snbc==1) then
      do j=2,jmm
      do i=2,imm
       stf(i,j)   = t(i,j,1,n,tau)*emp(i,j)
       wka(i,j,1) = wka(i,j,1) + stf(i,j)*rdz(1)*aidif
      enddo
      enddo
      else
      do j=2,jmm
      do i=2,imm
       stf(i,j)   = gamma_s*(bcs(i,j)-t(i,j,1,n,tau))*pbar(i,j)*dz(1)
       wka(i,j,1) = wka(i,j,1) + stf(i,j)*rdz(1)*aidif
      enddo
      enddo
      end if
      
!
      endif
!
!
!-----------------------------------------------------------------------
!     compute T&S at tau+1 time level. Doesn't include implicit diffusion
!-----------------------------------------------------------------------
      do k=1,km
      do j=2,jmm
      do i=2,imm
      wka(i,j,k) = t(i,j,k,n,taum) +                         &
                   wka(i,j,k)/pbar(i,j)*tmask(i,j,k)*dtts2
      enddo
      enddo
      enddo
!
!
!-----------------------------------------------------------------------
!     add the component due to implicit diffusion
!-----------------------------------------------------------------------
!
      if (gm90==1) then
      call invtri(wka,stf,wkb,aidif,dtts2,pbar,gravr,itn,tmask,rdz,rdzw,   &
                   imt,jmt,km,imm,jmm)
      end if
      if (implicitvmix==1) then
      call invtri(wka,stf,kappa_h,aidif,dtts2,pbar,gravr,itn,tmask,rdz,rdzw,   &
                   imt,jmt,km,imm,jmm)
      end if
!
!-----------------------------------------------------------------------
!     set T > -1.5 temporarily since no seaice model
!-----------------------------------------------------------------------
      if(n.eq.1)then
      do k=1,km
      do j=2,jmm
      do i=2,imm
      wka(i,j,k) = dmax1(tbice, wka(i,j,k))
      enddo
      enddo
      enddo
      endif
!
!
!-----------------------------------------------------------------------
!     T&S=wka  & Filter & periodic b.c
!-----------------------------------------------------------------------
!
      if(leapfrog_t) then
      if (asselin_t==1) then
      do k=1,km
      do j=2,jmm
      do i=2,imm
      t(i,j,k,n,taum) = t(i,j,k,n,tau)*aft2+aft1*(t(i,j,k,n,taum)+wka(i,j,k))
      enddo
      enddo
      enddo  
      else
      do k=1,km
      do j=2,jmm
      do i=2,imm
      t(i,j,k,n,taum) = t(i,j,k,n,tau)
      enddo
      enddo
      enddo
      end if
      endif
!
      do k=1,km
      do j=2,jmm
      do i=2,imm
      t(i,j,k,n,tau) = wka(i,j,k)
      enddo
      enddo
!
!      do j=2,jmm
!      t(1  ,j,k,n,tau) = t(imm,j,k,n,tau)
!      t(imt,j,k,n,tau) = t(2  ,j,k,n,tau)
!      t(1  ,j,k,n,taum)= t(imm,j,k,n,taum)
!      t(imt,j,k,n,taum)= t(2  ,j,k,n,taum)
!      enddo
      enddo
      
      call swap_array_real5d(t,imt,jmt,km,nt,2,west,east,north,south)
!
1000  continue
!
      if(.not.leapfrog_t) leapfrog_t = .true.
!
      return
      end
