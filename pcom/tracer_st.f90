!
!     =================
      subroutine tracer_st(imt,jmt,km,nt,imm,jmm,kmp1,dtts,t,w,ump,vmp,sdxt,rdxt,  &
                    rdxu,rdyt,rdy,rdz,rdzw,r1a,r1b,lat,cosu,tmask,itn,dz,dz0,z0,z,  &
                    pn,bct,bcs,ddd,emp,bcp,gamma_t,gamma_s,du,dv,ah,am,sma_c,acfl,kappa_h,  &
                    kappa_m,gravr,west,east,north,south,myid,implicitvmix,snbc,  &
                    tnbc,K1,K2,K3,adv_vetiso,adv_vntiso,adv_vbtiso,rhoi,e,slmxr,ahisop,  &
                    athkdf,fzisop,kref,rdz0,umn,vmn,wmn,pmn,pbt_st,rho,phib, &
                    gr_mass,dpo_adv,dpo_hdif,dpo_vdif,dpo_bc,dpo_imp,dpo_ice,  &
                    dpo_bar,dpo_bcp,din_adv,din_hdif,din_vdif,din_bc,din_imp,din_ice,  &
                    din_bar,p,pbar,mass_up,total_in,gr_mass2)
!     =================
!     compute potential temperature and salinity at tau+1 time level
!
      implicit none
      include 'pconst.h'
!      include 'isopyc.h'
      include 'mpif.h'
      integer imt,jmt,km,nt,imm,jmm,kmp1,i,j,k,n,kk
      integer implicitvmix,snbc,tnbc,bous
      real abc,aidif,dtts2
      real decibar,t0,s0,p0,unrdens,gsw_internal_energy
      real a(imt,jmt),b(imt,jmt),c(imt,jmt),wa(imt,jmt),wb(imt,jmt)
      real flxa(imt,jmt),flxb(imt,jmt),dt(imt,jmt),fx(imt,jmt),fy(imt,jmt)
      real wka(imt,jmt,km),wkb(imt,jmt,km),pbar(imt,jmt),stf(imt,jmt)
      
      real dtts,gravr,gamma_t,gamma_s,rdy
      real ah(imt,jmt,km),kappa_h(imt,jmt,km)
      real am(imt,jmt,km),kappa_m(imt,jmt,km)
      
      real t(imt,jmt,km,nt,2),w(imt,jmt,kmp1)
      real emp(imt,jmt)
      real pbt_st(imt,jmt,4)
      real ump(imt,jmt,km),vmp(imt,jmt,km)
      real rdxt(jmt),rdyt(jmt),r1a(jmt),r1b(jmt)
      real tmask(imt,jmt,km),rdz(km),rdzw(km),dz(km),pn(imt,jmt),sdxt(jmt)
      integer itn(imt,jmt)
      real bct(imt,jmt),bcs(imt,jmt),ddd(imt,jmt)
      real du(imt,jmt,km),dv(imt,jmt,km)
      
      real slmxr,fzisop(km),kref(km),rdz0(km),dz0(km),z0(km),cosu(jmt)
      real K1(imt,km,jmt,3:3),K2(imt,km,jmt,3:3),K3(imt,km,jmt,1:3)
      real adv_vetiso(imt,km,jmt),adv_vntiso(imt,km,jmt),adv_vbtiso(imt,0:km,jmt)
      real rhoi(imt,km,jmt,1:3),e(imt,kmp1,jmt,3),ahisop(imt,jmt,km),athkdf(imt,jmt,km)
      
      real gr_mass(imt,jmt,km),gr_h(imt,jmt,km),gr_mass2(imt,jmt,km)
      real umn(imt,jmt,km),vmn(imt,jmt,km),wmn(imt,jmt,km),pmn(imt,jmt),rho(imt,jmt,km)
      real psi(imt,jmt),phib(imt,jmt)
      integer west,east,north,south,myid
      
      real dt_am(imt,jmt,km),ds_am(imt,jmt,km),rdxu(jmt),lat(jmt),acfl(jmt),cofe,sma_c
      real ri,rhop_up,rhop_down,duvdz,drhopdz
      
      real dpo_adv (imt,jmt,km),din_adv (imt,jmt,km),ts_adv (imt,jmt,km,2)
      real dpo_hdif(imt,jmt,km),din_hdif(imt,jmt,km),ts_hdif(imt,jmt,km,2)
      real dpo_vdif(imt,jmt,km),din_vdif(imt,jmt,km),ts_vdif(imt,jmt,km,2)
      real din_bc  (imt,jmt)   ,dpo_bc  (imt,jmt)   ,ts_bc  (imt,jmt,2)
      real dpo_imp (imt,jmt,km),din_imp (imt,jmt,km),ts_imp (imt,jmt,km,2)
      real dpo_ice (imt,jmt,km),din_ice (imt,jmt,km),ts_ice (imt,jmt,km,2)
      real dpo_bar (imt,jmt,km),din_bar (imt,jmt,km)
      real dpo_bcp (imt,jmt,km)
      real mass_up(imt,jmt,km),ts_ini(imt,jmt,km,2)
      real bcp(imt,jmt),z(km),p(imt,jmt,km),total_in(imt,jmt,km)
      real rho_temp,dz_temp,ini_temp,wka_temp(imt,jmt,km)
!
      aidif = p5
      
      decibar = 1.0e-5
!
!-----------------------------------------------------------------------
!     calculate time-averaged pbt and mass advections
!-----------------------------------------------------------------------
      do k=1,km
         do j=1,jmt
            do i=1,imt
               wka(i,j,k) = ump(i,j,k)
               wkb(i,j,k) = vmp(i,j,k)
            end do
         end do
      end do
      
!
!---------------------------------------------------------------------
!     calculate vertical mass advection
!---------------------------------------------------------------------
      call upwelling(wka,wkb,w,rdxt,rdyt,tmask,dz,pn,imt,jmt,km,imm,jmm,kmp1,  &
                     west,east,north,south,snbc,emp)
!
!---------------------------------------------------------------------
!     calculate instantaneous umn vmn wmn pmn psi gr_mass po_en for output
!---------------------------------------------------------------------
      do j=1,jmt
         do i=1,imt
            pbar(i,j) = pbt_st(i,j,3)
         end do
      end do

      do k=1,km
        do j=1,jmt
          do i=1,imt
            umn(i,j,k)=wka(i,j,k)/pbar(i,j)
            vmn(i,j,k)=wkb(i,j,k)/(pbar(i,j)*cosu(j))
            wmn(i,j,k)=w(i,j,k)*rho(i,j,k)/(pbar(i,j)*grav)
          end do
        end do
      end do

      do j=1,jmt
         do i=1,imt
            pbar(i,j) = pbt_st(i,j,2)
         end do
      end do

      do j=1,jmt
         do i=1,imt
            do k=1,itn(i,j)
               do n=1,nt
                  ts_ini(i,j,k,n)=t(i,j,k,n,taum)
               end do
               p(i,j,k)=(pbar(i,j)*z(k) + bcp(i,j))*decibar
               t0=ts_ini(i,j,k,1)
               s0=ts_ini(i,j,k,2)
               p0=p(i,j,k)
               rho(i,j,k)=unrdens(t0,s0,p0)
            end do
         end do
      end do
      
      do j=1,jmt
      do i=1,imt
      if(itn(i,j).gt.0)then
        psi(i,j) = phib(i,j)/grav
        do k=1,itn(i,j)
          gr_h(i,j,k) = dz(k)*pbar(i,j)*rho(i,j,k)/grav
          psi(i,j) = psi(i,j)+gr_h(i,j,k)
          gr_mass(i,j,k) = dz(k)/grav
          gr_mass2(i,j,k)= dz(k)*pbar(i,j)/grav
        end do
      endif
      enddo
      enddo
      
      do j=1,jmt
        do i=1,imt
          pmn(i,j)=psi(i,j)
        end do
      end do
      
!---------------------------------------------------------------------
!     calculate new am by Smagorinsky Scheme
!---------------------------------------------------------------------
      do k=1,km
        do j=2,jmm
          do i=2,imm
            dt_am(i,j,k)=2*rdxu(j)*(umn(i+1,j,k)-umn(i-1,j,k))-  &
                      2*rdy*(vmn(i,j+1,k)-vmn(i,j-1,k))+vmn(i,j,k)*tan(lat(j)*torad)/radius
            ds_am(i,j,k)=2*rdxu(j)*(vmn(i+1,j,k)-vmn(i-1,j,k))+  &
                      2*rdy*(umn(i,j+1,k)-umn(i,j-1,k))-umn(i,j,k)*tan(lat(j)*torad)/radius
            cofe=(p5*sma_c*(c1/rdy+c1/rdxu(j))/pi)**2
            am(i,j,k)=cofe*sqrt(dt_am(i,j,k)**2+ds_am(i,j,k)**2)
            if (am(i,j,k).gt.acfl(j)) then
               am(i,j,k)=acfl(j)
            end if
          end do
        end do
      end do
      call swap_array_real3d(am,imt,jmt,km,west,east,north,south)

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
      call isopyc(imt,jmt,km,kmp1,imm,jmm,nt,itn,tmask,kref,fzisop,rdz0,dz0,z0,  &
                  rdxt,rdyt,rdy,rdzw,dz,pn,cosu,t,slmxr,athkdf,K1,K2,K3,  &
                  adv_vetiso,adv_vntiso,adv_vbtiso,rhoi,e,fx,fy,west,east,north, &
                  south)
!
      do k=1,km
      do j=1,jmt
      do i=1,imt
      wkb(i,j,k) = kappa_h(i,j,k) + ahisop(i,j,k)*k3(i,k,j,3)
      enddo
      enddo
      enddo
!
      dtts2 = dtts
!
!---------------------------------------------------------------------
!     solve for one tracer at a time
!---------------------------------------------------------------------
!     n = 1 => temperature
!     n = 2 => salinity
!
      do 1000 n=1,nt
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
      
100   continue
!
!-----------------------------------------------------------------------
!     diffusions
!-----------------------------------------------------------------------
!
!     compute the isopycnal/dipycnal mixing
!     xz and yz isopycnal diffusive flux are solved explicitly;
!     while zz component will be solved implicitly.
!
      call isoflux(wka,n,imt,jmt,imm,jmm,nt,km,kmp1,itn,gravr,rdz0,dz0,rdz,     &
                   rdzw,rdxt,rdyt,rdy,cosu,tmask,t,ah,K1,K2,K3,adv_vetiso,  &
                   adv_vntiso,adv_vbtiso,fx,fy)

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
      
      do j=2,jmm
      do i=2,imm
      wka(i,j,k) = wka(i,j,k) + rdz(k)*(flxa(i,j)-flxb(i,j))*aidif
      enddo
      enddo
      
300   continue

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
      if (tnbc==1) then
      do j=2,jmm
      do i=2,imm
       stf(i,j)   = ddd(i,j)
       wka(i,j,1) = wka(i,j,1) + stf(i,j)*rdz(1)*aidif
      enddo
      enddo
      else
      do j=2,jmm
      do i=2,imm
       stf(i,j)   = gamma_t*(bct(i,j)-t(i,j,1,n,tau))
       wka(i,j,1) = wka(i,j,1) + stf(i,j)*rdz(1)*aidif
      enddo
      enddo
      end if
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

      if implicitvmix==1) then
      call invtri(wka,stf,wkb,aidif,dtts2,pbar,gravr,itn,tmask,rdz,rdzw,   &
                   imt,jmt,km,imm,jmm)
      end if
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
!-----------------------------------------------------------------------
!     T&S=wka  & Filter & periodic b.c
!-----------------------------------------------------------------------
!
      do k=1,km
      do j=2,jmm
      do i=2,imm
      t(i,j,k,n,tau) = wka(i,j,k)
      enddo
      enddo
      enddo
      
!
1000  continue

      do j=1,jmt
         do i=1,imt
            pbt_st(i,j,1)=pbt_st(i,j,3)
            pbt_st(i,j,4)=pbt_st(i,j,2)
         end do
      end do
      
!      call swap_array_real5d(t,imt,jmt,km,nt,2,west,east,north,south)
!
      return
      end subroutine tracer_st
