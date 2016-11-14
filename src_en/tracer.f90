!
!     =================
      subroutine tracer(t,w,pmtp,pmtm,umm,vmm,ump,vmp,onbc,oncc,sdxt,rdxt,rdyt,    &
                  rdz,rdzw,r1a,r1b,tmask,itn,dz,pn,bct,bcs,ddd,gamma_t,gamma_s,  &
                  du,dv,ah,am,sma_c,acfl,kappa_h,gravr,leapfrog_t,c2dtts,dtts,  &
                  imt,jmt,km,nt,imm,jmm,kmp1,west,east,north,south,myid,asselin_t,  &
                  gm90,implicitvmix,snbc,tnbc,aft1,aft2,emp,dz0,z0,rdy,cosu,K1,   &
                  K2,K3,adv_vetiso,adv_vntiso,adv_vbtiso,rhoi,e,slmxr,ahisop,    &
                  athkdf,fzisop,kref,rdz0,unesco,umn,vmn,wmn,pmn,rho,phib,gr_mass, &
                  bottom_h,ncpux,ncpuy,mat_myid,z,nss,rdxu,lat,bcp,  &
                  energydiag,dpo_adv,dpo_hdif,dpo_vdif,dpo_bc,din_adv,din_hdif,  &
                  din_vdif,din_bc,dpo_imp,din_imp,dpo_ice,din_ice,dpo_ast,din_ast, &
                  dsa_ast,p,pbar,mass_up)
!     =================
!     compute potential temperature and salinity at tau+1 time level
!
      implicit none
      include 'pconst.h'
!      include 'isopyc.h'
      include 'mpif.h'
      integer imt,jmt,km,nt,imm,jmm,kmp1,i,j,k,n,kk
      integer asselin_t,implicitvmix,gm90,snbc,tnbc,unesco,bous,energydiag
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
      real bct(imt,jmt),bcs(imt,jmt),ddd(imt,jmt)
      real du(imt,jmt,km),dv(imt,jmt,km)
      logical leapfrog_t
      
      real slmxr,fzisop(km),kref(km),rdz0(km),dz0(km),z0(km),cosu(jmt)
      real K1(imt,km,jmt,3:3),K2(imt,km,jmt,3:3),K3(imt,km,jmt,1:3)
      real adv_vetiso(imt,km,jmt),adv_vntiso(imt,km,jmt),adv_vbtiso(imt,0:km,jmt)
      real rhoi(imt,km,jmt,1:3),e(imt,kmp1,jmt,3),ahisop(imt,jmt,km),athkdf(imt,jmt,km)
      
      real gr_mass(imt,jmt,km),gr_h(imt,jmt,km),bottom_h(imt,jmt)
      
      integer ncpux,ncpuy,nss
      integer mat_myid(ncpux+2,ncpuy)
      real decibar,t0,s0,p0,unrdens,gsw_internal_energy
      real am(imt,jmt,km)
      
      real umn(imt,jmt,km),vmn(imt,jmt,km),wmn(imt,jmt,km),pmn(imt,jmt),rho(imt,jmt,km)
      real psi(imt,jmt),phib(imt,jmt)
      integer west,east,north,south,myid
      
      real dt_am(imt,jmt,km),ds_am(imt,jmt,km),rdxu(jmt),lat(jmt),acfl(jmt),cofe,sma_c
      
      real dpo_adv (imt,jmt,km),din_adv (imt,jmt,km),ts_adv (imt,jmt,km,2)
      real dpo_hdif(imt,jmt,km),din_hdif(imt,jmt,km),ts_hdif(imt,jmt,km,2)
      real dpo_vdif(imt,jmt,km),din_vdif(imt,jmt,km),ts_vdif(imt,jmt,km,2)
      real din_bc  (imt,jmt)   ,dpo_bc  (imt,jmt)   ,ts_bc  (imt,jmt,2)
      real dpo_imp (imt,jmt,km),din_imp (imt,jmt,km),ts_imp (imt,jmt,km,2)
      real dpo_ice (imt,jmt,km),din_ice (imt,jmt,km),ts_ice (imt,jmt,km,2)
      real dpo_ast (imt,jmt,km),din_ast (imt,jmt,km),ts_ast1(imt,jmt,km,2)
      real ts_ast2(imt,jmt,km,2),dsa_ast(imt,jmt,km)
      real mass_up(imt,jmt,km),ts_ini(imt,jmt,km,2)
      real bcp(imt,jmt),z(km),p(imt,jmt,km)
      real rho_temp,dz_temp,ini_temp,wka_temp(imt,jmt,km)
!
      if ((gm90==1).or.(implicitvmix==1)) then
      aidif = p5
      else
      aidif = c1
      end if
      
      decibar = 1.0e-5
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
!---------------------------------------------------------------------
!     calculate instantaneous umn vmn wmn pmn psi gr_mass po_en for output
!---------------------------------------------------------------------
      do k=1,km
        do j=1,jmt
          do i=1,imt
            umn(i,j,k)=wka(i,j,k)/pbar(i,j)
            vmn(i,j,k)=wkb(i,j,k)/(pbar(i,j)*cosu(j))
            wmn(i,j,k)=w(i,j,k)*rho(i,j,k)/(pbar(i,j)*grav)
          end do
        end do
      end do
      
      do j=1,jmm
         do i=1,imm
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
          gr_mass(i,j,k) = dz(k)*pbar(i,j)/grav
        end do
      endif
      enddo
      enddo
      
      do j=1,jmt
        do i=1,imt
          pmn(i,j)=psi(i,j)
        end do
      end do
      
      if (energydiag.eq.1) then
         mass_up=c0
         do j=1,jmt
            do i=1,imt
               do k=1,itn(i,j)
                  do kk=1,k-1
                     mass_up(i,j,k)=mass_up(i,j,k)+gr_mass(i,j,kk)
                  end do
                  mass_up(i,j,k)=mass_up(i,j,k)+p5*gr_mass(i,j,k)
               end do
            end do
         end do
      end if

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
      wkb(i,j,k) = kappa_h(i,j,k) + ahisop(i,j,k)*k3(i,k,j,3)
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

      if (energydiag.eq.1) then
         do j=2,jmm
            do i=2,imm
               ts_adv(i,j,k,n)=ts_ini(i,j,k,n)+wka(i,j,k)/pbar(i,j)*tmask(i,j,k)*dtts2
            end do
         end do
      end if

!
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
      if (gm90==1) then
      
      if (energydiag.eq.1) then
      do j=2,jmm
      do i=2,imm
      do k=1,itn(i,j)
      wka_temp(i,j,k) = wka(i,j,k)
      end do
      end do
      end do
      end if
      
      call isoflux(wka,n,imt,jmt,imm,jmm,nt,km,kmp1,itn,gravr,rdz0,dz0,rdz,     &
                   rdzw,rdxt,rdyt,rdy,cosu,tmask,t,ah,K1,K2,K3,adv_vetiso,  &
                   adv_vntiso,adv_vbtiso,fx,fy)
      
      if (energydiag.eq.1) then
      do j=2,jmm
      do i=2,imm
      do k=1,itn(i,j)
      wka_temp(i,j,k) = wka(i,j,k)-wka_temp(i,j,k)
      end do
      end do
      end do
      end if
      
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

      if (energydiag.eq.1) then
      do j=2,jmm
      do i=2,imm
      wka_temp(i,j,k) = ah(i,j,k)*(sdxt(j)*(a(i+1,j)-a(i,j))+  &
                                r1a(j)*b(i,j+1)-r1b(j)*b(i,j))
      enddo
      enddo
      end if

200   continue

      end if
      
      if (energydiag.eq.1) then
      do j=2,jmm
         do i=2,imm
            do k=1,itn(i,j)
               ts_hdif(i,j,k,n)=ts_ini(i,j,k,n)+wka_temp(i,j,k)/pbar(i,j)*dtts2
            end do
         end do
      end do
      end if
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
      
      if (energydiag.eq.1) then
      do j=2,jmm
      do i=2,imm
      wka_temp(i,j,k) = rdz(k)*(flxa(i,j)-flxb(i,j))*aidif
      enddo
      enddo
      end if
      
300   continue

      if (energydiag.eq.1) then
      do j=2,jmm
         do i=2,imm
            do k=1,itn(i,j)
               ts_vdif(i,j,k,n)=ts_ini(i,j,k,n)+wka_temp(i,j,k)/pbar(i,j)*dtts2
            end do
         end do
      end do
      end if
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

      if (energydiag.eq.1) then
      do j=2,jmm
         do i=2,imm
            ts_bc(i,j,n)=ts_ini(i,j,1,n)+stf(i,j)*rdz(1)*aidif/pbar(i,j)*tmask(i,j,1)*dtts2
         end do
      end do
      end if

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

      if (energydiag.eq.1) then
      do j=2,jmm
      do i=2,imm
      do k=1,itn(i,j)
      wka_temp(i,j,k) = wka(i,j,k)
      end do
      end do
      end do
      end if

      if ((gm90==1).and.(implicitvmix==1)) then
      call invtri(wka,stf,wkb,aidif,dtts2,pbar,gravr,itn,tmask,rdz,rdzw,   &
                   imt,jmt,km,imm,jmm)
      end if
      if ((gm90==0).and.(implicitvmix==1)) then
      call invtri(wka,stf,kappa_h,aidif,dtts2,pbar,gravr,itn,tmask,rdz,rdzw,   &
                   imt,jmt,km,imm,jmm)
      end if
      
      if (energydiag.eq.1) then
      do j=2,jmm
         do i=2,imm
            do k=1,itn(i,j)
               ts_imp(i,j,k,n)=ts_ini(i,j,k,n)+wka(i,j,k)-wka_temp(i,j,k)
            end do
         end do
      end do
      end if
!
!-----------------------------------------------------------------------
!     set T > -1.5 temporarily since no seaice model
!-----------------------------------------------------------------------
      if (energydiag.eq.1) then
      do j=2,jmm
      do i=2,imm
      do k=1,itn(i,j)
      wka_temp(i,j,k) = wka(i,j,k)
      end do
      end do
      end do
      end if

      if(n.eq.1)then
      do k=1,km
      do j=2,jmm
      do i=2,imm
      wka(i,j,k) = dmax1(tbice, wka(i,j,k))
      enddo
      enddo
      enddo
      endif

      if (energydiag.eq.1) then
      do j=2,jmm
         do i=2,imm
            do k=1,itn(i,j)
               ts_ice(i,j,k,n)=ts_ini(i,j,k,n)+wka(i,j,k)-wka_temp(i,j,k)
            end do
         end do
      end do
      end if
!
!
!-----------------------------------------------------------------------
!     T&S=wka  & Filter & periodic b.c
!-----------------------------------------------------------------------
!
      if(leapfrog_t) then
      
      if (energydiag.eq.1) then
      do j=2,jmm
      do i=2,imm
      do k=1,itn(i,j)
      ts_ast1(i,j,k,n) = t(i,j,k,n,tau)*aft2+aft1*(t(i,j,k,n,taum)+wka(i,j,k))
      ts_ast2(i,j,k,n) = t(i,j,k,n,tau)
      enddo
      enddo
      enddo 
      end if
      
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
      else
      
      if (energydiag.eq.1) then
      do j=2,jmm
      do i=2,imm
      do k=1,itn(i,j)
      ts_ast1(i,j,k,n) = wka(i,j,k)
      ts_ast2(i,j,k,n) = wka(i,j,k)
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

      end do
!
1000  continue

      call swap_array_real5d(t,imt,jmt,km,nt,2,west,east,north,south)

!
      if(.not.leapfrog_t) leapfrog_t = .true.

      if (energydiag.eq.1) then
         do j=2,jmm
            do i=2,imm
               do k=1,itn(i,j)
                  t0=ts_ini(i,j,k,1)
                  s0=ts_ini(i,j,k,2)
                  p0=p(i,j,k)
                  ini_temp=gsw_internal_energy(s0,t0,p0)
                  
                  t0=ts_adv(i,j,k,1)
                  s0=ts_adv(i,j,k,2)
                  rho_temp=unrdens(t0,s0,p0)
                  dz_temp=dz(k)*pbar(i,j)*rho_temp/grav-gr_h(i,j,k)
                  dpo_adv(i,j,k)=dz_temp*mass_up(i,j,k)*grav
                  din_adv(i,j,k)=gsw_internal_energy(s0,t0,p0)-ini_temp
                  
                  t0=ts_hdif(i,j,k,1)
                  s0=ts_hdif(i,j,k,2)
                  rho_temp=unrdens(t0,s0,p0)
                  dz_temp=dz(k)*pbar(i,j)*rho_temp/grav-gr_h(i,j,k)
                  dpo_hdif(i,j,k)=dz_temp*mass_up(i,j,k)*grav
                  din_hdif(i,j,k)=gsw_internal_energy(s0,t0,p0)-ini_temp
                  
                  t0=ts_vdif(i,j,k,1)
                  s0=ts_vdif(i,j,k,2)
                  rho_temp=unrdens(t0,s0,p0)
                  dz_temp=dz(k)*pbar(i,j)*rho_temp/grav-gr_h(i,j,k)
                  dpo_vdif(i,j,k)=dz_temp*mass_up(i,j,k)*grav
                  din_vdif(i,j,k)=gsw_internal_energy(s0,t0,p0)-ini_temp
                  
                  t0=ts_imp(i,j,k,1)
                  s0=ts_imp(i,j,k,2)
                  rho_temp=unrdens(t0,s0,p0)
                  dz_temp=dz(k)*pbar(i,j)*rho_temp/grav-gr_h(i,j,k)
                  dpo_imp(i,j,k)=dz_temp*mass_up(i,j,k)*grav
                  din_imp(i,j,k)=gsw_internal_energy(s0,t0,p0)-ini_temp
                  
                  t0=ts_ice(i,j,k,1)
                  s0=ts_ice(i,j,k,2)
                  rho_temp=unrdens(t0,s0,p0)
                  dz_temp=dz(k)*pbar(i,j)*rho_temp/grav-gr_h(i,j,k)
                  dpo_ice(i,j,k)=dz_temp*mass_up(i,j,k)*grav
                  din_ice(i,j,k)=gsw_internal_energy(s0,t0,p0)-ini_temp
                  
                  t0=ts_ast1(i,j,k,1)
                  s0=ts_ast1(i,j,k,2)
                  rho_temp=unrdens(t0,s0,p0)
                  dz_temp=dz(k)*pbar(i,j)*rho_temp/grav
                  din_ast(i,j,k)=gsw_internal_energy(s0,t0,p0)
                  t0=ts_ast2(i,j,k,1)
                  s0=ts_ast2(i,j,k,2)
                  rho_temp=unrdens(t0,s0,p0)
                  dz_temp=dz_temp-dz(k)*pbar(i,j)*rho_temp/grav
                  dpo_ast(i,j,k)=dz_temp*mass_up(i,j,k)*grav
                  din_ast(i,j,k)=din_ast(i,j,k)-gsw_internal_energy(s0,t0,p0)
                  dsa_ast(i,j,k)=ts_ast1(i,j,k,2)-ts_ast2(i,j,k,2)
                  
                  if (k==1) then
                  t0=ts_bc(i,j,1)
                  s0=ts_bc(i,j,2)
                  rho_temp=unrdens(t0,s0,p0)
                  dz_temp=dz(k)*pbar(i,j)*rho_temp/grav-gr_h(i,j,k)
                  dpo_bc(i,j)=p5*dz_temp*gr_mass(i,j,k)*grav
                  din_bc(i,j)=gsw_internal_energy(s0,t0,p0)-ini_temp
                  end if
               end do
            end do
         end do
      end if

      return
      end
