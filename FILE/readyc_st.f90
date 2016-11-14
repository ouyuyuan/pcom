!
!     =================
      subroutine readyc_st(umask,tmask,ivn,pbt_st,du,dv,adv_u,adv_v,dub,  &
                    dvb,up,vp,cosu,rdxt,rdxu,rdyu,rdyt,sdxu,r1c,r1d,cv1,cv2,dz,  &
                    rdz,rdzw,rzu,pn,w,pax,pay,diffu,diffv,am,kappa_m,gravr,cdbot,  &
                    bcu,bcv,imt,jmt,km,imm,jmm,kmp1,west,east,north,  &
                    south,snbc,emp,t_stepu,energydiag,dke_bcf,dke_fri)
                        
!     =================
!     momentum advections, viscosities & atmpospheric pressure terms
!
      implicit none
      include 'pconst.h'
      include 'mpif.h'
!
      integer imt,jmt,km,imm,jmm,kmp1,i,j,k,snbc,energydiag
      real am(imt,jmt,km),kappa_m(imt,jmt,km),cdbot,gravr
      real ubar,vbar,abc,uvmag,t1,t2,t3
      real a(imt,jmt),b(imt,jmt),c(imt,jmt)
      real u(imt,jmt),v(imt,jmt),w(imt,jmt,kmp1)
      real ux(imt,jmt),vx(imt,jmt)
      real wua(imt,jmt),wub(imt,jmt)
      real wva(imt,jmt),wvb(imt,jmt)
      real tmask(imt,jmt,km),umask(imt,jmt,km)
      integer ivn(imt,jmt)
      real bcu(imt,jmt),bcv(imt,jmt),emp(imt,jmt)
      real pbt(imt,jmt,2),pbt_st(imt,jmt,4),spbt(imt,jmt)
      real du(imt,jmt,km),dv(imt,jmt,km),dub(imt,jmt),dvb(imt,jmt)
      real up(imt,jmt,km,2),vp(imt,jmt,km,2) 
      real z(km),dz(km),rdz(km),rdzw(km),pn(imt,jmt),rzu(imt,jmt)
      real rdxt(jmt),rdxu(jmt),rdyt(jmt),sdxu(jmt),rdyu(jmt)
      real pax(imt,jmt),pay(imt,jmt)
      real diffu(imt,jmt,km),diffv(imt,jmt,km)
      real cv1(jmt),cv2(jmt),r1c(jmt),r1d(jmt),cosu(jmt)
      
      real adv_u(imt,jmt,km,3),adv_v(imt,jmt,km,3)
      integer t_stepu
      
      integer west,east,north,south
      
      real dke_bcf(imt,jmt,2),dke_fri(imt,jmt,km,2)
!
!-----------------------------------------------------------------------
!     initialize time-averaged pbt (used in baroclinic Eq)
!-----------------------------------------------------------------------
!
      do j=1,jmt
         do i=1,imt
            pbt(i,j,tau)=pbt_st(i,j,1)
         end do
      end do
!
!-----------------------------------------------------------------------
!     calculate vertical mass advection (dz/dt * pbt)
!-----------------------------------------------------------------------
!
      do j=2,jmm
      do i=2,imm
      if(ivn(i,j).gt.0)then
      spbt(i,j) = p5*sqrt(pbt(i,j  ,tau) + pbt(i+1,j  ,tau) + &
                          pbt(i,j+1,tau) + pbt(i+1,j+1,tau))
      endif
      enddo
      enddo
      call swap_array_real2d(spbt,imt,jmt,west,east,north,south)
!
!     du,dv are used temporarily for mass advection
!
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
!
      call upwelling(du,dv,w,rdxt,rdyt,tmask,dz,pn,imt,jmt,km,imm,jmm,kmp1,  &
                     west,east,north,south,snbc,emp)
!
!     calculate vertical velocity on the surface of U cell
!
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
!
!-----------------------------------------------------------------------
!     calculate advections & horizontal viscosities
!-----------------------------------------------------------------------
!
      do 100 k=1,km
!
!     calculate horizontal advection velocities
!
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
!
!     U-advection
!
      a=0
      b=0
      c=0
      do j=2,jmm
      do i=2,imt
      ubar   = p5*(u(i,j)+u(i-1,j))
      a(i,j) = p5*(up(i,j,k,tau)+up(i-1,j,k,tau))*ubar
      b(i,j) = p5*(vp(i,j,k,tau)+vp(i-1,j,k,tau))*ubar
      c(i,j) = ubar
      enddo
      enddo
!
!      ux=0
!      vx=0
      do j=2,jmm
      do i=2,imm
      if(umask(i,j,k).gt.c0)then
      abc     = c(i+1,j)-c(i,j)
      ux(i,j) = rdxu(j)*((a(i+1,j)-a(i,j)) - p5*up(i,j,k,tau)*abc)
      vx(i,j) = rdxu(j)*((b(i+1,j)-b(i,j)) - p5*vp(i,j,k,tau)*abc)
      endif
      enddo
      enddo
!
!     V-advection
!
      do j=2,jmt
      do i=2,imm
      vbar   = p5*(v(i,j)*cosu(j)+v(i,j-1)*cosu(j-1))
      a(i,j) = p5*(up(i,j,k,tau)+up(i,j-1,k,tau))*vbar
      b(i,j) = p5*(vp(i,j,k,tau)+vp(i,j-1,k,tau))*vbar
      c(i,j) = vbar
      enddo
      enddo
!
      do j=2,jmm
      do i=2,imm
      if(umask(i,j,k).gt.c0)then
      abc     = c(i,j+1)-c(i,j)
      ux(i,j) = ux(i,j)+rdyu(j)*((a(i,j+1)-a(i,j))-p5*up(i,j,k,tau)*abc)
      vx(i,j) = vx(i,j)+rdyu(j)*((b(i,j+1)-b(i,j))-p5*vp(i,j,k,tau)*abc)
      endif
      enddo
      enddo
!
!
!     W-advection
!
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
!
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
!
      do j=2,jmm
      do i=2,imm
      if(umask(i,j,k).gt.c0)then
      abc     = w(i,j,k) - w(i,j,k+1)
      ux(i,j) = ux(i,j)+p5*rdz(k)*(wua(i,j)-wub(i,j)-up(i,j,k,tau)*abc)
      vx(i,j) = vx(i,j)+p5*rdz(k)*(wva(i,j)-wvb(i,j)-vp(i,j,k,tau)*abc)
      endif
      enddo
      enddo
!
!
!-----------------------------------------------------------------------
!     + advection + Pa + viscosities(t=taum)
!-----------------------------------------------------------------------
!
      do j=2,jmm
         do i=2,imm
            du(i,j,k) = - ux(i,j)
            dv(i,j,k) = - vx(i,j)
         end do
      end do

      if (t_stepu.eq.1) then 
        do j=2,jmm
           do i=2,imm
              adv_u(i,j,k,1)=du(i,j,k)
              adv_v(i,j,k,1)=dv(i,j,k)
              adv_u(i,j,k,2)=adv_u(i,j,k,1)
              adv_v(i,j,k,2)=adv_v(i,j,k,1)
              adv_u(i,j,k,3)=adv_u(i,j,k,2)
              adv_v(i,j,k,3)=adv_v(i,j,k,2)
           end do
        end do
      end if
      if (t_stepu.eq.2) then 
        do j=2,jmm
           do i=2,imm
              adv_u(i,j,k,2)=adv_u(i,j,k,1)
              adv_v(i,j,k,2)=adv_v(i,j,k,1)
              adv_u(i,j,k,1)=du(i,j,k)
              adv_v(i,j,k,1)=dv(i,j,k)
              adv_u(i,j,k,3)=adv_u(i,j,k,2)
              adv_v(i,j,k,3)=adv_v(i,j,k,2)
           end do
        end do
      end if
      if (t_stepu.gt.2) then 
        do j=2,jmm
           do i=2,imm
              adv_u(i,j,k,3)=adv_u(i,j,k,2)
              adv_v(i,j,k,3)=adv_v(i,j,k,2)
              adv_u(i,j,k,2)=adv_u(i,j,k,1)
              adv_v(i,j,k,2)=adv_v(i,j,k,1)
              adv_u(i,j,k,1)=du(i,j,k)
              adv_v(i,j,k,1)=dv(i,j,k)
           end do
        end do
      end if

      do j=2,jmm
         do i=2,imm
            du(i,j,k) = - ux(i,j) - pax(i,j)*spbt(i,j)
            dv(i,j,k) = - vx(i,j) - pay(i,j)*spbt(i,j)
         end do
      end do
!
!
!---------------------------------------------------------------------
!     calculate horizontal viscosities in tau time level
!---------------------------------------------------------------------
!
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
!
      do j=2,jmm
      do i=2,imm
      if(umask(i,j,k).gt.c0)then
      diffu(i,j,k) = am(i,j,k)*( (sdxu(j)*(ux(i+1,j)-ux(i,j))+r1c(j)*a(i,j+1)  &
                           -r1d(j)*a(i,j))/spbt(i,j)                    &
                   + cv1(j)*up(i,j,k,tau)                               &
                   - cv2(j)*(v(i+1,j)-v(i-1,j))*spbt(i,j) )
!
      diffv(i,j,k) = am(i,j,k)*( (sdxu(j)*(vx(i+1,j)-vx(i,j))+r1c(j)*b(i,j+1)  &
                           -r1d(j)*b(i,j))/spbt(i,j)                    &
                   + cv1(j)*vp(i,j,k,tau)                               &
                   + cv2(j)*(u(i+1,j)-u(i-1,j))*spbt(i,j) )
      endif
      enddo
      enddo
!
100   continue
!
!
!---------------------------------------------------------------------
!     NOTE: u(i,j) & v(i,j) = u&v at bottom (k=ivn(i,j))
!---------------------------------------------------------------------
!
      if (energydiag==1) then
      do j=2,jmm
         do i=2,imm
            k=1
            if (umask(i,j,k)==1) then
               dke_bcf(i,j,1)=rdz(k)*bcu(i,j)*gravr/spbt(i,j)*rrho_0
               dke_bcf(i,j,2)=rdz(k)*bcv(i,j)*gravr/spbt(i,j)*rrho_0
            end if
         end do
      end do
      end if
!
!---------------------------------------------------------------------
!     calculate vertical viscosities in tau time level
!---------------------------------------------------------------------
!
!     2-D working array

      do 200 k=1,km
!
      if(k.eq.1)then  
!       windstress bcu&bcv
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
!
!------BUG ??-------------
      do j=2,jmm
      do i=2,imm
      if(umask(i,j,k).gt.c0)then
      u(i,j) = up(i,j,k,tau)/spbt(i,j)
      v(i,j) = vp(i,j,k,tau)/spbt(i,j)
      else
      u(i,j) = c0
      v(i,j) = c0
      endif
      enddo
      enddo
!-------------------------
      do j=2,jmm
      do i=2,imm
      if(k.eq.ivn(i,j))then
!       bottom drag (assume rho=c1 for both PCOM and BCOM)
        uvmag    = sqrt(u(i,j)**2+v(i,j)**2)
        abc      = grav/spbt(i,j)*cdbot*uvmag
        wub(i,j) = u(i,j)*abc
        wvb(i,j) = v(i,j)*abc
      else if(k.lt.ivn(i,j))then
        a(i,j)=kappa_m(i,j,k)*gravr/spbt(i,j)**2
        wub(i,j) = a(i,j)*(up(i,j,k,tau)-up(i,j,k+1,tau))*rdzw(k)
        wvb(i,j) = a(i,j)*(vp(i,j,k,tau)-vp(i,j,k+1,tau))*rdzw(k)
      else
        wub(i,j) = c0
        wvb(i,j) = c0
      endif
      enddo
      enddo
!
!
      do j=2,jmm
      do i=2,imm
      diffu(i,j,k) = diffu(i,j,k) + rdz(k)*(wua(i,j)-wub(i,j))
      diffv(i,j,k) = diffv(i,j,k) + rdz(k)*(wva(i,j)-wvb(i,j))
      enddo
      enddo
      
!
!
!-----------------------------------------------------------------------
!     + diffusion (at tau time level) on euler forward time step
!-----------------------------------------------------------------------
!
      do j=2,jmm
      do i=2,imm
      du(i,j,k) = du(i,j,k) + diffu(i,j,k)
      dv(i,j,k) = dv(i,j,k) + diffv(i,j,k)
      enddo
      enddo
!
200   continue

      if (energydiag==1) then
      do j=2,jmm
         do i=2,imm
            do k=2,ivn(i,j)
               dke_fri(i,j,k,1)=diffu(i,j,k)
               dke_fri(i,j,k,2)=diffv(i,j,k)
            end do
            k=1
            if (umask(i,j,k)==1) then
               dke_fri(i,j,k,1)=diffu(i,j,k)-dke_bcf(i,j,1)
               dke_fri(i,j,k,2)=diffv(i,j,k)-dke_bcf(i,j,2)
            end if
         end do
      end do
      end if

!
!-----------------------------------------------------------------------
!     vertical integration of du&dv
!-----------------------------------------------------------------------
!
      call vinteg(du,dub,ivn,dz,rzu,imt,jmt,km)
      call vinteg(dv,dvb,ivn,dz,rzu,imt,jmt,km)
!
      return
      end subroutine readyc_st
