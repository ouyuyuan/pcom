!
!     =================
      subroutine bclinc_st(phib,spbt,pbt,pbt_st,rhodp,rho,rdxt,rdy,ump,vmp,upb,  &
                     vpb,up,vp,du,dv,adv_u,adv_v,diffu,diffv,epla,eplb,epea,epeb,  &
                     pax,pay,ff,umask,ivn,itn,z,dz,rzu,onbb,dtuv,  &
                     cosu,imt,jmt,km,imm,jmm,west,east,north,south,  &
                     dke_pre,dke_adv,dke_ape,dke_bar,dke_bcf,dke_fri,dke_cor,myid)
!     =================
!
!     compute up & vp at tau+1 time level
!
      implicit none
      include 'pconst.h'
      include 'mpif.h'
!
      integer imt,jmt,km,imm,jmm,i,j,k
      integer ivn(imt,jmt),itn(imt,jmt)
      real t1,t2
      real dtuv,onbb,z(km),dz(km),rzu(imt,jmt),cosu(jmt)
      real a(imt,jmt),b(imt,jmt),c(imt,jmt),px(imt,jmt,km),py(imt,jmt,km)
      real etax(imt,jmt),etay(imt,jmt),pbar(imt,jmt)
      real pbt_st(imt,jmt,4),pbt(imt,jmt,2),spbt(imt,jmt),phib(imt,jmt)
      real rhodp(imt,jmt,km),rho(imt,jmt,km)
      real rdxt(jmt),rdy,ff(jmt),umask(imt,jmt,km)
      real ump(imt,jmt,km),vmp(imt,jmt,km),upb(imt,jmt,2),vpb(imt,jmt,2)
      real up(imt,jmt,km,2),vp(imt,jmt,km,2),du(imt,jmt,km),dv(imt,jmt,km)
      real epea(jmt),epeb(jmt),epla(jmt),eplb(jmt)
      
      real adv_u(imt,jmt,km,3),adv_v(imt,jmt,km,3),advu(imt,jmt,km),advv(imt,jmt,km)
      real diffu(imt,jmt,km),diffv(imt,jmt,km)
      real c_1,c_2,c_3
      real pax(imt,jmt),pay(imt,jmt)
      
      integer myid,west,east,north,south
      
      real dke_pre(imt,jmt,km,2),dke_adv(imt,jmt,km,2),dke_bar(imt,jmt,km,2)
      real dke_bcf(imt,jmt,2),dke_ape(imt,jmt,km,2),dke_fri(imt,jmt,km,2)
      real uv_ini(imt,jmt,km,2),dke_cor(imt,jmt,km,2)
      real ka,kb,kc,kd,ke,kf,kg,kh,temp

      c_1=23/12
      c_2=-16/12
      c_3=5/12
!
!-----------------------------------------------------------------------
!     calculate time-averaged pbt & sqrt(pbt)
!-----------------------------------------------------------------------

      do j=1,jmt
         do i=1,imt
            pbar(i,j) = pbt_st(i,j,2)
         end do
      end do
!
      do j=2,jmm
      do i=2,imm
      if(ivn(i,j).gt.0)then
      spbt(i,j) = p5*sqrt(pbar(i,j  ) + pbar(i+1,j  ) +  &
                          pbar(i,j+1) + pbar(i+1,j+1))
      endif
      enddo
      enddo
      call swap_array_real2d(spbt,imt,jmt,west,east,north,south)
!
!
      px=0
      py=0
      
!
!-----------------------------------------------------------------------
!     calculate pressure gradients  * spbt
!-----------------------------------------------------------------------
!
!     a  = geopotential
!     b  = pressure
!     px = (a)x + (rho)xbar*(b)x
!     py = (a)y + (rho)ybar*(b)y
!
      do k=1,km
      do j=1,jmt
      do i=1,imt
      a(i,j) = phib(i,j) + pbar(i,j)*rhodp(i,j,k)
      b(i,j) = pbar(i,j)*z(k)
      enddo
      enddo
!
      do j=2,jmm
      do i=2,imm
      px(i,j,k) = ( (rho(i,j,k)+rho(i+1,j,k))*p5*(b(i+1,j)-b(i,j)) + &
                  (a(i+1,j)-a(i,j)) )*rdxt(j)
      py(i,j,k) = ( (rho(i,j,k)+rho(i,j+1,k))*p5*(b(i,j+1)-b(i,j)) + &
                  (a(i,j+1)-a(i,j)) )*rdy
      enddo
      enddo
      end do
      
      call swap_ns_real3d(px,imt,jmt,km,north,south)
      call swap_ew_real3d(py,imt,jmt,km,west,east)

      do 100 k=1,km
      do j=2,jmm
      do i=2,imm
      etax(i,j) = (px(i,j,k)+px(i,j+1,k))*p5*spbt(i,j)
      etay(i,j) = (py(i,j,k)+py(i+1,j,k))*p5*spbt(i,j)
      enddo
      enddo
!
!
!---------------------------------------------------------------------
!     advection + diffusion + pressure gradients + coriolis
!---------------------------------------------------------------------
      do j=2,jmm
         do i=2,imm
            advu(i,j,k)=c_1*adv_u(i,j,k,1)+c_2*adv_u(i,j,k,2)+c_3*adv_u(i,j,k,3)
            advv(i,j,k)=c_1*adv_v(i,j,k,1)+c_2*adv_v(i,j,k,2)+c_3*adv_v(i,j,k,3)
         end do
      end do


      do j=2,jmm
      do i=2,imm
      a(i,j) = advu(i,j,k) + diffu(i,j,k) + ff(j)*vp(i,j,k,tau) - etax(i,j) - pax(i,j)*spbt(i,j)
      b(i,j) = advv(i,j,k) + diffv(i,j,k) - ff(j)*up(i,j,k,tau) - etay(i,j) - pay(i,j)*spbt(i,j)
      enddo
      enddo
!
!
!-----------------------------------------------------------------------
!     coriolis adjustment
!-----------------------------------------------------------------------
      do j=2,jmm
         do i=2,imm
            du(i,j,k) = (epea(j)*a(i,j) + epeb(j)*b(i,j))*umask(i,j,k)
            dv(i,j,k) = (epea(j)*b(i,j) - epeb(j)*a(i,j))*umask(i,j,k)
         end do
      end do
      
100   continue 
!
!-----------------------------------------------------------------------
!     compute up & vp at tau+1 time level
!-----------------------------------------------------------------------
!
      do k=1,km
      do j=2,jmm
      do i=2,imm
      if(umask(i,j,k).gt.c0)then
      up(i,j,k,tau) = up(i,j,k,taum) + du(i,j,k)*dtuv
      vp(i,j,k,tau) = vp(i,j,k,taum) + dv(i,j,k)*dtuv
      endif
      enddo
      enddo
      enddo
!
!@@@  interaction between barotropic and baroclinic modes
!
      call vinteg(up(1,1,1,tau),a,ivn,dz,rzu,imt,jmt,km)
      call vinteg(vp(1,1,1,tau),b,ivn,dz,rzu,imt,jmt,km)
!
      do k=1,km
      do j=2,jmm
      do i=2,imm
      if(umask(i,j,k).gt.c0)then
      up(i,j,k,tau) = up(i,j,k,tau) - a(i,j) + upb(i,j,tau)
      vp(i,j,k,tau) = vp(i,j,k,tau) - b(i,j) + vpb(i,j,tau)
      up(i,j,k,taum) = up(i,j,k,tau)
      vp(i,j,k,taum) = vp(i,j,k,tau)
      endif
      enddo
      enddo
      enddo
      
      call swap_array_real4d(up,imt,jmt,km,2,west,east,north,south)
      call swap_array_real4d(vp,imt,jmt,km,2,west,east,north,south)
!
!
!-----------------------------------------------------------------------
!     constract time-averaged mass advections
!-----------------------------------------------------------------------
      do j=1,jmt
         do i=1,imt
            pbt(i,j,tau)=pbt_st(i,j,3)
         end do
      end do
      
      
      do j=2,jmm
      do i=2,imm
      spbt(i,j) = p5*sqrt(pbt(i,j  ,tau) + pbt(i+1,j  ,tau) + &
                          pbt(i,j+1,tau) + pbt(i+1,j+1,tau))
      enddo
      enddo
      call swap_array_real2d(spbt,imt,jmt,west,east,north,south)
!
      do k=1,km
      do j=2,jmm
      do i=2,imm
      if(umask(i,j,k).gt.c0)then
      ump(i,j,k) = up(i,j,k,tau)*spbt(i,j)
      vmp(i,j,k) = vp(i,j,k,tau)*spbt(i,j)*cosu(j)
      endif
      enddo
      enddo
      enddo
      call swap_array_real3d(ump,imt,jmt,km,west,east,north,south)
      call swap_array_real3d(vmp,imt,jmt,km,west,east,north,south)

      return
      end subroutine bclinc_st
