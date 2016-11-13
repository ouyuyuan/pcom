!
!     =================
      subroutine bclinc(pmup,pmum,phib,spbt,pbt,rhodp,rho,rdxt,rdy,ump,vmp,upb,    &
                        vpb,up,vp,du,dv,epla,eplb,epea,epeb,ff,umask,ivn,z,dz,rzu, &
                        onbb,leapfrog_c,dtuv,c2dtuv,afc1,afc2,cosu,imt,jmt,km,imm, &
                        jmm,west,east,north,south,asselin_c,boussinesq,itn)
!     =================
!
!     compute up & vp at tau+1 time level
!
      implicit none
      include 'pconst.h'
      include 'mpif.h'
!
      integer imt,jmt,km,imm,jmm,i,j,k,asselin_c,boussinesq
      integer ivn(imt,jmt),itn(imt,jmt)
      logical leapfrog_c
      real t1,t2
      real c2dtuv,dtuv,afc1,afc2,onbb,z(km),dz(km),rzu(imt,jmt),cosu(jmt)
      real a(imt,jmt),b(imt,jmt),c(imt,jmt),px(imt,jmt,km),py(imt,jmt,km)
      real etax(imt,jmt),etay(imt,jmt),pbar(imt,jmt)
      real pmup(imt,jmt),pmum(imt,jmt),pbt(imt,jmt,2),spbt(imt,jmt),phib(imt,jmt)
      real rhodp(imt,jmt,km),rho(imt,jmt,km)
      real rdxt(jmt),rdy,ff(jmt),umask(imt,jmt,km)
      real ump(imt,jmt,km),vmp(imt,jmt,km),upb(imt,jmt,2),vpb(imt,jmt,2)
      real up(imt,jmt,km,2),vp(imt,jmt,km,2),du(imt,jmt,km),dv(imt,jmt,km)
      real epea(jmt),epeb(jmt),epla(jmt),eplb(jmt)
      
      integer west,east,north,south
!
!-----------------------------------------------------------------------
!     calculate time-averaged pbt & sqrt(pbt)
!-----------------------------------------------------------------------
      do j=1,jmt
      do i=1,imt
      pmup(i,j) = pmup(i,j)*onbb
      enddo
      enddo
!
      if(leapfrog_c)then
        do j=1,jmt
        do i=1,imt
        pbar(i,j) = (pmup(i,j)+pmum(i,j))*p5
        enddo
        enddo
      else
        do j=1,jmt
        do i=1,imt
        pbar(i,j) = pmup(i,j)
        enddo
        enddo
      endif
!
      do j=2,jmm
      do i=2,imm
      if(ivn(i,j).gt.0)then
      spbt(i,j) = p5*sqrt(pbar(i,j  ) + pbar(i+1,j  ) +  &
                          pbar(i,j+1) + pbar(i+1,j+1))
      endif
!      spbt(1  ,j) = spbt(imm,j)
!      spbt(imt,j) = spbt(2  ,j)
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
      if (boussinesq==1) then
!
!     a  = geopotential
!     b  = pressure
!     px = (rho)xbar*(a)x + (b)x
!     py = (rho)ybar*(a)y + (b)y
!
      do k=1,km
      do j=1,jmt
      do i=1,imt
      a(i,j) = phib(i,j)-pbar(i,j)*(z(k)+phib(i,j))
      b(i,j) = pbar(i,j)*rhodp(i,j,k)
      enddo
      enddo
      do j=2,jmm
      do i=2,imm
      px(i,j,k) = ( (rho(i,j,k)+rho(i+1,j,k))*p5*(a(i+1,j)-a(i,j)) + &
                  (b(i+1,j)-b(i,j)) )*rdxt(j)
      py(i,j,k) = ( (rho(i,j,k)+rho(i,j+1,k))*p5*(a(i,j+1)-a(i,j)) + &
                  (b(i,j+1)-b(i,j)) )*rdy
      enddo
      enddo
      end do
      
      else
      
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
      end if
      
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
      a(i,j) = du(i,j,k) + ff(j)*vp(i,j,k,tau) - etax(i,j)
      b(i,j) = dv(i,j,k) - ff(j)*up(i,j,k,tau) - etay(i,j)
      enddo
      enddo
!
!
!-----------------------------------------------------------------------
!     coriolis adjustment
!-----------------------------------------------------------------------
      if(leapfrog_c)then
        do j=2,jmm
        do i=2,imm
        du(i,j,k) = (epla(j)*a(i,j) + eplb(j)*b(i,j))*umask(i,j,k)
        dv(i,j,k) = (epla(j)*b(i,j) - eplb(j)*a(i,j))*umask(i,j,k)
        enddo
        enddo
      else
        do j=2,jmm
        do i=2,imm
        du(i,j,k) = (epea(j)*a(i,j) + epeb(j)*b(i,j))*umask(i,j,k)
        dv(i,j,k) = (epea(j)*b(i,j) - epeb(j)*a(i,j))*umask(i,j,k)
        enddo
        enddo
      endif
!
100   continue 
!
!
!-----------------------------------------------------------------------
!     compute up & vp at tau+1 time level
!-----------------------------------------------------------------------
!
      if(leapfrog_c) go to 120
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
      endif
      enddo
      enddo
      enddo
!
      leapfrog_c = .true.
!
      go to 150
!
120   continue 
!
      do k=1,km
      do j=2,jmm
      do i=2,imm
      if(umask(i,j,k).gt.c0)then
      du(i,j,k) = up(i,j,k,taum) + du(i,j,k)*c2dtuv
      dv(i,j,k) = vp(i,j,k,taum) + dv(i,j,k)*c2dtuv
      endif
      enddo
      enddo
      enddo
!
!@@@  interaction between barotropic and baroclinic modes
!
      call vinteg(du,a,ivn,dz,rzu,imt,jmt,km)
      call vinteg(dv,b,ivn,dz,rzu,imt,jmt,km)
!
      if (asselin_c==1) then
      do k=1,km
      do j=2,jmm
      do i=2,imm
      if(umask(i,j,k).gt.c0)then
      t1             = du(i,j,k) - a(i,j) + upb(i,j,tau)
      t2             = dv(i,j,k) - b(i,j) + vpb(i,j,tau)
      up(i,j,k,taum) = afc2*up(i,j,k,tau)+afc1*(up(i,j,k,taum)+t1)
      vp(i,j,k,taum) = afc2*vp(i,j,k,tau)+afc1*(vp(i,j,k,taum)+t2)
      up(i,j,k,tau)  = t1
      vp(i,j,k,tau)  = t2
      endif
      enddo
      enddo
      enddo
      else
      do k=1,km
      do j=2,jmm
      do i=2,imm
      if(umask(i,j,k).gt.c0)then
      t1             = du(i,j,k) - a(i,j) + upb(i,j,tau)
      t2             = dv(i,j,k) - b(i,j) + vpb(i,j,tau)
      up(i,j,k,taum) = up(i,j,k,tau)
      vp(i,j,k,taum) = vp(i,j,k,tau)
      up(i,j,k,tau)  = t1
      vp(i,j,k,tau)  = t2
      endif
      enddo
      enddo
      enddo
      end if
!
150   continue
!
!
!      do k=1,km
!      do j=2,jmm
!      up(1  ,j,k,tau) = up(imm,j,k,tau)
!      up(imt,j,k,tau) = up(2  ,j,k,tau)
!      vp(1  ,j,k,tau) = vp(imm,j,k,tau)
!      vp(imt,j,k,tau) = vp(2  ,j,k,tau)
!      up(1  ,j,k,taum) = up(imm,j,k,taum)
!      up(imt,j,k,taum) = up(2  ,j,k,taum)
!      vp(1  ,j,k,taum) = vp(imm,j,k,taum)
!      vp(imt,j,k,taum) = vp(2  ,j,k,taum)
!      enddo
!      enddo
      call swap_array_real4d(up,imt,jmt,km,2,west,east,north,south)
      call swap_array_real4d(vp,imt,jmt,km,2,west,east,north,south)
!
!
!-----------------------------------------------------------------------
!     constract time-averaged mass advections
!-----------------------------------------------------------------------
      do j=2,jmm
      do i=2,imm
      spbt(i,j) = p5*sqrt(pbt(i,j  ,tau) + pbt(i+1,j  ,tau) + &
                         pbt(i,j+1,tau) + pbt(i+1,j+1,tau))
      enddo
!      spbt(1  ,j) = spbt(imm,j)
!      spbt(imt,j) = spbt(2  ,j)
      enddo
      call swap_array_real2d(spbt,imt,jmt,west,east,north,south)
!
      do k=1,km
      do j=2,jmm
      do i=2,imm
      if(umask(i,j,k).gt.c0)then
      ump(i,j,k) = ump(i,j,k) + up(i,j,k,tau)*spbt(i,j)
      vmp(i,j,k) = vmp(i,j,k) + vp(i,j,k,tau)*spbt(i,j)*cosu(j)
      endif
      enddo
!      ump(1  ,j,k) = ump(imm,j,k)
!      ump(imt,j,k) = ump(2  ,j,k)
!      vmp(1  ,j,k) = vmp(imm,j,k)
!      vmp(imt,j,k) = vmp(2  ,j,k)
      enddo
      enddo
      call swap_array_real3d(ump,imt,jmt,km,west,east,north,south)
      call swap_array_real3d(vmp,imt,jmt,km,west,east,north,south)

      return
      end
