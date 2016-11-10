c
c     =================
      subroutine bclinc
c     =================
c
c     compute up & vp at tau+1 time level
c
      implicit none
      include 'param.h'
      include 'pconst.h'
      include 'scalar.h'
      include 'grdvar.h'
      include 'prog.h'
c
      real pbar,etax,etay,px,py,a,b,c,t1,t2
      common /works/ a(imt,jmt),b(imt,jmt),c(imt,jmt)
      common /works/ px(imt,jmt),py(imt,jmt)
      common /works/ etax(imt,jmt),etay(imt,jmt)
      common /works/ pbar(imt,jmt)
c
c
c-----------------------------------------------------------------------
c     calculate time-averaged pbt & sqrt(pbt)
c-----------------------------------------------------------------------
      do j=1,jmt
      do i=1,imt
      pmup(i,j) = pmup(i,j)*onbb
      enddo
      enddo
c
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
c
      do j=2,jmm
      do i=2,imm
      if(ivn(i,j).gt.0)then
      spbt(i,j) = p5*sqrt(pbar(i,j  ) + pbar(i+1,j  ) +
     &                    pbar(i,j+1) + pbar(i+1,j+1))
      endif
      spbt(1  ,j) = spbt(imm,j)
      spbt(imt,j) = spbt(2  ,j)
      enddo
      enddo
c
c
      do 100 k=1,km
c
c-----------------------------------------------------------------------
c     calculate pressure gradients  * spbt
c-----------------------------------------------------------------------
#ifdef boussinesq
c
c     a  = geopotential
c     b  = pressure
c     px = (rho)xbar*(a)x + (b)x
c     py = (rho)ybar*(a)y + (b)y
c
      do j=2,jmm
      do i=1,imt
      a(i,j) = phib(i,j)-pbar(i,j)*(z(k)+phib(i,j))
      b(i,j) = pbar(i,j)*rhodp(i,j,k)
      enddo
      enddo
      do j=2,jmm
      do i=2,imm
      px(i,j) = ( (rho(i,j,k)+rho(i+1,j,k))*p5*(a(i+1,j)-a(i,j)) +
     &            (b(i+1,j)-b(i,j)) )*rdxt(j)
      py(i,j) = ( (rho(i,j,k)+rho(i,j+1,k))*p5*(a(i,j+1)-a(i,j)) +
     &            (b(i,j+1)-b(i,j)) )*rdy
      enddo
      enddo
#else
c
c     a  = geopotential
c     b  = pressure
c     px = (a)x + (rho)xbar*(b)x
c     py = (a)y + (rho)ybar*(b)y
c
      do j=2,jmm
      do i=1,imt
      a(i,j) = phib(i,j) + pbar(i,j)*rhodp(i,j,k)
      b(i,j) = pbar(i,j)*z(k)
      enddo
      enddo
c
      do j=2,jmm
      do i=2,imm
      px(i,j) = ( (rho(i,j,k)+rho(i+1,j,k))*p5*(b(i+1,j)-b(i,j)) +
     &            (a(i+1,j)-a(i,j)) )*rdxt(j)
      py(i,j) = ( (rho(i,j,k)+rho(i,j+1,k))*p5*(b(i,j+1)-b(i,j)) +
     &            (a(i,j+1)-a(i,j)) )*rdy
      enddo
      enddo
#endif
c
      do j=2,jmm
      do i=2,imm
      etax(i,j) = (px(i,j)+px(i,j+1))*p5*spbt(i,j)
      etay(i,j) = (py(i,j)+py(i+1,j))*p5*spbt(i,j)
      enddo
      enddo
c
c
c---------------------------------------------------------------------
c     advection + diffusion + pressure gradients + coriolis
c---------------------------------------------------------------------
      do j=2,jmm
      do i=2,imm
      a(i,j) = du(i,j,k) + ff(j)*vp(i,j,k,tau) - etax(i,j)
      b(i,j) = dv(i,j,k) - ff(j)*up(i,j,k,tau) - etay(i,j)
      enddo
      enddo
c
c
c-----------------------------------------------------------------------
c     coriolis adjustment
c-----------------------------------------------------------------------
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
c
100   continue 
c
c
c-----------------------------------------------------------------------
c     compute up & vp at tau+1 time level
c-----------------------------------------------------------------------
c
      if(leapfrog_c) go to 120
c
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
c
c@@@  interaction between barotropic and baroclinic modes
c
      call vinteg(up(1,1,1,tau),a)
      call vinteg(vp(1,1,1,tau),b)
c
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
c
      leapfrog_c = .true.
c
      go to 150
c
120   continue 
c
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
c
c@@@  interaction between barotropic and baroclinic modes
c
      call vinteg(du,a)
      call vinteg(dv,b)
c
      do k=1,km
      do j=2,jmm
      do i=2,imm
      if(umask(i,j,k).gt.c0)then
      t1             = du(i,j,k) - a(i,j) + upb(i,j,tau)
      t2             = dv(i,j,k) - b(i,j) + vpb(i,j,tau)
#ifdef asselin_c
      up(i,j,k,taum) = afc2*up(i,j,k,tau)+afc1*(up(i,j,k,taum)+t1)
      vp(i,j,k,taum) = afc2*vp(i,j,k,tau)+afc1*(vp(i,j,k,taum)+t2)
#else
      up(i,j,k,taum) = up(i,j,k,tau)
      vp(i,j,k,taum) = vp(i,j,k,tau)
#endif
      up(i,j,k,tau)  = t1
      vp(i,j,k,tau)  = t2
      endif
      enddo
      enddo
      enddo
c
150   continue
c
c
      do k=1,km
      do j=2,jmm
      up(1  ,j,k,tau) = up(imm,j,k,tau)
      up(imt,j,k,tau) = up(2  ,j,k,tau)
      vp(1  ,j,k,tau) = vp(imm,j,k,tau)
      vp(imt,j,k,tau) = vp(2  ,j,k,tau)
      up(1  ,j,k,taum) = up(imm,j,k,taum)
      up(imt,j,k,taum) = up(2  ,j,k,taum)
      vp(1  ,j,k,taum) = vp(imm,j,k,taum)
      vp(imt,j,k,taum) = vp(2  ,j,k,taum)
      enddo
      enddo
c
c
c-----------------------------------------------------------------------
c     constract time-averaged mass advections
c-----------------------------------------------------------------------
      do j=2,jmm
      do i=2,imm
      spbt(i,j) = p5*sqrt(pbt(i,j  ,tau) + pbt(i+1,j  ,tau) +
     &                    pbt(i,j+1,tau) + pbt(i+1,j+1,tau))
      enddo
      spbt(1  ,j) = spbt(imm,j)
      spbt(imt,j) = spbt(2  ,j)
      enddo
c
      do k=1,km
      do j=2,jmm
      do i=2,imm
      if(umask(i,j,k).gt.c0)then
      ump(i,j,k) = ump(i,j,k) + up(i,j,k,tau)*spbt(i,j)
      vmp(i,j,k) = vmp(i,j,k) + vp(i,j,k,tau)*spbt(i,j)*cosu(j)
      endif
      enddo
      ump(1  ,j,k) = ump(imm,j,k)
      ump(imt,j,k) = ump(2  ,j,k)
      vmp(1  ,j,k) = vmp(imm,j,k)
      vmp(imt,j,k) = vmp(2  ,j,k)
      enddo
      enddo
c
c
      return
      end
