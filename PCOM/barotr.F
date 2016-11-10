c
c     =================
      subroutine barotr
c     =================
c     compute pbt, upb & vpb at "tau+1" time level
c
      implicit none
      include 'param.h'
      include 'pconst.h'
      include 'scalar.h'
      include 'grdvar.h'
      include 'prog.h'
      include 'cvbc.h'
c
      integer mode_b
      real a,b,dp,adu,adv,etax,etay,du2,dv2,t1,t2,t3
c
      common /works/ a(imt,jmt),b(imt,jmt),dp(imt,jmt)
      common /works/ adu(imt,jmt),adv(imt,jmt)
      common /works/ etax(imt,jmt),etay(imt,jmt)
      common /works/ du2(imt,jmt),dv2(imt,jmt)
c
c
      do 1000 mode_b = 1,nbb
c
c-----------------------------------------------------------------------
c     compute the "artificial" horizontal viscosity
c-----------------------------------------------------------------------
c
      do j=2,jmm
      do i=2,imm
      if(ivn(i,j).gt.0)then
      adu(i,j) = am*( r1c(j)*(upb(i,j+1,taum)-upb(i,j,taum)) -
     &                r1d(j)*(upb(i,j,taum)-upb(i,j-1,taum)) +
     &                sdxu(j)*(upb(i+1,j,taum)-c2*upb(i,j,taum)+
     &                         upb(i-1,j,taum)) )*c2
      adv(i,j) = am*( r1c(j)*(vpb(i,j+1,taum)-vpb(i,j,taum)) -
     &                r1d(j)*(vpb(i,j,taum)-vpb(i,j-1,taum)) +
     &                sdxu(j)*(vpb(i+1,j,taum)-c2*vpb(i,j,taum)+
     &                         vpb(i-1,j,taum)) )*c2
      endif
      enddo
      enddo
c
      if(mode_b.eq.1) then
        do j=2,jmm
        do i=2,imm
        dub(i,j) = dub(i,j) - adu(i,j)
        dvb(i,j) = dvb(i,j) - adv(i,j)
        enddo
        enddo
      end if
c
c
1001  continue
c
c-----------------------------------------------------------------------
c     square root of pbt at the U cell
c-----------------------------------------------------------------------
      do j=2,jmm
      do i=2,imm
      if(ivn(i,j).gt.0)then
      spbt(i,j) = p5*sqrt(pbt(i,j  ,tau) + pbt(i+1,j  ,tau) +
     &                    pbt(i,j+1,tau) + pbt(i+1,j+1,tau))
      endif
      spbt(1  ,j) = spbt(imm,j)
      spbt(imt,j) = spbt(2  ,j)
      enddo
      enddo
c
c
c-----------------------------------------------------------------------
c     calculate the pressure gradients due to change of free surface
c-----------------------------------------------------------------------
#ifdef boussinesq
      do j=2,jmm
      do i=1,imt
      a(i,j) = phib(i,j)*(c1-pbt(i,j,tau))
      enddo
      enddo
#endif
c
c
      do j=2,jmm
      do i=2,imm
      if(ivn(i,j).gt.0)then
c
      etax(i,j) = spbt(i,j)*p5*(
#ifdef boussinesq
     &        (a(i+1,j  )-a(i,j  ))*pdxn(i,j)
     &       +(a(i+1,j+1)-a(i,j+1))*pdxs(i,j+1) +
#else
     &       + phibx(i,j) +
#endif
     &        (pbt(i+1,j  ,tau)-pbt(i,j  ,tau))*pbxn(i,j)
     &       +(pbt(i+1,j+1,tau)-pbt(i,j+1,tau))*pbxs(i,j+1) +
     &        (pbt(i+1,j  ,tau)+pbt(i,j  ,tau))*pcxn(i,j)
     &       +(pbt(i+1,j+1,tau)+pbt(i,j+1,tau))*pcxs(i,j+1) )
c
c
      etay(i,j) = spbt(i,j)*p5*(
#ifdef boussinesq
     &        (a(i  ,j+1)-a(i  ,j))*pdye(i,j)
     &       +(a(i+1,j+1)-a(i+1,j))*pdyw(i+1,j) +
#else
     &       + phiby(i,j) +
#endif
     &        (pbt(i  ,j+1,tau)-pbt(i,  j,tau))*pbye(i  ,j)
     &       +(pbt(i+1,j+1,tau)-pbt(i+1,j,tau))*pbyw(i+1,j) +
     &        (pbt(i  ,j+1,tau)+pbt(i,  j,tau))*pcye(i  ,j)
     &       +(pbt(i+1,j+1,tau)+pbt(i+1,j,tau))*pcyw(i+1,j) )
c
      endif
      enddo
      enddo
c
c-----------------------------------------------------------------------
c     advection + viscosiy + pressure gradient + coriolis
c-----------------------------------------------------------------------
c
      do j=1,jmt
      do i=1,imt
      if(ivn(i,j).gt.0)then
        a(i,j) = dub(i,j) + adu(i,j) + ff(j)*vpb(i,j,tau) - etax(i,j)
        b(i,j) = dvb(i,j) + adv(i,j) - ff(j)*upb(i,j,tau) - etay(i,j)
      else
        a(i,j) = c0
        b(i,j) = c0
      endif
      enddo
      enddo
c
c
c-----------------------------------------------------------------------
c     coriolis adjustment
c-----------------------------------------------------------------------
c
      if(leapfrog_b)then
        do j=1,jmt
        do i=1,imt
        du2(i,j) = ebla(j)*a(i,j) + eblb(j)*b(i,j)
        dv2(i,j) = ebla(j)*b(i,j) - eblb(j)*a(i,j)
        enddo
        enddo
      else
        do j=1,jmt
        do i=1,imt
        du2(i,j) = ebea(j)*a(i,j) + ebeb(j)*b(i,j)
        dv2(i,j) = ebea(j)*b(i,j) - ebeb(j)*a(i,j)
        enddo
        enddo
      endif
c
c
c-----------------------------------------------------------------------
c     calculate the change of Pbt
c-----------------------------------------------------------------------
c
      do j=1,jmt
      do i=1,imt
      a(i,j) = zu(i,j)*spbt(i,j)*upb(i,j,tau)
      b(i,j) = zu(i,j)*spbt(i,j)*vpb(i,j,tau)*cosu(j)
      enddo
      enddo
c
c
      do j=2,jmm
      do i=2,imm
      if(itn(i,j).gt.0)then
       dp(i,j) = ( - p5*(rdxt(j)*(a(i,j)+a(i,j-1)-a(i-1,j)-a(i-1,j-1))
     &                  +rdyt(j)*(b(i,j)+b(i-1,j)-b(i,j-1)-b(i-1,j-1)))
#ifdef snbc
     &             - emp(i,j)
#endif
     &             )/pn(i,j)
      else
       dp(i,j) = c0
      endif
      enddo
      enddo
c
c
#ifdef smth
      print*, ' check I&J boundary for smooth!!!!!!'
      stop 111222
      call smths(du2,umask, 2, 5)	!! from 68S-74S
      call smths(dv2,umask, 2, 5)
      call smths(du2,umask,73,76)	!! from 68N-74N
      call smths(dv2,umask,73,76)
#endif
c
c
c
c-----------------------------------------------------------------------
c     compute pbt, upb & vpb at tau+1 time level
c-----------------------------------------------------------------------
c
      if(leapfrog_b) go to 120
c
       do j=2,jmm
       do i=2,imm
       upb(i,j,tau) = upb(i,j,taum) + du2(i,j)*dtsf
       vpb(i,j,tau) = vpb(i,j,taum) + dv2(i,j)*dtsf
       pbt(i,j,tau) = pbt(i,j,taum) + dp(i,j) *dtsf
       enddo
       enddo
c
       if(euler_back) then
        euler_back = .false.
        go to 1001
       endif
c
       leapfrog_b = .true.
c
       go to 150
c
c
120   continue
c
c
      do j=2,jmm
      do i=2,imm
      t1            = pbt(i,j,taum) + dp (i,j)*c2dtsf
      t2            = upb(i,j,taum) + du2(i,j)*c2dtsf
      t3            = vpb(i,j,taum) + dv2(i,j)*c2dtsf
#ifdef asselin_b
      pbt(i,j,taum) = afb2*pbt(i,j,tau)+afb1*(pbt(i,j,taum)+t1)
      upb(i,j,taum) = afb2*upb(i,j,tau)+afb1*(upb(i,j,taum)+t2)
      vpb(i,j,taum) = afb2*vpb(i,j,tau)+afb1*(vpb(i,j,taum)+t3)
#else
      pbt(i,j,taum) = pbt(i,j,tau)
      upb(i,j,taum) = upb(i,j,tau)
      vpb(i,j,taum) = vpb(i,j,tau)
#endif
      pbt(i,j,tau)  = t1
      upb(i,j,tau)  = t2
      vpb(i,j,tau)  = t3
      enddo
      enddo
c
150   continue
c
c
      do j=2,jmm
      upb(1  ,j,tau) = upb(imm,j,tau)
      upb(imt,j,tau) = upb(2  ,j,tau)
      vpb(1  ,j,tau) = vpb(imm,j,tau)
      vpb(imt,j,tau) = vpb(2  ,j,tau)
      pbt(1  ,j,tau) = pbt(imm,j,tau)
      pbt(imt,j,tau) = pbt(2  ,j,tau)
      upb(1  ,j,taum) = upb(imm,j,taum)
      upb(imt,j,taum) = upb(2  ,j,taum)
      vpb(1  ,j,taum) = vpb(imm,j,taum)
      vpb(imt,j,taum) = vpb(2  ,j,taum)
      pbt(1  ,j,taum) = pbt(imm,j,taum)
      pbt(imt,j,taum) = pbt(2  ,j,taum)
      enddo
c
c
c-----------------------------------------------------------------------
c     constract time-averaged pbt for solving baroclinc and tracer Eqs.
c-----------------------------------------------------------------------
      do j=2,jmm
      do i=1,imt
      pmup(i,j) = pmup(i,j) + pbt(i,j,tau)
      pmtp(i,j) = pmtp(i,j) + pbt(i,j,tau)
      enddo
      enddo
c
c
1000  continue
c
      return
      end
