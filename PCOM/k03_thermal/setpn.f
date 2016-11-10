c
c     ================
      subroutine setpn
c     ================
c     calculate pn, z & dz
c
c     PCOM  pn = constant bottom pressure   (dy/cm2)
c     BCOM  pn = grav * thickness           (cm*cm/s2)
c     z = eta; dz=1/z
c
      implicit none
      include 'param.h'
      include 'pconst.h'
      include 'scalar.h'
      include 'grdvar.h'
c
      real p0,dens,pres
      real pre(kmp1),denz(km)
      real t30(km),s30(km)
c
c---------------------------------------------------------------------
c     input global averaged TS stratification
c---------------------------------------------------------------------
      open(82,file='ts30.data',status='old',form='formatted')
      read(82,*) t30,s30
      close(82)
c
c
c---------------------------------------------------------------------
c     calculate depth of each layer (in cm)
c---------------------------------------------------------------------
      dz0(1) = z0(1)*c2
      do k=2,km
      dz0(k) = (z0(k)-z0(k-1))*c2 - dz0(k-1)
      enddo
c
c---------------------------------------------------------------------
c     calculate pressure at the bottom of each layer
c---------------------------------------------------------------------
c     set an initial density
      do k=1,km
      denz(k) = c1
      enddo
c
      pres = c0
c
25    pre(1)   = c0
      do k=1,km
      pre(k+1) = pre(k) + denz(k)*grav*dz0(k)
      enddo
c
      do k=1,km
      p0      = (pre(k)+pre(k+1))*p5*decibar
      denz(k) = dens(t30(k),s30(k),p0)
      enddo
c
      if(abs(pre(kmp1)-pres).gt.c1em4)then
       pres = pre(kmp1)
       goto 25
      endif
c
c---------------------------------------------------------------------
c     calculate z & dz
c---------------------------------------------------------------------
#ifdef boussinesq
      do k=1,km
       z(k)  =  z0(k) * grav
       dz(k) = dz0(k) * grav
      enddo
#else
      do k=1,km
        z(k) = (pre(k+1)+pre(k))*p5
       dz(k) =  pre(k+1)-pre(k)
      enddo
#endif
c
c---------------------------------------------------------------------
c     calculate pn
c---------------------------------------------------------------------
#ifdef boussinesq
      do j=1,jmt
      do i=1,imt
      if(itn(i,j).eq.0) then
       pn(i,j) = c1
      else
       pn(i,j) = abs(phib(i,j))
      endif
      enddo
      enddo
#else
      do j=1,jmt
      do i=1,imt
      if(itn(i,j).eq.0) then
       pn(i,j) = c1
      else
       pn(i,j) = pre(itn(i,j)+1)
      endif
      enddo
      enddo
#endif
c
ctmp-------------------------------------------------------------
c     print*, 'pn(90,40) & phib(90,40) =', pn(90,40),phib(90,40)
ctmp-------------------------------------------------------------
      return
      end
