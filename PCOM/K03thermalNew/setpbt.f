c     =================
      subroutine setpbt
c     =================
      implicit none
      include 'param.h'
      include 'pconst.h'
      include 'scalar.h'
      include 'grdvar.h'
      include 'prog.h'
c
      real t0,s0,p0,dens,prea,preb,pre(kmp1)
c
      do k=1,km
      do j=1,jmt
      do i=1,imt
      rho(i,j,k)   = c1
      pbt(i,j,tau) = c1
      enddo
      enddo
      enddo
c
#ifdef boussinesq
      do k=1,km
      do j=1,jmt
      do i=1,imt
      fixp(i,j,k) = c0
      enddo
      enddo
      enddo
#endif
c
c-----------------------------------------------------------------------
c     initialize pbt & rho
c-----------------------------------------------------------------------
      do 100 j=1,jmt
      do 100 i=1,imt
      if(itn(i,j).eq.0) goto 100
c
      pre(1) = c0
      do k=1,itn(i,j)
      pre(k+1) = pre(k) + dz(k)
      enddo
c
      do k=1,itn(i,j)
      t0         = pt(i,j,k)
      s0         = ps(i,j,k)
      p0         = (pre(k)+pre(k+1))*p5*decibar
      rho(i,j,k) = dens(t0,s0,p0)
      enddo
c
      pbt(i,j,tau) = c1
c
c
#ifdef boussinesq
      do k=1,itn(i,j)
      fixp(i,j,k) = (pre(k)+pre(k+1))*p5*decibar
      enddo
#endif
c
100   continue
c
c
#ifdef boussinesq
      do j=1,jmt
      do i=1,imt
      pbt(i,j,tau) = c1
      enddo
      enddo
#endif
c
ctmp-------------------------------------------------------------
c     print*, 'pbt(90,40,1) =', pbt(90,40,tau)
ctmp-------------------------------------------------------------
      return
      end
