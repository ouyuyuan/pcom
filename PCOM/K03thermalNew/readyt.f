c
c     =================
      subroutine readyt
c     =================
c     calculate the variables depended on stratification
c
      implicit none
      include 'param.h'
      include 'pconst.h'
      include 'scalar.h'
      include 'grdvar.h'
      include 'prog.h'
c
      real a,t1,t2,t3,rdens
c
      common /works/ a(imt,jmt)
c
c
c-----------------------------------------------------------------------
c     calculate density/rho_0 (BCOM) or 1/density (PCOM)
c-----------------------------------------------------------------------
      call density
c
c
c-----------------------------------------------------------------------
c     initialize the time-averaged pbt & up vp
c-----------------------------------------------------------------------
      do j=1,jmt
      do i=1,imt
      pmtm(i,j) = pmtp(i,j)
      pmtp(i,j) = pbt (i,j,tau)
      enddo
      enddo
c
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
      do j=1,jmt
      do i=1,imt
      umm(i,j,k) = ump(i,j,k)
      vmm(i,j,k) = vmp(i,j,k)
      ump(i,j,k) = up(i,j,k,tau)*spbt(i,j)
      vmp(i,j,k) = vp(i,j,k,tau)*spbt(i,j)*cosu(j)
      enddo
      enddo
      enddo
c
c
c-----------------------------------------------------------------------
c     calculate pressure gradients related to stratification
c-----------------------------------------------------------------------
      call prsgrd
c
c
      return
      end
