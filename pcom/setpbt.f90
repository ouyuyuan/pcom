!     =================
!BOP
!
! !MODULE: setpbt.f90
! !DESCRIPTION: \input{sections/code-setpbt}
!
! !INTERFACE:
!
      subroutine setpbt(rho,pbt,itn,pt,ps,dz,decibar,imt,jmt,km,kmp1,nt,  &
                        fixp)
!EOP
!-------------------------------------------------------------------------------
!     =================
      implicit none
      include 'pconst.h'
      include 'mpif.h'
!
      integer imt,jmt,km,kmp1,nt,i,j,k
      real t0,s0,p0,dens,undens,prea,preb,decibar,pre(kmp1),dz(km),fixp(imt,jmt,km)
      integer itn(imt,jmt)
      real rho(imt,jmt,km),pbt(imt,jmt,2),pt(imt,jmt,km),ps(imt,jmt,km)
!
      do k=1,km
      do j=1,jmt
      do i=1,imt
      rho(i,j,k)   = c1
      pbt(i,j,tau) = c1
      enddo
      enddo
      enddo
      
!
!-----------------------------------------------------------------------
!     initialize pbt & rho
!-----------------------------------------------------------------------

      do j=1,jmt
      do i=1,imt
      if(itn(i,j).eq.0) goto 100
!
      pre(1) = c0
      do k=1,itn(i,j)
      pre(k+1) = pre(k) + dz(k)
      enddo
!
      do k=1,itn(i,j)
      t0         = pt(i,j,k)
      s0         = ps(i,j,k)
      p0         = (pre(k)+pre(k+1))*p5*decibar
      rho(i,j,k) = undens(t0,s0,p0)
      enddo
!
      pbt(i,j,tau) = c1
      
100   continue
      enddo
      enddo

      return
      end
