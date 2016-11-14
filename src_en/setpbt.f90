!     =================
      subroutine setpbt(rho,pbt,itn,pt,ps,dz,decibar,imt,jmt,km,kmp1,nt,unesco,  &
                        boussinesq,fixp)
!     =================
      implicit none
      include 'pconst.h'
      include 'mpif.h'
!
      integer imt,jmt,km,kmp1,nt,i,j,k,unesco,boussinesq
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
      
      if (boussinesq==1) then
      do k=1,km
      do j=1,jmt
      do i=1,imt
      fixp(i,j,k) = c0
      enddo
      enddo
      enddo
      end if
!
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
      if (unesco==1) then
      do k=1,itn(i,j)
      t0         = pt(i,j,k)
      s0         = ps(i,j,k)
      p0         = (pre(k)+pre(k+1))*p5*decibar
      rho(i,j,k) = undens(t0,s0,p0)
      enddo
      else
      do k=1,itn(i,j)
      t0         = pt(i,j,k)
      s0         = ps(i,j,k)
      p0         = (pre(k)+pre(k+1))*p5*decibar
      rho(i,j,k) = dens(t0,s0,p0)
      enddo
      end if
!
      pbt(i,j,tau) = c1
      
      if (boussinesq==1) then
      do k=1,itn(i,j)
      fixp(i,j,k) = (pre(k)+pre(k+1))*p5*decibar
      enddo
      end if
      
100   continue
      enddo
      enddo

      if (boussinesq==1) then
      do j=1,jmt
      do i=1,imt
      pbt(i,j,tau) = c1
      enddo
      enddo
      end if

      return
      end
