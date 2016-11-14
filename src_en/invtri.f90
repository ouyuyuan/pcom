! 2006-1-10 #ifdef gm90
!     =================
      subroutine invtri (wk,topbc,dcb,aidif,dtts2,pbar,gravr,itn,tmask,rdz,rdzw,   &
                         imt,jmt,km,imm,jmm)
!     =================
      implicit none
      include 'pconst.h'
      include 'mpif.h'
!
      integer imt,jmt,km,imm,jmm,i,j,k,m
      real a8(km),b8(km),c8(km),d8(km),e8(0:km),f8(0:km),g8
      real wk(imt,jmt,km),topbc(imt,jmt),dcb(imt,jmt,km),pbar(imt,jmt)
      real aidif,dtts2,grp,topbct,gravr
      integer itn(imt,jmt)
      real tmask(imt,jmt,km),rdz(km),rdzw(km)
!
      do 100 j=2,jmm
      do 100 i=2,imm
      if(itn(i,j)<=1)goto 100
!
      grp    = gravr/pbar(i,j)
      topbct = topbc(i,j)/gravr
!
      do k=2,itn(i,j)
        a8(k)   = dcb(i,j,k-1)*rdzw(k-1)*grp*rdz(k)*dtts2*aidif
        c8(k)   = dcb(i,j,k  )*rdzw(k  )*grp*rdz(k)*dtts2*aidif
        b8(k)   = c1 + a8(k) + c8(k)
        d8(k)   = wk(i,j,k)
        e8(k-1) = c0
        f8(k-1) = c0
      enddo
!
!     b. c. at top
      k = 1
      a8(k)   =                    grp*rdz(k)*dtts2*aidif
      c8(k)   = dcb(i,j,k)*rdzw(k)*grp*rdz(k)*dtts2*aidif
      b8(k)   = c1 + c8(k)
      d8(k)   = wk(i,j,k)
      e8(k-1) = c0
      f8(k-1) = c0
!
!     b. c. at bottom
      m = itn(i,j)
      b8(m) = c1 + a8(m)
      c8(m) = grp*rdz(m)*dtts2*aidif
      e8(m) = c0
      f8(m) = c0
!
!     now invert
      do k=itn(i,j),1,-1
        g8      = c1/(b8(k)-c8(k)*e8(k))
        e8(k-1) = a8(k)*g8
        f8(k-1) = (d8(k)+c8(k)*f8(k))*g8        
      enddo
!
!     b.c. at surface
      wk(i,j,1) = (e8(0)*topbct + f8(0))*tmask(i,j,1)
      do k=2,itn(i,j)
      wk(i,j,k)=(e8(k-1)*wk(i,j,k-1)+f8(k-1))*tmask(i,j,k)
      enddo
!
100   continue
      return
      end
