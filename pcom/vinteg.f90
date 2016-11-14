!     ==========================
      subroutine vinteg(wk3,wk2,ivn,dz,rzu,imt,jmt,km)
!     ==========================
!     calculate vertical integration at "u" grid
!
      implicit none
      include 'pconst.h'
      include 'mpif.h'
      
      integer imt,jmt,km,i,j,k
      real wk3(imt,jmt,km),wk2(imt,jmt)
      integer ivn(imt,jmt)
      real dz(km),rzu(imt,jmt)
!
      wk2=c0
      do j=2,jmt-1
      do i=2,imt-1
!      wk2(i,j)=c0
      do k=1,ivn(i,j)
      wk2(i,j)=wk2(i,j)+dz(k)*wk3(i,j,k)
      enddo
      wk2(i,j)=wk2(i,j)*rzu(i,j)
      enddo
      enddo
!
      return
      end
!
!
!
!     ==============================
      subroutine vinteg_ns(wk,an,as,ivn,dz,rzu,imt,jmt,km,imm,jmm)
!     ==============================
!     calculate vertical integration at "u" grid
!
      implicit none
      include 'pconst.h'
      include 'mpif.h'
!
      integer imt,jmt,km,imm,jmm,i,j,k,n,m
      real wk(imt,jmt,km),an(imt,jmt),as(imt,jmt)
      real tmax,tmin
      integer ivn(imt,jmt)
      real dz(km),rzu(imt,jmt)
      
!
      do j=2,jmm
      do i=2,imm
       m = min(ivn(i,j),ivn(i,j-1))
       n = max(ivn(i,j),ivn(i,j-1))
       tmin = c0
       do k=1,m
       tmin = tmin + wk(i,j,k)*dz(k)
       enddo
       tmax = tmin
       do k=m+1,n
       tmax = tmax + wk(i,j,k)*dz(k)
       enddo
       if(ivn(i,j).eq.m) then
        an(i,j) = tmin*rzu(i,j)
        as(i,j) = tmax*rzu(i,j-1)
       else
        an(i,j) = tmax*rzu(i,j)
        as(i,j) = tmin*rzu(i,j-1)
       endif
      enddo
      enddo
      return
      end
!
!
!
!     ==============================
      subroutine vinteg_ew(wk,ae,aw,ivn,dz,rzu,imt,jmt,km,imm,jmm)
!     ==============================
!     calculate vertical integration at "u" grid
!
      implicit none
      include 'pconst.h'
      include 'mpif.h'
!
      integer imt,jmt,km,imm,jmm,i,j,k,n,m
      real wk(imt,jmt,km),ae(imt,jmt),aw(imt,jmt)
      real tmax,tmin
      integer ivn(imt,jmt)
      real dz(km),rzu(imt,jmt)
!
      do j=2,jmm
      do i=2,imm
       m = min(ivn(i,j),ivn(i-1,j))
       n = max(ivn(i,j),ivn(i-1,j))
       tmin = c0
       do k=1,m
       tmin = tmin + wk(i,j,k)*dz(k)
       enddo
       tmax = tmin
       do k=m+1,n
       tmax = tmax + wk(i,j,k)*dz(k)
       enddo
       if(ivn(i,j).eq.m) then
        ae(i,j) = tmin*rzu(i,j)
        aw(i,j) = tmax*rzu(i-1,j)
       else
        ae(i,j) = tmax*rzu(i,j)
        aw(i,j) = tmin*rzu(i-1,j)
       endif
      enddo
      enddo
      return
      end
