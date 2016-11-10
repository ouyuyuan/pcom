c     ==========================
      subroutine vinteg(wk3,wk2)
c     ==========================
c     calculate vertical integration at "u" grid
c
      implicit none
      include 'param.h'
      include 'pconst.h'
      include 'grdvar.h'
      real wk3(imt,jmt,km),wk2(imt,jmt)
c
      do j=1,jmt
      do i=1,imt
      wk2(i,j)=c0
      do k=1,ivn(i,j)
      wk2(i,j)=wk2(i,j)+dz(k)*wk3(i,j,k)
      enddo
      wk2(i,j)=wk2(i,j)*rzu(i,j)
      enddo
      enddo
c
      return
      end
c
c
c
c     ==============================
      subroutine vinteg_ns(wk,an,as)
c     ==============================
c     calculate vertical integration at "u" grid
c
      implicit none
      include 'param.h'
      include 'pconst.h'
      include 'grdvar.h'
c
      real wk(imt,jmt,km),an(imt,jmt),as(imt,jmt)
      real tmax,tmin
c
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
c
c
c
c     ==============================
      subroutine vinteg_ew(wk,ae,aw)
c     ==============================
c     calculate vertical integration at "u" grid
c
      implicit none
      include 'param.h'
      include 'pconst.h'
      include 'grdvar.h'
c
      real wk(imt,jmt,km),ae(imt,jmt),aw(imt,jmt)
      real tmax,tmin
c
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
