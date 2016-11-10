c
c     =================
      subroutine interp
c     =================
c
c     interpolates monthly mean forcing fileds
c
      implicit none
      include 'param.h'
      include 'pconst.h'
      include 'cvbc.h'
      include 'calendar.h'
c
      integer  pt1,pt2
      real     abc,factor
c
      abc = float(day-daymd(mth))
c
c
      if(abc.le.0.0) then
        pt1    = mth -1
        if(pt1.eq.0 ) pt1 = 12
        pt2    = mth
        factor = abc/float(daypm(pt1)) + c1
      else
        pt1    = mth
        pt2    = mod(mth,12)+1
        factor = abc/float(daypm(pt1))
      end if
c
c
      do j=1,jmt
      do i=1,imt
      bcu(i,j) = (bcf(i,j,pt2,1)-bcf(i,j,pt1,1))*factor + bcf(i,j,pt1,1)
      bcv(i,j) = (bcf(i,j,pt2,2)-bcf(i,j,pt1,2))*factor + bcf(i,j,pt1,2)
      bct(i,j) = (bcf(i,j,pt2,3)-bcf(i,j,pt1,3))*factor + bcf(i,j,pt1,3)
      bcp(i,j) = (bcf(i,j,pt2,4)-bcf(i,j,pt1,4))*factor + bcf(i,j,pt1,4)
      bcs(i,j) = (bcf(i,j,pt2,5)-bcf(i,j,pt1,5))*factor + bcf(i,j,pt1,5)
      emp(i,j) = (bcf(i,j,pt2,6)-bcf(i,j,pt1,6))*factor + bcf(i,j,pt1,6)
      ddd(i,j) = (bcf(i,j,pt2,7)-bcf(i,j,pt1,7))*factor + bcf(i,j,pt1,7)
      enddo
      enddo
c
c---------------------------------------------------------------------
c     ddd:  newtonia restoring coefficient for heat flux
c---------------------------------------------------------------------
c     ddd => grav*ddd/cp (1.0e-1=unit change)
c
      do j=1,jmt
      do i=1,imt
      ddd(i,j) = grav*ddd(i,j)/3901.0*1.0e-1
      enddo
      enddo
c
      return
      end
