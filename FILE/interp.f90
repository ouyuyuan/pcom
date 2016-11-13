!
!     =================
      subroutine interp(day,daymd,daypm,mth,bcf,bcu,bcv,bct,bcp,bcs,emp,ddd,imt,jmt,  &
                        bcfrec,simt,sjmt,myid,ncpux,ncpuy,mat_myid)
!     =================
!
!     interpolates monthly mean forcing fileds
!
      implicit none
      include 'pconst.h'
      include 'mpif.h'
!
      integer  pt1,pt2
      real     abc,factor
      integer imt,jmt,i,j,k,simt,sjmt
      integer   daypm(12),daymd(12)
      integer   mth,day
      integer bcfrec
      real bcf(imt,jmt,12,7)
      real bcu(imt,jmt),bcv(imt,jmt),bct(imt,jmt)
      real bcp(imt,jmt),bcs(imt,jmt),emp(imt,jmt),ddd(imt,jmt)
      
      integer myid,ncpux,ncpuy
      integer mat_myid(ncpux+2,ncpuy)
      real sbcf(simt,sjmt,7),bcf2(imt,jmt,7)
!
!      abc = real(day-daymd(mth))
!
!
!      if(abc.le.0.0) then
!        pt1    = mth -1
!        if(pt1.eq.0) pt1 = 12
!        pt2    = mth
!        factor = abc/real(daypm(pt1)) + c1
!      else
!        pt1    = mth
!        pt2    = mod(mth,12)+1
!        factor = abc/real(daypm(pt1))
!      end if
!
!
!      do j=1,jmt
!      do i=1,imt
!      bcu(i,j) = (bcf(i,j,pt2,1)-bcf(i,j,pt1,1))*factor + bcf(i,j,pt1,1)
!      bcv(i,j) = (bcf(i,j,pt2,2)-bcf(i,j,pt1,2))*factor + bcf(i,j,pt1,2)
!      bct(i,j) = (bcf(i,j,pt2,3)-bcf(i,j,pt1,3))*factor + bcf(i,j,pt1,3)
!      bcp(i,j) = (bcf(i,j,pt2,4)-bcf(i,j,pt1,4))*factor + bcf(i,j,pt1,4)
!      bcs(i,j) = (bcf(i,j,pt2,5)-bcf(i,j,pt1,5))*factor + bcf(i,j,pt1,5)
!      emp(i,j) = (bcf(i,j,pt2,6)-bcf(i,j,pt1,6))*factor + bcf(i,j,pt1,6)
!      ddd(i,j) = (bcf(i,j,pt2,7)-bcf(i,j,pt1,7))*factor + bcf(i,j,pt1,7)
!      enddo
!      enddo
!
!------new bcf input------------------------------------
      if (myid==0) then
!      if ((mth.eq.1).and.(day.eq.1)) then
!      bcfrec=1
!      end if
      read(80,rec=bcfrec) sbcf
      bcfrec=bcfrec+1
      end if
      call div_array_real3d(sbcf,bcf2,mat_myid,ncpux,ncpuy,simt,sjmt,  &
                            7,imt,jmt,myid)
      do j=1,jmt
      do i=1,imt
      bcu(i,j) = bcf2(i,j,1)
      bcv(i,j) = bcf2(i,j,2)
      bct(i,j) = bcf2(i,j,3)
      bcp(i,j) = bcf2(i,j,4)
      bcs(i,j) = bcf2(i,j,5)
      emp(i,j) = bcf2(i,j,6)
      ddd(i,j) = bcf2(i,j,7)
      enddo
      enddo
!-------------------------------------------------------
!---------------------------------------------------------------------
!     ddd:  newtonia restoring coefficient for heat flux
!---------------------------------------------------------------------
!     ddd => grav*ddd/cp (1.0e-1=unit change)
!
      do j=1,jmt
      do i=1,imt
      ddd(i,j) = grav*ddd(i,j)/3901.0*1.0e-1
      enddo
      enddo
!
      return
      end
