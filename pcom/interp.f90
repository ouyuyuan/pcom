!
!     =================
!BOP
!
! !MODULE: interp.f90
! !DESCRIPTION: \input{sections/code-interp}
!
! !INTERFACE:
!
      subroutine interp(day,daymd,daypm,mth,bcf,bcu,bcv,bct,bcp,bcs,emp,ddd,imt,jmt,  &
                        simt,sjmt,myid,ncpux,ncpuy,mat_myid,monloop,yearloop,monlong)
!EOP
!-------------------------------------------------------------------------------
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
      integer monloop,yearloop,monlong
      real bcf(imt,jmt,12,7)
      real bcu(imt,jmt),bcv(imt,jmt),bct(imt,jmt)
      real bcp(imt,jmt),bcs(imt,jmt),emp(imt,jmt),ddd(imt,jmt)
      
      integer myid,ncpux,ncpuy
      integer mat_myid(ncpux+2,ncpuy)

      if (monloop==1) then
      abc = real(day-daymd(mth))
      if(abc.le.0.0) then
        pt1    = mth -1
        if(pt1.eq.0 ) pt1 = 12
        pt2    = mth
        factor = abc/real(daypm(pt1)) + c1
      else
        pt1    = mth
        pt2    = mod(mth,12)+1
        factor = abc/real(daypm(pt1))
      end if
      
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
      end if
      
      if (yearloop==1) then
      do j=1,jmt
      do i=1,imt
      bcu(i,j) = bcf(i,j,1,1)
      bcv(i,j) = bcf(i,j,1,2)
      bct(i,j) = bcf(i,j,1,3)
      bcp(i,j) = bcf(i,j,1,4)
      bcs(i,j) = bcf(i,j,1,5)
      emp(i,j) = bcf(i,j,1,6)
      ddd(i,j) = bcf(i,j,1,7)
      enddo
      enddo
      end if
      
      if (monlong==1) then
      abc = real(day-daymd(mth))
      if(abc.le.0.0) then
        pt1    = 1
        pt2    = 2
        factor = abc/real(daypm(pt1)) + c1
      else
        pt1    = 2
        pt2    = 3
        factor = abc/real(daypm(pt1))
      end if
      
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
      
      end if
      
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
      end subroutine interp
      
!
!     =================
      subroutine dismonbcf(month,runlen_mon,bcf,imt,jmt,simt,sjmt,myid,ncpux,ncpuy,mat_myid)
!     =================
!
!     distribute real monthly mean bcf run into 1 (past) 2 (now) 3 (furture) month
!
      implicit none
      include 'pconst.h'
      include 'mpif.h'
!    
      integer month,imt,jmt,simt,sjmt
      real runlen_mon,missvalue
      integer i,j,k,j2,ts
      
      real sbcu2d(simt-2,sjmt),sbcv2d(simt-2,sjmt),sbct2d(simt-2,sjmt)
      real sbcp2d(simt-2,sjmt),sbcs2d(simt-2,sjmt),semp2d(simt-2,sjmt)
      real sddd2d(simt-2,sjmt)
      real sbcf(simt,sjmt,12,7)
      
      real bcf(imt,jmt,12,7)
      
      integer myid,ncpux,ncpuy
      integer mat_myid(ncpux+2,ncpuy)
      
      
      if (myid==0) then
      
      do k=1,3
      
      ts=month+k-2
      if (ts.eq.0) ts = 1
      if (ts.eq.int(runlen_mon)) ts = int(runlen_mon)
      
      call netcdf_re_var2d_s(bcfname,bcuname,sbcu2d,simt-2,sjmt,missvalue,ts)
      call netcdf_re_var2d_s(bcfname,bcvname,sbcv2d,simt-2,sjmt,missvalue,ts)
      call netcdf_re_var2d_s(bcfname,bctname,sbct2d,simt-2,sjmt,missvalue,ts)
      call netcdf_re_var2d_s(bcfname,bcpname,sbcp2d,simt-2,sjmt,missvalue,ts)
      call netcdf_re_var2d_s(bcfname,bcsname,sbcs2d,simt-2,sjmt,missvalue,ts)
      call netcdf_re_var2d_s(bcfname,empname,semp2d,simt-2,sjmt,missvalue,ts)
      call netcdf_re_var2d_s(bcfname,dddname,sddd2d,simt-2,sjmt,missvalue,ts)
      
        do j=1,sjmt
          do i=1,simt-2
            sbcf(i+1,j,k,1)=sbcu2d(i,j)
            sbcf(i+1,j,k,2)=sbcv2d(i,j)
            sbcf(i+1,j,k,3)=sbct2d(i,j)
            sbcf(i+1,j,k,4)=sbcp2d(i,j)
            sbcf(i+1,j,k,5)=sbcs2d(i,j)
            sbcf(i+1,j,k,6)=semp2d(i,j)
            sbcf(i+1,j,k,7)=sddd2d(i,j)
          end do
          sbcf(1,j,k,1)=sbcu2d(simt-2,j)
          sbcf(1,j,k,2)=sbcv2d(simt-2,j)
          sbcf(1,j,k,3)=sbct2d(simt-2,j)
          sbcf(1,j,k,4)=sbcp2d(simt-2,j)
          sbcf(1,j,k,5)=sbcs2d(simt-2,j)
          sbcf(1,j,k,6)=semp2d(simt-2,j)
          sbcf(1,j,k,7)=sddd2d(simt-2,j)
          sbcf(simt,j,k,1)=sbcu2d(1,j)
          sbcf(simt,j,k,2)=sbcv2d(1,j)
          sbcf(simt,j,k,3)=sbct2d(1,j)
          sbcf(simt,j,k,4)=sbcp2d(1,j)
          sbcf(simt,j,k,5)=sbcs2d(1,j)
          sbcf(simt,j,k,6)=semp2d(1,j)
          sbcf(simt,j,k,7)=sddd2d(1,j)
        end do
      end do
      
      do j2=1,7
        do k=1,3
          do j=1,sjmt
            do i=1,simt
              if (sbcf(i,j,k,j2)==missvalue) then
                sbcf(i,j,k,j2)=0
              end if
            end do
          end do
        end do
      end do
      end if !(myid==0)

      call div_array_real4d(sbcf,bcf,mat_myid,ncpux,ncpuy,simt,sjmt,  &
                            12,7,imt,jmt,myid)

      return
      end subroutine dismonbcf

