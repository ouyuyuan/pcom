!     =============================================
      subroutine smth3(x,mask,jstn,jedn,jsts,jeds,imt,jmt,km,fcof)
!     =============================================
!     1-d zonal smoother
!
      implicit none
!      include 'param.h'
      include 'pconst.h'
!
      integer i,j,k,n,imt,jmt,km
      real x(imt,jmt,km,2),xs(imt),masks(imt),mask(imt,jmt,km),fcof
      integer jstn,jedn,jsts,jeds
!
      do n=1,2
         do k=1,km
!
            do j=jstn,jedn
               do i=1,imt
                  xs(i)=x(i,j,k,n)
                  masks(i)=mask(i,j,k)
               end do
               do i=2,imt-1
                  if ((masks(i-1)==1).and.(masks(i+1)==1)) then
                     x(i,j,k,n)=((c1-c2*fcof)*xs(i)+fcof*(xs(i-1)+xs(i+1)))*masks(i)
                  end if
                  if ((masks(i-1)==1).and.(masks(i+1)==0)) then
                     x(i,j,k,n)=((c1-fcof)*xs(i)+fcof*xs(i-1))*masks(i)
                  end if
                  if ((masks(i-1)==0).and.(masks(i+1)==1)) then
                     x(i,j,k,n)=((c1-fcof)*xs(i)+fcof*xs(i+1))*masks(i)
                  end if
                  if ((masks(i-1)==0).and.(masks(i+1)==0)) then
                     x(i,j,k,n)=xs(i)*masks(i)
                  end if
               end do
            end do
!
            do j=jsts,jeds
               do i=1,imt
                  xs(i)=x(i,j,k,n)
                  masks(i)=mask(i,j,k)
               end do
               do i=2,imt-1
                  if ((masks(i-1)==1).and.(masks(i+1)==1)) then
                     x(i,j,k,n)=((c1-c2*fcof)*xs(i)+fcof*(xs(i-1)+xs(i+1)))*masks(i)
                  end if
                  if ((masks(i-1)==1).and.(masks(i+1)==0)) then
                     x(i,j,k,n)=((c1-fcof)*xs(i)+fcof*xs(i-1))*masks(i)
                  end if
                  if ((masks(i-1)==0).and.(masks(i+1)==1)) then
                     x(i,j,k,n)=((c1-fcof)*xs(i)+fcof*xs(i+1))*masks(i)
                  end if
                  if ((masks(i-1)==0).and.(masks(i+1)==0)) then
                     x(i,j,k,n)=xs(i)*masks(i)
                  end if
               end do
            end do
!
         end do
      end do
!
      return
      end subroutine smth3
!
!     =================================
      subroutine smth2(x,mask,jstn,jedn,jsts,jeds,imt,jmt,km,fcof)
!     =================================
!     1-d zonal smoother
!
      implicit none
!      include 'param.h'
      include 'pconst.h'
!
      integer imt,jmt,km,i,j,n
      real x(imt,jmt,2),xs(imt),masks(imt),mask(imt,jmt,km),fcof
      integer jstn,jedn,jsts,jeds
!
      do n=1,2
         do j=jstn,jedn
            do i=1,imt
               xs(i)=x(i,j,n)
               masks(i)=mask(i,j,1)
            end do
            do i=2,imt-1
               if ((masks(i-1)==1).and.(masks(i+1)==1)) then
                  x(i,j,n)=((c1-c2*fcof)*xs(i)+fcof*(xs(i-1)+xs(i+1)))*masks(i)
               end if
               if ((masks(i-1)==1).and.(masks(i+1)==0)) then
                  x(i,j,n)=((c1-fcof)*xs(i)+fcof*xs(i-1))*masks(i)
               end if
               if ((masks(i-1)==0).and.(masks(i+1)==1)) then
                  x(i,j,n)=((c1-fcof)*xs(i)+fcof*xs(i+1))*masks(i)
               end if
               if ((masks(i-1)==0).and.(masks(i+1)==0)) then
                  x(i,j,n)=xs(i)*masks(i)
               end if
            end do
         end do
!
         do j=jsts,jeds
            do i=1,imt
               xs(i)=x(i,j,n)
               masks(i)=mask(i,j,1)
            end do
            do i=2,imt-1
               if ((masks(i-1)==1).and.(masks(i+1)==1)) then
                  x(i,j,n)=((c1-c2*fcof)*xs(i)+fcof*(xs(i-1)+xs(i+1)))*masks(i)
               end if
               if ((masks(i-1)==1).and.(masks(i+1)==0)) then
                  x(i,j,n)=((c1-fcof)*xs(i)+fcof*xs(i-1))*masks(i)
               end if
               if ((masks(i-1)==0).and.(masks(i+1)==1)) then
                  x(i,j,n)=((c1-fcof)*xs(i)+fcof*xs(i+1))*masks(i)
               end if
               if ((masks(i-1)==0).and.(masks(i+1)==0)) then
                  x(i,j,n)=xs(i)*masks(i)
               end if
            end do
         end do
!
      end do
!
      return
      end subroutine smth2
!
!
!     =================================
      subroutine smths(x,mask,jst,jed,imt,jmt,km,fcof)
!     =================================
!     1-d zonal smoother
!
      implicit none
!      include 'param.h'
      include 'pconst.h'
!
      integer i,j,imt,jmt,km,tnum
      real x(imt,jmt),xs(imt),mask(imt,jmt,km),masks(imt),fcof
      integer jst,jed
!
      do j=jst,jed
         do i=1,imt
            xs(i)=x(i,j)
            masks(i)=mask(i,j,1)
         end do
         do i=2,imt-1
            if ((masks(i-1)==1).and.(masks(i+1)==1)) then
               x(i,j)=((c1-c2*fcof)*xs(i)+fcof*(xs(i-1)+xs(i+1)))*masks(i)
            end if
            if ((masks(i-1)==1).and.(masks(i+1)==0)) then
               x(i,j)=((c1-fcof)*xs(i)+fcof*xs(i-1))*masks(i)
            end if
            if ((masks(i-1)==0).and.(masks(i+1)==1)) then
               x(i,j)=((c1-fcof)*xs(i)+fcof*xs(i+1))*masks(i)
            end if
            if ((masks(i-1)==0).and.(masks(i+1)==0)) then
               x(i,j)=xs(i)*masks(i)
            end if
         end do
      end do
!
      return
      end subroutine smths
