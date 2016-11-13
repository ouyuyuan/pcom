!     =============================================
      subroutine smth3(x,mask,jstn,jedn,jsts,jeds,imt,jmt,km)
!     =============================================
!     1-d zonal smoother
!
      implicit none
!      include 'param.h'
      include 'pconst.h'
!
      integer i,j,k,n,imt,jmt,km
      real x(imt,jmt,km,2),xs(imt),mask(imt,jmt,km)
      integer jstn,jedn,jsts,jeds
!
      do n=1,2
      do k=1,km
!
      do j=jstn,jedn
      do i=1,imt
      xs(i)=x(i,j,k,n)*mask(i,j,k)
      enddo
      do i=2,imt-1
      x(i,j,k,n)=(p5*xs(i)+p25*(xs(i-1)+xs(i+1)))*mask(i,j,k)
      enddo
      enddo
!
      do j=jsts,jeds
      do i=1,imt
      xs(i)=x(i,j,k,n)*mask(i,j,k)
      enddo
      do i=2,imt-1
      x(i,j,k,n)=(p5*xs(i)+p25*(xs(i-1)+xs(i+1)))*mask(i,j,k)
      enddo
      enddo
!
      enddo
      enddo
!
      return
      end
!
!     =================================
      subroutine smth2(x,mask,jstn,jedn,jsts,jeds,imt,jmt,km)
!     =================================
!     1-d zonal smoother
!
      implicit none
!      include 'param.h'
      include 'pconst.h'
!
      integer imt,jmt,km,i,j,n
      real x(imt,jmt,2),xs(imt),mask(imt,jmt,km)
      integer jstn,jedn,jsts,jeds
!
      do n=1,2
      do j=jstn,jedn
      do i=1,imt
      xs(i)=x(i,j,n)*mask(i,j,1)
      enddo
      do i=2,imt-1
      x(i,j,n)=(p5*xs(i)+p25*(xs(i-1)+xs(i+1)))*mask(i,j,1)
      enddo
      enddo
!
      do j=jsts,jeds
      do i=1,imt
      xs(i)=x(i,j,n)*mask(i,j,1)
      enddo
      do i=2,imt-1
      x(i,j,n)=(p5*xs(i)+p25*(xs(i-1)+xs(i+1)))*mask(i,j,1)
      enddo
      enddo
!
      enddo
!
      return
      end
!
!
!     =================================
      subroutine smths(x,mask,jst,jed,imt,jmt,km)
!     =================================
!     1-d zonal smoother
!
      implicit none
!      include 'param.h'
      include 'pconst.h'
!
      integer i,j,imt,jmt,km
      real x(imt,jmt),xs(imt),mask(imt,jmt,km)
      integer jst,jed
!
      do j=jst,jed
      do i=1,imt
      xs(i)=x(i,j)*mask(i,j,1)
      enddo
      do i=2,imt-1
      x(i,j)=(p5*xs(i)+p25*(xs(i-1)+xs(i+1)))*mask(i,j,1)
      enddo
      enddo
!
      return
      end
