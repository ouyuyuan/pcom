c     =============================================
      subroutine smth3(x,mask,jstn,jedn,jsts,jeds)
c     =============================================
c     1-d zonal smoother
c
      implicit none
      include 'param.h'
      include 'pconst.h'
c
      real x(imt,jmt,km,2),xs(imt),mask(imt,jmt,km)
      integer jstn,jedn,jsts,jeds
c
      do n=1,2
      do k=1,km
c
      do j=jstn,jedn
      do i=2,imm
      xs(i)=x(i,j,k,n)*mask(i,j,k)
      enddo
      xs(1)=xs(imm)
      xs(imt)=xs(2)
      do i=2,imm
      x(i,j,k,n)=(p5*xs(i)+p25*(xs(i-1)+xs(i+1)))*mask(i,j,k)
      enddo
      x(1  ,j,k,n) = x(imm,j,k,n)
      x(imt,j,k,n) = x(2  ,j,k,n)
      enddo
c
      do j=jsts,jeds
      do i=2,imm
      xs(i)=x(i,j,k,n)*mask(i,j,k)
      enddo
      xs(1)=xs(imm)
      xs(imt)=xs(2)
      do i=2,imm
      x(i,j,k,n)=(p5*xs(i)+p25*(xs(i-1)+xs(i+1)))*mask(i,j,k)
      enddo
      x(1  ,j,k,n) = x(imm,j,k,n)
      x(imt,j,k,n) = x(2  ,j,k,n)
      enddo
c
      enddo
      enddo
c
      return
      end
c
c     =================================
      subroutine smth2(x,mask,jstn,jedn,jsts,jeds)
c     =================================
c     1-d zonal smoother
c
      implicit none
      include 'param.h'
      include 'pconst.h'
c
      real x(imt,jmt,2),xs(imt),mask(imt,jmt,km)
      integer jstn,jedn,jsts,jeds
c
      do n=1,2
      do j=jstn,jedn
      do i=2,imm
      xs(i)=x(i,j,n)*mask(i,j,1)
      enddo
      xs(1)=xs(imm)
      xs(imt)=xs(2)
      do i=2,imm
      x(i,j,n)=(p5*xs(i)+p25*(xs(i-1)+xs(i+1)))*mask(i,j,1)
      enddo
      x(1  ,j,n) = x(imm,j,n)
      x(imt,j,n) = x(2  ,j,n)
      enddo
c
      do j=jsts,jeds
      do i=2,imm
      xs(i)=x(i,j,n)*mask(i,j,1)
      enddo
      xs(1)=xs(imm)
      xs(imt)=xs(2)
      do i=2,imm
      x(i,j,n)=(p5*xs(i)+p25*(xs(i-1)+xs(i+1)))*mask(i,j,1)
      enddo
      x(1  ,j,n) = x(imm,j,n)
      x(imt,j,n) = x(2  ,j,n)
      enddo
c
      enddo
c
      return
      end
c
c
c     =================================
      subroutine smths(x,mask,jst,jed)
c     =================================
c     1-d zonal smoother
c
      implicit none
      include 'param.h'
      include 'pconst.h'
c
      real x(imt,jmt),xs(imt),mask(imt,jmt,km)
      integer jst,jed
c
      do j=jst,jed
      do i=2,imm
      xs(i)=x(i,j)*mask(i,j,1)
      enddo
      xs(1)=xs(imm)
      xs(imt)=xs(2)
      do i=2,imm
      x(i,j)=(p5*xs(i)+p25*(xs(i-1)+xs(i+1)))*mask(i,j,1)
      enddo
      x(1  ,j) = x(imm,j)
      x(imt,j) = x(2  ,j)
      enddo
c
      return
      end
