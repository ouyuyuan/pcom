!
!     ================
      subroutine setpn(dz0,z0,dz,z,pn,imt,jmt,km,kmp1,itn,decibar,unesco,boussinesq,   &
                       phib,myid,ncpux,ncpuy,mat_myid,west,east,north,south)
!     ================
!     calculate pn, z & dz
!
!     PCOM  pn = constant bottom pressure   (dy/cm2)
!     BCOM  pn = grav * thickness           (cm*cm/s2)
!     z = eta; dz=1/z
!
      implicit none
      include 'pconst.h'
      include 'mpif.h'
!
      integer imt,jmt,km,kmp1,i,j,k,itn(imt,jmt),unesco,boussinesq
      real p0,dens,undens,pres,decibar,missvalue
      real pre(kmp1),denz(km)
      real t30(km),s30(km)
      real z0(km),dz0(km),z(km),dz(km),pn(imt,jmt),phib(imt,jmt)
      
      integer myid,ncpux,ncpuy,west,east,north,south
      integer mat_myid(ncpux+2,ncpuy)
!
!---------------------------------------------------------------------
!     input global averaged TS stratification
!---------------------------------------------------------------------
      if (myid==0) then
      call netcdf_read_var1d(ncname,t30name,t30,km,missvalue)
      call netcdf_read_var1d(ncname,s30name,s30,km,missvalue)
      end if
      call dis_var_real1d(t30,km,mat_myid,ncpux,ncpuy,myid)
      call dis_var_real1d(s30,km,mat_myid,ncpux,ncpuy,myid)
!
!
!---------------------------------------------------------------------
!     calculate depth of each layer (in cm)
!---------------------------------------------------------------------
      dz0(1) = z0(1)*c2
      do k=2,km
      dz0(k) = (z0(k)-z0(k-1))*c2 - dz0(k-1)
      enddo
!
!---------------------------------------------------------------------
!     calculate pressure at the bottom of each layer
!---------------------------------------------------------------------
!     set an initial density
      do k=1,km
      denz(k) = c1
      enddo
!
      pres = c0
!
25    pre(1)   = c0
      do k=1,km
      pre(k+1) = pre(k) + denz(k)*grav*dz0(k)
      enddo
!
      if (unesco==1) then
      do k=1,km
      p0      = (pre(k)+pre(k+1))*p5*decibar
      denz(k) = undens(t30(k),s30(k),p0)
      enddo
      else
      do k=1,km
      p0      = (pre(k)+pre(k+1))*p5*decibar
      denz(k) = dens(t30(k),s30(k),p0)
      enddo
      end if
!
      if(abs(pre(kmp1)-pres).gt.c1em4)then
       pres = pre(kmp1)
       goto 25
      endif
!
!---------------------------------------------------------------------
!     calculate z & dz
!---------------------------------------------------------------------
      if (boussinesq==1) then
      do k=1,km
        z(k) = z0(k) * grav
       dz(k) = dz0(k) * grav
      enddo
      else
      do k=1,km
        z(k) = (pre(k+1)+pre(k))*p5
       dz(k) =  pre(k+1)-pre(k)
      enddo
      end if
!
!---------------------------------------------------------------------
!     calculate pn
!---------------------------------------------------------------------
      if (boussinesq==1) then
      do j=1,jmt
      do i=1,imt
      if(itn(i,j).eq.0) then
       pn(i,j) = c1
      else
       pn(i,j) = abs(phib(i,j))
      endif
      enddo
      enddo
      
      else
      
      do j=1,jmt
      do i=1,imt
      if(itn(i,j).eq.0) then
       pn(i,j) = c1
      else
       pn(i,j) = pre(itn(i,j)+1)
      endif
      enddo
      enddo
      end if
      
      call swap_array_real2d(pn,imt,jmt,west,east,north,south)

      return
      end
