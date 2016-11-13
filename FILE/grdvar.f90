!
!     =================
      subroutine grdvar(z,z0,dz0,pn,itn,ivn,tmask,umask,phib,dz,rdz,rdzw,zu,rzu,   &
                        jstn,jedn,jeds,jsts,lat,cost,cosu,ff,rdxt,rdxu,rdyt,rdyu,  &
                        sdxt,sdxu,r1a,r1b,r1c,r1d,cv1,cv2,dxdyt,dxdyu,area,rdy,    &
                        decibar,imt,jmt,km,dlam,dphi,phis,imm,jmm,kmp1,kmm1,unesco, &
                        myid,ncpux,ncpuy,west,east,north,south,mat_myid,simt,sjmt,  &
                        smth_start_nlat,smth_start_slat,boussinesq)
!     =================
!
!     set resolution, t/u mask and the j,k-depended parameters
!
      implicit none
      include 'pconst.h'
      include 'mpif.h'
!
      integer imt,jmt,km,imm,jmm,kmp1,kmm1,i,j,k,n,m,unesco,boussinesq
      real dlam,dphi,phis
      real rdx,rdy,t1,area,decibar,smth_start_nlat,smth_start_slat
      real lat(jmt),cost(jmt),cosu(jmt),ff(jmt)
      real rdxt(jmt),rdxu(jmt),rdyt(jmt),rdyu(jmt)
      real dxdyt(jmt),dxdyu(jmt)
      real cv1(jmt),cv2(jmt),sdxt(jmt),sdxu(jmt),r1a(jmt),r1b(jmt),r1c(jmt),r1d(jmt)

      real z0(km),dz(km),rdz(km),rdzw(km),dz0(km),z(km),pn(imt,jmt)
      real phib(imt,jmt),rzu(imt,jmt),zu(imt,jmt)
      real tmask(imt,jmt,km),umask(imt,jmt,km)
      integer itn(imt,jmt),ivn(imt,jmt)
      integer jstn,jedn,jsts,jeds

      integer myid,ncpux,ncpuy,west,east,north,south,simt,sjmt
      integer mat_myid(ncpux+2,ncpuy),sitn(simt,sjmt)
      real slat(sjmt),scost(sjmt),scosu(sjmt),sff(sjmt)
      real srdxt(sjmt),srdxu(sjmt),srdyt(sjmt),srdyu(sjmt)
      real sdxdyt(sjmt),sdxdyu(sjmt)
      real scv1(sjmt),scv2(sjmt),ssdxt(sjmt),ssdxu(sjmt),sr1a(sjmt)
      real sr1b(sjmt),sr1c(sjmt),sr1d(sjmt)
      real stmask(simt,sjmt,km)

!
!     ---------------------------------------------------------------
!     model's topography
!     ---------------------------------------------------------------
      if (myid==0) then
      open(72,file='topog.data',form='unformatted',status='old')
      read(72) sitn
      close(72)
      end if
      
      call div_array_int2d(sitn,itn,mat_myid,ncpux,ncpuy,simt,sjmt,  &
                           imt,jmt,myid)
!
!     ---------------------------------------------------------------
!     set up T/U mask
!     ---------------------------------------------------------------
      do k=1,km
      do j=1,jmt
      do i=1,imt
      if(k.le.itn(i,j))then
        tmask(i,j,k)=c1
      else
        tmask(i,j,k)=c0
      endif
      enddo
      enddo
      enddo
!
      do k=1,km
      do j=1,jmt
      do i=1,imt
      umask(i,j,k)=c0
      enddo
      enddo
      do j=2,jmm
      do i=2,imm
      umask(i,j,k)=tmask(i,j,k)*tmask(i+1,j,k)  &
                  *tmask(i,j+1,k)*tmask(i+1,j+1,k)
      enddo
!      umask(1  ,j,k) = umask(imm,j,k)
!      umask(imt,j,k) = umask(2  ,j,k)
      enddo
      enddo
      call swap_array_real3d(umask,imt,jmt,km,west,east,north,south)
!
!
      do j=1,jmt
      do i=1,imt
      ivn(i,j)=0
      do k=1,km
      if(umask(i,j,k).gt.c0) ivn(i,j) = ivn(i,j) + 1
      enddo
      enddo
      enddo
!
!     ---------------------------------------------------------------
!     geopotential at bottom
!     ---------------------------------------------------------------
      do j=1,jmt
      do i=1,imt
        t1 = c0
        do k=1,itn(i,j)
        t1 = -t1-z0(k)*c2
        enddo
        phib(i,j) = t1 * grav
      enddo
      enddo
!
!     ---------------------------------------------------------------
!     calculate pn & z & dz
!     ---------------------------------------------------------------
!     pn = constant bottom pressure for PCOM; -1*phib for BCOM
!     z  = eta at the center of T grids
!     dz = thinkness of T grids
!
      call setpn(dz0,z0,dz,z,pn,imt,jmt,km,kmp1,itn,decibar,unesco,boussinesq,   &
                 phib,myid,ncpux,ncpuy,mat_myid,west,east,north,south)
!
!
!     ---------------------------------------------------------------
!     rdz  = 1/dz
!     zu   = pn at "u" cell
!     rzu  = 1/zu
!     rdzw = 1/(layer thickness at "w" cell (cm))
!     ---------------------------------------------------------------
      do k=1,km
      rdz(k)=c1/dz(k)
      enddo
!
      do j=1,jmt
      do i=1,imt
        zu (i,j) = c0
        rzu(i,j) = c0
        do k=1,ivn(i,j)
        zu(i,j) = zu(i,j) + dz(k)
        enddo
        if(zu(i,j).gt.c0) rzu(i,j) = c1/zu(i,j)
      enddo
      enddo
      call swap_array_real2d(zu,imt,jmt,west,east,north,south)
      call swap_array_real2d(rzu,imt,jmt,west,east,north,south)
!
      do k=1,kmm1
      rdzw(k) = c1/(z0(k+1) - z0(k))
      enddo
      rdzw(km) = rdzw(kmm1)
!
!     ---------------------------------------------------------------
!     calculate latitude for starting filter
!     ---------------------------------------------------------------
!      do j=1,ncpuy
!      do i=2,ncpux+1
!      if (mat_myid(i,j)==myid) then
!      i_id=i
!      j_id=j
!      end if
!      end do
!      end do
!      J-1=j-2+(jmt-2)*(j_id-1)
      call gath_array_real3d(stmask,tmask,mat_myid,ncpux,ncpuy,simt,sjmt,  &
                             km,imt,jmt,myid)
      
      rdx     = c1/(radius*dlam*torad)
      rdy     = c1/(radius*dphi*torad)
      
      if (myid==0) then
!
!
!     ---------------------------------------------------------------
!     calculate parameters which are the function of grid
!     ---------------------------------------------------------------
!
      do j=1,sjmt
      slat(j) = phis + dphi*(j-1)
      enddo
!
      do j=1,sjmt
      scost(j) = cos((slat(j)-dphi*p5)*torad)
      scosu(j) = cos(slat(j)*torad)
      sff  (j) = sin(slat(j)*torad)*c2*omega
      enddo
!
      
!
      do j=1,sjmt
      srdxt(j) = rdx / scost(j)
      srdxu(j) = rdx / scosu(j)
      srdyt(j) = rdy / scost(j)
      srdyu(j) = rdy / scosu(j)
      enddo
!
      do j=1,sjmt
      ssdxt(j) = srdxt(j)*srdxt(j)
      ssdxu(j) = srdxu(j)*srdxu(j)*p5
      enddo
!
      do j=2,sjmt
      sr1a(j)   = scosu(j)  /scost(j)*(rdy**2)
      sr1b(j)   = scosu(j-1)/scost(j)*(rdy**2)
      enddo
      sr1a(1) = sr1a(2)
      sr1b(1) = sr1b(2)
!
      do j=1,sjmt-1
      sr1c(j)   = p5*scost(j+1)/scosu(j)*(rdy**2)
      sr1d(j)   = p5*scost(j)  /scosu(j)*(rdy**2)
      enddo
      sr1c(sjmt) = sr1c(sjmt-1)
      sr1d(sjmt) = sr1d(sjmt-1)
!
      do j=1,sjmt
      t1     = sqrt(c1-scosu(j)**2)/scosu(j)
      scv1(j) = (c1-t1*t1)/(radius**2)
      scv2(j) = t1*srdxu(j)/radius
      enddo
!
!
      do j=1,sjmt
      sdxdyt(j) = c1/rdy/srdxt(j)
      sdxdyu(j) = c1/rdy/srdxu(j)
      enddo
!
      area = c0
      do j=1,sjmt
      do i=2,simt-1
      area = area + sdxdyt(j)*stmask(i,j,1)
      enddo
      enddo
      end if
      
      call dis_var_real(area,mat_myid,ncpux,ncpuy,myid)
      call div_array_real1d(slat,lat,mat_myid,ncpux,ncpuy,sjmt,jmt,myid)
      call div_array_real1d(scost,cost,mat_myid,ncpux,ncpuy,sjmt,jmt,myid)
      call div_array_real1d(scosu,cosu,mat_myid,ncpux,ncpuy,sjmt,jmt,myid)
      call div_array_real1d(sff,ff,mat_myid,ncpux,ncpuy,sjmt,jmt,myid)
      call div_array_real1d(srdxt,rdxt,mat_myid,ncpux,ncpuy,sjmt,jmt,myid)
      call div_array_real1d(srdxu,rdxu,mat_myid,ncpux,ncpuy,sjmt,jmt,myid)
      call div_array_real1d(srdyt,rdyt,mat_myid,ncpux,ncpuy,sjmt,jmt,myid)
      call div_array_real1d(srdyu,rdyu,mat_myid,ncpux,ncpuy,sjmt,jmt,myid)
      call div_array_real1d(srdyt,rdyt,mat_myid,ncpux,ncpuy,sjmt,jmt,myid)
      call div_array_real1d(ssdxt,sdxt,mat_myid,ncpux,ncpuy,sjmt,jmt,myid)
      call div_array_real1d(ssdxu,sdxu,mat_myid,ncpux,ncpuy,sjmt,jmt,myid)
      call div_array_real1d(sr1a,r1a,mat_myid,ncpux,ncpuy,sjmt,jmt,myid)
      call div_array_real1d(sr1b,r1b,mat_myid,ncpux,ncpuy,sjmt,jmt,myid)
      call div_array_real1d(sr1c,r1c,mat_myid,ncpux,ncpuy,sjmt,jmt,myid)
      call div_array_real1d(sr1d,r1d,mat_myid,ncpux,ncpuy,sjmt,jmt,myid)
      call div_array_real1d(scv1,cv1,mat_myid,ncpux,ncpuy,sjmt,jmt,myid)
      call div_array_real1d(scv2,cv2,mat_myid,ncpux,ncpuy,sjmt,jmt,myid)
      call div_array_real1d(sdxdyt,dxdyt,mat_myid,ncpux,ncpuy,sjmt,jmt,myid)
      call div_array_real1d(sdxdyu,dxdyu,mat_myid,ncpux,ncpuy,sjmt,jmt,myid)

      do j=2,jmt
      jstn = j
      if(lat(j).ge.smth_start_nlat) go to 10
      end do
10    jedn = jmt-1
!
      do j=1,jmt-1
      jeds = j
      if(lat(j).ge.smth_start_slat) go to 20
      end do
20    jsts = 2
      
      return
      end
