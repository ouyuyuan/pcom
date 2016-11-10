c
c     =================
      subroutine grdvar
c     =================
c
c     set resolution, t/u mask and the j,k-depended parameters
c
      implicit none
      include 'param.h'
      include 'pconst.h'
      include 'scalar.h'
      include 'grdvar.h'
c
      real rdx,t1
c
c     ---------------------------------------------------------------
c     horizontal resolution (-1S to 61N) 2x2
c     ---------------------------------------------------------------
      real dlam,dphi,phis
      parameter (dlam=2.0)
      parameter (dphi=2.0)
      parameter (phis=0.0)
c
c     ---------------------------------------------------------------
c     vertical resolution  (depth at the center of T grids, in cm)
c     ---------------------------------------------------------------
      data z0/
     &     1500.00,  4574.87,  7873.52, 11542.41, 15723.96, 20554.94,
     &    26165.04, 32675.40, 40197.28, 48830.88, 58664.21, 69772.14,
     &    82215.58, 96040.78,111278.89,127945.55,146040.78,165548.91,
     &   186438.81,208664.22,232164.22,256863.95,282675.41,309498.38,
     &   337221.59,365723.94,394875.72,424540.16,454574.84,484833.32/
c
c     ---------------------------------------------------------------
c     model's topography
c     ---------------------------------------------------------------
      open(72,file='topog.data',form='unformatted',status='old')
      read(72) itn
      close(72)
c
c     ---------------------------------------------------------------
c     set up T/U mask
c     ---------------------------------------------------------------
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
c
      do k=1,km
      do j=1,jmt
      do i=1,imt
      umask(i,j,k)=c0
      enddo
      enddo
      do j=2,jmm
      do i=2,imm
      umask(i,j,k)=tmask(i,j,k)*tmask(i+1,j,k)
     &            *tmask(i,j+1,k)*tmask(i+1,j+1,k)
      enddo
ccc   umask(1  ,j,k) = umask(imm,j,k)
ccc   umask(imt,j,k) = umask(2  ,j,k)
      enddo
      enddo
c
c
      do j=1,jmt
      do i=1,imt
      ivn(i,j)=0
      do k=1,km
      if(umask(i,j,k).gt.c0) ivn(i,j) = ivn(i,j) + 1
      enddo
      enddo
      enddo
c
c     ---------------------------------------------------------------
c     geopotential at bottom
c     ---------------------------------------------------------------
      do j=1,jmt
      do i=1,imt
        t1 = c0
        do k=1,itn(i,j)
        t1 = -t1-z0(k)*c2
        enddo
        phib(i,j) = t1 * grav
      enddo
      enddo
c
c     ---------------------------------------------------------------
c     calculate pn & z & dz
c     ---------------------------------------------------------------
c     pn = constant bottom pressure for PCOM; -1*phib for BCOM
c     z  = eta at the center of T grids
c     dz = thinkness of T grids
c
      call setpn
c
c
c     ---------------------------------------------------------------
c     rdz  = 1/dz
c     zu   = pn at "u" cell
c     rzu  = 1/zu
c     rdzw = 1/(layer thickness at "w" cell (cm))
c     ---------------------------------------------------------------
      do k=1,km
      rdz(k)=c1/dz(k)
      enddo
c
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
c
      do k=1,kmm1
      rdzw(k) = c1/(z0(k+1) - z0(k))
      enddo
      rdzw(km) = rdzw(kmm1)
c
c     ---------------------------------------------------------------
c     calculate latitude for starting filter
c     ---------------------------------------------------------------
      do j=1,jmt
      jstn = j
      if((phis+dphi*(j-1)).ge.65) go to 10
      enddo
10    jedn = jmm
c
      do j=1,jmt
      jeds = j
      if((phis+dphi*(j-1)).ge.-65) then
        jeds=j
        go to 20
      endif
      enddo
20    jsts = 2
c
c
c     ---------------------------------------------------------------
c     calculate parameters which are the function of grid
c     ---------------------------------------------------------------
c
      do j=1,jmt
      lat(j) = phis + dphi*(j-1)
      enddo
c
      do j=1,jmt
      cost(j) = cos((lat(j)-dphi*p5)*torad)
      cosu(j) = cos(lat(j)*torad)
      ff  (j) = sin(lat(j)*torad)*c2*omega
      enddo
c
      rdx     = c1/(radius*dlam*torad)
      rdy     = c1/(radius*dphi*torad)
c
      do j=1,jmt
      rdxt(j) = rdx / cost(j)
      rdxu(j) = rdx / cosu(j)
      rdyt(j) = rdy / cost(j)
      rdyu(j) = rdy / cosu(j)
      enddo
c
      do j=1,jmt
      sdxt(j) = rdxt(j)*rdxt(j)
      sdxu(j) = rdxu(j)*rdxu(j)*p5
      enddo
c
      do j=2,jmt
      r1a(j)   = cosu(j)  /cost(j)*(rdy**2)
      r1b(j)   = cosu(j-1)/cost(j)*(rdy**2)
      enddo
      r1a(1) = r1a(2)
      r1b(1) = r1b(2)
c
      do j=1,jmm
      r1c(j)   = p5*cost(j+1)/cosu(j)*(rdy**2)
      r1d(j)   = p5*cost(j)  /cosu(j)*(rdy**2)
      enddo
      r1c(jmt) = r1c(jmm)
      r1d(jmt) = r1d(jmm)
c
      do j=1,jmt
      t1     = sqrt(c1-cosu(j)**2)/cosu(j)
      cv1(j) = (c1-t1*t1)/(radius**2)
      cv2(j) = t1*rdxu(j)/radius
      enddo
c
c
      do j=1,jmt
      dxdyt(j) = c1/rdy/rdxt(j)
      dxdyu(j) = c1/rdy/rdxu(j)
      enddo
c
      area = c0
      do j=1,jmt
      do i=2,imm
      area = area + dxdyt(j)*tmask(i,j,1)
      enddo
      enddo
c
      return
      end
