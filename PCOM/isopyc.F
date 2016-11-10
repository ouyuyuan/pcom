#ifdef gm90
c
c     =================
      subroutine isopyi
c     =================
c
c=======================================================================
c
c       Initialization for isopycnal mixing scheme
c
c       Redi/Cox version + Gent_McWilliams version
c
cifdef isopycxixspatialvar
c       dciso1 = isopycnal tracer diffusivity coeffs modified based
c              on the slopes of the isopycnal surfaces on the east face
c              of "T" cells.
c       dciso2 = isopycnal tracer diffusivity coeffs modified based
c              on the slopes of the isopycnal surfaces on the north face
c              of "T" cells.
c       dslope = half length of the interval in which "ahisop" changes
c              with a steep slope from about 0.9*"ahisop" to about
c              0.1*"ahisop"
c       slopec = slope at which "ahisop" is equal to half of its 
c              original  value
cendif
c
c       dptlim = depth limits for the reference pressure levels (in cm). 
c
c=======================================================================
c
      implicit none
      include 'param.h'
      include 'pconst.h'
      include 'grdvar.h'
      include 'isopyc.h'
c
      real dptmid,t1,t2
c
c-----------------------------------------------------------------------
c     USER INPUT
c     initialize variables (all mixing units are cm**2/sec.)
c-----------------------------------------------------------------------
c
      slmxr  = 100.0
      ahisop = 1.e7
c
c     define the isopycnal thickness diffusion coefficient
c
      athkdf = 1.0e7
c
c     reference pressure level intervals are defined (see "isopyc.h").
c     "dptlim" must have "nrpl+1" elements (see "param.h"). also,
c     "dptlim" are in cm.
c
c
c     REMARK: the first and the last elements of "dptlim" must be the 
c             depth at the top (0cm) and the maximum bottom depth,
c             respectively. Also, the elements of "dptlim" must be in
c             increasing order.
c
c     dptlim(1) = 0.0e2
c     dptlim(2) = 1000.0e2
c     dptlim(3) = 2000.0e2
c     dptlim(4) = 3000.0e2
c     dptlim(5) = 4000.0e2
c     dptlim(6) = 5000.0e2
c
c
c-----------------------------------------------------------------------
c     determine the isopycnal reference pressure levels for the "t"
c     grid point levels, using the depths at the "t" grid points as the
c     reference depth (pressure)
c-----------------------------------------------------------------------
c
c     do k=1,km
c       do m=2,nrpl+1
c         if (z0(k).gt.dptlim(m-1) .and. z0(k).le.dptlim(m)) then
c           kisrpl(k) = m-1
c           go to 101
c         endif
c       enddo
101     continue
c
c       if (kisrpl(k) .lt. 1 .or. kisrpl(k) .gt. nrpl) then
c         write (6,9100) kisrpl(k), k
c         stop 9999
c       endif
c     enddo
9100  format (/,' =>Error: kisrpl is ',i3,' at k ',1x,i3,' in isopyc.F')
c
c-----------------------------------------------------------------------
c     the indices used in isopycnal mixing indicating the location of
c     the reference pressure levels in the 20-level table of polynomial
c     expansion variables are computed
c
c     REMARK: because the polynomial expansion coefficients are
c             functions of the reference potential temperature and
c             salinity profiles, at the reference pressure level
c             the corresponding potential temperature and salinity
c             values will be used.
c-----------------------------------------------------------------------
c
c     do m=1,nrpl
c       krplin(m) = 0
c     enddo
c
c     do m=2,nrpl+1
c       dptmid = 0.5*(dptlim(m-1)+dptlim(m))
c       if (dptmid .le. z0(1)) then
c         krplin(m-1) = 1
c       elseif (dptmid .gt. z0(km)) then
c         krplin(m-1) = km
c       elseif (dptmid.gt.z0(1) .and. dptmid.le.z0(km)) then
c         do k=2,km
c           if (z0(k) .ge. dptmid) then
c             t1 = z0(k)-dptmid
c             t2 = dptmid-z0(k-1)
c             if (t1 .gt. t2) then
c               krplin(m-1) = k-1
c             else
c               krplin(m-1) = k
c             endif
c             go to 102 
c           endif
c         enddo
102       continue
c       endif
c       if (krplin(m-1) .lt. 1 .or. krplin(m-1) .gt. km) then
c         write (6,9110) krplin(m-1),m-1
c       endif
c     enddo
c
ccc   write (6,96) (kisrpl(k),k=1,km)
ccc   write (6,97) (krplin(m),m=1,nrpl)
c
c-----------------------------------------------------------------------
c     the isopycnal diffusion coefficient may be a function of depth. in 
c     the default configuration, the isopycnal diffusion coefficient is 
c     a constant: "fzisop", which multiplies "ahisop", is set to unity. 
c     if "ahisop" varies in the vertical, "fzisop" should contain this 
c     variation. the value of "ahisop" should be adjusted accordingly.
c-----------------------------------------------------------------------
c
      do k=1,km
        fzisop(k) = c1
      enddo
c
ccc   write (6,98)
ccc   write (6,99) (fzisop(k),k=1,km)
c
  96  format (/,' isopycnal reference pressure levels (kisrpl) = ',
     &        20(1x,i4))
  97  format (/,' reference pressure level indices (krplin) = ',
     &        20(1x,i4))
  98  format (/,' vertical variation of "ahisop" (fzisop) = ')
  99  format (5(1x,e12.6))
9110  format (/,' =>Error: krplin is ',i3,' at m ',1x,i3)
c
c
c-----------------------------------------------------------------------
c     set reference depth for calculation of rhoi  (m)
c-----------------------------------------------------------------------
      do k=1,km
      kref(k) = z0(k)*0.01
      enddo
      do k=1,km
      rdz0(k) = c1/dz0(k)
      enddo
c
c
c-----------------------------------------------------------------------
c     initialize arrays
c-----------------------------------------------------------------------
c
      do j=1,jmt
	do k=1,km
          do i=1,imt
            K1(i,k,j,3)       = c0
            K3(i,k,j,1)       = c0
            K3(i,k,j,2)       = c0
            K3(i,k,j,3)       = c0
            K2(i,k,j,3)       = c0
            adv_vetiso(i,k,j) = c0
            adv_vntiso(i,k,j) = c0
          enddo
        enddo
      enddo
c
      do j=1,jmt
	do k=1,km
          do i=1,imt
	    rhoi(i,k,j,xup) = c1
	    rhoi(i,k,j,xmd) = c1
	    rhoi(i,k,j,xlo) = c1
	  enddo
	enddo
      enddo
c
      do m=1,3
        do j=1,jmt
	  do k=1,km+1
            do i=1,imt
	      e(i,k,j,m) = c0
	    enddo
	  enddo
        enddo
      enddo
c
      do j=1,jmt
        do k=0,km
          do i=1,imt
            adv_vbtiso(i,k,j) = c0
	  enddo
        enddo
      enddo
c
      return
      end
c
c     =================
      subroutine isopyc
c     =================
c
c=======================================================================
c
c     Compute the isopycnal mixing tensor components and the
c     isopycnal advection velocities which parameterize the effect
c     of eddies on the isopycnals.
c
c
c     Mixing tensor "K" is ...
c
c          | 1.0            K1(,,,2)          K1(,,,3) | 
c          |                                           |
c     K =  | K2(,,,1)        1.0              K2(,,,3) |
c          |                                           |
c          | K3(,,,1)       K3(,,,2)          K3(,,,3) |
c
c     where K1(,,,2) and K2(,,,1) are set to 0.0 (neglected)
c
c
c     output:
c       rhoi = density at tau-1 referenced to pressure levels
c       K1   = tensor components (1,2), (1,3) centered on east face
c              of "T" cells
c       K2   = tensor components (2,1), (2,3) centered on north face
c              of "T" cells
c       K3   = tensor components (3,1), (3,2), (3,3) centered on
c              bottom face of "T" cells  
cifdef gent_mcwilliams
c       adv_vetiso = isopycnal advective vel on east face of "T" cell
c       adv_vntiso = isopycnal advective vel on north face of "T" cell
c               (Note: this includes the cosine factor as in "adv_vnt")
c       adv_vbtiso = isopycnal advective vel on bottom face of "T" cell       
cendif
c
c=======================================================================
c
c
      implicit none
      include 'param.h'
      include 'pconst.h'
      include 'grdvar.h'
      include 'prog.h'
      include 'isopyc.h'
c
      real    t0,s0,dens
      integer mm1,mm2
c
c-----------------------------------------------------------------------
c     compute normalized densities for each isopycnal reference pressure
c     level using a 3rd order polynomial fit to the equation of state.
c     for each isopycnal reference pressure level, the same reference
c     potential temperature, reference salinity and expansion coeff
c     values are used at all of the vertical levels.
c
c     Note: this density is used for the mixing tensor in both the 
c     Redi/Cox and Gent/McWilliams options
c-----------------------------------------------------------------------
c
cxjin calculate densities using unesco(1981) state equation
c
      do k=1,km
      do j=2,jmm
      do i=2,imm
      if(tmask(i,j,k).gt.c0)then
      mm1             = max(k-1,1)
      mm2             = min(k+1,km)
      t0              = t(i,j,k,1,tau)
      s0              = t(i,j,k,2,tau)
      rhoi(i,k,j,xup) = dens(t0,s0,kref(mm1))
      rhoi(i,k,j,xmd) = dens(t0,s0,kref(k  ))
      rhoi(i,k,j,xlo) = dens(t0,s0,kref(mm2))
      endif
      enddo
      enddo
      enddo
c
c
      do j=2,jmm
        call setbcx (rhoi(1,1,j,xup), imt, km)
        call setbcx (rhoi(1,1,j,xmd), imt, km)
        call setbcx (rhoi(1,1,j,xlo), imt, km)
      enddo
c
c-----------------------------------------------------------------------
c     evaluate K2(,,3) centered on the northern face of "T" cells
c-----------------------------------------------------------------------
c
      call k2_3
c
c-----------------------------------------------------------------------
c     evaluate K1(,,3) centered on eastern face of "T" cells
c-----------------------------------------------------------------------
c
      call k1_3
c
c-----------------------------------------------------------------------
c     evaluate K3(,,1..3) centered on bottom face of "T" cells
c-----------------------------------------------------------------------
c
      call k3_123
c
c
c-----------------------------------------------------------------------
c     compute isopycnal advective velocities for tracers
c-----------------------------------------------------------------------
c
      call isoadv
c
c
      return
      end
c
c
c     ===============
      subroutine k1_3
c     ===============
c
c=======================================================================
c     compute "K1(,,3)" at the center of the eastern face of "T" cells
c     use "c1e10" to keep the exponents in range.
c=======================================================================
c
      implicit none
      include 'param.h'
      include 'pconst.h'
      include 'grdvar.h'
      include 'isopyc.h'
c
      real c1e10,eps,fxd,fxe,fxc,fxa,fxb,chkslp,olmask
c
c-----------------------------------------------------------------------
c     set local constants 
c-----------------------------------------------------------------------
      c1e10 = 1.0e10
      eps   = 1.0e-25
c
c-----------------------------------------------------------------------
c     d(rho_barx_barz)/dz centered on eastern face of "t" cells
c     Note: values involving ocean surface and ocean bottom are
c           estimated afterwards using a linear extrapolation
c-----------------------------------------------------------------------
c
      do j=2,jmm
        do k=2,km-1
          fxd = c1e10*p25*rdz0(k)
          do i=1,imm
            e(i,k,j,3) = fxd*(rhoi(i  ,k-1,j,xlo) - rhoi(i  ,k+1,j,xup)
     &                       +rhoi(i+1,k-1,j,xlo) - rhoi(i+1,k+1,j,xup))
          enddo
        enddo
      enddo
c
c
c-----------------------------------------------------------------------
c     linearly extrapolate densities to ocean surface for calculation
c     of d(rho_barx_barz)/dz involving level 1.
c
c     REMARK: requires min(kmt(i,jrow)) = 2 cells in ocean.
c-----------------------------------------------------------------------
c
      k   = 1
      fxd = c1e10*rdz0(k)
      fxe = dz0(k) + p5*dz0(k+1)
      do j=2,jmm
        do i=1,imm
          fxa        = p5*(rhoi(i,k+1,j,xup) + rhoi(i+1,k+1,j,xup))
          fxb        = p5*(rhoi(i,k  ,j,xmd) + rhoi(i+1,k  ,j,xmd))
          fxc        = rdzw(k)*(fxb*fxe - fxa*p5*dz0(k))
          e(i,k,j,3) = fxd*(fxc - p5*(fxa+fxb))
        enddo
      enddo
c
c-----------------------------------------------------------------------
c     linearly extrapolate densities to ocean bottom for calculation
c     of d(rho_barx_barz)/dz involving bottom level.
c-----------------------------------------------------------------------
c
      do j=2,jmm
        do i=1,imm
          e(i,km,j,3) = c0
        enddo
      enddo
c
      do j=2,jmm
        do i=1,imm
          k = min(itn(i,j),itn(i+1,j))
          if (k .ne. 0) then
            fxe        = dz0(k) + p5*dz0(k-1)
            fxa        = p5*(rhoi(i,k-1,j,xlo) + rhoi(i+1,k-1,j,xlo))
            fxb        = p5*(rhoi(i,k  ,j,xmd) + rhoi(i+1,k  ,j,xmd))
            fxc        = rdzw(k-1)*(fxb*fxe - fxa*p5*dz0(k))
            e(i,k,j,3) = rdz0(k)*c1e10*(p5*(fxa+fxb) - fxc)
          endif
        enddo
      enddo
c
c
c-----------------------------------------------------------------------
c     "e(,,,1)" = d(rho)/dx centered on east face of "T" cells
c     "e(,,,2)" = d(rho_barx_bary)/dy centered on east face of "T" cells
c-----------------------------------------------------------------------
c
      do j=2,jmm
        do k=1,km
          do i=1,imm
            e(i,k,j,1) = tmask(i,j,k)*tmask(i+1,j,k)*rdxt(j)*c1e10*
     &                   (rhoi(i+1,k,j,xmd) - rhoi(i,k,j,xmd))
            e(i,k,j,2) = p25*rdy*c1e10*(
     &                    rhoi(i  ,k,j+1,xmd) - rhoi(i  ,k,j-1,xmd)
     &                  + rhoi(i+1,k,j+1,xmd) - rhoi(i+1,k,j-1,xmd))
           enddo
        enddo
      enddo
c
c-----------------------------------------------------------------------
c     if any one of the 4 neighboring corner grid points is a land point,
c     set "e(i,k,j,2)" to zero. note that "e(i,k,j,2)" will be used
c     only in the slope check.
c-----------------------------------------------------------------------
c
      do j=2,jmm
        do k=1,km
          do i=1,imm
            olmask = tmask(i,j-1,k)*tmask(i,j+1,k)*tmask(i+1,j-1,k)
     &   	    *tmask(i+1,j+1,k)
            if (olmask .eq. c0)  e(i,k,j,2) = c0
           enddo
        enddo
      enddo
c
c-----------------------------------------------------------------------
c     impose zonal boundary conditions at "i"=1 and "imt"
c-----------------------------------------------------------------------
      do j=2,jmm
        call setbcx (e(1,1,j,1), imt, km)
        call setbcx (e(1,1,j,2), imt, km)
        call setbcx (e(1,1,j,3), imt, km)
      enddo
c
c
c-----------------------------------------------------------------------
c     compute "K1", using "slmxr" to limit vertical slope of isopycnal
c     to guard against numerical instabilities. 
c-----------------------------------------------------------------------
c
cjxz-------------------------------------------------------------------
c# ifdef isopycmixspatialvar
c          do i=1,imt
c            fxa = c0
c            fxb = sign(1.,e(i,k,j,3))/(abs(e(i,k,j,3))+eps)
c            slope = fxb*sqrt(e(i,k,j,1)**2+e(i,k,j,2)**2)
c            if (slope .le. c0 .and. slope .ge. (-c1/slmxr)) then
c              fxa = p5*(c1+tanh((slope+slopec)/dslope))
c            endif
c            dciso1(i,k,j) = fxa
c            K1(i,k,j,3)  = -fxb*e(i,k,j,1)*dciso1(i,k,j)
c          enddo
c# else # endif
cjxz-------------------------------------------------------------------
      do j=2,jmm
        do k=1,km
          do i=1,imt
cjxz        chkslp = -sqrt(e(i,k,j,1)**2+e(i,k,j,2)**2)*slmxr*dtxsqr(k)
            chkslp = -sqrt(e(i,k,j,1)**2+e(i,k,j,2)**2)*slmxr
            if (e(i,k,j,3) .gt. chkslp)  e(i,k,j,3) = chkslp
          enddo
          do i=1,imt
            K1(i,k,j,3) = (-e(i,k,j,1)*e(i,k,j,3)*fzisop(k))
     &                     /(e(i,k,j,3)**2+eps)
          enddo
        enddo
      enddo
      return
      end
c
c
c     ===============
      subroutine k2_3
c     ===============
c     
c=======================================================================
c     compute "K2(,,3)" at the center of the northern face of "T" cells
c     use "c1e10" to keep the exponents in range.
c=======================================================================
c
      implicit none
      include 'param.h'
      include 'pconst.h'
      include 'grdvar.h'
      include 'isopyc.h'
c
      real c1e10,eps,fxd,fxe,fxc,fxa,fxb,chkslp,olmask
c
c-----------------------------------------------------------------------
c     set local constants 
c-----------------------------------------------------------------------
      c1e10 = 1.0e10
      eps   = 1.0e-25
c
c-----------------------------------------------------------------------
c     d(rho_bary_barz)/dz centered on northern face of "T" cells
c     Note: values involving ocean surface and ocean bottom are
c           estimated afterwards using a linear extrapolation
c-----------------------------------------------------------------------
c
      do j=2,jmm
        do k=2,km-1
          fxd = c1e10*p25*rdz0(k)
          do i=2,imm
            e(i,k,j,3) = fxd*(rhoi(i,k-1,j  ,xlo) - rhoi(i,k+1,j  ,xup)
     &                       +rhoi(i,k-1,j+1,xlo) - rhoi(i,k+1,j+1,xup))
          enddo
        enddo
      enddo
c
c-----------------------------------------------------------------------
c     linearly extrapolate densities to ocean surface for calculation
c     of d(rho_bary_barz)/dz involving level 1.
c-----------------------------------------------------------------------
c
      k   = 1
      fxd = c1e10*rdz0(k) 
      fxe = dz0(k) + p5*dz0(k+1)
      do j=2,jmm
        do i=2,imm
          fxa        = p5*(rhoi(i,k+1,j,xup) + rhoi(i,k+1,j+1,xup))
          fxb        = p5*(rhoi(i,k  ,j,xmd) + rhoi(i,k  ,j+1,xmd))
          fxc        = rdzw(k)*(fxb*fxe - fxa*p5*dz0(k))
          e(i,k,j,3) = fxd*(fxc - p5*(fxa+fxb))
        enddo
      enddo
c
c-----------------------------------------------------------------------
c     linearly extrapolate densities to ocean bottom for calculation
c     of d(rho_bary_barz)/dz involving bottom level.
c-----------------------------------------------------------------------
c
      do j=2,jmm
        do i=2,imm
          e(i,km,j,3) = c0
        enddo
      enddo
c
      do j=2,jmm
        do i=2,imm
          k = min(itn(i,j),itn(i,j+1))
          if (k .ne. 0) then
            fxe        = dz0(k) + p5*dz0(k-1)
            fxa        = p5*(rhoi(i,k-1,j,xlo) + rhoi(i,k-1,j+1,xlo))
            fxb        = p5*(rhoi(i,k  ,j,xmd) + rhoi(i,k  ,j+1,xmd))
            fxc        = rdzw(k-1)*(fxb*fxe - fxa*p5*dz0(k))
            e(i,k,j,3) = rdz0(k)*c1e10*(p5*(fxa+fxb) - fxc)
          endif
        enddo
      enddo
c
c-----------------------------------------------------------------------
c     "e(,,,1)" = d(rho_barx_bary)/dx centered on north face of "T" cells
c     "e(,,,2)" = d(rho)/dy on north face of "T" cells
c-----------------------------------------------------------------------
c
      do j=2,jmm
        do k=1,km
          do i=2,imm
            e(i,k,j,1) = rdxt(j)*p25*c1e10*(
     &                    rhoi(i+1,k,j+1,xmd) - rhoi(i-1,k,j+1,xmd)
     &                  + rhoi(i+1,k,j  ,xmd) - rhoi(i-1,k,j  ,xmd))
            e(i,k,j,2) = tmask(i,j,k)*tmask(i,j+1,k)*rdy*c1e10
     &                  *(rhoi(i,k,j+1,xmd) - rhoi(i,k,j,xmd))   
           enddo
        enddo
      enddo
c
c-----------------------------------------------------------------------
c     if any one of the 4 neighboring corner grid points is a land point,
c     set "e(i,k,j,1)" to zero. note that "e(i,k,j,1)" will be used
c     only in the slope check.
c-----------------------------------------------------------------------
c
      do j=2,jmm
        do k=1,km
          do i=2,imm
            olmask = tmask(i-1,j+1,k)*tmask(i+1,j+1,k)*tmask(i-1,j,k)
     &              *tmask(i+1,j,k)
          if (olmask .eq. c0)  e(i,k,j,1) = c0
           enddo
        enddo
      enddo
c
c-----------------------------------------------------------------------
c     impose zonal boundary conditions at "i"=1 and "imt"
c-----------------------------------------------------------------------
c
      do j=2,jmm
        call setbcx (e(1,1,j,1), imt, km)
        call setbcx (e(1,1,j,2), imt, km)
        call setbcx (e(1,1,j,3), imt, km)
      enddo
c
c
c-----------------------------------------------------------------------
c     compute "K2", using "slmxr" to limit vertical slope of isopycnal
c     to guard against numerical instabilities. 
c-----------------------------------------------------------------------
c
c# ifdef isopycmixspatialvar
c          do i=1,imt
c            fxa = c0
c            fxb = sign(1.,e(i,k,j,3))/(abs(e(i,k,j,3))+eps)
c            slope = fxb*sqrt(e(i,k,j,1)**2+e(i,k,j,2)**2)
c            if (slope .le. c0 .and. slope .ge. (-c1/slmxr)) then
c              fxa = p5*(c1+tanh((slope+slopec)/dslope))
c            endif
c            dciso2(i,k,j) = fxa
c            K2(i,k,j,3)  = -fxb*e(i,k,j,2)*dciso2(i,k,j)
c          enddo
c# else # endif
cjxz-------------------------------------------------------------------
      do j=2,jmm
        do k=1,km
          do i=1,imt
cjxz        chkslp = -sqrt(e(i,k,j,1)**2+e(i,k,j,2)**2)*slmxr*dtxsqr(k)
            chkslp = -sqrt(e(i,k,j,1)**2+e(i,k,j,2)**2)*slmxr
            if (e(i,k,j,3) .gt. chkslp)  e(i,k,j,3) = chkslp
          enddo
          do i=1,imt
            K2(i,k,j,3) = (-e(i,k,j,2)*e(i,k,j,3)*fzisop(k))
     &                     /(e(i,k,j,3)**2+eps)
          enddo
        enddo
      enddo
      return
      end
c
c     =================
      subroutine k3_123
c     =================
c
c=======================================================================
c     compute K2(,,,1:3) at the center of the bottom face of "T" cells
c     use "c1e10" to keep the exponents in range.
c=======================================================================
c
      implicit none
      include 'param.h'
      include 'pconst.h'
      include 'grdvar.h'
      include 'isopyc.h'
c
      real c1e10,eps,chkslp,ahfctr,fxa
c
c-----------------------------------------------------------------------
c     set local constants 
c-----------------------------------------------------------------------
      c1e10 = 1.0e10
      eps   = 1.0e-25
c
c
      do j=2,jmm
        do k=2,km
          do i=2,imm
            e(i,k-1,j,1) = rdxt(j)*p25*c1e10
     &            *(tmask(i-1,j,k-1)*tmask(i,j,k-1)*(rhoi(i,k-1,j,xlo)
     &            -rhoi(i-1,k-1,j,xlo))
     &            +tmask(i,j,k-1)*tmask(i+1,j,k-1)*(rhoi(i+1,k-1,j,xlo)
     &            -rhoi(i,k-1,j,xlo))
     &            +tmask(i-1,j,k)*tmask(i,j,k)*(rhoi(i,k,j,xmd)
     &            -rhoi(i-1,k,j,xmd))
     &            +tmask(i,j,k)*tmask(i+1,j,k)*(rhoi(i+1,k,j,xmd)
     &            -rhoi(i,k,j,xmd)))
c
            e(i,k-1,j,2) = rdy*p25*c1e10
     &            *(tmask(i,j-1,k-1)*tmask(i,j,k-1)*(rhoi(i,k-1,j,xlo)
     &            -rhoi(i,k-1,j-1,xlo))
     &            +tmask(i,j,k-1)*tmask(i,j+1,k-1)*(rhoi(i,k-1,j+1,xlo)
     &            -rhoi(i,k-1,j,xlo))
     &            +tmask(i,j-1,k)*tmask(i,j,k)*(rhoi(i,k,j,xmd)
     &            -rhoi(i,k,j-1,xmd))
     &            +tmask(i,j,k)*tmask(i,j+1,k)*(rhoi(i,k,j+1,xmd)
     &            -rhoi(i,k,j,xmd)))
c
            e(i,k-1,j,3) = rdzw(k-1)*tmask(i,j,k-1)*tmask(i,j,k)*c1e10
     &               *(rhoi(i,k-1,j,xlo) - rhoi(i,k,j,xmd))
c
           enddo
        enddo
        k = km
        do i=2,imm
	  e(i,k,j,1) = c0
	  e(i,k,j,2) = c0
	  e(i,k,j,3) = c0
        enddo
      enddo
c
c-----------------------------------------------------------------------
c     compute "K3", using "slmxr" to limit vertical slope of isopycnal
c     to guard against numerical instabilities.  
c-----------------------------------------------------------------------
c
c# ifdef isopycmixspatialvar
c          do i=2,imt-1
c            fxa = c0
c            fxb = sign(1.,e(i,k,j,3))/(abs(e(i,k,j,3))+eps)
c            slope = fxb*sqrt(e(i,k,j,1)**2+e(i,k,j,2)**2)
c            if (slope .le. c0 .and. slope .ge. (-c1/slmxr)) then
c              fxa = p5*(c1+tanh((slope+slopec)/dslope))
c            endif
c            fxc = fxb*fxa
c            K3(i,k,j,1) = -fxc*e(i,k,j,1)
c            K3(i,k,j,2) = -fxc*e(i,k,j,2)
c            K3(i,k,j,3) = fxb*fxb*(e(i,k,j,1)**2+e(i,k,j,2)**2)
c     &                    *fxa
c          enddo
c# else # endif
      do j=2,jmm
        do k=1,km
          do i=2,imt-1
cjxz        chkslp = -sqrt(e(i,k,j,1)**2+e(i,k,j,2)**2)*slmxr*dtxsqr(k)
            chkslp = -sqrt(e(i,k,j,1)**2+e(i,k,j,2)**2)*slmxr
            if (e(i,k,j,3) .gt. chkslp)  e(i,k,j,3) = chkslp
          enddo
          fxa = p5*(fzisop(min(km,k+1))+fzisop(k))
          do i=2,imt-1
            ahfctr = fxa/(e(i,k,j,3)**2+eps)
            K3(i,k,j,1) = -e(i,k,j,3)*e(i,k,j,1)*ahfctr
            K3(i,k,j,2) = -e(i,k,j,3)*e(i,k,j,2)*ahfctr
            K3(i,k,j,3) = (e(i,k,j,1)**2+e(i,k,j,2)**2)*ahfctr
          enddo
        enddo
      enddo
c
c-----------------------------------------------------------------------
c     impose zonal boundary conditions at "i"=1 and "imt"
c-----------------------------------------------------------------------
c
      do j=2,jmm
        call setbcx (K3(1,1,j,1), imt, km)
        call setbcx (K3(1,1,j,2), imt, km)
        call setbcx (K3(1,1,j,3), imt, km)
      enddo
      return
      end               
c
c
c
c     ===================================
      subroutine setbcx (a, imt, jmtorkm)
c     ===================================
      implicit none
      integer imt,jmtorkm,k
      real a(imt,jmtorkm)
      do k=1,jmtorkm
        a(1,k)   = a(imt-1,k)
        a(imt,k) = a(2,k)
      enddo
      return
      end
c
c     =================
      subroutine isoadv
c     =================
c
c=======================================================================
c     compute isopycnal transport velocities.
c=======================================================================
c
      implicit none
      include 'param.h'
      include 'pconst.h'
      include 'grdvar.h'
      include 'isopyc.h'
c
      integer jstrt
      real    fxa
c
      real fx,fy
      common /wisop/ fx(imt,jmt),fy(imt,jmt)
c
c-----------------------------------------------------------------------
c     compute the meridional component of the isopycnal mixing velocity
c     at the center of the northern face of the "t" cells.
c-----------------------------------------------------------------------
c
      do j=2,jmm
        do k=2,km-1
          fxa = -p5*rdz0(k)*athkdf*cosu(j)
          do i=1,imt
            adv_vntiso(i,k,j) = fxa*tmask(i,j,k)*tmask(i,j+1,k)*(
     &                          K2(i,k-1,j,3) - K2(i,k+1,j,3))
          enddo
        enddo
      enddo
c
c     consider the top and bottom levels. "K2" is assumed to be zero
c     at the ocean top and bottom.
c
      k = 1
      fxa = -p5*rdz0(k)*athkdf
      do j=2,jmm
        do i=1,imt
          adv_vntiso(i,k,j) = -fxa*tmask(i,j,k)*tmask(i,j+1,k)*cosu(j)
     &                        *(K2(i,k,j,3) + K2(i,k+1,j,3))
        enddo
      enddo
c
      do j=2,jmm
        do i=1,imt
          adv_vntiso(i,km,j) = c0
        enddo
      enddo
c
      do j=2,jmm
        do i=1,imt
          k = min(itn(i,j),itn(i,j+1))
          if (k .ne. 0) then
            adv_vntiso(i,k,j) = -p5*rdz0(k)*athkdf*cosu(j)*tmask(i,j,k)
     &                   *tmask(i,j+1,k)*(K2(i,k,j,3) + K2(i,k-1,j,3))
          endif
        enddo
      enddo
c
c
cxjin------------------------------
      do j=2,jmm
      do k=1,km
      do i=1,imt
        adv_vntiso(i,k,j) = adv_vntiso(i,k,j) * fy(i,j+1)
      enddo
      enddo
      enddo
cxjin------------------------------
c
c-----------------------------------------------------------------------
c     compute the zonal component of the isopycnal mixing velocity
c     at the center of the eastern face of the "t" grid box.
c-----------------------------------------------------------------------
c
cjxz  jstrt = max(js,jsmw)
      jstrt = 2
c
      do j=jstrt,jmm
        do k=2,km-1
          fxa = -p5*rdz0(k)*athkdf
          do i=1,imm
            adv_vetiso(i,k,j) = fxa*tmask(i,j,k)*tmask(i+1,j,k)
     &                          *(K1(i,k-1,j,3) - K1(i,k+1,j,3))
          enddo
        enddo
      enddo
c
c     consider the top and bottom levels. "K1" is assumed to be zero
c     at the ocean top and bottom.
c
      k = 1
      fxa = -p5*rdz0(k)*athkdf
      do j=jstrt,jmm
        do i=1,imm
          adv_vetiso(i,k,j) = -fxa*tmask(i,j,k)*tmask(i+1,j,k)
     &                        *(K1(i,k,j,3)+K1(i,k+1,j,3))
        enddo
      enddo
c
      do j=jstrt,jmm
        do i=1,imm
          adv_vetiso(i,km,j) = c0
        enddo
      enddo
c
      do j=jstrt,jmm
        do i=1,imm
          k = min(itn(i,j),itn(i+1,j))
          if (k .ne. 0) then
            adv_vetiso(i,k,j) = -p5*rdz0(k)*athkdf*tmask(i,j,k)
     &                     *tmask(i+1,j,k)*(K1(i,k,j,3)+K1(i,k-1,j,3))
          endif
        enddo
      enddo
c
cxjin------------------------------
      do j=2,jmm
      do k=1,km
      do i=1,imm
        adv_vetiso(i,k,j) = adv_vetiso(i,k,j) * fx(i+1,j)
      enddo
      enddo
      enddo
cxjin------------------------------
c
c----------------------------------------------------------------------
c     set the boundary conditions
c----------------------------------------------------------------------
      do j=jstrt,jmm
        call setbcx (adv_vetiso(1,1,j), imt, km)
      enddo 
c
c----------------------------------------------------------------------
c     compute the vertical component of the isopycnal mixing velocity
c     at the center of the bottom face of the "t" cells, using the
c     continuity equation for the isopycnal mixing velocities
c-----------------------------------------------------------------------
c
c      do j=jstrt,jmm
c        do i=1,imt
c          adv_vbtiso(i,0,j) = c0
c        enddo
c      enddo
c
c      do j=jstrt,jmm
c        do k=1,km-1
c          do i=2,imt
c            adv_vbtiso(i,k,j) = dz0(k)*(
c     &      (adv_vetiso(i,k,j) - adv_vetiso(i-1,k,j))*rdxt(j) + 
c     &      (adv_vntiso(i,k,j) - adv_vntiso(i,k,j-1))*rdyt(j)) 
c          enddo
c        enddo
c      enddo
c
c      do j=jstrt,jmm
c        do k=1,km-1
c          do i=2,imt
c            adv_vbtiso(i,k,j) = adv_vbtiso(i,k,j) + adv_vbtiso(i,k-1,j)
c          enddo
c        enddo
c      enddo
cc
c      do j=jstrt,jmm
c        do i=2,imt
c          adv_vbtiso(i,itn(i,j),j) = c0
c        enddo
c      enddo
c
c-----------------------------------------------------------------------
c     set the boundary conditions
c-----------------------------------------------------------------------
c      do j=jstrt,jmm
c        call setbcx (adv_vbtiso(1,0,j), imt, km+1)
c      enddo
c
c
cxjin-------------------------------------------------
      call isow (adv_vetiso,adv_vntiso,adv_vbtiso)
cxjin-------------------------------------------------
c
c
      return
      end
c
c
c
c     =============================
      subroutine isoflux (tf,mtrace)
c     =============================
c     isopycnal diffusive tracer fluxes are computed.
c
      implicit none
      include 'param.h'
      include 'pconst.h'
      include 'grdvar.h'
      include 'isopyc.h'
      include 'prog.h'
      include 'scalar.h'
c
      real  fx(imt,jmt),fy(imt,jmt)
c
      real  temp(imt,jmt,km),tf(imt,jmt,km)
      real  wki(imt),wkj(jmt),wkk(kmp1)
      real  fk31(0:km),fk32(0:km)
      integer mtrace
      real    fxa,fxb,fxc,fxe
c
      wki(imt)  = c0
      wkj(1)    = c0
      wkj(jmt)  = c0
      wkk(1)    = c0
      wkk(kmp1) = c0
c
      do i=1,imt
      do j=1,jmt
      do k=1,km
      temp(i,j,k)=c0
      enddo
      enddo
      enddo
c
      m = mtrace
c
c-----------------------------------------------------------------------
c     first compute the vertical tracer flux "temp" at the northern
c     face of "t" cells.
c-----------------------------------------------------------------------
c
      do k=2,km-1
        do j=1,jmm
          do i=1,imt
            temp(i,j,k)=p25*rdz0(k)*
     &         (t(i,j+1,k-1,m,taum)-t(i,j+1,k+1,m,taum)
     &         +t(i,j  ,k-1,m,taum)-t(i,j,  k+1,m,taum))
           enddo
        enddo
      enddo
c
c-----------------------------------------------------------------------
c     now consider the top level, assuming that the surface tracer 
c     values are the same as the ones at "k"=1
c-----------------------------------------------------------------------
c
      k = 1
      do j=1,jmm
        do i=1,imt
          temp(i,j,k) = p25*rdz0(k)*
     &         (t(i,j+1,k,m,taum)-t(i,j+1,k+1,m,taum)
     &         +t(i,j  ,k,m,taum)-t(i,j,  k+1,m,taum))
        enddo
      enddo
c
c-----------------------------------------------------------------------
c     finally, consider the bottom level. the extrapolative estimator
c     is used to compute the tracer values at the ocean bottom.
c-----------------------------------------------------------------------
c
      do j=1,jmt
        do i=1,imt
          temp(i,j,km) = 0.0
        enddo
      enddo
c
      do j=1,jmm
        do i=1,imt
          k = min(itn(i,j),itn(i,j+1))
          if (k .ne. 0) then
            fxe = dz0(k) + p5*dz0(k-1)
            fxa = p5*(t(i,j+1,k-1,m,taum) + t(i,j,k-1,m,taum))
            fxb = p5*(t(i,j+1,k  ,m,taum) + t(i,j,k  ,m,taum))
            fxc = rdzw(k-1)*(fxb*fxe-fxa*p5*dz0(k))
            temp(i,j,k) = rdz0(k)*(p5*(fxa+fxb) - fxc)
          endif
        enddo
      enddo
c
c-----------------------------------------------------------------------
c     compute of meridional tracer flux
c     add in the effects of the along isopycnal diffusion computed
c     using "K2" component of the tensor and apply land/sea masks
c-----------------------------------------------------------------------
c
c     fzisop(k)=1.0
c
      do 100 k=1,km
      do 100 i=2,imm
        do j=1,jmm
          wkj(j) = ( rdy*(t(i,j+1,k,m,taum)-t(i,j,k,m,taum))
     &             + K2(i,k,j,3)*temp(i,j,k) )*
     &             cosu(j)*tmask(i,j,k)*tmask(i,j+1,k)*fy(i,j+1)*p5
        enddo
        do j=2,jmm
          tf(i,j,k) = tf(i,j,k) + ahisop*rdyt(j)*(wkj(j)-wkj(j-1))
        enddo
100   continue
c
c
c-----------------------------------------------------------------------
c     compute the vertical tracer flux "temp" at the eastern
c     face of "t" cells.
c-----------------------------------------------------------------------
c
      do k=2,km-1
        do j=2,jmm
          do i=1,imm
            temp(i,j,k)=p25*rdz0(k)*
     &          (t(i+1,j,k-1,m,taum)-t(i+1,j,k+1,m,taum)
     &          +t(i  ,j,k-1,m,taum)-t(i  ,j,k+1,m,taum))
          enddo
        enddo
      enddo
c
c-----------------------------------------------------------------------
c     now consider the top level, assuming that the surface tracer 
c     values are the same as the ones at "k"=1
c-----------------------------------------------------------------------
c
      k = 1
      do j=2,jmm
        do i=1,imm
          temp(i,j,k)=p25*rdz0(k)*
     &          (t(i+1,j,k,m,taum)-t(i+1,j,k+1,m,taum)
     &          +t(i  ,j,k,m,taum)-t(i  ,j,k+1,m,taum))
        enddo
      enddo
c
c-----------------------------------------------------------------------
c     finally, consider the bottom level. the extrapolative estimator
c     is used to compute the tracer values at the ocean bottom.
c-----------------------------------------------------------------------
c
      do j=2,jmm
        do i=1,imm
          temp(i,j,km) = 0.0
        enddo
      enddo
c
      do j=2,jmm
        do i=1,imm
          k = min(itn(i,j),itn(i+1,j))
          if (k .ne. 0) then
            fxe          = dz0(k) + p5*dz0(k-1)
            fxa          = p5*(t(i,j,k-1,m,taum)+t(i+1,j,k-1,m,taum))
            fxb          = p5*(t(i,j,k,m,taum)+t(i+1,j,k,m,taum))
            fxc          = rdzw(k-1)*(fxb*fxe - fxa*p5*dz0(k))
            temp(i,j,k)  = rdz0(k)*(p5*(fxa+fxb)-fxc)
          endif
        enddo
      enddo
c
c-----------------------------------------------------------------------
c     compute of zonal tracer flux 
c     add in the effects of the along isopycnal diffusion computed
c     using "K1" component of the tensor and apply land/sea masks      
c-----------------------------------------------------------------------
c      
c     fzisop(k)=1.0
c     dxtr=const
c
      do 120 k=1,km
      do 120 j=2,jmm
        do i=1,imm
          wki(i) = ( rdxt(j)*(t(i+1,j,k,m,taum)-t(i,j,k,m,taum))
     &             + K1(i,k,j,3)*temp(i,j,k) )*
     &             tmask(i+1,j,k)*tmask(i,j,k)*fx(i+1,j)*p5
        enddo
        do i=2,imm
          tf(i,j,k) = tf(i,j,k) + ahisop*rdxt(j)*(wki(i)-wki(i-1))
        enddo
120   continue
c
c
c-----------------------------------------------------------------------
c     compute the vertical tracer flux "diff_fbiso" containing the K31
c     and K32 components which are to be solved explicitly. The K33
c     component will be treated semi-implicitly
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c     at ocean surface the flux is set to zero to reflect the no tracer 
c     flux condition. Same condition is also imposed at ocean bottom.
c-----------------------------------------------------------------------
      do k=0,km
      fk31(k)=0.0
      fk32(k)=0.0
      enddo
c
      do 140 j=2,jmm
      do 140 i=2,imm
c
        do k=2, min(itn(i,j)+1,km)
          fk31(k-1) = ahisop*p25*rdxt(j)*K3(i,k-1,j,1)*tmask(i,j,k)*(
     & tmask(i-1,j,k  )*(t(i,j,k  ,m,taum)-t(i-1,j,k  ,m,taum))+
     & tmask(i-1,j,k-1)*(t(i,j,k-1,m,taum)-t(i-1,j,k-1,m,taum))+
     & tmask(i+1,j,k  )*(t(i+1,j,k,m,taum)-t(i,j  ,k  ,m,taum))+
     & tmask(i+1,j,k-1)*(t(i+1,j,k-1,m,taum)-t(i,j,k-1,m,taum)) )
        enddo
c
        do k=2, min(itn(i,j)+1,km)
          fk32(k-1) = ahisop*p25*rdy*K3(i,k-1,j,2)*tmask(i,j,k)*(
     & tmask(i,j-1,k  )*(t(i,j,k  ,m,taum)-t(i,j-1,k  ,m,taum))+
     & tmask(i,j-1,k-1)*(t(i,j,k-1,m,taum)-t(i,j-1,k-1,m,taum))+
     & tmask(i,j+1,k  )*(t(i,j+1,k,m,taum)-t(i,j  ,k  ,m,taum))+
     & tmask(i,j+1,k-1)*(t(i,j+1,k-1,m,taum)-t(i,j,k-1,m,taum)) )
        enddo
c
        do k=1,km
          tf(i,j,k) = tf(i,j,k) + gravr*rdz(k)*(
     &                (fk31(k-1)-fk31(k))+(fk32(k-1)-fk32(k)) )
        enddo
140   continue
c
c
cifdef gent_mcwilliams
c
c-----------------------------------------------------------------------
c     compute the meridional component of the isopycnal velocity mixing 
c-----------------------------------------------------------------------
c
c
      do k=1,km
        do i=2,imm
          do j=2,jmm
            wkj(j)    = adv_vntiso(i,k,j)*
     &           (t(i,j+1,k,m,taum)-t(i,j,k,m,taum))
          enddo
          do j=2,jmm
            tf(i,j,k) = tf(i,j,k) - rdyt(j)*(wkj(j)+wkj(j-1))*p5
          enddo
        enddo
      enddo
c
c
c-----------------------------------------------------------------------
c     compute the zonal component of the isopycnal velocity mixing 
c-----------------------------------------------------------------------
c
      do k=1,km
        do j=2,jmm
          do i=1,imm
            wki(i)    = adv_vetiso(i,k,j)*
     &        (t(i+1,j,k,m,taum)-t(i,j,k,m,taum))
          enddo
          do i=2,imm
            tf(i,j,k) = tf(i,j,k) - rdxt(j)*(wki(i)+wki(i-1))*p5
          enddo
        enddo
      enddo
c
c-----------------------------------------------------------------------
c     compute the vertical component of the isopycnal velocity mixing 
c-----------------------------------------------------------------------
c
      do j=2,jmm
        do i=2,imm
          do k=2, min(itn(i,j)+1,km)
            wkk(k)  = adv_vbtiso(i,k-1,j)*
     &           (t(i,j,k-1,m,taum)-t(i,j,k,m,taum))
          enddo
        do k=1,itn(i,j)
          tf(i,j,k) = tf(i,j,k) - (wkk(k)+wkk(k+1))*rdz(k)*p5
          enddo
        enddo
      enddo
c
      return
      end
#endif
c
c
c
c
c     ======================
      subroutine isow(u,v,w)
c     ======================
c     calculate vertical mass advection: Ps*dz/dt
c     upward is positive
c
      implicit none
      include 'param.h'
      include 'pconst.h'
      include 'grdvar.h'
      include 'scalar.h'
c
      real u(imt,km,jmt),v(imt,km,jmt),w(imt,0:km,jmt)
      real dp(imt,jmt)
c
      do j=1,jmt
      do i=1,imt
      w(i,0,j) = c0
      dp(i,j)  = c0
      enddo
      enddo
c
      do j=2,jmm
      do k=1,km
      do i=2,imm
      w(i,k,j) = dz(k)*((u(i,k,j)-u(i-1,k,j))*rdxt(j) +
     &                  (v(i,k,j)-v(i,k,j-1))*rdyt(j))
      dp(i,j)  = dp(i,j) - w(i,k,j)
      enddo
      enddo
      enddo
c
      do j=2,jmm
      do k=1,km-1
      do i=2,imm
      w(i,k,j)=(w(i,k-1,j)+dz(k)*dp(i,j)/pn(i,j)+w(i,k,j))
     &         *tmask(i,j,k+1)
      enddo
      enddo
      enddo
c
      do j=2,jmm
        do i=2,imt
          w(i,itn(i,j),j) = c0
        enddo
      enddo
c
c
      do k=1,km
      do j=2,jmm
      w(1  ,k,j) = w(imm,k,j)
      w(imt,k,j) = w(2  ,k,j)
      enddo
      enddo
c
      return
      end
