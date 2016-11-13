!     =================
      subroutine isopyi(imt,jmt,km,kmp1,slmxr,ahisop,athkdf,fzisop,kref,rdz0,   &
                        dz0,z0,K1,K2,K3,adv_vetiso,adv_vntiso,rhoi,e,adv_vbtiso)
!     =================
!
!=======================================================================
!
!       Initialization for isopycnal mixing scheme
!
!       Redi/Cox version + Gent_McWilliams version
!
!ifdef isopycxixspatialvar
!       dciso1 = isopycnal tracer diffusivity coeffs modified based
!              on the slopes of the isopycnal surfaces on the east face
!              of "T" cells.
!       dciso2 = isopycnal tracer diffusivity coeffs modified based
!              on the slopes of the isopycnal surfaces on the north face
!              of "T" cells.
!       dslope = half length of the interval in which "ahisop" changes
!              with a steep slope from about 0.9*"ahisop" to about
!              0.1*"ahisop"
!       slopec = slope at which "ahisop" is equal to half of its 
!              original  value
!endif
!
!       dptlim = depth limits for the reference pressure levels (in cm). 
!
!=======================================================================
!
      implicit none
!      include 'param.h'
      include 'pconst.h'
!      include 'grdvar.h'
!      include 'isopyc.h'
      include 'mpif.h'
!
      integer xup,xmd,xlo
      parameter(xup=1,xmd=2,xlo=3)
      real dptmid,t1,t2
      
      integer i,j,k,m,imt,jmt,km,kmp1
      real slmxr,ahisop,athkdf,fzisop(km),kref(km),rdz0(km),dz0(km),z0(km)
      real K1(imt,km,jmt,3:3),K2(imt,km,jmt,3:3),K3(imt,km,jmt,1:3)
      real adv_vetiso(imt,km,jmt),adv_vntiso(imt,km,jmt)
      real rhoi(imt,km,jmt,xup:xlo),e(imt,kmp1,jmt,3),adv_vbtiso(imt,0:km,jmt)
!
!-----------------------------------------------------------------------
!     USER INPUT
!     initialize variables (all mixing units are cm**2/sec.)
!-----------------------------------------------------------------------
!
      slmxr  = 100.0
      ahisop = 1.e7
!
!     define the isopycnal thickness diffusion coefficient
!
      athkdf = 1.0e7
!
!     reference pressure level intervals are defined (see "isopyc.h").
!     "dptlim" must have "nrpl+1" elements (see "param.h"). also,
!     "dptlim" are in cm.
!
!
!     REMARK: the first and the last elements of "dptlim" must be the 
!             depth at the top (0cm) and the maximum bottom depth,
!             respectively. Also, the elements of "dptlim" must be in
!             increasing order.
!
!     dptlim(1) = 0.0e2
!     dptlim(2) = 1000.0e2
!     dptlim(3) = 2000.0e2
!     dptlim(4) = 3000.0e2
!     dptlim(5) = 4000.0e2
!     dptlim(6) = 5000.0e2
!
!
!-----------------------------------------------------------------------
!     determine the isopycnal reference pressure levels for the "t"
!     grid point levels, using the depths at the "t" grid points as the
!     reference depth (pressure)
!-----------------------------------------------------------------------
!
!     do k=1,km
!       do m=2,nrpl+1
!         if (z0(k).gt.dptlim(m-1) .and. z0(k).le.dptlim(m)) then
!           kisrpl(k) = m-1
!           go to 101
!         endif
!       enddo
101     continue
!
!       if (kisrpl(k) .lt. 1 .or. kisrpl(k) .gt. nrpl) then
!         write (6,9100) kisrpl(k), k
!         stop 9999
!       endif
!     enddo
9100  format (/,' =>Error: kisrpl is ',i3,' at k ',1x,i3,' in isopyc.F')
!
!-----------------------------------------------------------------------
!     the indices used in isopycnal mixing indicating the location of
!     the reference pressure levels in the 20-level table of polynomial
!     expansion variables are computed
!
!     REMARK: because the polynomial expansion coefficients are
!             functions of the reference potential temperature and
!             salinity profiles, at the reference pressure level
!             the corresponding potential temperature and salinity
!             values will be used.
!-----------------------------------------------------------------------
!
!     do m=1,nrpl
!       krplin(m) = 0
!     enddo
!
!     do m=2,nrpl+1
!       dptmid = 0.5*(dptlim(m-1)+dptlim(m))
!       if (dptmid .le. z0(1)) then
!         krplin(m-1) = 1
!       elseif (dptmid .gt. z0(km)) then
!         krplin(m-1) = km
!       elseif (dptmid.gt.z0(1) .and. dptmid.le.z0(km)) then
!         do k=2,km
!           if (z0(k) .ge. dptmid) then
!             t1 = z0(k)-dptmid
!             t2 = dptmid-z0(k-1)
!             if (t1 .gt. t2) then
!               krplin(m-1) = k-1
!             else
!               krplin(m-1) = k
!             endif
!             go to 102 
!           endif
!         enddo
102       continue
!       endif
!       if (krplin(m-1) .lt. 1 .or. krplin(m-1) .gt. km) then
!         write (6,9110) krplin(m-1),m-1
!       endif
!     enddo
!
!cc   write (6,96) (kisrpl(k),k=1,km)
!cc   write (6,97) (krplin(m),m=1,nrpl)
!
!-----------------------------------------------------------------------
!     the isopycnal diffusion coefficient may be a function of depth. in 
!     the default configuration, the isopycnal diffusion coefficient is 
!     a constant: "fzisop", which multiplies "ahisop", is set to unity. 
!     if "ahisop" varies in the vertical, "fzisop" should contain this 
!     variation. the value of "ahisop" should be adjusted accordingly.
!-----------------------------------------------------------------------
!
      do k=1,km
        fzisop(k) = c1
      enddo
!
!cc   write (6,98)
!cc   write (6,99) (fzisop(k),k=1,km)
!
  96  format (/,' isopycnal reference pressure levels (kisrpl) = ', &
              20(1x,i4))
  97  format (/,' reference pressure level indices (krplin) = ', &
              20(1x,i4))
  98  format (/,' vertical variation of "ahisop" (fzisop) = ')
  99  format (5(1x,e12.6))
9110  format (/,' =>Error: krplin is ',i3,' at m ',1x,i3)
!
!
!-----------------------------------------------------------------------
!     set reference depth for calculation of rhoi  (m)
!-----------------------------------------------------------------------
      do k=1,km
      kref(k) = z0(k)*0.01
      enddo
      do k=1,km
      rdz0(k) = c1/dz0(k)
      enddo
!
!
!-----------------------------------------------------------------------
!     initialize arrays
!-----------------------------------------------------------------------
!
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
!
      do j=1,jmt
	    do k=1,km
          do i=1,imt
	        rhoi(i,k,j,xup) = c1
	        rhoi(i,k,j,xmd) = c1
	        rhoi(i,k,j,xlo) = c1
	      enddo
	    enddo
      enddo
!
      do m=1,3
        do j=1,jmt
	      do k=1,km+1
            do i=1,imt
	          e(i,k,j,m) = c0
	        enddo
	      enddo
        enddo
      enddo
!
      do j=1,jmt
        do k=0,km
          do i=1,imt
            adv_vbtiso(i,k,j) = c0
	      enddo
        enddo
      enddo
!
      return
      end
!
!     =================
      subroutine isopyc(imt,jmt,km,kmp1,imm,jmm,nt,itn,tmask,kref,fzisop,rdz0,dz0,z0,  &
                        rdxt,rdyt,rdy,rdzw,dz,pn,cosu,t,slmxr,athkdf,K1,K2,K3,  &
                        adv_vetiso,adv_vntiso,adv_vbtiso,rhoi,e,fx,fy,west,east,north, &
                        south,unesco)
!     =================
!
!=======================================================================
!
!     Compute the isopycnal mixing tensor components and the
!     isopycnal advection velocities which parameterize the effect
!     of eddies on the isopycnals.
!
!
!     Mixing tensor "K" is ...
!
!          | 1.0            K1(,,,2)          K1(,,,3) | 
!          |                                           |
!     K =  | K2(,,,1)        1.0              K2(,,,3) |
!          |                                           |
!          | K3(,,,1)       K3(,,,2)          K3(,,,3) |
!
!     where K1(,,,2) and K2(,,,1) are set to 0.0 (neglected)
!
!
!     output:
!       rhoi = density at tau-1 referenced to pressure levels
!       K1   = tensor components (1,2), (1,3) centered on east face
!              of "T" cells
!       K2   = tensor components (2,1), (2,3) centered on north face
!              of "T" cells
!       K3   = tensor components (3,1), (3,2), (3,3) centered on
!              bottom face of "T" cells  
!ifdef gent_mcwilliams
!       adv_vetiso = isopycnal advective vel on east face of "T" cell
!       adv_vntiso = isopycnal advective vel on north face of "T" cell
!               (Note: this includes the cosine factor as in "adv_vnt")
!       adv_vbtiso = isopycnal advective vel on bottom face of "T" cell       
!endif
!
!=======================================================================
!
!
      implicit none
!      include 'param.h'
      include 'pconst.h'
!      include 'grdvar.h'
!      include 'prog.h'
!      include 'isopyc.h'
      include 'mpif.h'
!
      integer xup,xmd,xlo
      parameter(xup=1,xmd=2,xlo=3)
      
      integer i,j,k,imt,jmt,km,kmp1,imm,jmm,nt,unesco
      real    t0,s0,dens,undens
      integer mm1,mm2
      
      integer itn(imt,jmt)
      real slmxr,athkdf,rdy
      real fx(imt,jmt),fy(imt,jmt)
      real tmask(imt,jmt,km),rhoi(imt,km,jmt,xup:xlo),kref(km),fzisop(km)
      real e(imt,kmp1,jmt,3),adv_vbtiso(imt,0:km,jmt)
      real rdz0(km),dz0(km),rdxt(jmt),rdzw(km),rdyt(jmt),dz(km),pn(imt,jmt)
      real cosu(jmt),z0(km),t(imt,jmt,km,nt,2)
      real K1(imt,km,jmt,3:3),K2(imt,km,jmt,3:3),K3(imt,km,jmt,1:3)
      real adv_vetiso(imt,km,jmt),adv_vntiso(imt,km,jmt)
      
      integer west,east,north,south
!
!-----------------------------------------------------------------------
!     compute normalized densities for each isopycnal reference pressure
!     level using a 3rd order polynomial fit to the equation of state.
!     for each isopycnal reference pressure level, the same reference
!     potential temperature, reference salinity and expansion coeff
!     values are used at all of the vertical levels.
!
!     Note: this density is used for the mixing tensor in both the 
!     Redi/Cox and Gent/McWilliams options
!-----------------------------------------------------------------------
!
!xjin calculate densities using unesco(1981) state equation
!
      if (unesco==1) then
      do k=1,km
      do j=1,jmt
      do i=1,imt
      if(tmask(i,j,k).gt.c0)then
      mm1             = max(k-1,1)
      mm2             = min(k+1,km)
      t0              = t(i,j,k,1,tau)
      s0              = t(i,j,k,2,tau)
      rhoi(i,k,j,xup) = undens(t0,s0,kref(mm1))
      rhoi(i,k,j,xmd) = undens(t0,s0,kref(k  ))
      rhoi(i,k,j,xlo) = undens(t0,s0,kref(mm2))
      endif
      enddo
      enddo
      enddo
      else
      do k=1,km
      do j=1,jmt
      do i=1,imt
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
      end if
!
!
!      do j=1,jmt
!        call setbcx (rhoi(1,1,j,xup), imt, km)
!        call setbcx (rhoi(1,1,j,xmd), imt, km)
!        call setbcx (rhoi(1,1,j,xlo), imt, km)
!      enddo
!      call swap_isopy_real4d(rhoi,imt,jmt,km,3,west,east,north,south)
!
!-----------------------------------------------------------------------
!     evaluate K2(,,3) centered on the northern face of "T" cells
!-----------------------------------------------------------------------
!
      call k2_3(imt,jmt,imm,jmm,km,kmp1,slmxr,itn,e,rhoi,rdz0,dz0,tmask,   &
                rdxt,rdy,rdzw,fzisop,K2,west,east,north,south)
!
!-----------------------------------------------------------------------
!     evaluate K1(,,3) centered on eastern face of "T" cells
!-----------------------------------------------------------------------
!
      call k1_3(imt,jmt,imm,jmm,km,kmp1,slmxr,itn,e,rhoi,rdz0,dz0,tmask,   &
                rdxt,rdy,rdzw,fzisop,K1,west,east,north,south)
!
!-----------------------------------------------------------------------
!     evaluate K3(,,1..3) centered on bottom face of "T" cells
!-----------------------------------------------------------------------
!
      call k3_123(imt,jmt,imm,jmm,km,kmp1,slmxr,e,rhoi,rdz0,dz0,tmask,   &
                  rdxt,rdy,rdzw,fzisop,K3)
!
!
!-----------------------------------------------------------------------
!     compute isopycnal advective velocities for tracers
!-----------------------------------------------------------------------
!
      call isoadv(imt,jmt,imm,jmm,km,athkdf,itn,rdxt,rdyt,dz,pn,fx,fy,rdz0, &
                  cosu,tmask,adv_vetiso,adv_vntiso,adv_vbtiso,K1,K2,west,east,  &
                  north,south)
!
      return
      end
!
!
!     ===============
      subroutine k1_3(imt,jmt,imm,jmm,km,kmp1,slmxr,itn,e,rhoi,rdz0,dz0,tmask,   &
                      rdxt,rdy,rdzw,fzisop,K1,west,east,north,south)
!     ===============
!
!=======================================================================
!     compute "K1(,,3)" at the center of the eastern face of "T" cells
!     use "c1e10" to keep the exponents in range.
!=======================================================================
!
      implicit none
!      include 'param.h'
      include 'pconst.h'
!      include 'grdvar.h'
!      include 'isopyc.h'
      include 'mpif.h'
!
      integer xup,xmd,xlo
      parameter(xup=1,xmd=2,xlo=3)
      real c1e10,eps,fxd,fxe,fxc,fxa,fxb,chkslp,olmask

      integer i,j,k,imt,jmt,imm,jmm,km,kmp1
      integer itn(imt,jmt)
      integer west,east,north,south
      real slmxr,rdy
      real e(imt,kmp1,jmt,3),rhoi(imt,km,jmt,xup:xlo),rdz0(km),dz0(km)
      real tmask(imt,jmt,km),rdxt(jmt),rdzw(km),fzisop(km)
      real K1(imt,km,jmt,3:3)
!
!-----------------------------------------------------------------------
!     set local constants 
!-----------------------------------------------------------------------
      c1e10 = 1.0e10
      eps   = 1.0e-25
!
!-----------------------------------------------------------------------
!     d(rho_barx_barz)/dz centered on eastern face of "t" cells
!     Note: values involving ocean surface and ocean bottom are
!           estimated afterwards using a linear extrapolation
!-----------------------------------------------------------------------
!
      do j=2,jmm
        do k=2,km-1
          fxd = c1e10*p25*rdz0(k)
          do i=1,imm
            e(i,k,j,3) = fxd*(rhoi(i  ,k-1,j,xlo) - rhoi(i  ,k+1,j,xup)  &
                             +rhoi(i+1,k-1,j,xlo) - rhoi(i+1,k+1,j,xup))
          enddo
        enddo
      enddo
!
!
!-----------------------------------------------------------------------
!     linearly extrapolate densities to ocean surface for calculation
!     of d(rho_barx_barz)/dz involving level 1.
!
!     REMARK: requires min(kmt(i,jrow)) = 2 cells in ocean.
!-----------------------------------------------------------------------
!
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
!
!-----------------------------------------------------------------------
!     linearly extrapolate densities to ocean bottom for calculation
!     of d(rho_barx_barz)/dz involving bottom level.
!-----------------------------------------------------------------------
!
      do j=2,jmm
        do i=1,imm
          e(i,km,j,3) = c0
        enddo
      enddo
!
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
!
!
!-----------------------------------------------------------------------
!     "e(,,,1)" = d(rho)/dx centered on east face of "T" cells
!     "e(,,,2)" = d(rho_barx_bary)/dy centered on east face of "T" cells
!-----------------------------------------------------------------------
!
      do j=2,jmm
        do k=1,km
          do i=1,imm
            e(i,k,j,1) = tmask(i,j,k)*tmask(i+1,j,k)*rdxt(j)*c1e10* &
                         (rhoi(i+1,k,j,xmd) - rhoi(i,k,j,xmd))
            e(i,k,j,2) = p25*rdy*c1e10*(  &
                          rhoi(i  ,k,j+1,xmd) - rhoi(i  ,k,j-1,xmd) &
                        + rhoi(i+1,k,j+1,xmd) - rhoi(i+1,k,j-1,xmd))
           enddo
        enddo
      enddo
!
!-----------------------------------------------------------------------
!     if any one of the 4 neighboring corner grid points is a land point,
!     set "e(i,k,j,2)" to zero. note that "e(i,k,j,2)" will be used
!     only in the slope check.
!-----------------------------------------------------------------------
!
      do j=2,jmm
        do k=1,km
          do i=1,imm
            olmask = tmask(i,j-1,k)*tmask(i,j+1,k)*tmask(i+1,j-1,k) &
         	    *tmask(i+1,j+1,k)
            if (olmask .eq. c0)  e(i,k,j,2) = c0
           enddo
        enddo
      enddo
!
!-----------------------------------------------------------------------
!     impose zonal boundary conditions at "i"=1 and "imt"
!-----------------------------------------------------------------------
!      do j=2,jmm
!        call setbcx (e(1,1,j,1), imt, km)
!        call setbcx (e(1,1,j,2), imt, km)
!        call setbcx (e(1,1,j,3), imt, km)
!      enddo
!      call swap_isopy_real4d(e,imt,jmt,kmp1,3,west,east,north,south)
!
!
!-----------------------------------------------------------------------
!     compute "K1", using "slmxr" to limit vertical slope of isopycnal
!     to guard against numerical instabilities. 
!-----------------------------------------------------------------------
!
!jxz-------------------------------------------------------------------
!# ifdef isopycmixspatialvar
!          do i=1,imt
!            fxa = c0
!            fxb = sign(1.,e(i,k,j,3))/(abs(e(i,k,j,3))+eps)
!            slope = fxb*sqrt(e(i,k,j,1)**2+e(i,k,j,2)**2)
!            if (slope .le. c0 .and. slope .ge. (-c1/slmxr)) then
!              fxa = p5*(c1+tanh((slope+slopec)/dslope))
!            endif
!            dciso1(i,k,j) = fxa
!            K1(i,k,j,3)  = -fxb*e(i,k,j,1)*dciso1(i,k,j)
!          enddo
!# else # endif
!jxz-------------------------------------------------------------------
      do j=2,jmm
        do k=1,km
          do i=1,imt
!jxz        chkslp = -sqrt(e(i,k,j,1)**2+e(i,k,j,2)**2)*slmxr*dtxsqr(k)
            chkslp = -sqrt(e(i,k,j,1)**2+e(i,k,j,2)**2)*slmxr
            if (e(i,k,j,3) .gt. chkslp)  e(i,k,j,3) = chkslp
          enddo
          do i=1,imt
            K1(i,k,j,3) = (-e(i,k,j,1)*e(i,k,j,3)*fzisop(k)) &
                           /(e(i,k,j,3)**2+eps)
          enddo
        enddo
      enddo
      return
      end
!
!
!     ===============
      subroutine k2_3(imt,jmt,imm,jmm,km,kmp1,slmxr,itn,e,rhoi,rdz0,dz0,tmask,   &
                      rdxt,rdy,rdzw,fzisop,K2,west,east,north,south)
!     ===============
!     
!=======================================================================
!     compute "K2(,,3)" at the center of the northern face of "T" cells
!     use "c1e10" to keep the exponents in range.
!=======================================================================
!
      implicit none
!      include 'param.h'
      include 'pconst.h'
!      include 'grdvar.h'
!      include 'isopyc.h'
      include 'mpif.h'
!
      integer xup,xmd,xlo
      parameter(xup=1,xmd=2,xlo=3)
      real c1e10,eps,fxd,fxe,fxc,fxa,fxb,chkslp,olmask
      
      integer i,j,k,imt,jmt,km,imm,jmm,kmp1
      integer itn(imt,jmt)
      integer west,east,north,south
      real slmxr,rdy
      real e(imt,kmp1,jmt,3),rhoi(imt,km,jmt,xup:xlo),rdz0(km),dz0(km)
      real tmask(imt,jmt,km),rdxt(jmt),rdzw(km),fzisop(km)
      real K2(imt,km,jmt,3:3)
!
!-----------------------------------------------------------------------
!     set local constants 
!-----------------------------------------------------------------------
      c1e10 = 1.0e10
      eps   = 1.0e-25
!
!-----------------------------------------------------------------------
!     d(rho_bary_barz)/dz centered on northern face of "T" cells
!     Note: values involving ocean surface and ocean bottom are
!           estimated afterwards using a linear extrapolation
!-----------------------------------------------------------------------
!
      do j=2,jmm
        do k=2,km-1
          fxd = c1e10*p25*rdz0(k)
          do i=2,imm
            e(i,k,j,3) = fxd*(rhoi(i,k-1,j  ,xlo) - rhoi(i,k+1,j  ,xup) &
                             +rhoi(i,k-1,j+1,xlo) - rhoi(i,k+1,j+1,xup))
          enddo
        enddo
      enddo
!
!-----------------------------------------------------------------------
!     linearly extrapolate densities to ocean surface for calculation
!     of d(rho_bary_barz)/dz involving level 1.
!-----------------------------------------------------------------------
!
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
!
!-----------------------------------------------------------------------
!     linearly extrapolate densities to ocean bottom for calculation
!     of d(rho_bary_barz)/dz involving bottom level.
!-----------------------------------------------------------------------
!
      do j=2,jmm
        do i=2,imm
          e(i,km,j,3) = c0
        enddo
      enddo
!
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
!
!-----------------------------------------------------------------------
!     "e(,,,1)" = d(rho_barx_bary)/dx centered on north face of "T" cells
!     "e(,,,2)" = d(rho)/dy on north face of "T" cells
!-----------------------------------------------------------------------
!
      do j=2,jmm
        do k=1,km
          do i=2,imm
            e(i,k,j,1) = rdxt(j)*p25*c1e10*( &
                         rhoi(i+1,k,j+1,xmd) - rhoi(i-1,k,j+1,xmd) &
                       + rhoi(i+1,k,j  ,xmd) - rhoi(i-1,k,j  ,xmd))
            e(i,k,j,2) = tmask(i,j,k)*tmask(i,j+1,k)*rdy*c1e10 &
                        *(rhoi(i,k,j+1,xmd) - rhoi(i,k,j,xmd))   
           enddo
        enddo
      enddo
!
!-----------------------------------------------------------------------
!     if any one of the 4 neighboring corner grid points is a land point,
!     set "e(i,k,j,1)" to zero. note that "e(i,k,j,1)" will be used
!     only in the slope check.
!-----------------------------------------------------------------------
!
      do j=2,jmm
        do k=1,km
          do i=2,imm
            olmask = tmask(i-1,j+1,k)*tmask(i+1,j+1,k)*tmask(i-1,j,k) &
                    *tmask(i+1,j,k)
          if (olmask .eq. c0)  e(i,k,j,1) = c0
           enddo
        enddo
      enddo
!
!-----------------------------------------------------------------------
!     impose zonal boundary conditions at "i"=1 and "imt"
!-----------------------------------------------------------------------
!
!      do j=2,jmm
!        call setbcx (e(1,1,j,1), imt, km)
!        call setbcx (e(1,1,j,2), imt, km)
!        call setbcx (e(1,1,j,3), imt, km)
!      enddo
!      call swap_isopy_real4d(e,imt,jmt,kmp1,3,west,east,north,south)
!
!
!-----------------------------------------------------------------------
!     compute "K2", using "slmxr" to limit vertical slope of isopycnal
!     to guard against numerical instabilities. 
!-----------------------------------------------------------------------
!
!# ifdef isopycmixspatialvar
!          do i=1,imt
!            fxa = c0
!            fxb = sign(1.,e(i,k,j,3))/(abs(e(i,k,j,3))+eps)
!            slope = fxb*sqrt(e(i,k,j,1)**2+e(i,k,j,2)**2)
!            if (slope .le. c0 .and. slope .ge. (-c1/slmxr)) then
!              fxa = p5*(c1+tanh((slope+slopec)/dslope))
!            endif
!            dciso2(i,k,j) = fxa
!            K2(i,k,j,3)  = -fxb*e(i,k,j,2)*dciso2(i,k,j)
!          enddo
!# else # endif
!jxz-------------------------------------------------------------------
      do j=2,jmm
        do k=1,km
          do i=1,imt
!jxz        chkslp = -sqrt(e(i,k,j,1)**2+e(i,k,j,2)**2)*slmxr*dtxsqr(k)
            chkslp = -sqrt(e(i,k,j,1)**2+e(i,k,j,2)**2)*slmxr
            if (e(i,k,j,3) .gt. chkslp)  e(i,k,j,3) = chkslp
          enddo
          do i=1,imt
            K2(i,k,j,3) = (-e(i,k,j,2)*e(i,k,j,3)*fzisop(k)) &
                          /(e(i,k,j,3)**2+eps)
          enddo
        enddo
      enddo
      return
      end
!
!     =================
      subroutine k3_123(imt,jmt,imm,jmm,km,kmp1,slmxr,e,rhoi,rdz0,dz0,tmask,   &
                        rdxt,rdy,rdzw,fzisop,K3)
!     =================
!
!=======================================================================
!     compute K2(,,,1:3) at the center of the bottom face of "T" cells
!     use "c1e10" to keep the exponents in range.
!=======================================================================
!
      implicit none
!      include 'param.h'
      include 'pconst.h'
!      include 'grdvar.h'
!      include 'isopyc.h'
      include 'mpif.h'
!
      integer xup,xmd,xlo
      parameter(xup=1,xmd=2,xlo=3)
      real c1e10,eps,chkslp,ahfctr,fxa
      
      integer i,j,k,imt,jmt,km,imm,jmm,kmp1
      real slmxr,rdy
      real e(imt,kmp1,jmt,3),rhoi(imt,km,jmt,xup:xlo),rdz0(km),dz0(km)
      real tmask(imt,jmt,km),rdxt(jmt),rdzw(km),fzisop(km)
      real K3(imt,km,jmt,1:3)
!
!-----------------------------------------------------------------------
!     set local constants 
!-----------------------------------------------------------------------
      c1e10 = 1.0e10
      eps   = 1.0e-25
!
!
      do j=2,jmm
        do k=2,km
          do i=2,imm
            e(i,k-1,j,1) = rdxt(j)*p25*c1e10                            &
                  *(tmask(i-1,j,k-1)*tmask(i,j,k-1)*(rhoi(i,k-1,j,xlo)  &
                  -rhoi(i-1,k-1,j,xlo))                                 &
                  +tmask(i,j,k-1)*tmask(i+1,j,k-1)*(rhoi(i+1,k-1,j,xlo) &
                  -rhoi(i,k-1,j,xlo))                                   &
                  +tmask(i-1,j,k)*tmask(i,j,k)*(rhoi(i,k,j,xmd)         &
                  -rhoi(i-1,k,j,xmd))                                   &
                  +tmask(i,j,k)*tmask(i+1,j,k)*(rhoi(i+1,k,j,xmd)       &
                  -rhoi(i,k,j,xmd)))
!
            e(i,k-1,j,2) = rdy*p25*c1e10                                &
                  *(tmask(i,j-1,k-1)*tmask(i,j,k-1)*(rhoi(i,k-1,j,xlo)  &
                  -rhoi(i,k-1,j-1,xlo))                                 &
                  +tmask(i,j,k-1)*tmask(i,j+1,k-1)*(rhoi(i,k-1,j+1,xlo) &
                  -rhoi(i,k-1,j,xlo))                                   &
                  +tmask(i,j-1,k)*tmask(i,j,k)*(rhoi(i,k,j,xmd)         &
                  -rhoi(i,k,j-1,xmd))                                   &
                  +tmask(i,j,k)*tmask(i,j+1,k)*(rhoi(i,k,j+1,xmd)       &
                  -rhoi(i,k,j,xmd)))
!
            e(i,k-1,j,3) = rdzw(k-1)*tmask(i,j,k-1)*tmask(i,j,k)*c1e10  &
                     *(rhoi(i,k-1,j,xlo) - rhoi(i,k,j,xmd))
!
           enddo
        enddo
        k = km
        do i=2,imm
	  e(i,k,j,1) = c0
	  e(i,k,j,2) = c0
	  e(i,k,j,3) = c0
        enddo
      enddo
!
!-----------------------------------------------------------------------
!     compute "K3", using "slmxr" to limit vertical slope of isopycnal
!     to guard against numerical instabilities.  
!-----------------------------------------------------------------------
!
!# ifdef isopycmixspatialvar
!          do i=2,imt-1
!            fxa = c0
!            fxb = sign(1.,e(i,k,j,3))/(abs(e(i,k,j,3))+eps)
!            slope = fxb*sqrt(e(i,k,j,1)**2+e(i,k,j,2)**2)
!            if (slope .le. c0 .and. slope .ge. (-c1/slmxr)) then
!              fxa = p5*(c1+tanh((slope+slopec)/dslope))
!            endif
!            fxc = fxb*fxa
!            K3(i,k,j,1) = -fxc*e(i,k,j,1)
!            K3(i,k,j,2) = -fxc*e(i,k,j,2)
!            K3(i,k,j,3) = fxb*fxb*(e(i,k,j,1)**2+e(i,k,j,2)**2) &
!                          *fxa
!          enddo
!# else # endif
      do j=2,jmm
        do k=1,km
          do i=2,imt-1
!jxz        chkslp = -sqrt(e(i,k,j,1)**2+e(i,k,j,2)**2)*slmxr*dtxsqr(k)
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
!
!-----------------------------------------------------------------------
!     impose zonal boundary conditions at "i"=1 and "imt"
!-----------------------------------------------------------------------
!
!      do j=2,jmm
!        call setbcx (K3(1,1,j,1), imt, km)
!        call setbcx (K3(1,1,j,2), imt, km)
!        call setbcx (K3(1,1,j,3), imt, km)
!      enddo
      return
      end               
!
!
!
!     ===================================
      subroutine setbcx (a, imt, jmtorkm)
!     ===================================
      implicit none
      integer imt,jmtorkm,k
      real a(imt,jmtorkm)
      do k=1,jmtorkm
        a(1,k)   = a(imt-1,k)
        a(imt,k) = a(2,k)
      enddo
      return
      end
!
!     =================
      subroutine isoadv(imt,jmt,imm,jmm,km,athkdf,itn,rdxt,rdyt,dz,pn,fx,fy,rdz0,  &
                        cosu,tmask,adv_vetiso,adv_vntiso,adv_vbtiso,K1,K2,west,    &
                        east,north,south)
!     =================
!
!=======================================================================
!     compute isopycnal transport velocities.
!=======================================================================
!
      implicit none
!      include 'param.h'
      include 'pconst.h'
!      include 'grdvar.h'
!      include 'isopyc.h'
      include 'mpif.h'
!
      integer i,j,k,imt,jmt,imm,jmm,km
      integer jstrt
      real    fxa
!
      real athkdf
      integer itn(imt,jmt),rdxt(jmt),rdyt(jmt),dz(km),pn(imt,jmt)
      real fx(imt,jmt),fy(imt,jmt)
      real rdz0(km),cosu(jmt),tmask(imt,jmt,km)
      real adv_vetiso(imt,km,jmt),adv_vntiso(imt,km,jmt),adv_vbtiso(imt,0:km,jmt)
      real K1(imt,km,jmt,3:3),K2(imt,km,jmt,3:3)
      
      integer west,east,north,south
      
!      common /wisop/ fx(imt,jmt),fy(imt,jmt)
!
!-----------------------------------------------------------------------
!     compute the meridional component of the isopycnal mixing velocity
!     at the center of the northern face of the "t" cells.
!-----------------------------------------------------------------------
!
      do j=2,jmm
        do k=2,km-1
          fxa = -p5*rdz0(k)*athkdf*cosu(j)
          do i=1,imt
            adv_vntiso(i,k,j) = fxa*tmask(i,j,k)*tmask(i,j+1,k)*( &
                                K2(i,k-1,j,3) - K2(i,k+1,j,3))
          enddo
        enddo
      enddo
!
!     consider the top and bottom levels. "K2" is assumed to be zero
!     at the ocean top and bottom.
!
      k = 1
      fxa = -p5*rdz0(k)*athkdf
      do j=2,jmm
        do i=1,imt
          adv_vntiso(i,k,j) = -fxa*tmask(i,j,k)*tmask(i,j+1,k)*cosu(j) &
                              *(K2(i,k,j,3) + K2(i,k+1,j,3))
        enddo
      enddo
!
      do j=2,jmm
        do i=1,imt
          adv_vntiso(i,km,j) = c0
        enddo
      enddo
!
      do j=2,jmm
        do i=1,imt
          k = min(itn(i,j),itn(i,j+1))
          if (k .ne. 0) then
            adv_vntiso(i,k,j) = -p5*rdz0(k)*athkdf*cosu(j)*tmask(i,j,k) &
                         *tmask(i,j+1,k)*(K2(i,k,j,3) + K2(i,k-1,j,3))
          endif
        enddo
      enddo
!
!
!xjin------------------------------
      do j=2,jmm
      do k=1,km
      do i=1,imt
        adv_vntiso(i,k,j) = adv_vntiso(i,k,j) * fy(i,j+1)
      enddo
      enddo
      enddo
!xjin------------------------------
!
!-----------------------------------------------------------------------
!     compute the zonal component of the isopycnal mixing velocity
!     at the center of the eastern face of the "t" grid box.
!-----------------------------------------------------------------------
!
!jxz  jstrt = max(js,jsmw)
      jstrt = 2
!
      do j=jstrt,jmm
        do k=2,km-1
          fxa = -p5*rdz0(k)*athkdf
          do i=1,imm
            adv_vetiso(i,k,j) = fxa*tmask(i,j,k)*tmask(i+1,j,k)  &
                                *(K1(i,k-1,j,3) - K1(i,k+1,j,3))
          enddo
        enddo
      enddo
!
!     consider the top and bottom levels. "K1" is assumed to be zero
!     at the ocean top and bottom.
!
      k = 1
      fxa = -p5*rdz0(k)*athkdf
      do j=jstrt,jmm
        do i=1,imm
          adv_vetiso(i,k,j) = -fxa*tmask(i,j,k)*tmask(i+1,j,k) &
                              *(K1(i,k,j,3)+K1(i,k+1,j,3))
        enddo
      enddo
!
      do j=jstrt,jmm
        do i=1,imm
          adv_vetiso(i,km,j) = c0
        enddo
      enddo
!
      do j=jstrt,jmm
        do i=1,imm
          k = min(itn(i,j),itn(i+1,j))
          if (k .ne. 0) then
            adv_vetiso(i,k,j) = -p5*rdz0(k)*athkdf*tmask(i,j,k)  &
                           *tmask(i+1,j,k)*(K1(i,k,j,3)+K1(i,k-1,j,3))
          endif
        enddo
      enddo
!
!xjin------------------------------
      do j=2,jmm
      do k=1,km
      do i=1,imm
        adv_vetiso(i,k,j) = adv_vetiso(i,k,j) * fx(i+1,j)
      enddo
      enddo
      enddo
!xjin------------------------------
!
!----------------------------------------------------------------------
!     set the boundary conditions
!----------------------------------------------------------------------
!      do j=jstrt,jmm
!        call setbcx (adv_vetiso(1,1,j), imt, km)
!      enddo 
!      call swap_isopy_real3d(adv_vetiso,imt,jmt,km,west,east,north,south)
!
!----------------------------------------------------------------------
!     compute the vertical component of the isopycnal mixing velocity
!     at the center of the bottom face of the "t" cells, using the
!     continuity equation for the isopycnal mixing velocities
!-----------------------------------------------------------------------
!
!      do j=jstrt,jmm
!        do i=1,imt
!          adv_vbtiso(i,0,j) = c0
!        enddo
!      enddo
!
!      do j=jstrt,jmm
!        do k=1,km-1
!          do i=2,imt
!            adv_vbtiso(i,k,j) = dz0(k)*(  &
!            (adv_vetiso(i,k,j) - adv_vetiso(i-1,k,j))*rdxt(j) + &
!            (adv_vntiso(i,k,j) - adv_vntiso(i,k,j-1))*rdyt(j)) &
!          enddo
!        enddo
!      enddo
!
!      do j=jstrt,jmm
!        do k=1,km-1
!          do i=2,imt
!            adv_vbtiso(i,k,j) = adv_vbtiso(i,k,j) + adv_vbtiso(i,k-1,j)
!          enddo
!        enddo
!      enddo
!c
!      do j=jstrt,jmm
!        do i=2,imt
!          adv_vbtiso(i,itn(i,j),j) = c0
!        enddo
!      enddo
!
!-----------------------------------------------------------------------
!     set the boundary conditions
!-----------------------------------------------------------------------
!      do j=jstrt,jmm
!        call setbcx (adv_vbtiso(1,0,j), imt, km+1)
!      enddo
!
!
!xjin-------------------------------------------------
      call isow(adv_vetiso,adv_vntiso,adv_vbtiso,dz,rdxt,rdyt,pn,tmask,itn,  &
                imt,jmt,km,imm,jmm,west,east,north,south)
!xjin-------------------------------------------------
!
!
      return
      end
!
!
!
!     =============================
      subroutine isoflux(tf,mtrace,imt,jmt,imm,jmm,nt,km,kmp1,itn,gravr,rdz0,dz0,rdz,  &
                         rdzw,rdxt,rdyt,rdy,cosu,tmask,t,ahisop,K1,K2,K3,adv_vetiso, &
                         adv_vntiso,adv_vbtiso,fx,fy)
!     =============================
!     isopycnal diffusive tracer fluxes are computed.
!
      implicit none
!      include 'param.h'
      include 'pconst.h'
!      include 'grdvar.h'
!      include 'isopyc.h'
!      include 'prog.h'
!      include 'scalar.h'
      include 'mpif.h'
!
      integer i,j,k,imt,jmt,km,m,imm,jmm,nt,kmp1
      real  fx(imt,jmt),fy(imt,jmt)
!
      real  temp(imt,jmt,km),tf(imt,jmt,km)
      real  wki(imt),wkj(jmt),wkk(kmp1)
      real  fk31(0:km),fk32(0:km)
      integer mtrace
      real    fxa,fxb,fxc,fxe
      
      integer itn(imt,jmt)
      real t(imt,jmt,km,nt,2)
      real ahisop(imt,jmt,km),K1(imt,km,jmt,3:3),K2(imt,km,jmt,3:3),K3(imt,km,jmt,1:3)
      real gravr,rdz0(km),dz0(km),rdz(km),rdzw(km),rdxt(jmt),rdyt(jmt),cosu(jmt),rdy
      real tmask(imt,jmt,km)
      real adv_vetiso(imt,km,jmt),adv_vntiso(imt,km,jmt),adv_vbtiso(imt,0:km,jmt)
!
      wki(imt)  = c0
      wkj(1)    = c0
      wkj(jmt)  = c0
      wkk(1)    = c0
      wkk(kmp1) = c0
!
      do i=1,imt
      do j=1,jmt
      do k=1,km
      temp(i,j,k)=c0
      enddo
      enddo
      enddo
!
      m = mtrace
!
!-----------------------------------------------------------------------
!     first compute the vertical tracer flux "temp" at the northern
!     face of "t" cells.
!-----------------------------------------------------------------------
!
      do k=2,km-1
        do j=1,jmm
          do i=1,imt
            temp(i,j,k)=p25*rdz0(k)*                       &
               (t(i,j+1,k-1,m,taum)-t(i,j+1,k+1,m,taum)    &
               +t(i,j  ,k-1,m,taum)-t(i,j,  k+1,m,taum))   
           enddo
        enddo
      enddo
!
!-----------------------------------------------------------------------
!     now consider the top level, assuming that the surface tracer 
!     values are the same as the ones at "k"=1
!-----------------------------------------------------------------------
!
      k = 1
      do j=1,jmm
        do i=1,imt
          temp(i,j,k) = p25*rdz0(k)*                     &
               (t(i,j+1,k,m,taum)-t(i,j+1,k+1,m,taum)    &
               +t(i,j  ,k,m,taum)-t(i,j,  k+1,m,taum))
        enddo
      enddo
!
!-----------------------------------------------------------------------
!     finally, consider the bottom level. the extrapolative estimator
!     is used to compute the tracer values at the ocean bottom.
!-----------------------------------------------------------------------
!
      do j=1,jmt
        do i=1,imt
          temp(i,j,km) = 0.0
        enddo
      enddo
!
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
!
!-----------------------------------------------------------------------
!     compute of meridional tracer flux
!     add in the effects of the along isopycnal diffusion computed
!     using "K2" component of the tensor and apply land/sea masks
!-----------------------------------------------------------------------
!
!     fzisop(k)=1.0
!
      do 100 k=1,km
      do 100 i=2,imm
        do j=1,jmm
          wkj(j) = ( rdy*(t(i,j+1,k,m,taum)-t(i,j,k,m,taum)) &
                   + K2(i,k,j,3)*temp(i,j,k) )*              &
                   cosu(j)*tmask(i,j,k)*tmask(i,j+1,k)*fy(i,j+1)*p5
        enddo
        do j=2,jmm
          tf(i,j,k) = tf(i,j,k) + ahisop(i,j,k)*rdyt(j)*(wkj(j)-wkj(j-1))
        enddo
100   continue
!
!
!-----------------------------------------------------------------------
!     compute the vertical tracer flux "temp" at the eastern
!     face of "t" cells.
!-----------------------------------------------------------------------
!
      do k=2,km-1
        do j=2,jmm
          do i=1,imm
            temp(i,j,k)=p25*rdz0(k)*                        &
                (t(i+1,j,k-1,m,taum)-t(i+1,j,k+1,m,taum)    &
                +t(i  ,j,k-1,m,taum)-t(i  ,j,k+1,m,taum))
          enddo
        enddo
      enddo
!
!-----------------------------------------------------------------------
!     now consider the top level, assuming that the surface tracer 
!     values are the same as the ones at "k"=1
!-----------------------------------------------------------------------
!
      k = 1
      do j=2,jmm
        do i=1,imm
          temp(i,j,k)=p25*rdz0(k)*                         &
                (t(i+1,j,k,m,taum)-t(i+1,j,k+1,m,taum)     &
                +t(i  ,j,k,m,taum)-t(i  ,j,k+1,m,taum))
        enddo
      enddo
!
!-----------------------------------------------------------------------
!     finally, consider the bottom level. the extrapolative estimator
!     is used to compute the tracer values at the ocean bottom.
!-----------------------------------------------------------------------
!
      do j=2,jmm
        do i=1,imm
          temp(i,j,km) = 0.0
        enddo
      enddo
!
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
!
!-----------------------------------------------------------------------
!     compute of zonal tracer flux 
!     add in the effects of the along isopycnal diffusion computed
!     using "K1" component of the tensor and apply land/sea masks      
!-----------------------------------------------------------------------
!      
!     fzisop(k)=1.0
!     dxtr=const
!
      do 120 k=1,km
      do 120 j=2,jmm
        do i=1,imm
          wki(i) = ( rdxt(j)*(t(i+1,j,k,m,taum)-t(i,j,k,m,taum))  &
                   + K1(i,k,j,3)*temp(i,j,k) )*                   &
                   tmask(i+1,j,k)*tmask(i,j,k)*fx(i+1,j)*p5
        enddo
        do i=2,imm
          tf(i,j,k) = tf(i,j,k) + ahisop(i,j,k)*rdxt(j)*(wki(i)-wki(i-1))
        enddo
120   continue
!
!
!-----------------------------------------------------------------------
!     compute the vertical tracer flux "diff_fbiso" containing the K31
!     and K32 components which are to be solved explicitly. The K33
!     component will be treated semi-implicitly
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!     at ocean surface the flux is set to zero to reflect the no tracer 
!     flux condition. Same condition is also imposed at ocean bottom.
!-----------------------------------------------------------------------
      do k=0,km
      fk31(k)=0.0
      fk32(k)=0.0
      enddo
!
      do 140 j=2,jmm
      do 140 i=2,imm
!
        do k=2, min(itn(i,j)+1,km)
          fk31(k-1) = ahisop(i,j,k)*p25*rdxt(j)*K3(i,k-1,j,1)*tmask(i,j,k)*(  &
       tmask(i-1,j,k  )*(t(i,j,k  ,m,taum)-t(i-1,j,k  ,m,taum))+       &
       tmask(i-1,j,k-1)*(t(i,j,k-1,m,taum)-t(i-1,j,k-1,m,taum))+       &
       tmask(i+1,j,k  )*(t(i+1,j,k,m,taum)-t(i,j  ,k  ,m,taum))+       &
       tmask(i+1,j,k-1)*(t(i+1,j,k-1,m,taum)-t(i,j,k-1,m,taum)) )
        enddo
!
        do k=2, min(itn(i,j)+1,km)
          fk32(k-1) = ahisop(i,j,k)*p25*rdy*K3(i,k-1,j,2)*tmask(i,j,k)*(      &
       tmask(i,j-1,k  )*(t(i,j,k  ,m,taum)-t(i,j-1,k  ,m,taum))+       &
       tmask(i,j-1,k-1)*(t(i,j,k-1,m,taum)-t(i,j-1,k-1,m,taum))+       &
       tmask(i,j+1,k  )*(t(i,j+1,k,m,taum)-t(i,j  ,k  ,m,taum))+       &
       tmask(i,j+1,k-1)*(t(i,j+1,k-1,m,taum)-t(i,j,k-1,m,taum)) )
        enddo
!
        do k=1,km
          tf(i,j,k) = tf(i,j,k) + gravr*rdz(k)*(                       &
                      (fk31(k-1)-fk31(k))+(fk32(k-1)-fk32(k)) )
        enddo
140   continue
!
!
!ifdef gent_mcwilliams
!
!-----------------------------------------------------------------------
!     compute the meridional component of the isopycnal velocity mixing 
!-----------------------------------------------------------------------
!
!
      do k=1,km
        do i=2,imm
          do j=2,jmm
            wkj(j)    = adv_vntiso(i,k,j)*              &
                 (t(i,j+1,k,m,taum)-t(i,j,k,m,taum))
          enddo
          do j=2,jmm
            tf(i,j,k) = tf(i,j,k) - rdyt(j)*(wkj(j)+wkj(j-1))*p5
          enddo
        enddo
      enddo
!
!
!-----------------------------------------------------------------------
!     compute the zonal component of the isopycnal velocity mixing 
!-----------------------------------------------------------------------
!
      do k=1,km
        do j=2,jmm
          do i=1,imm
            wki(i)    = adv_vetiso(i,k,j)*           &
              (t(i+1,j,k,m,taum)-t(i,j,k,m,taum))
          enddo
          do i=2,imm
            tf(i,j,k) = tf(i,j,k) - rdxt(j)*(wki(i)+wki(i-1))*p5
          enddo
        enddo
      enddo
!
!-----------------------------------------------------------------------
!     compute the vertical component of the isopycnal velocity mixing 
!-----------------------------------------------------------------------
!
      do j=2,jmm
        do i=2,imm
          do k=2, min(itn(i,j)+1,km)
            wkk(k)  = adv_vbtiso(i,k-1,j)*                &
                 (t(i,j,k-1,m,taum)-t(i,j,k,m,taum))
          enddo
        do k=1,itn(i,j)
          tf(i,j,k) = tf(i,j,k) - (wkk(k)+wkk(k+1))*rdz(k)*p5
          enddo
        enddo
      enddo
!
      return
      end
!
!     ======================
      subroutine isow(u,v,w,dz,rdxt,rdyt,pn,tmask,itn,imt,jmt,km,imm,jmm,west,east,  &
                      north,south)
!     ======================
!     calculate vertical mass advection: Ps*dz/dt
!     upward is positive
!
      implicit none
      include 'pconst.h'
      include 'mpif.h'
!
      integer imt,jmt,km,imm,jmm,i,j,k
      integer itn(imt,jmt)
      real tmask(imt,jmt,km),pn(imt,jmt)
      real u(imt,km,jmt),v(imt,km,jmt),w(imt,0:km,jmt)
      real dp(imt,jmt),rdxt(jmt),rdyt(jmt),dz(km)
      
      integer west,east,north,south
!
      do j=1,jmt
      do i=1,imt
      w(i,0,j) = c0
      dp(i,j)  = c0
      enddo
      enddo
!
      do j=2,jmm
      do k=1,km
      do i=2,imm
      w(i,k,j) = dz(k)*((u(i,k,j)-u(i-1,k,j))*rdxt(j) +  &
                        (v(i,k,j)-v(i,k,j-1))*rdyt(j))
      dp(i,j)  = dp(i,j) - w(i,k,j)
      enddo
      enddo
      enddo
!
      do j=2,jmm
      do k=1,km-1
      do i=2,imm
      w(i,k,j)=(w(i,k-1,j)+dz(k)*dp(i,j)/pn(i,j)+w(i,k,j))  &
               *tmask(i,j,k+1)
      enddo
      enddo
      enddo
!
      do j=2,jmm
        do i=2,imt
          w(i,itn(i,j),j) = c0
        enddo
      enddo
!
!
!      do k=1,km
!      do j=2,jmm
!      w(1  ,k,j) = w(imm,k,j)
!      w(imt,k,j) = w(2  ,k,j)
!      enddo
!      enddo
!      call swap_isopyw_real3d(w,imt,jmt,km,west,east,north,south)
!
      return
      end
