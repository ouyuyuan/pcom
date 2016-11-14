!====================== include file "prog.h" ==========================
!
!     variables used for prognostic and diagnostic equations
!
!     tau  = tau   time level
!     taum = tau-1 time level
!
      integer tau,taum
      parameter (tau=1,taum=2)
!
!-----------------------------------------------------------------------
!     arrays for prognostic equations:
!-----------------------------------------------------------------------
!
!     t   = tracer   (nt=1 for temperature, 2 for salinity)
!     up  = u * sqrt (pbt)                                 (cm/s * spbt)
!     vp  = v * sqrt (pbt)                                 (cm/s * spbt)
!     pbt = pressure difference between bottom and surface (dynes/cm**2)
!     upb = barotropic component of up                     (cm/s * spbt)
!     vpb = barotropic component of vp                     (cm/s * spbt)
!
!     note: temperature is potential temperature in degrees Celsius and
!     salinity is in "model units", the deviation from 0.035 g of salt/g
!     of water, assuming a water density of 1 g/cm**3
!
      real t,up,vp,pbt,upb,vpb
!
      common /vprog/ t(imt,jmt,km,nt,2)
      common /vprog/ up(imt,jmt,km,2),vp(imt,jmt,km,2)
      common /vprog/ pbt(imt,jmt,2)
      common /vprog/ upb(imt,jmt,2),vpb(imt,jmt,2)
!
!
!-----------------------------------------------------------------------
!     arrays for diagnostic equations:
!-----------------------------------------------------------------------
!
!     spbt  = sqrt(pbt)
!     rho   = density
!     w     = pbt * dz/dt at "t" cell or dz/dt at "u" cell
!
      real spbt,rho,w
!
      common /vdiag/ spbt(imt,jmt)
      common /vdiag/ rho(imt,jmt,km)
      common /vdiag/ w(imt,jmt,kmp1)
!
!
!-----------------------------------------------------------------------
!     woring arrays
!-----------------------------------------------------------------------
!     du    = advection + diffusion + pressure gradients for up
!     dv    = advection + diffusion + pressure dradients for vp
!     dub   = vertical integration of du
!     dvb   = vertical integration of dv
!
      real du,dv,dub,dvb
      common /vwork/ du(imt,jmt,km),dv(imt,jmt,km)
      common /vwork/ dub(imt,jmt),dvb(imt,jmt)
!
!
!     diffu = viscosity term for up eq.
!     diffv = viscosity term for vp eq.
!     pt    = {partial rho} over {partial reference temperature}
!     ps    = {partial rho} over {partial reference salinity}
!
      real diffu,diffv,pt,ps
      common /vwork/ diffu(imt,jmt,km),diffv(imt,jmt,km)
      common /vwork/ pt(imt,jmt,km),ps(imt,jmt,km)
!
!
!     pmup  = pbt at tau+1/2 time level (baroclinic eq)
!     pmum  = pbt at tau-1/2 time level  ''
!     pmtp  = pbt at tau+1/2 time level (tracer eq)
!     pmtm  = pbt at tau-1/2 time level  ''
!     ump   = u   at tau+1/2 time level  ''
!     umm   = u   at tau-1/2 time level  ''
!     vmp   = v   at tau+1/2 time level  ''
!     vmm   = v   at tau-1/2 time level  ''
!
      real pmup,pmum,pmtp,pmtm,ump,umm,vmp,vmm
      common /vwork/ pmup(imt,jmt),pmum(imt,jmt)
      common /vwork/ pmtp(imt,jmt),pmtm(imt,jmt)
      common /vwork/ ump(imt,jmt,km),umm(imt,jmt,km)
      common /vwork/ vmp(imt,jmt,km),vmm(imt,jmt,km)
!
!
!-----------------------------------------------------------------------
!     arrays rlated to pressure gradents
!-----------------------------------------------------------------------
      real rhodp,pax,pay
      real pbxn,pbxs,pcxn,pcxs,pdxn,pdxs
      real pbye,pbyw,pcye,pcyw,pdye,pdyw
      real phibx,phiby
!
      common /vprsg/ rhodp(imt,jmt,km)
      common /vprsg/ pax(imt,jmt),pay(imt,jmt)
      common /vprsg/ pbxn(imt,jmt),pbxs(imt,jmt)
      common /vprsg/ pcxn(imt,jmt),pcxs(imt,jmt)
      common /vprsg/ pdxn(imt,jmt),pdxs(imt,jmt)
      common /vprsg/ pbye(imt,jmt),pbyw(imt,jmt)
      common /vprsg/ pcye(imt,jmt),pcyw(imt,jmt)
      common /vprsg/ pdye(imt,jmt),pdyw(imt,jmt)
      common /vprsg/ phibx(imt,jmt),phiby(imt,jmt)
!
!#ifdef boussinesq
      real           fixp
      common /vprsg/ fixp(imt,jmt,km)
!#endif
