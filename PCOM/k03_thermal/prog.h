c====================== include file "prog.h" ==========================
c
c     variables used for prognostic and diagnostic equations
c
c     tau  = tau   time level
c     taum = tau-1 time level
c
      integer tau,taum
      parameter (tau=1,taum=2)
c
c-----------------------------------------------------------------------
c     arrays for prognostic equations:
c-----------------------------------------------------------------------
c
c     t   = tracer   (nt=1 for temperature, 2 for salinity)
c     up  = u * sqrt (pbt)                                 (cm/s * spbt)
c     vp  = v * sqrt (pbt)                                 (cm/s * spbt)
c     pbt = pressure difference between bottom and surface (dynes/cm**2)
c     upb = barotropic component of up                     (cm/s * spbt)
c     vpb = barotropic component of vp                     (cm/s * spbt)
c
c     note: temperature is potential temperature in degrees Celsius and
c     salinity is in "model units", the deviation from 0.035 g of salt/g
c     of water, assuming a water density of 1 g/cm**3
c
      real t,up,vp,pbt,upb,vpb
c
      common /vprog/ t(imt,jmt,km,nt,2)
      common /vprog/ up(imt,jmt,km,2),vp(imt,jmt,km,2)
      common /vprog/ pbt(imt,jmt,2)
      common /vprog/ upb(imt,jmt,2),vpb(imt,jmt,2)
c
c
c-----------------------------------------------------------------------
c     arrays for diagnostic equations:
c-----------------------------------------------------------------------
c
c     spbt  = sqrt(pbt)
c     rho   = density
c     w     = pbt * dz/dt at "t" cell or dz/dt at "u" cell
c
      real spbt,rho,w
c
      common /vdiag/ spbt(imt,jmt)
      common /vdiag/ rho(imt,jmt,km)
      common /vdiag/ w(imt,jmt,kmp1)
c
c
c-----------------------------------------------------------------------
c     woring arrays
c-----------------------------------------------------------------------
c     du    = advection + diffusion + pressure gradients for up
c     dv    = advection + diffusion + pressure dradients for vp
c     dub   = vertical integration of du
c     dvb   = vertical integration of dv
c
      real du,dv,dub,dvb
      common /vwork/ du(imt,jmt,km),dv(imt,jmt,km)
      common /vwork/ dub(imt,jmt),dvb(imt,jmt)
c
c
c     diffu = viscosity term for up eq.
c     diffv = viscosity term for vp eq.
c     pt    = {partial rho} over {partial reference temperature}
c     ps    = {partial rho} over {partial reference salinity}
c
      real diffu,diffv,pt,ps
      common /vwork/ diffu(imt,jmt,km),diffv(imt,jmt,km)
      common /vwork/ pt(imt,jmt,km),ps(imt,jmt,km)
c
c
c     pmup  = pbt at tau+1/2 time level (baroclinic eq)
c     pmum  = pbt at tau-1/2 time level  ''
c     pmtp  = pbt at tau+1/2 time level (tracer eq)
c     pmtm  = pbt at tau-1/2 time level  ''
c     ump   = u   at tau+1/2 time level  ''
c     umm   = u   at tau-1/2 time level  ''
c     vmp   = v   at tau+1/2 time level  ''
c     vmm   = v   at tau-1/2 time level  ''
c
      real pmup,pmum,pmtp,pmtm,ump,umm,vmp,vmm
      common /vwork/ pmup(imt,jmt),pmum(imt,jmt)
      common /vwork/ pmtp(imt,jmt),pmtm(imt,jmt)
      common /vwork/ ump(imt,jmt,km),umm(imt,jmt,km)
      common /vwork/ vmp(imt,jmt,km),vmm(imt,jmt,km)
c
c
c-----------------------------------------------------------------------
c     arrays rlated to pressure gradents
c-----------------------------------------------------------------------
      real rhodp,pax,pay
      real pbxn,pbxs,pcxn,pcxs,pdxn,pdxs
      real pbye,pbyw,pcye,pcyw,pdye,pdyw
      real phibx,phiby
c
      common /vprsg/ rhodp(imt,jmt,km)
      common /vprsg/ pax(imt,jmt),pay(imt,jmt)
      common /vprsg/ pbxn(imt,jmt),pbxs(imt,jmt)
      common /vprsg/ pcxn(imt,jmt),pcxs(imt,jmt)
      common /vprsg/ pdxn(imt,jmt),pdxs(imt,jmt)
      common /vprsg/ pbye(imt,jmt),pbyw(imt,jmt)
      common /vprsg/ pcye(imt,jmt),pcyw(imt,jmt)
      common /vprsg/ pdye(imt,jmt),pdyw(imt,jmt)
      common /vprsg/ phibx(imt,jmt),phiby(imt,jmt)
c
c#ifdef boussinesq
      real           fixp
      common /vprsg/ fixp(imt,jmt,km)
c#endif
