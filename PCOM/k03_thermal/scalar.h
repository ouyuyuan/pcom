c====================== include file "scalar.h" ========================
c
c     ==================
c     scalar quantities:
c     ==================
c
c     dtts    = time step for tracers (in seconds)
c     dtuv    = time step for solving baroclinic eq (in seconds)
c     dtsf    = time step for solving barotropic eq (in seconds)
c     c2dtts  = 2*dtts
c     c2dtuv  = 2*dtuv
c     c2dtsf  = 2*dtsf
c     nss     = number of time steps for tracer eq
c     ncc     = number of time steps for baroclinic eq
c     nbb     = number of time steps for barotropic eq
c
c     af#     = asselin temporal filter parameters
c     jst/jed = starting/ending latitude for filter
c     decibar = unit factor for density calculation
c     delta#  = delta #, for calculation of {partial rho}/{partial #}
c
c     am      = constant lateral viscosity coeff for momentum
c     ah      = constant lateral diffusion coeff for tracers
c     kappa_m = constant vertical viscosity coefficient (cm**2/sec)
c     kappa_h = constant vertical diffusion coefficient (cm**2/sec)
c     gamma_t = parameter for calculation of surface heat flux
c     gamma_s = parameter for surface salinity b.c.
c     cdbot   = parameter for bottom drag
c
c     ==================
c     control variables:
c     ==================
c
c     runlen     = integration length (in months)
c
c     restrt     = (true,false) indicates that this run is a
c                  (restart, start from initial conditions)
c
c     euler_back = (false,true) on a (eular forward, backward) time step
c                  in "adv-1"
c
c     leapfrog_b = (false,true) on a (eular backward, normal leapfrog)
c                  time step in "adv-1"
c
c     leapfrog_c = (false,true) on a (eular forward, normal leapfrog)
c                  time step in "adv-2"
c
c     leapfrog_t = (false,true) on a (eular forward, normal leapfrog)
c                  time step in "tracer"
c
c     io_tsuvp   = interval for annual/monthly mean output (in year)
c     io_restr   = interval for saving data for restarting (in year)
c
      integer nss,ncc,nbb,jstn,jedn,jsts,jeds
      real    onbb,oncc,onbc
      real    dtts,dtuv,dtsf,c2dtts,c2dtuv,c2dtsf
      real    afb1,afc1,aft1,afb2,afc2,aft2
      real    decibar,deltap,rdeltap,deltat,deltas,rdeltat,rdeltas
      real am,ah,kappa_m,kappa_h(imt,jmt,km),gamma_t,gamma_s,cdbot,gravr
c
      common /scalar/ nss,ncc,nbb,onbb,oncc,onbc
      common /scalar/ dtts,dtuv,dtsf,c2dtts,c2dtuv,c2dtsf
      common /scalar/ afb1,afc1,aft1,afb2,afc2,aft2,jstn,jedn,jsts,jeds
      common /scalar/ decibar
      common /scalar/ deltap,rdeltap,deltat,deltas,rdeltat,rdeltas
      common /scalar/ am,ah,kappa_m,kappa_h,gamma_t,gamma_s,cdbot
      common /scalar/ gravr
c
      real    runlen
      integer io_tsuvp,io_restr
      logical euler_back,leapfrog_b,leapfrog_c,leapfrog_t
      logical restrt
c
      common /switcr/ runlen,restrt
      common /switcr/ euler_back,leapfrog_b,leapfrog_c,leapfrog_t
      common /switcr/ io_tsuvp,io_restr
