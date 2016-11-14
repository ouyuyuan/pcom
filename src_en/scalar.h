!====================== include file "scalar.h" ========================
!
!     ==================
!     scalar quantities:
!     ==================
!
!     dtts    = time step for tracers (in seconds)
!     dtuv    = time step for solving baroclinic eq (in seconds)
!     dtsf    = time step for solving barotropic eq (in seconds)
!     c2dtts  = 2*dtts
!     c2dtuv  = 2*dtuv
!     c2dtsf  = 2*dtsf
!     nss     = number of time steps for tracer eq
!     ncc     = number of time steps for baroclinic eq
!     nbb     = number of time steps for barotropic eq
!
!     af#     = asselin temporal filter parameters
!     jst/jed = starting/ending latitude for filter
!     decibar = unit factor for density calculation
!     delta#  = delta #, for calculation of {partial rho}/{partial #}
!
!     am      = constant lateral viscosity coeff for momentum
!     ah      = constant lateral diffusion coeff for tracers
!     kappa_m = constant vertical viscosity coefficient (cm**2/sec)
!     kappa_h = constant vertical diffusion coefficient (cm**2/sec)
!     gamma_t = parameter for calculation of surface heat flux
!     gamma_s = parameter for surface salinity b.c.
!     cdbot   = parameter for bottom drag
!
!     ==================
!     control variables:
!     ==================
!
!     runlen     = integration length (in months)
!
!     restrt     = (true,false) indicates that this run is a
!                  (restart, start from initial conditions)
!
!     euler_back = (false,true) on a (eular forward, backward) time step
!                  in "adv-1"
!
!     leapfrog_b = (false,true) on a (eular backward, normal leapfrog)
!                  time step in "adv-1"
!
!     leapfrog_c = (false,true) on a (eular forward, normal leapfrog)
!                  time step in "adv-2"
!
!     leapfrog_t = (false,true) on a (eular forward, normal leapfrog)
!                  time step in "tracer"
!
!     io_tsuvp   = interval for annual/monthly mean output (in year)
!     io_restr   = interval for saving data for restarting (in year)
!
      integer nss,ncc,nbb,jstn,jedn,jsts,jeds
      real    onbb,oncc,onbc
      real    dtts,dtuv,dtsf,c2dtts,c2dtuv,c2dtsf
      real    afb1,afc1,aft1,afb2,afc2,aft2
      real    decibar,deltap,rdeltap,deltat,deltas,rdeltat,rdeltas
      real    am,ah,kappa_m,gamma_t,gamma_s,cdbot,gravr
!      real   kappa_h(imt,jmt,km)
!      common /scalar/ nss,ncc,nbb,onbb,oncc,onbc
!      common /scalar/ dtts,dtuv,dtsf,c2dtts,c2dtuv,c2dtsf
!      common /scalar/ afb1,afc1,aft1,afb2,afc2,aft2,jstn,jedn,jsts,jeds
!      common /scalar/ decibar
!      common /scalar/ deltap,rdeltap,deltat,deltas,rdeltat,rdeltas
!      common /scalar/ am,ah,kappa_m,gamma_t,gamma_s,cdbot
!      common /scalar/ kappa_h
!      common /scalar/ gravr
!
      real    runlen
      integer io_tsuvp,io_restr
      logical euler_back,leapfrog_b,leapfrog_c,leapfrog_t
      logical restrt
!
!      common /switcr/ runlen,restrt
!      common /switcr/ euler_back,leapfrog_b,leapfrog_c,leapfrog_t
!      common /switcr/ io_tsuvp,io_restr
