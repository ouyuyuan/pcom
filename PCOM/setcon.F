c     =================
      subroutine setcon
c     =================
c     set scalar quantities
c
      implicit none
      include 'param.h'
      include 'pconst.h'
      include 'scalar.h'
      include 'calendar.h'
c
c---------------------------------------------------------------------
c     calendar
c---------------------------------------------------------------------
      data monname /'jan','feb','mar','apr','may','jun',
     &              'jul','aug','sep','oct','nov','dec'/
c     data daymd   /16,14,16,15,16,15,16,16,15,16,15,16/
c     data daypm   /31,28,31,30,31,30,31,31,30,31,30,31/
      data daymd   /15,15,15,15,15,15,15,15,15,15,15,15/
      data daypm   /30,30,30,30,30,30,30,30,30,30,30,30/
c
c
c---------------------------------------------------------------------
c     physical parameters
c---------------------------------------------------------------------
c     am      = constant lateral viscosity coeff for momentum
c     ah      = constant lateral diffusion coeff for tracers
c     kappa_m = constant vertical viscosity coefficient (cm**2/sec) 
c     kappa_h = constant vertical diffusion coefficient (cm**2/sec)
c     cdbot   = parameter for bottom drag
c
      am      = 3.0e8
      ah      = 1.0e7
      kappa_m = 1.0
      kappa_h = 0.3
      cdbot   = 2.6e-3
c
#ifdef boussinesq
      gravr = grav
#else
      gravr = grav*rho_0
#endif
c
c     change the unit of pressure from dynes/cm**2 to decibars
c     1 dynes/cm**2 = 0.1 N/m**2 = 1.0e-3 mbar = 1.0e-5 decibars
c
      decibar = 1.0e-5
c
c     deltax is used to compute {partial rho}/{partial x}, x=p,t,s
      deltap  = 1.0e-2
      deltat  = 1.0e-2
      deltas  = 1.0e-3
      rdeltap = p5/deltap
      rdeltat = p5/deltat
      rdeltas = p5/deltas
c
c
c---------------------------------------------------------------------
c     newtonia restoring coefficient for heat flux
c---------------------------------------------------------------------
c     gamma_t = grav * ddd/cp
c     ddd     = 40W/m**2/K
c     cp      = cpecific heat capacity of sea water, 3901 J/kg/K
c     1.0e-1  = unit change
c
      gamma_t = grav*40.0/3901.0*1.0e-1
c
c
c---------------------------------------------------------------------
c     newtonia restoring time scale for salinity
c---------------------------------------------------------------------
      gamma_s = r120*secday
c
c
c---------------------------------------------------------------------
c     asselin temporal filter parameter
c---------------------------------------------------------------------
      afb1 = 0.25
      afc1 = 0.25
      aft1 = 0.25
c
c
c---------------------------------------------------------------------
c     time-step of integration
c---------------------------------------------------------------------
c     dtsf:   time step for adv_1                 (input unit: min)
c     dtuv:   time step for adv_2                 (input unit: hour)
c     dtts:   time step for tracer                (input unit: hour)
c
      dtsf = 2.0
      dtuv = 1.0
      dtts = 24.0
c
c
c---------------------------------------------------------------------
c     parameters for model's output
c---------------------------------------------------------------------
c     io_tsuvp  = interval for output time average variables (unit year)
c     io_restr  = interval for output restr.data             (unit year)
c
      io_tsuvp  = 1
      io_restr  = 1
c
c
      return
      end
