!====================== include file "param.h" =========================
!
!     main parameter file which sets ocean characteristics:
!
!                                                                       
!     imt = number of grid points in longitudinal direction          
!     jmt = number of grid points in latitudinal  direction           
!     km  = number of the sigma levels
!     nt  = number of tracers (temperature, salinity, ...)
!                                                                       
      integer imt,jmt,km,nt,imm,jmm,kmp1,kmm1,i,j,k,n,m
      real dlam,dphi,phis
      integer simt,sjmt,bcfrec
      integer snbc,gm90,implicitvmix,asselin_b,asselin_c,asselin_t,smth,unesco
      integer boussinesq
      real slmxr,ahisop,athkdf,smth_start_nlat,smth_start_slat
!
!      parameter (imt=362,jmt=141,km=60)
!      parameter (nt=2)
!      parameter (imm=imt-1,jmm=jmt-1,kmp1=km+1,kmm1=km-1)
!      parameter (dlam=1.0,dphi=1.0,phis=-70.0)

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
      integer fam,fah,fkm,fkh
      real    am_c,ah_c,km_c,kh_c,gamma_t,gamma_s,cdbot,gravr
      real    runlen
      integer io_tsuvp,io_restr
      logical euler_back,leapfrog_b,leapfrog_c,leapfrog_t
      logical restrt
      
!                       calendar specification arrays
!
!-----------------------------------------------------------------------
!
!     monname = character names of months
!     daypm   = array of month lengths in days
!     daymd   = array of the middle date in monthes
!     month   = cumulative monthes of the integration
!     year    = the year of model
!     mth     = calendar month (from 1 to 12)
!     day     = calendar day (from 1 to daypm(mth))
!
!-----------------------------------------------------------------------
      integer   daypm(12),daymd(12)
      integer   month,year,mth,day
      character monname(12)*3
!

      data monname /'jan','feb','mar','apr','may','jun', &
                    'jul','aug','sep','oct','nov','dec'/
      data daymd   /15,15,15,15,15,15,15,15,15,15,15,15/
      data daypm   /30,30,30,30,30,30,30,30,30,30,30,30/