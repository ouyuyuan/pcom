c    +                                                               +
c    +===============================================================+
c    +             PCOM and BCOM in eta-coordinates                  +
c    +===============================================================+
c    +                                                               +
      implicit none
      include 'param.h'
      include 'pconst.h'
      include 'scalar.h'
      include 'calendar.h'
      include 'prog.h'
      include 'grdvar.h'
c
      integer  mode_c,mode_t
c
      namelist /contrl/ runlen,restrt
      namelist /tsteps/ dtts,dtuv,dtsf
      namelist /mixing/ am,ah,kappa_m,kappa_h,cdbot
      namelist /filter/ afb1,afc1,aft1
      namelist /io/     io_tsuvp,io_restr
c
c     set scalar quantities
      call setcon
c
c
c     set resolution, t/u mask and the j-depended parameters
      call grdvar
c
c
c     read in model's control parameters
      open(15,file='namelist')
      read(15,contrl)
      read(15,tsteps)
      read(15,mixing)
      read(15,filter)
      read(15,io)
      close(15)
c
c
c     initialization
      call inirun
c
#ifdef gm90
      call isopyi
#endif
c
c-----------------------------------------------------------------------
c     monthly cycle
c-----------------------------------------------------------------------
c
10    year = 1 + (month-1)/12
      mth  = month - (month-1)/12*12
c
c
c     set euler forward/backward scheme at the beginning of every month
c
      leapfrog_t = .false.
      leapfrog_c = .false.
      leapfrog_b = .false.
      euler_back = .true.
c    
c
c-----------------------------------------------------------------------
c     daily cycle
c-----------------------------------------------------------------------
c
      do 20 day=1,daypm(mth)
c
c
c     compute coefficients for calculation of potential density anomaly
      call rho_ref
c
c
c     interpolate the observed monthly mean data
      call interp
c
c
c-----------------------------------------------------------------------
c     for thermal cycle
c-----------------------------------------------------------------------
      do 30 mode_t = 1,nss
c
c
c     calculate  baroclinic pressure and the relavant variables
      call readyt
c
c
c-----------------------------------------------------------------------
c     for baroclinic & barotropic cycle
c-----------------------------------------------------------------------
      do 40 mode_c = 1,ncc
c
c     calculate momentum advection, diffusion & their vertical integrals
      call readyc
c
c
c     prediction of barotropic mode
      call barotr
c
c
c     prediction of baroclinic mode
      call bclinc
c
#ifdef smth
      print*, ' check I&J boundary for smooth!!!!!!'
      stop 111111
      call smth2(upb,umask,jstn,jedn,jsts,jeds)
      call smth2(vpb,umask,jstn,jedn,jsts,jeds)
#endif
c
40    continue
c
c
c     prediction of temperature and salinity
      call tracer
c
c
c     compute sea ice
c     call seaice
c
c
c     convective adjustment if unstable stratification ocurs
      call convect
c
c
30    continue
c
c
c#ifdef accelerate
c      if(day.ne.daypm(mth)) goto 50
c#endif
#ifdef smth
      print*, ' check I&J boundary for smooth!!!!!!'
      stop 111333
      call smth3(up,umask,jstn,jedn,jsts,jeds)
      call smth3(vp,umask,jstn,jedn,jsts,jeds)
#endif
c
c50    continue
c
c     output
      call diag
c
20    continue
c
      month  = month  + 1
      runlen = runlen - 1
c
      if(runlen.gt.0) goto 10
c
      end
