c======================= include file "calendar.h"======================
c
c                       calendar specification arrays
c
c-----------------------------------------------------------------------
c
c     monname = character names of months
c     daypm   = array of month lengths in days
c     daymd   = array of the middle date in monthes
c     month   = cumulative monthes of the integration
c     year    = the year of model
c     mth     = calendar month (from 1 to 12)
c     day     = calendar day (from 1 to daypm(mth))
c
c-----------------------------------------------------------------------
      integer   daypm,daymd
      integer   month,year,mth,day
      character monname(12)*3
c
      common /calend/ monname
      common /calend/ daypm(12), daymd(12)
      common /calend/ month,year,mth,day


c====================== include file "cvbc.h" ==========================
c
c     vertical boundary condition variables:
c
c     bcf = monthly mean surface forcing fileds. use linear
c           interpolation to produce:
c
c     bcu = bcf(1) :   sea surface zonal windstres       (dynes/cm**2)
c     bcv = bcf(2) :   sea surface meridional windstres  (dynes/cm**2)
c     bct = bcf(3) :   sea surface air temperature       (celsius)
c     bcp = bcf(4) :   sea surface air presure           (dynes/cm**2)
c     bcs = bcf(5) :   sea surface salinity              (model unit)
c     emp = bcf(6) :   rate of evaporation minus precipitation (cm/s)
c     ddd = bcf(7) :   coefficient for calculation of HF (w/m2/c)
c
      real*4 bcf,bcu,bcv,bct,bcp,bcs,emp,ddd
c
      common /cvbc/ bcf(imt,jmt,12,7)
      common /cvbc/ bcu(imt,jmt),bcv(imt,jmt)
      common /cvbc/ bct(imt,jmt),bcp(imt,jmt),bcs(imt,jmt)
      common /cvbc/ emp(imt,jmt),ddd(imt,jmt)


c====================== include file "diag.h" ==========================
c
c     tmn    = monthly/annual averaged potential temperature (c)
c     smn    = monthly/annual averaged salinity            (model unit)
c     pmn    = monthly/annual averaged pbt                 (dynes/cm**2)
c     umn    = monthly/annual averaged up
c     vmn    = monthly/annual averaged vp
c
      real tmn,smn,pmn,umn,vmn
c
      common /cdiag/ tmn(imt,jmt,km),smn(imt,jmt,km)
      common /cdiag/ umn(imt,jmt,km),vmn(imt,jmt,km)
      common /cdiag/ pmn(imt,jmt)

c====================== include file "grdvar.h" ========================
c
c     variables which are functions of the grid
c
c     z0    = depth at the center of "t" and "u" grid cells
c     dz0   = layer-thickness in z-coordinate
c     z     = i.e. "eat",  at the center of "t" and "u" grid cells
c     dz    = thickness of "t" and "u" grid cells
c     rdz   = reciprocal of dz
c     rdzw  = reciprocal of thickness of "w" cell (in cm)
c     phib  = potential high at bottom
c     pn    = constant bottom pressure (PCOM); (-1)*phib (BCOM)
c     zu    = int(dz) of "u" grid cell
c     rzu   = reciprocal of zu
c
      real z0,dz0,z,dz,rdz,rdzw,phib,pn,rzu,zu
      common /grdv/ z0(km),dz0(km),z(km),dz(km),rdz(km),rdzw(km)
      common /grdv/ phib(imt,jmt),pn(imt,jmt),rzu(imt,jmt),zu(imt,jmt)
c
c
c     tmask = tracer cell land/sea mask   = (0.0, 1.0) on (land, sea)
c     umask = velocity cell land/sea mask = (0.0, 1.0) on (land, sea)
      real tmask,umask
      integer itn,ivn
      common /grdv/ tmask(imt,jmt,km),umask(imt,jmt,km)
      common /grdv/ itn(imt,jmt),ivn(imt,jmt)
c
c
c     cost  = cosine of "t" grid point latitude
c     cosu  = cosine of "u,v" grid point latitude
c     ff    = 2*omega*sine(j)
c     rdxt  = reciprocal of longitudinal width of "t" grid box (in cm)
c     rdxu  = reciprocal of longitudinal width of "u" grid box (in cm)
c     rdy   = reciprocal of latitudinal height
c     rdyt  = reciprocal of (latitudinal height * cost)
c     rdyu  = reciprocal of (latitudinal height * cosu)
c     dxdyt = dxt * dy
c     dxdyu = dxu * dy
c     area  = global ocean area
c     cv1   = square of sinu/cosu/a
c     cv2   = sinu/cosu * rdxu/radius
c     sdxt  = 1/2 square of rdxt
c     sdxu  = 1/2 square of rdxu
c     r1#   = coefficients for calculation of diffusion
c     ep#   = coefficients for semi-implicitly handle of Coriolis term
c
      real lat,cost,cosu,ff
      real rdxt,rdxu,rdyt,rdyu,rdy
      real dxdyt,dxdyu,area
      real cv1,cv2,sdxt,sdxu,r1a,r1b,r1c,r1d
      real ebea,ebeb,ebla,eblb,epea,epeb,epla,eplb
c
      common /grdv/ lat(jmt),cost(jmt),cosu(jmt),ff(jmt)
      common /grdv/ rdxt(jmt),rdxu(jmt),rdyt(jmt),rdyu(jmt),rdy
      common /grdv/ dxdyt(jmt),dxdyu(jmt),area
      common /grdv/ cv1(jmt),cv2(jmt),sdxt(jmt),sdxu(jmt)
      common /grdv/ r1a(jmt),r1b(jmt),r1c(jmt),r1d(jmt)
      common /grdv/ ebea(jmt),ebeb(jmt),ebla(jmt),eblb(jmt)
      common /grdv/ epea(jmt),epeb(jmt),epla(jmt),eplb(jmt)
c


c=========================isopyc.f=========================================
c     isopycnal diffusion variables:
c
c     ahisop = isopycnal tracer diffusivity (cm**2/sec) 
c     athkdf = isopycnal thickness diffusivity (cm**2/sec)
c     dptlim = depth limits for the reference pressure levels (in cm). 
c              the mid point of the two neighboring "dptlim" elements is
c              used as the reference pressure (depth) for that interval,
c              i.e.,
c
c             reference pressure level        reference depth
c              ----------------------- -----------------------------
c                          1           0.5*(dptlim(1)+dptlim(2))
c                          2           0.5*(dptlim(2)+dptlim(3))
c                          .           .
c                        nrpl          0.5*(dptlim(nrpl)+dptlim(nrpl+1))
c
c              REMARK: the first and the last elements of "dptlim" must
c                      be the depth at the top (0m) and the maximum bottom
c                      depth, respectively. also, the elements of "dptlim"
c                      must be in increasing order.
c     e      = scratch array 
c     K1     = (1,3) component of the isopycnal mixing tensor computed
c              at the center of the eastern face of the "t" grid cell
c     K2     = (2,3) component of the isopycnal mixing tensor computed
c              at the center of the northern face of the "t" grid cell
c     K3     = (3,.) components of the isopycnal mixing tensor 
c              computed at the top face of the "t" grid cell
c               (,,1) --> (3,1) component
c               (,,2) --> (3,2) component
c               (,,3) --> (3,3) component 
c              REMARK: (1,1) and (2,2) components of the tensor are
c                      assumed to be unity. also, (1,2) and (2,1)
c                      components are assumed to be zero.
c     fzisop = function containing the vertical variation of the isopycnal
c              diffusion coefficient. "fzisop" multiplies "ahisop".
c     kisrpl = isopycnal reference pressure levels for the "t" grid 
c              point levels computed based on the depth (pressure) at the
c              "t" grid points
c     krplin = indices indicating the location of the reference pressure 
c              depth in the 20-level table of polynomial expansion
c              variables
c     slmxr  = reciprocal of maximum slope of isopycnals allowed in mixing
c              scheme to prevent excessively large vertical mixing that
c              could create numerical instabilities. furthermore, the
c              form of the isopycnal diffusion tensor incorporates the
c              assumption that horizontal density gradients are much
c              smaller than vertical gradients. "slmxr" is also used to
c              satisfy this assumption. a value of 100 for "slmxr"
c              translates to a slope of 1:100.
c
c     rhoi    = potential density at "t" cell centers
c
c# ifdef isopycmixspatialvar
c     dciso1 = isopycnal tracer diffusivity coefficients modified based
c              on the slopes of the isopycnal surfaces on the east face
c              of "T" cells.
c     dciso2 = isopycnal tracer diffusivity coefficients modified based
c              on the slopes of the isopycnal surfaces on the north face
c              of "T" cells.
c     dslope = half length of the interval in which "ahisop" changes
c              with a steep slope from about 0.9*"ahisop" to about
c              0.1*"ahisop"
c     slopec = slope at which "ahisop" is equal to half of its original
c              value
c
c              REMARK: 0 <= "slopec", "dslope" <= 1/"slmxr".
c              REMARK: because the vertical gradient of density must be
c                      less than zero for isopycnal mixing, 1/"slmxr" is
c                      actually a negative maximum slope. this fact is
c                      taken into account in "isopyc.F". consequently,
c                      "slmxr" must be a positive number in "blkdta.F".
c# endif
c# ifdef gent_mcwilliams
c     adv_vetiso = zonal isopycnal mixing velocity computed at the 
c                  center of the eastern face of the "t" cells
c     adv_vntiso = meridional isopycnal mixing velocity computed at
c                  the center of the northern face of the "t" cells
c                  (Note: this includes the cosine as in "adv_vnt")
c     adv_vbtiso = vertical isopycnal mixing velocity computed at the
c                  center of the top face of the "t" cells
c     adv_fbiso  = "adv_vbtiso" * (tracer) evaluated at the center of
c                  the bottom face of the "t" cells
c# endif
c
cmove nrpl = number of reference presure levels used in isopycnal
c     integer nrpl
c     parameter (nrpl=5)

      integer xup,xmd,xlo
      parameter(xup=1,xmd=2,xlo=3)
c
      real    rhoi,e
      real    K1,K2,K3
      real    ahisop, athkdf, fzisop, slmxr
      real    adv_vetiso
      real    adv_vntiso
      real    adv_vbtiso
      real    dciso1, dciso2
      real    dslope, slopec
cxjin
      real    kref,rdz0
c
      common /cisop/ rhoi(imt,km,jmt,xup:xlo) 
      common /cisop/ e(imt,kmp1,jmt,3)
      common /cisop/ K1(imt,km,jmt,3:3), K2(imt,km,jmt,3:3) 
      common /cisop/ K3(imt,km,jmt,1:3)   
      common /cisop/ ahisop, athkdf, fzisop(km), slmxr 
cmove common /cisopi/ kisrpl(km), krplin(nrpl),dptlim(nrpl+1)
      common /cisop/ adv_vetiso(imt,km,jmt)
      common /cisop/ adv_vntiso(imt,km,jmt)
      common /cisop/ adv_vbtiso(imt,0:km,jmt)   
      common /cisop/ dciso1(imt,km,jmt), dciso2(imt,km,jmt)
      common /cisop/ dslope, slopec
      common /cisop/ kref(km),rdz0(km)


c====================== include file "param.h" =========================
c
c     main parameter file which sets ocean characteristics:
c
c                                                                       
c     imt = number of grid points in longitudinal direction          
c     jmt = number of grid points in latitudinal  direction           
c     km  = number of the sigma levels
c     nt  = number of tracers (temperature, salinity, ...)
c                                                                       
      integer imt,jmt,km,nt,imm,jmm,kmp1,kmm1,i,j,k,n,m
c
      parameter (imt=32,jmt=32,km=30)
      parameter (nt=2)
      parameter (imm=imt-1,jmm=jmt-1,kmp1=km+1,kmm1=km-1)


c====================== include file "pconst.h" ========================
c
c
c     rules for parameter constants
c
c     use prefix of "c" for whole real numbers (eg: c57 for 57.0)
c     use "m" after prefix to designate negative values (minus sign)
c       (eg: cm7 for -7.0)
c     use prefix of "p" for non repeating fractions (eg: p5 for 0.5)
c     use prefix of "r" for reciprocals (eg: r3 for 1/3.0)
c     combine use of prefix above and "e" for scientific notation, with
c       (eg: c5e4 for 5.0e4, c1em10 for 1.0e-10)
c
      real c0,c1,c2,c3,c4,c5,c6,c25,c35
      parameter (c0=0.0,c1=1.0,c2=2.0,c3=3.0,c4=4.0)
      parameter (c5=5.0,c6=6.0,c25=25.0,c35=35.0)
c
      real p125,p25,p5
      parameter (p125=0.125,p25=0.25, p5=0.5)
c
      real c60,c1440,c3600
      parameter (c60=60.0, c1440=1440.0, c3600=3600.0)
c
      real c1e3,c1em4,c1em12
      parameter (c1e3=1.0e3,c1em4=1.0e-4,c1em12=1.0e-12)
c
      real r120,r150,r180,r365
      parameter (r120=c1/120.0,r150=c1/150.0)
      parameter (r180=c1/180.0,r365=c1/365.0)
c
      real secday
      parameter (secday=c1/(c60*c1440))
c
      real rho_0,rrho_0
      parameter (rho_0=1.029,rrho_0=c1/rho_0)
c
      real tbice
      parameter (tbice=-1.5)
c
c
      real pi,torad,radius,omega,grav
      parameter (pi     = 3.141592653589793)
      parameter (torad  = pi*r180    )
      parameter (radius = 6370.0d5   )
      parameter (omega  = pi/43082.0 )
      parameter (grav   = 980.0)
c
c     grav = earth's gravitational acceleration (cm/sec**2)
c


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
      real    am,ah,kappa_m,kappa_h,gamma_t,gamma_s,cdbot,gravr
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


