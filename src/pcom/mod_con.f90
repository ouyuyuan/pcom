
! Description: physical constants and tunable coefficients
!
!      Author: OU Yuyuan <ouyuyuan@lasg.iap.ac.cn>
!     Created: 2015-10-03 07:33:17 BJT
! Last Change: 2016-03-29 15:43:07 BJT

module mod_con

  use mod_kind, only: wp
  implicit none
  public

  ! physical 
  real (kind=wp), parameter :: &
    rho0 = 1029.0, & !    rho0 = 1035, & ! kg/m^3, reference density, \citep[P.47]{Gill1982book}
    rrho0 = 1.0/rho0, &
    g = 9.8,     & ! m/s^2, gravity accerleration, or 9.7963, P.46 of Griffies2004 
    a = 6370e3,  & ! 6371e3 m, mean radius of Earth \citep{Griffies2008}, \citep[P.597]{Gill1982book]
    pi = 3.141592653589793, &
    omega = pi / 43082.0, & ! \citep[P.43]{Griffies2004}
    torad = pi / 180.0

  ! tunable
  real (kind=wp), parameter :: &
    gamma_t = g * 40.0 / 3901.0, & ! kg/(m*s^3), see below
    gamma_s = 1.0 / (120 * 24*60*60), & ! 1/s, restore per 120 days (3 months)
    gamma_b = 0.2, & ! predict-corrector algorithm coef. in barotropic integration
    deltat = 1.0e-2, & ! Cel. Deg, temperature increament in calc. prho
    deltas = 1.0e-3, & ! psu, salinity increament in calc. prho
    cdbot = 2.6,& ! kg/m^3, bottom drag coefficient
    tice = - 1.5,  & ! Celcius degree, average frozon point of sea water
    km_c = 1.0e-5,    & ! m^2/s, background vertical momentum viscosity coefficient
    am_c = 1.0e3,  & ! m^2/s, initial horizontal momentum viscosity, 
    smag_c = 0.56, & ! , Smagorinsky scheme coef.
    ah_c = 1.0e3     ! m^2/s, horizontal tracer diffusity, 

    ! gamma_t = g * ddd/cp, newtonia restoring coefficient for heat flux
    !   where ddd = 40W/m**2/K, is the surface heat flux
    !         cp  = 3901 J/kg/K, is cpecific heat capacity of sea water 

end module mod_con!{{{1
!-------------------------------------------------------{{{1
! vim:fdm=marker:fdl=0:
! vim:foldtext=getline(v\:foldstart).'...'.(v\:foldend-v\:foldstart):
