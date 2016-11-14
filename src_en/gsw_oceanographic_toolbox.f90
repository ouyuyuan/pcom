!==========================================================================
! Gibbs SeaWater (GSW) Oceanographic Toolbox of TEOS-10 version 3.0 (Fortran)
!==========================================================================
!
! This is a subset of functions contained in the Gibbs SeaWater (GSW) 
! Oceanographic Toolbox of TEOS-10 (version 3.0).
! 
! salinity and temperature conversions
! gsw_sa_from_sp          - Absolute Salinity from Practical Salinity
! gsw_sstar_from_sp       - Preformed Salinity from Practical Salinity
! gsw_ct_from_t           - Conservative Temperature from in-situ temperature
!
! gsw_deltasa_from_sp     - Absolute Salinity Anomaly from Practical Salinity
! gsw_sr_from_sp          - Reference Salinity from Practical Salinity
! gsw_sp_from_sr          - Practical Salinity from Reference Salinity
! gsw_sp_from_sa          - Practical Salinity from Absolute Salinity
! gsw_sstar_from_sa       - Preformed Salinity from Absolute Salinity
! gsw_sp_from_sstar       - Practical Salinity from Preformed Salinity
! gsw_sa_from_sstar       - Absolute Salinity from Preformed Salinity
! gsw_ct_from_pt          - Conservative Temperature from potential temperature
! gsw_pt_from_ct          - potential temperature from Conservative Temperature
! gsw_pt0_from_t          - potential temperature with reference pressure of 0 dbar
! gsw_pt_from_t           - potential temperature 
!
! density and enthalpy, based on the 48-term expression for density
! gsw_rho                 - in-situ density from CT, and potential density
! gsw_alpha               - thermal expansion coefficient with respect to CT
! gsw_beta                - saline contraction coefficient at constant CT
! gsw_specvol             - specific volume
! gsw_specvol_anom        - specific volume anomaly
! gsw_sound_speed         - sound speed
! gsw_internal_energy     - internal energy
! gsw_enthalpy            - enthalpy
! gsw_dynamic_enthalpy    - dynamic enthalpy
!
! freezing temperatures
! gsw_ct_freezing         - Conservative Temperature freezing temperature of seawater
! gsw_t_freezing          - in-situ temperature freezing temperature of seawater
!
! isobaric melting enthalpy and isobaric evaporation enthalpy
! gsw_latentheat_melting  - latent heat of melting
! gsw_latentheat_evap_ct  - latent heat of evaporation
! gsw_latentheat_evap_t   - latent heat of evaporation
!
! basic thermodynamic properties in terms of in-situ t, based on the exact Gibbs function
! gsw_rho_t_exact         - in-situ density
! gsw_pot_rho_t_exact     - potential density
! gsw_alpha_wrt_t_exact   - thermal expansion coefficient with respect to in-situ temperature
! gsw_beta_const_t_exact  - saline contraction coefficient at constant in-situ temperature
! gsw_specvol_t_exact     - specific volume
! gsw_sound_speed_t_exact - sound speed
! gsw_kappa_t_exact       - isentropic compressibility
! gsw_enthalpy_t_exact    - enthalpy
! gsw_entropy_t_exact     - entropy
! gsw_cp_t_exact          - isobaric heat capacity
!
! Library functions of the GSW toolbox
! gsw_gibbs               - the TEOS-10 Gibbs function and its derivatives
! gsw_saar                - Absolute Salinity Anomaly Ratio (excluding the Baltic Sea)
! gsw_delta_sa_ref        - Absolute Salinity Anomaly ref. value (excluding the Baltic Sea)
! gsw_fdelta              - ratio of Absolute to Preformed Salinity, minus 1
! gsw_sa_from_sp_baltic   - Absolute Salinity Anomaly from Practical Salinity in the Baltic Sea
! gsw_sp_from_sa_baltic   - Practical Salinity from Absolute Salinity in the Baltic Sea
! gsw_entropy_part        - entropy minus the terms that are a function of only SA
! gsw_entropy_part_zerop  - entropy_part evaluated at 0 dbar
! gsw_gibbs_pt0_pt0       - gibbs(0,2,0,SA,t,0)
!
!
! Version 1.0 written by David Jackett
! Modified by Paul Barker (version 3.0)
!
! For help with this Oceanographic Toolbox email: help@teos-10.org
!
! This software is available from http://www.teos-10.org
!
!==========================================================================

!--------------------------------------------------------------------------
! salinity and temperature conversions
!--------------------------------------------------------------------------

!==========================================================================
function gsw_sa_from_sp(sp,p,long,lat)       
!==========================================================================

! Calculates Absolute Salinity, SA, from Practical Salinity, SP
!
! sp     : Practical Salinity                              [unitless]
! p      : sea pressure                                    [dbar]
! long   : longitude                                       [DEG E]     
! lat    : latitude                                        [DEG N]
!
! gsw_sa_from_sp   : Absolute Salinity                     [g/kg]

implicit none
integer, parameter :: r14 = selected_real_kind(14,30)

real (r14) :: sp, long, lat, p, gsw_sa_from_sp, gsw_saar, saar
real (r14) :: gsw_sa_baltic, gsw_sa_from_sp_baltic

saar = gsw_saar(p,long,lat)

gsw_sa_from_sp = (35.16504d0/35.d0)*sp*(1.d0 + saar)

gsw_sa_baltic = gsw_sa_from_sp_baltic(sp,long,lat)

if (gsw_sa_baltic.lt.1d10) then
   gsw_sa_from_sp = gsw_sa_baltic
end if

if (saar.eq.9d15) then
   gsw_sa_from_sp = 9d15
end if

return
end function

!--------------------------------------------------------------------------

!==========================================================================
function gsw_sstar_from_sp(sp,p,long,lat) 
!==========================================================================

! Calculates Preformed Salinity, Sstar, from Practical Salinity, SP. 
!
! sp     : Practical Salinity                              [unitless]
! p      : sea pressure                                    [dbar]
! long   : longitude                                       [deg E]     
! lat    : latitude                                        [deg N]
!
! gsw_sstar_from_sp  : Preformed Salinity                  [g/kg]

implicit none

integer, parameter :: r14 = selected_real_kind(14,30)

real (r14) :: sp, long, lat, p, gsw_saar, gsw_sa_from_sp_baltic
real (r14) :: saar, gsw_sstar_from_sp, sstar_baltic

saar = gsw_saar(p,long,lat)

gsw_sstar_from_sp = (35.16504d0/35.d0)*sp*(1 - 0.35d0*saar);

!In the Baltic Sea, Sstar = SA.
sstar_baltic = gsw_sa_from_sp_baltic(sp,long,lat);

if (sstar_baltic.lt.1d10) then
    gsw_sstar_from_sp = sstar_baltic;
end if

if (saar.eq.9d15) then
    gsw_sstar_from_sp = 9d15
end if

return
end function

!--------------------------------------------------------------------------

!==========================================================================
function gsw_ct_from_t(sa,t,p) 
!==========================================================================
   
! Calculates Conservative Temperature from in-situ temperature
!
! sa     : Absolute Salinity                               [g/kg]
! t      : in-situ temperature                             [deg C]
! p      : sea pressure                                    [dbar]
!
! gsw_ct_from_t : Conservative Temperature                 [deg C]

implicit none

integer, parameter :: r14 = selected_real_kind(14,30)

real (r14) :: sa, t, p, pt0, gsw_pt0_from_t, gsw_ct_from_pt,  gsw_ct_from_t

pt0 = gsw_pt0_from_t(sa,t,p)
gsw_ct_from_t = gsw_ct_from_pt(sa,pt0)

return
end function

!--------------------------------------------------------------------------

!==========================================================================
function gsw_deltasa_from_sp(sp,p,long,lat) 
!==========================================================================

! Calculates Absolute Salinity Anomaly, deltaSA, from Practical Salinity, SP. 
!
! sp     : Practical Salinity                              [unitless]
! p      : sea pressure                                    [dbar]
! long   : longitude                                       [deg E]     
! lat    : latitude                                        [deg N]
!
! gsw_deltasa_from_sp : Absolute Salinty Anomaly           [g/kg]

implicit none

integer, parameter :: r14 = selected_real_kind(14,30)

real (r14) :: sp, long, lat, p, gsw_sa_from_sp, gsw_sr_from_sp
real (r14) :: gsw_deltasa_from_sp

gsw_deltasa_from_sp = gsw_sa_from_sp(sp,p,long,lat) - gsw_sr_from_sp(sp)

if (gsw_deltasa_from_sp.gt.1d10) then
    gsw_deltasa_from_sp = 9d15
end if

return
end function

!--------------------------------------------------------------------------

!==========================================================================
function gsw_sr_from_sp(sp) 
!==========================================================================

! Calculates Reference Salinity, SR, from Practical Salinity, SP. 
!
! sp     : Practical Salinity                              [unitless]
!
! gsw_sr_from_sp : Reference Salinity                      [g/kg]

implicit none

integer, parameter :: r14 = selected_real_kind(14,30)

real (r14) :: sp, gsw_sr_from_sp

gsw_sr_from_sp = 1.004715428571429d0*sp;

if (gsw_sr_from_sp.ge.1.d10) then
    gsw_sr_from_sp = 9.d15
end if

return
end function

!--------------------------------------------------------------------------

!==========================================================================
function gsw_sp_from_sr(sr)  
!==========================================================================

! Calculates Practical Salinity, sp, from Reference Salinity, sr. 
!
! sr     : Reference Salinity                              [g/kg]
!
! gsw_sp_from_sr  : Practical Salinity                     [unitless]

implicit none

integer, parameter :: r14 = selected_real_kind(14,30)

real (r14) :: sr, gsw_sp_from_sr

gsw_sp_from_sr = 0.995306702338459d0*sr;

if (gsw_sp_from_sr.gt.1d10) then
    gsw_sp_from_sr = 9d15
end if

return
end function

!--------------------------------------------------------------------------

!==========================================================================
function gsw_sp_from_sa(sa,p,long,lat) 
!==========================================================================

! Calculates Practical salinity, sp, from Absolute salinity, sa  
!
! sa     : Absolute Salinity                               [g/kg]
! p      : sea pressure                                    [dbar]
! long   : longitude                                       [DEG E]     
! lat    : latitude                                        [DEG N]
!
! gsw_sp_from_sa      : Practical Salinity                 [unitless]

implicit none

integer, parameter :: r14 = selected_real_kind(14,30)

real (r14) :: sa, long, lat, p, gsw_sp_from_sa, gsw_saar, saar
real (r14) :: gsw_sp_baltic, gsw_sp_from_sa_baltic

saar = gsw_saar(p,long,lat)

gsw_sp_from_sa = (35.d0/35.16504d0)*sa/(1d0 + saar)

gsw_sp_baltic = gsw_sp_from_sa_baltic(sa,long,lat);

if (gsw_sp_baltic.lt.1d10) then
   gsw_sp_from_sa = gsw_sp_baltic
end if

if (saar.eq.9d15) then
   gsw_sp_from_sa = 9d15
end if

return
end function

!--------------------------------------------------------------------------

!==========================================================================
function gsw_sstar_from_sa(sa,p,long,lat) 
!==========================================================================

! Calculates Preformed Salinity, Sstar, from Absolute Salinity, SA. 
!
! sa     : Absolute Salinity                               [g/kg]
! p      : sea pressure                                    [dbar]
! long   : longitude                                       [deg E]     
! lat    : latitude                                        [deg N]
!
! gsw_sstar_from_sa : Preformed Salinity                   [g/kg]

implicit none

integer, parameter :: r14 = selected_real_kind(14,30)

real (r14) :: sa, long, lat, p, gsw_saar, gsw_sp_from_sa_baltic
real (r14) :: saar, gsw_sstar_from_sa

saar = gsw_saar(p,long,lat)

gsw_sstar_from_sa = sa*(1d0 - 0.35d0*saar)/(1d0 + saar)

! In the Baltic Sea, Sstar = sa, and note that gsw_saar returns zero
! for saar in the Baltic.

if (saar.eq.9d15) then
    gsw_sstar_from_sa = 9d15
end if

return
end function

!--------------------------------------------------------------------------

!==========================================================================
function gsw_sa_from_sstar(sstar,p,long,lat)  
!==========================================================================

! Calculates Absolute Salinity, SA, from Preformed Salinity, Sstar.
!
! Sstar  : Preformed Salinity                              [g/kg]
! p      : sea pressure                                    [dbar]
! long   : longitude                                       [deg E]     
! lat    : latitude                                        [deg N]
!
! gsw_sa_from_sstar   : Absolute Salinity                  [g/kg]

implicit none

integer, parameter :: r14 = selected_real_kind(14,30)

real (r14) :: sa, long, lat, p, gsw_saar, gsw_sp_from_sa_baltic
real (r14) :: saar, gsw_sa_from_sstar, sstar

saar = gsw_saar(p,long,lat)

gsw_sa_from_sstar = sstar*(1d0 + saar)/(1d0 - 0.35d0*saar)

! In the Baltic Sea, Sstar = SA, and note that gsw_saar returns zero
! for SAAR in the Baltic.

if (saar.eq.9d15) then
    gsw_sa_from_sstar = 9d15
end if

return
end function

!--------------------------------------------------------------------------

!==========================================================================
function gsw_sp_from_sstar(sstar,p,long,lat)  
!==========================================================================

! Calculates Practical Salinity, SP, from Preformed Salinity, Sstar. 
!
! sstar  : Preformed Salinity                              [g/kg]
! p      : sea pressure                                    [dbar]
! long   : longitude                                       [deg E]     
! lat    : latitude                                        [deg N]
!
! gsw_sp_from_Sstar : Preformed Salinity                   [g/kg]

implicit none

integer, parameter :: r14 = selected_real_kind(14,30)

real (r14) :: long, lat, p, gsw_saar, gsw_sp_from_sa_baltic
real (r14) :: saar, gsw_sp_from_sstar, sp_baltic, Sstar

saar = gsw_saar(p,long,lat)

gsw_sp_from_sstar = (35.d0/35.16504d0)*Sstar/(1 - 0.35d0*saar);

!In the Baltic Sea, SA = Sstar.
sp_baltic = gsw_sp_from_sa_baltic(sstar,long,lat);

if (sp_baltic.lt.1d10) then
    gsw_sp_from_sstar = sp_baltic;
end if

if (saar.eq.9d15) then
    gsw_sp_from_sstar = 9d15
end if

return
end function

!--------------------------------------------------------------------------

!==========================================================================
function gsw_t_from_ct(sa,ct,p)  
!==========================================================================

! Calculates in-situ temperature from Conservative Temperature of seawater  
!
! sa      : Absolute Salinity                              [g/kg]
! ct      : Conservative Temperature                       [deg C]
!
! gsw_t_from_ct : in-situ temperature                      [deg C]

implicit none

integer, parameter :: r14 = selected_real_kind(14,30)

real (r14) :: sa, ct, p, pt0, p0, gsw_t_from_ct, gsw_pt_from_ct, gsw_pt_from_t

p0 = 0d0
pt0 = gsw_pt_from_ct(sa,ct)
gsw_t_from_ct = gsw_pt_from_t(sa,pt0,p0,p)

return
end function

!--------------------------------------------------------------------------

!==========================================================================
function gsw_ct_from_pt(sa,pt) 
!==========================================================================

! Calculates Conservative Temperature from potential temperature of seawater  
!
! sa      : Absolute Salinity                              [g/kg]
! pt      : potential temperature with                     [deg C]
!           reference pressure of 0 dbar
!
! gsw_ct_from_pt : Conservative Temperature                [deg C]

implicit none

integer, parameter :: r14 = selected_real_kind(14,30)

real (r14) :: sa, pt, p, pot_enthalpy, gsw_ct_from_pt, sfac 
real (r14) :: x2, x, y, cp0

sfac = 0.0248826675584615d0 

x2 = sfac*sa
x = sqrt(x2)
y = pt*0.025d0        ! normalize for F03 and F08

pot_enthalpy =  61.01362420681071d0 + y*(168776.46138048015d0 + &
               y*(-2735.2785605119625d0 + y*(2574.2164453821433d0 + &
               y*(-1536.6644434977543d0 + y*(545.7340497931629d0 + &
               (-50.91091728474331d0 - 18.30489878927802d0*y)*y))))) + &
               x2*(268.5520265845071d0 + y*(-12019.028203559312d0 + &
               y*(3734.858026725145d0 + y*(-2046.7671145057618d0 + &
               y*(465.28655623826234d0 + (-0.6370820302376359d0 - &
               10.650848542359153d0*y)*y)))) + &
               x*(937.2099110620707d0 + y*(588.1802812170108d0 + &
               y*(248.39476522971285d0 + (-3.871557904936333d0 - &
               2.6268019854268356d0*y)*y)) + &
               x*(-1687.914374187449d0 + x*(246.9598888781377d0 + &
               x*(123.59576582457964d0 - 48.5891069025409d0*x)) + &
               y*(936.3206544460336d0 + &
               y*(-942.7827304544439d0 + y*(369.4389437509002d0 + &
               (-33.83664947895248d0 - 9.987880382780322d0*y)*y))))))

cp0 = 3991.86795711963d0

gsw_ct_from_pt = pot_enthalpy/cp0

end

!--------------------------------------------------------------------------

!==========================================================================
function gsw_pt_from_t(sa,t,p,p_ref) 
!==========================================================================
   
! Calculates potential temperature of seawater from in-situ temperature 
!
! sa     : Absolute Salinity                               [g/kg]
! t      : in-situ temperature                             [deg C]
! p      : sea pressure                                    [dbar]
! p_ref  : reference sea pressure                          [dbar]
!
! gsw_pt_from_t : potential temperature                    [deg C]

implicit none

integer, parameter :: r14 = selected_real_kind(14,30)

integer n0, n2, n, no_iter
real (r14) :: sa, t, p, p_ref, s1, gsw_entropy_t_exact, gsw_gibbs, gsw_pt_from_t
real (r14) :: pt, pt_old, de_dt, dentropy, dentropy_dt, sso, cp0
real (r14) :: gsw_entropy_part, true_entropy_part, ptm

n0 = 0
n2 = 2

cp0 = 3991.86795711963d0
sso = 35.16504d0

s1 = sa*35d0/sso

pt = t + (p-p_ref)*( 8.65483913395442d-6  - &
               s1 *  1.41636299744881d-6  - &
        (p+p_ref) *  7.38286467135737d-9  + &
               t  *(-8.38241357039698d-6  + &
               s1 *  2.83933368585534d-8  + &
               t  *  1.77803965218656d-8  + &
        (p+p_ref) *  1.71155619208233d-10))

dentropy_dt = cp0/((273.15d0 + pt)*(1d0-0.05d0*(1d0 - sa/sso)));

true_entropy_part = gsw_entropy_part(sa,t,p)

do no_iter = 1,2
    pt_old = pt
    dentropy = gsw_entropy_part(sa,pt_old,p_ref) - true_entropy_part
    pt = pt_old - dentropy/dentropy_dt 
    ptm = 0.5d0*(pt + pt_old)
    dentropy_dt = -gsw_gibbs(n0,n2,n0,sa,ptm,p_ref)
    pt = pt_old - dentropy/dentropy_dt
end do

gsw_pt_from_t = pt

end

!--------------------------------------------------------------------------

!==========================================================================
function gsw_pt0_from_t(sa,t,p) 
!==========================================================================
   
! Calculates potential temperature with reference pressure, p_ref = 0 dbar. 
!
! sa     : Absolute Salinity                               [g/kg]
! t      : in-situ temperature                             [deg C]
! p      : sea pressure                                    [dbar]
!
! gsw_pt0_from_t : potential temperature, p_ref = 0        [deg C]

implicit none

integer, parameter :: r14 = selected_real_kind(14,30)

integer n0, n2, n, no_iter
real (r14) :: sa, t, p, s1, gsw_entropy_t_exact, gsw_gibbs_pt0_pt0, gsw_pt0_from_t
real (r14) :: pt0, pt0_old, de_dt, dentropy, dentropy_dt, sso, cp0
real (r14) :: gsw_entropy_part_zerop, true_entropy_part, pt0m, gsw_entropy_part

n0 = 0
n2 = 2

cp0 = 3991.86795711963d0
sso = 35.16504d0

s1 = sa*35d0/sso

pt0 = t + p*( 8.65483913395442d-6  - &
        s1 *  1.41636299744881d-6  - &
         p *  7.38286467135737d-9  + &
         t *(-8.38241357039698d-6  + &
        s1 *  2.83933368585534d-8  + &
         t *  1.77803965218656d-8  + &
         p *  1.71155619208233d-10))

dentropy_dt = cp0/((273.15d0 + pt0)*(1d0-0.05d0*(1d0 - sa/sso)))

true_entropy_part = gsw_entropy_part(sa,t,p)

do no_iter = 1,2
    pt0_old = pt0
    dentropy = gsw_entropy_part_zerop(sa,pt0_old) - true_entropy_part
    pt0 = pt0_old - dentropy/dentropy_dt 
    pt0m = 0.5d0*(pt0 + pt0_old)
    dentropy_dt = -gsw_gibbs_pt0_pt0(SA,pt0m)
    pt0 = pt0_old - dentropy/dentropy_dt
end do

gsw_pt0_from_t = pt0

end

!--------------------------------------------------------------------------

!==========================================================================
function gsw_pt_from_ct(sa,ct) 
!==========================================================================

! potential temperature of seawater from conservative temperature
!
! sa     : Absolute Salinity                               [g/kg]
! ct     : Conservative Temperature                        [deg C]
! p      : sea pressure                                    [dbar]
!
! gsw_pt_from_ct : potential temperature with              [deg C]
!                  reference pressure of  0 dbar

implicit none

integer, parameter :: r14 = selected_real_kind(14,30)

integer n0, n2, nloops, n
real (r14) :: sa, ct, s1, p0, gsw_pt_from_ct, gsw_ct_from_pt, gsw_gibbs, cp0 
real (r14) :: a0, a1, a2, a3, a4, a5, b0, b1, b2, b3
real (r14) :: a5ct, b3ct, ct_factor, pt_num, pt_den, ct_diff
real (r14) :: ct0, pt, pt_old, ptm, dct, dct_dpt, gsw_gibbs_pt0_pt0

cp0 = 3991.86795711963d0    

n0 = 0
n2 = 2

s1 = sa*35.d0/35.16504d0
p0 = 0.d0

a0 = -1.446013646344788d-2;    
a1 = -3.305308995852924d-3;    
a2 =  1.062415929128982d-4;     
a3 =  9.477566673794488d-1;     
a4 =  2.166591947736613d-3
a5 =  3.828842955039902d-3

b0 =  1.000000000000000d0
b1 =  6.506097115635800d-4
b2 =  3.830289486850898d-3
b3 =  1.247811760368034d-6

a5ct = a5*ct
b3ct = b3*ct

ct_factor = (a3 + a4*s1 + a5ct)
pt_num = a0 + s1*(a1 + a2*s1) + ct*ct_factor
pt_den = b0 + b1*s1 + ct*(b2 + b3ct)
pt = (pt_num)/(pt_den)

dct_dpt = (pt_den)/(ct_factor + a5ct - (b2 + b3ct + b3ct)*pt);

! Start the 1.5 iterations through the modified Newton-Rapshon iterative,
! method, which is also know as the Newton-McDougall method. 

ct_diff = gsw_ct_from_pt(sa,pt) - ct
pt_old = pt
pt = pt_old - (ct_diff)/dct_dpt
ptm = 0.5d0*(pt + pt_old)

dct_dpt = -(ptm + 273.15d0)*gsw_gibbs_pt0_pt0(sa,ptm)/cp0

pt = pt_old - (ct_diff)/dct_dpt
ct_diff = gsw_ct_from_pt(sa,pt) - ct
pt_old = pt
gsw_pt_from_ct = pt_old - (ct_diff)/dct_dpt

return 
end function

!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
! density and enthalpy, based on the 48-term expression for density
!--------------------------------------------------------------------------

!==========================================================================
function gsw_rho(sa,ct,p) 
!==========================================================================

!  Calculates in-situ density from Absolute Salinity and Conservative 
!  Temperature, using the computationally-efficient 48-term expression for
!  density in terms of SA, CT and p (McDougall et al., 2011).
!
! sa     : Absolute Salinity                               [g/kg]
! ct     : Conservative Temperature                        [deg C]
! p      : sea pressure                                    [dbar]
! 
! gsw_rho  : in-situ density (48 term equation)

implicit none

integer, parameter :: r14 = selected_real_kind(14,30)

real (r14), parameter :: v01 =  9.998420897506056d2, v02 =  2.839940833161907d0
real (r14), parameter :: v03 = -3.147759265588511d-2, v04 =  1.181805545074306d-3
real (r14), parameter :: v05 = -6.698001071123802d0, v06 = -2.986498947203215d-2
real (r14), parameter :: v07 =  2.327859407479162d-4, v08 = -3.988822378968490d-2
real (r14), parameter :: v09 =  5.095422573880500d-4, v10 = -1.426984671633621d-5
real (r14), parameter :: v11 =  1.645039373682922d-7, v12 = -2.233269627352527d-2
real (r14), parameter :: v13 = -3.436090079851880d-4, v14 =  3.726050720345733d-6
real (r14), parameter :: v15 = -1.806789763745328d-4, v16 =  6.876837219536232d-7
real (r14), parameter :: v17 = -3.087032500374211d-7, v18 = -1.988366587925593d-8
real (r14), parameter :: v19 = -1.061519070296458d-11, v20 =  1.550932729220080d-10
real (r14), parameter :: v21 =  1.0d0, v22 =  2.775927747785646d-3, v23 = -2.349607444135925d-5
real (r14), parameter :: v24 =  1.119513357486743d-6, v25 =  6.743689325042773d-10
real (r14), parameter :: v26 = -7.521448093615448d-3, v27 = -2.764306979894411d-5
real (r14), parameter :: v28 =  1.262937315098546d-7, v29 =  9.527875081696435d-10
real (r14), parameter :: v30 = -1.811147201949891d-11, v31 = -3.303308871386421d-5
real (r14), parameter :: v32 =  3.801564588876298d-7, v33 = -7.672876869259043d-9
real (r14), parameter :: v34 = -4.634182341116144d-11, v35 =  2.681097235569143d-12
real (r14), parameter :: v36 =  5.419326551148740d-6, v37 = -2.742185394906099d-5
real (r14), parameter :: v38 = -3.212746477974189d-7, v39 =  3.191413910561627d-9
real (r14), parameter :: v40 = -1.931012931541776d-12, v41 = -1.105097577149576d-7
real (r14), parameter :: v42 =  6.211426728363857d-10, v43 = -1.119011592875110d-10
real (r14), parameter :: v44 = -1.941660213148725d-11, v45 = -1.864826425365600d-14
real (r14), parameter :: v46 =  1.119522344879478d-14, v47 = -1.200507748551599d-15
real (r14), parameter :: v48 =  6.057902487546866d-17 

real (r14) :: sa, ct, p, sqrtsa, v_hat_denominator, v_hat_numerator, gsw_rho

sqrtsa = sqrt(sa)

v_hat_denominator = v01 + ct*(v02 + ct*(v03 + v04*ct))  &
             + sa*(v05 + ct*(v06 + v07*ct) &
         + sqrtsa*(v08 + ct*(v09 + ct*(v10 + v11*ct)))) &
              + p*(v12 + ct*(v13 + v14*ct) + sa*(v15 + v16*ct) &
              + p*(v17 + ct*(v18 + v19*ct) + v20*sa))

v_hat_numerator = v21 + ct*(v22 + ct*(v23 + ct*(v24 + v25*ct))) &
           + sa*(v26 + ct*(v27 + ct*(v28 + ct*(v29 + v30*ct))) + v36*sa &
       + sqrtsa*(v31 + ct*(v32 + ct*(v33 + ct*(v34 + v35*ct)))))  &
            + p*(v37 + ct*(v38 + ct*(v39 + v40*ct))  &
           + sa*(v41 + v42*ct) &
            + p*(v43 + ct*(v44 + v45*ct + v46*sa) &
            + p*(v47 + v48*ct)))

gsw_rho = v_hat_denominator/v_hat_numerator

return
end function

!--------------------------------------------------------------------------

!==========================================================================
function gsw_alpha(sa,ct,p)  
!==========================================================================

!  Calculates the thermal expansion coefficient of seawater with respect to 
!  Conservative Temperature using the computationally-efficient 48-term 
!  expression for density in terms of SA, CT and p (McDougall et al., 2011)
!
! sa     : Absolute Salinity                               [g/kg]
! ct     : Conservative Temperature                        [deg C]
! p      : sea pressure                                    [dbar]
! 
! gsw_alpha : thermal expansion coefficient of seawater (48 term equation)

implicit none

integer, parameter :: r14 = selected_real_kind(14,30)

real (r14), parameter :: v01 =  9.998420897506056d2, v02 =  2.839940833161907d0
real (r14), parameter :: v03 = -3.147759265588511d-2, v04 =  1.181805545074306d-3
real (r14), parameter :: v05 = -6.698001071123802d0, v06 = -2.986498947203215d-2
real (r14), parameter :: v07 =  2.327859407479162d-4, v08 = -3.988822378968490d-2
real (r14), parameter :: v09 =  5.095422573880500d-4, v10 = -1.426984671633621d-5
real (r14), parameter :: v11 =  1.645039373682922d-7, v12 = -2.233269627352527d-2
real (r14), parameter :: v13 = -3.436090079851880d-4, v14 =  3.726050720345733d-6
real (r14), parameter :: v15 = -1.806789763745328d-4, v16 =  6.876837219536232d-7
real (r14), parameter :: v17 = -3.087032500374211d-7, v18 = -1.988366587925593d-8
real (r14), parameter :: v19 = -1.061519070296458d-11, v20 =  1.550932729220080d-10
real (r14), parameter :: v21 =  1.0d0, v22 = 2.775927747785646d-3, v23 = -2.349607444135925d-5
real (r14), parameter :: v24 =  1.119513357486743d-6, v25 =  6.743689325042773d-10
real (r14), parameter :: v26 = -7.521448093615448d-3, v27 = -2.764306979894411d-5
real (r14), parameter :: v28 =  1.262937315098546d-7, v29 =  9.527875081696435d-10
real (r14), parameter :: v30 = -1.811147201949891d-11, v31 = -3.303308871386421d-5
real (r14), parameter :: v32 =  3.801564588876298d-7, v33 = -7.672876869259043d-9
real (r14), parameter :: v34 = -4.634182341116144d-11, v35 =  2.681097235569143d-12
real (r14), parameter :: v36 =  5.419326551148740d-6, v37 = -2.742185394906099d-5
real (r14), parameter :: v38 = -3.212746477974189d-7, v39 =  3.191413910561627d-9
real (r14), parameter :: v40 = -1.931012931541776d-12, v41 = -1.105097577149576d-7
real (r14), parameter :: v42 =  6.211426728363857d-10, v43 = -1.119011592875110d-10
real (r14), parameter :: v44 = -1.941660213148725d-11, v45 = -1.864826425365600d-14
real (r14), parameter :: v46 =  1.119522344879478d-14, v47 = -1.200507748551599d-15
real (r14), parameter :: v48 =  6.057902487546866d-17 
real (r14), parameter :: a01 =  2.839940833161907d0, a02 = -6.295518531177023d-2
real (r14), parameter :: a03 =  3.545416635222918d-3, a04 = -2.986498947203215d-2
real (r14), parameter :: a05 =  4.655718814958324d-4, a06 =  5.095422573880500d-4
real (r14), parameter :: a07 = -2.853969343267241d-5, a08 =  4.935118121048767d-7
real (r14), parameter :: a09 = -3.436090079851880d-4, a10 =  7.452101440691467d-6
real (r14), parameter :: a11 =  6.876837219536232d-7, a12 = -1.988366587925593d-8
real (r14), parameter :: a13 = -2.123038140592916d-11, a14 =  2.775927747785646d-3
real (r14), parameter :: a15 = -4.699214888271850d-5, a16 =  3.358540072460230d-6
real (r14), parameter :: a17 =  2.697475730017109d-9, a18 = -2.764306979894411d-5
real (r14), parameter :: a19 =  2.525874630197091d-7, a20 =  2.858362524508931d-9
real (r14), parameter :: a21 = -7.244588807799565d-11, a22 =  3.801564588876298d-7
real (r14), parameter :: a23 = -1.534575373851809d-8, a24 = -1.390254702334843d-10
real (r14), parameter :: a25 =  1.072438894227657d-11, a26 = -3.212746477974189d-7
real (r14), parameter :: a27 =  6.382827821123254d-9, a28 = -5.793038794625329d-12
real (r14), parameter :: a29 =  6.211426728363857d-10, a30 = -1.941660213148725d-11
real (r14), parameter :: a31 = -3.729652850731201d-14, a32 =  1.119522344879478d-14
real (r14), parameter :: a33 =  6.057902487546866d-17

real (r14) :: sa, ct, p, sqrtsa, v_hat_denominator, v_hat_numerator, spec_vol
real (r14) :: dvhatden_dct, dvhatnum_dct, gsw_alpha

sqrtsa = sqrt(sa)

v_hat_denominator = v01 + ct*(v02 + ct*(v03 + v04*ct))  &
             + sa*(v05 + ct*(v06 + v07*ct) &
         + sqrtsa*(v08 + ct*(v09 + ct*(v10 + v11*ct)))) &
              + p*(v12 + ct*(v13 + v14*ct) + sa*(v15 + v16*ct) &
              + p*(v17 + ct*(v18 + v19*ct) + v20*sa))

v_hat_numerator = v21 + ct*(v22 + ct*(v23 + ct*(v24 + v25*ct))) &
           + sa*(v26 + ct*(v27 + ct*(v28 + ct*(v29 + v30*ct))) + v36*sa &
       + sqrtsa*(v31 + ct*(v32 + ct*(v33 + ct*(v34 + v35*ct)))))  &
            + p*(v37 + ct*(v38 + ct*(v39 + v40*ct))  &
           + sa*(v41 + v42*ct) &
            + p*(v43 + ct*(v44 + v45*ct + v46*sa) &
            + p*(v47 + v48*ct)))
       
spec_vol = v_hat_numerator/v_hat_denominator

dvhatden_dct = a01 + ct*(a02 + a03*ct) &
        + sa*(a04 + a05*ct &
    + sqrtsa*(a06 + ct*(a07 + a08*ct))) &
         + p*(a09 + a10*ct + a11*sa &
         + p*(a12 + a13*ct))

dvhatnum_dct = a14 + ct*(a15 + ct*(a16 + a17*ct)) &
        + sa*(a18 + ct*(a19 + ct*(a20 + a21*ct)) &
    + sqrtsa*(a22 + ct*(a23 + ct*(a24 + a25*ct)))) &
         + p*(a26 + ct*(a27 + a28*ct) + a29*sa &
         + p*(a30 + a31*ct + a32*sa + a33*p))
 
gsw_alpha = (dvhatnum_dct - dvhatden_dct*spec_vol)/v_hat_numerator

return
end function

!--------------------------------------------------------------------------

!==========================================================================
function gsw_beta(sa,ct,p) 
!==========================================================================

!  Calculates the saline (i.e. haline) contraction coefficient of seawater  
!  at constant Conservative Temperature using the computationally-efficient
!  48-term expression for density in terms of SA, CT and p 
!  (McDougall et al., 2011).
!
! sa     : Absolute Salinity                               [g/kg]
! ct     : Conservative Temperature                        [deg C]
! p      : sea pressure                                    [dbar]
! 
! gsw_beta : saline contraction coefficient of seawater (48 term equation)

implicit none

integer, parameter :: r14 = selected_real_kind(14,30)

real (r14), parameter :: v01 =  9.998420897506056d2, v02 =  2.839940833161907d0
real (r14), parameter :: v03 = -3.147759265588511d-2, v04 =  1.181805545074306d-3
real (r14), parameter :: v05 = -6.698001071123802d0, v06 = -2.986498947203215d-2
real (r14), parameter :: v07 =  2.327859407479162d-4, v08 = -3.988822378968490d-2
real (r14), parameter :: v09 =  5.095422573880500d-4, v10 = -1.426984671633621d-5
real (r14), parameter :: v11 =  1.645039373682922d-7, v12 = -2.233269627352527d-2
real (r14), parameter :: v13 = -3.436090079851880d-4, v14 =  3.726050720345733d-6
real (r14), parameter :: v15 = -1.806789763745328d-4, v16 =  6.876837219536232d-7
real (r14), parameter :: v17 = -3.087032500374211d-7, v18 = -1.988366587925593d-8
real (r14), parameter :: v19 = -1.061519070296458d-11, v20 =  1.550932729220080d-10
real (r14), parameter :: v21 =  1.0d0, v22 =  2.775927747785646d-3, v23 = -2.349607444135925d-5
real (r14), parameter :: v24 =  1.119513357486743d-6, v25 =  6.743689325042773d-10
real (r14), parameter :: v26 = -7.521448093615448d-3, v27 = -2.764306979894411d-5
real (r14), parameter :: v28 =  1.262937315098546d-7, v29 =  9.527875081696435d-10
real (r14), parameter :: v30 = -1.811147201949891d-11, v31 = -3.303308871386421d-5
real (r14), parameter :: v32 =  3.801564588876298d-7, v33 = -7.672876869259043d-9
real (r14), parameter :: v34 = -4.634182341116144d-11, v35 =  2.681097235569143d-12
real (r14), parameter :: v36 =  5.419326551148740d-6, v37 = -2.742185394906099d-5
real (r14), parameter :: v38 = -3.212746477974189d-7, v39 =  3.191413910561627d-9
real (r14), parameter :: v40 = -1.931012931541776d-12, v41 = -1.105097577149576d-7
real (r14), parameter :: v42 =  6.211426728363857d-10, v43 = -1.119011592875110d-10
real (r14), parameter :: v44 = -1.941660213148725d-11, v45 = -1.864826425365600d-14
real (r14), parameter :: v46 =  1.119522344879478d-14, v47 = -1.200507748551599d-15
real (r14), parameter :: v48 =  6.057902487546866d-17 
real (r14), parameter :: b01 = -6.698001071123802d0, b02 = -2.986498947203215d-2
real (r14), parameter :: b03 =  2.327859407479162d-4, b04 = -5.983233568452735d-2
real (r14), parameter :: b05 =  7.643133860820750d-4, b06 = -2.140477007450431d-5
real (r14), parameter :: b07 =  2.467559060524383d-7, b08 = -1.806789763745328d-4
real (r14), parameter :: b09 =  6.876837219536232d-7, b10 =  1.550932729220080d-10
real (r14), parameter :: b11 = -7.521448093615448d-3, b12 = -2.764306979894411d-5
real (r14), parameter :: b13 =  1.262937315098546d-7, b14 =  9.527875081696435d-10
real (r14), parameter :: b15 = -1.811147201949891d-11, b16 = -4.954963307079632d-5
real (r14), parameter :: b17 =  5.702346883314446d-7, b18 = -1.150931530388857d-8
real (r14), parameter :: b19 = -6.951273511674217d-11, b20 =  4.021645853353715d-12
real (r14), parameter :: b21 =  1.083865310229748d-5, b22 = -1.105097577149576d-7
real (r14), parameter :: b23 =  6.211426728363857d-10, b24 =  1.119522344879478d-14

real (r14) :: sa, ct, p, sqrtsa, v_hat_denominator, v_hat_numerator, gsw_beta
real (r14) :: spec_vol, dvhatden_dsa, dvhatnum_dsa 

sqrtsa = sqrt(sa)

v_hat_denominator = v01 + ct*(v02 + ct*(v03 + v04*ct))  &
             + sa*(v05 + ct*(v06 + v07*ct) &
         + sqrtsa*(v08 + ct*(v09 + ct*(v10 + v11*ct)))) &
              + p*(v12 + ct*(v13 + v14*ct) + sa*(v15 + v16*ct) &
              + p*(v17 + ct*(v18 + v19*ct) + v20*sa))

v_hat_numerator = v21 + ct*(v22 + ct*(v23 + ct*(v24 + v25*ct))) &
           + sa*(v26 + ct*(v27 + ct*(v28 + ct*(v29 + v30*ct))) + v36*sa &
       + sqrtsa*(v31 + ct*(v32 + ct*(v33 + ct*(v34 + v35*ct)))))  &
            + p*(v37 + ct*(v38 + ct*(v39 + v40*ct))  &
           + sa*(v41 + v42*ct) &
            + p*(v43 + ct*(v44 + v45*ct + v46*sa) &
            + p*(v47 + v48*ct)))

spec_vol = v_hat_numerator/v_hat_denominator

dvhatden_dsa = b01 + ct*(b02 + b03*ct) &
     + sqrtsa*(b04 + ct*(b05 + ct*(b06 + b07*ct))) &
          + p*(b08 + b09*ct + b10*p) 

dvhatnum_dsa = b11 + ct*(b12 + ct*(b13 + ct*(b14 + b15*ct))) &
     + sqrtsa*(b16 + ct*(b17 + ct*(b18 + ct*(b19 + b20*ct)))) + b21*sa &
          + p*(b22 + ct*(b23 + b24*p))

gsw_beta = (dvhatden_dsa*spec_vol - dvhatnum_dsa)/v_hat_numerator

return
end function

!--------------------------------------------------------------------------

!==========================================================================
function gsw_specvol(sa,ct,p) 
!==========================================================================

!  Calculates specific volume of seawater using the computationally-
!  efficient 48-term expression for density in terms of SA, CT and p
!  (McDougall et al., 2011)
!
! sa     : Absolute Salinity                               [g/kg]
! ct     : Conservative Temperature                        [deg C]
! p      : sea pressure                                    [dbar]
! 
! gsw_specvol  :  specific volume of seawater (48 term equation)

implicit none

integer, parameter :: r14 = selected_real_kind(14,30)

real (r14) :: sa, ct, p, gsw_specvol, gsw_rho

gsw_specvol = 1d0/gsw_rho(sa,ct,p)

return
end function

!--------------------------------------------------------------------------

!==========================================================================
function gsw_specvol_anom(sa,ct,p)  
!==========================================================================

!  Calculates specific volume anomaly of seawater using the computationally
!  efficient 48-term expression for density in terms of SA, CT and p
!  (McDougall et al., 2011)
!
! sa     : Absolute Salinity                               [g/kg]
! ct     : Conservative Temperature                        [deg C]
! p      : sea pressure                                    [dbar]
! 
! gsw_specvol_anom  :  specific volume anomaly of seawater (48 term equation)

implicit none

integer, parameter :: r14 = selected_real_kind(14,30)

real (r14) :: sa, ct, p, gsw_specvol, gsw_rho, gsw_specvol_sso_0_p
real (r14) :: gsw_specvol_anom

gsw_specvol_anom = gsw_specvol(sa,ct,p) - gsw_specvol_sso_0_p(p)

return
end function

!--------------------------------------------------------------------------

!==========================================================================
function gsw_sound_speed(sa,ct,p) 
!==========================================================================

!  Calculates sound speed of seawater using the computationally
!  efficient 48-term expression for density in terms of SA, CT and p
!  (McDougall et al., 2011)
!
! sa     : Absolute Salinity                               [g/kg]
! ct     : Conservative Temperature                        [deg C]
! p      : sea pressure                                    [dbar]
! 
! gsw_sound_speed  :  sound speed of seawater (48 term equation)

implicit none

integer, parameter :: r14 = selected_real_kind(14,30)

real (r14), parameter :: v01 =  9.998420897506056d2,   v02 =  2.839940833161907d0
real (r14), parameter :: v03 = -3.147759265588511d-2,  v04 =  1.181805545074306d-3
real (r14), parameter :: v05 = -6.698001071123802d0,   v06 = -2.986498947203215d-2
real (r14), parameter :: v07 =  2.327859407479162d-4,  v08 = -3.988822378968490d-2
real (r14), parameter :: v09 =  5.095422573880500d-4,  v10 = -1.426984671633621d-5
real (r14), parameter :: v11 =  1.645039373682922d-7,  v12 = -2.233269627352527d-2
real (r14), parameter :: v13 = -3.436090079851880d-4,  v14 =  3.726050720345733d-6
real (r14), parameter :: v15 = -1.806789763745328d-4,  v16 =  6.876837219536232d-7
real (r14), parameter :: v17 = -3.087032500374211d-7,  v18 = -1.988366587925593d-8
real (r14), parameter :: v19 = -1.061519070296458d-11, v20 =  1.550932729220080d-10
real (r14), parameter :: v21 =  1.0d0, v22 = 2.775927747785646d-3, v23 = -2.349607444135925d-5
real (r14), parameter :: v24 =  1.119513357486743d-6,  v25 =  6.743689325042773d-10
real (r14), parameter :: v26 = -7.521448093615448d-3,  v27 = -2.764306979894411d-5
real (r14), parameter :: v28 =  1.262937315098546d-7,  v29 =  9.527875081696435d-10
real (r14), parameter :: v30 = -1.811147201949891d-11, v31 = -3.303308871386421d-5
real (r14), parameter :: v32 =  3.801564588876298d-7,  v33 = -7.672876869259043d-9
real (r14), parameter :: v34 = -4.634182341116144d-11, v35 =  2.681097235569143d-12
real (r14), parameter :: v36 =  5.419326551148740d-6,  v37 = -2.742185394906099d-5
real (r14), parameter :: v38 = -3.212746477974189d-7,  v39 =  3.191413910561627d-9
real (r14), parameter :: v40 = -1.931012931541776d-12, v41 = -1.105097577149576d-7
real (r14), parameter :: v42 =  6.211426728363857d-10, v43 = -1.119011592875110d-10
real (r14), parameter :: v44 = -1.941660213148725d-11, v45 = -1.864826425365600d-14
real (r14), parameter :: v46 =  1.119522344879478d-14, v47 = -1.200507748551599d-15
real (r14), parameter :: v48 =  6.057902487546866d-17 
real (r14), parameter :: c01 = -2.233269627352527d-2,  c02 = -3.436090079851880d-4
real (r14), parameter :: c03 =  3.726050720345733d-6,  c04 = -1.806789763745328d-4
real (r14), parameter :: c05 =  6.876837219536232d-7,  c06 = -6.174065000748422d-7
real (r14), parameter :: c07 = -3.976733175851186d-8,  c08 = -2.123038140592916d-11
real (r14), parameter :: c09 =  3.101865458440160d-10, c10 = -2.742185394906099d-5
real (r14), parameter :: c11 = -3.212746477974189d-7,  c12 =  3.191413910561627d-9
real (r14), parameter :: c13 = -1.931012931541776d-12, c14 = -1.105097577149576d-7
real (r14), parameter :: c15 =  6.211426728363857d-10, c16 = -2.238023185750219d-10
real (r14), parameter :: c17 = -3.883320426297450d-11, c18 = -3.729652850731201d-14
real (r14), parameter :: c19 =  2.239044689758956d-14, c20 = -3.601523245654798d-15
real (r14), parameter :: c21 =  1.817370746264060d-16

real (r14) :: sa, ct, p, gsw_sound_speed, v_hat_denominator, v_hat_numerator
real (r14) :: sqrtsa, dvden_dp, dvnum_dp, dp_drho

sqrtsa = sqrt(sa)

v_hat_denominator = v01 + ct*(v02 + ct*(v03 + v04*ct))  &
             + sa*(v05 + ct*(v06 + v07*ct) &
         + sqrtsa*(v08 + ct*(v09 + ct*(v10 + v11*ct)))) &
              + p*(v12 + ct*(v13 + v14*ct) + sa*(v15 + v16*ct) &
              + p*(v17 + ct*(v18 + v19*ct) + v20*sa))

v_hat_numerator = v21 + ct*(v22 + ct*(v23 + ct*(v24 + v25*ct))) &
           + sa*(v26 + ct*(v27 + ct*(v28 + ct*(v29 + v30*ct))) + v36*sa &
       + sqrtsa*(v31 + ct*(v32 + ct*(v33 + ct*(v34 + v35*ct)))))  &
            + p*(v37 + ct*(v38 + ct*(v39 + v40*ct))  &
           + sa*(v41 + v42*ct) &
            + p*(v43 + ct*(v44 + v45*ct + v46*sa) &
            + p*(v47 + v48*ct)))

dvden_dp = c01 + ct*(c02 + c03*ct) &
    + sa*(c04 + c05*ct) &
    + p*(c06 + ct*(c07 + c08*ct) + c09*sa)

dvnum_dp = c10 + ct*(c11 + ct*(c12 + c13*ct)) &
    + sa*(c14 + c15*ct) &
    + p*(c16 + ct*(c17 + c18*ct + c19*SA) &
    + p*(c20 + c21*ct))

dp_drho = (v_hat_numerator*v_hat_numerator)/(dvden_dp*v_hat_numerator - dvnum_dp*v_hat_denominator)
    
gsw_sound_speed = 100*sqrt(dp_drho)

return
end function

!--------------------------------------------------------------------------

!==========================================================================
function gsw_internal_energy(sa,ct,p)  
!==========================================================================

!  Calculates internal energy of seawater using the computationally
!  efficient 48-term expression for density in terms of SA, CT and p
!  (McDougall et al., 2011)
!
! sa     : Absolute Salinity                               [g/kg]
! ct     : Conservative Temperature                        [deg C]
! p      : sea pressure                                    [dbar]
! 
! gsw_internal_energy  :  internal_energy of seawater (48 term equation)

implicit none

integer, parameter :: r14 = selected_real_kind(14,30)

real (r14), parameter :: p0 = 101325d0, db2pa = 1d4

real (r14) :: sa, ct, p, gsw_internal_energy, gsw_enthalpy, gsw_specvol

gsw_internal_energy = gsw_enthalpy(sa,ct,p) - (p0 + db2pa*p)*gsw_specvol(sa,ct,p)

return
end function

!--------------------------------------------------------------------------

!==========================================================================
function gsw_enthalpy(sa,ct,p)  
!==========================================================================

! Calculates specific enthalpy of seawater using the computationally-
! efficient 48-term expression for density in terms of SA, CT and p
! (McDougall et al., 2011)
!
! sa     : Absolute Salinity                               [g/kg]
! ct     : Conservative Temperature                        [deg C]
! p      : sea pressure                                    [dbar]
! 
! gsw_enthalpy  :  specific enthalpy of seawater (48 term equation)

implicit none

integer, parameter :: r14 = selected_real_kind(14,30)

real (r14), parameter :: v01 =  9.998420897506056d+2, v02 =  2.839940833161907d0
real (r14), parameter :: v03 = -3.147759265588511d-2, v04 =  1.181805545074306d-3
real (r14), parameter :: v05 = -6.698001071123802d0, v06 = -2.986498947203215d-2
real (r14), parameter :: v07 =  2.327859407479162d-4, v08 = -3.988822378968490d-2
real (r14), parameter :: v09 =  5.095422573880500d-4, v10 = -1.426984671633621d-5
real (r14), parameter :: v11 =  1.645039373682922d-7, v12 = -2.233269627352527d-2
real (r14), parameter :: v13 = -3.436090079851880d-4, v14 =  3.726050720345733d-6
real (r14), parameter :: v15 = -1.806789763745328d-4, v16 =  6.876837219536232d-7
real (r14), parameter :: v17 = -3.087032500374211d-7, v18 = -1.988366587925593d-8
real (r14), parameter :: v19 = -1.061519070296458d-11, v20 =  1.550932729220080d-10
real (r14), parameter :: v21 =  1.0d0, v22 =  2.775927747785646d-3, v23 = -2.349607444135925d-5
real (r14), parameter :: v24 =  1.119513357486743d-6, v25 =  6.743689325042773d-10
real (r14), parameter :: v26 = -7.521448093615448d-3, v27 = -2.764306979894411d-5
real (r14), parameter :: v28 =  1.262937315098546d-7, v29 =  9.527875081696435d-10
real (r14), parameter :: v30 = -1.811147201949891d-11, v31 = -3.303308871386421d-5
real (r14), parameter :: v32 =  3.801564588876298d-7, v33 = -7.672876869259043d-9
real (r14), parameter :: v34 = -4.634182341116144d-11, v35 =  2.681097235569143d-12
real (r14), parameter :: v36 =  5.419326551148740d-6, v37 = -2.742185394906099d-5
real (r14), parameter :: v38 = -3.212746477974189d-7, v39 =  3.191413910561627d-9
real (r14), parameter :: v40 = -1.931012931541776d-12, v41 = -1.105097577149576d-7
real (r14), parameter :: v42 =  6.211426728363857d-10, v43 = -1.119011592875110d-10
real (r14), parameter :: v44 = -1.941660213148725d-11, v45 = -1.864826425365600d-14
real (r14), parameter :: v46 =  1.119522344879478d-14, v47 = -1.200507748551599d-15
real (r14), parameter :: v48 =  6.057902487546866d-17 

real (r14) :: sa, ct, p, db2pa, cp0, sqrtsa, a0, a1, a2, a3, b0, b1, b2, b1sq
real (r14) :: sqrt_disc, ca, cb, cn, cm, part, gsw_enthalpy

db2pa = 1d4                             ! factor to convert from dbar to Pa
cp0 = 3991.86795711963d0           ! from Eqn. (3.3.3) of IOC et al. (2010)

sqrtsa = sqrt(sa)

a0 = v21 + ct*(v22 + ct*(v23 + ct*(v24 + v25*ct))) &
       + sa*(v26 + ct*(v27 + ct*(v28 + ct*(v29 + v30*ct))) + v36*sa &
   + sqrtsa*(v31 + ct*(v32 + ct*(v33 + ct*(v34 + v35*ct)))))
 
a1 = v37 + ct*(v38 + ct*(v39 + v40*ct)) + sa*(v41 + v42*ct)

a2 = v43 + ct*(v44 + v45*ct + v46*sa)

a3 = v47 + v48*ct

b0 = v01 + ct*(v02 + ct*(v03 + v04*ct))  &
         + sa*(v05 + ct*(v06 + v07*ct) &
     + sqrtsa*(v08 + ct*(v09 + ct*(v10 + v11*ct))))
 
b1 = 0.5d0*(v12 + ct*(v13 + v14*ct) + sa*(v15 + v16*ct))

b2 = v17 + ct*(v18 + v19*ct) + v20*sa

b1sq = b1*b1
sqrt_disc = sqrt(b1sq - b0*b2)

cn = a0 + (2d0*a3*b0*b1/b2 - a2*b0)/b2

cm = a1 + (4d0*a3*b1sq/b2 - a3*b0 - 2*a2*b1)/b2

ca = b1 - sqrt_disc
cb = b1 + sqrt_disc

part = (cn*b2 - cm*b1)/(b2*(cb - ca))

gsw_enthalpy = cp0*ct + &
           db2Pa*(p*(a2 - 2d0*a3*b1/b2 + 0.5d0*a3*p)/b2 + &
           (cm/(2d0*b2))*log(1 + p*(2d0*b1 + b2*p)/b0) + &
           part*log(1d0 + (b2*p*(cb - ca))/(ca*(cb + b2*p))))

return
end function

!--------------------------------------------------------------------------

!==========================================================================
function gsw_dynamic_enthalpy(sa,ct,p) 
!==========================================================================

! Calculates dynamic enthalpy of seawater using the computationally
! efficient 48-term expression for density in terms of SA, CT and p
! (McDougall et al., 2011)
!
! sa     : Absolute Salinity                               [g/kg]
! ct     : Conservative Temperature                        [deg C]
! p      : sea pressure                                    [dbar]
! 
! gsw_dynamic_enthalpy  :  dynamic enthalpy of seawater (48 term equation)

implicit none

integer, parameter :: r14 = selected_real_kind(14,30)

real (r14), parameter :: v01 =  9.998420897506056d+2, v02 =  2.839940833161907d0
real (r14), parameter :: v03 = -3.147759265588511d-2, v04 =  1.181805545074306d-3
real (r14), parameter :: v05 = -6.698001071123802d0, v06 = -2.986498947203215d-2
real (r14), parameter :: v07 =  2.327859407479162d-4, v08 = -3.988822378968490d-2
real (r14), parameter :: v09 =  5.095422573880500d-4, v10 = -1.426984671633621d-5
real (r14), parameter :: v11 =  1.645039373682922d-7, v12 = -2.233269627352527d-2
real (r14), parameter :: v13 = -3.436090079851880d-4, v14 =  3.726050720345733d-6
real (r14), parameter :: v15 = -1.806789763745328d-4, v16 =  6.876837219536232d-7
real (r14), parameter :: v17 = -3.087032500374211d-7, v18 = -1.988366587925593d-8
real (r14), parameter :: v19 = -1.061519070296458d-11, v20 =  1.550932729220080d-10
real (r14), parameter :: v21 =  1.0d0, v22 =  2.775927747785646d-3, v23 = -2.349607444135925d-5
real (r14), parameter :: v24 =  1.119513357486743d-6, v25 =  6.743689325042773d-10
real (r14), parameter :: v26 = -7.521448093615448d-3, v27 = -2.764306979894411d-5
real (r14), parameter :: v28 =  1.262937315098546d-7, v29 =  9.527875081696435d-10
real (r14), parameter :: v30 = -1.811147201949891d-11, v31 = -3.303308871386421d-5
real (r14), parameter :: v32 =  3.801564588876298d-7, v33 = -7.672876869259043d-9
real (r14), parameter :: v34 = -4.634182341116144d-11, v35 =  2.681097235569143d-12
real (r14), parameter :: v36 =  5.419326551148740d-6, v37 = -2.742185394906099d-5
real (r14), parameter :: v38 = -3.212746477974189d-7, v39 =  3.191413910561627d-9
real (r14), parameter :: v40 = -1.931012931541776d-12, v41 = -1.105097577149576d-7
real (r14), parameter :: v42 =  6.211426728363857d-10, v43 = -1.119011592875110d-10
real (r14), parameter :: v44 = -1.941660213148725d-11, v45 = -1.864826425365600d-14
real (r14), parameter :: v46 =  1.119522344879478d-14, v47 = -1.200507748551599d-15
real (r14), parameter :: v48 =  6.057902487546866d-17 

real (r14) :: sa, ct, p, db2pa, sqrtsa, a0, a1, a2, a3, b0, b1, b2, b1sq
real (r14) :: sqrt_disc, ca, cb, cn, cm, part, gsw_dynamic_enthalpy

db2pa = 1d4                             ! factor to convert from dbar to Pa

sqrtsa = sqrt(sa)

a0 = v21 + ct*(v22 + ct*(v23 + ct*(v24 + v25*ct))) &
       + sa*(v26 + ct*(v27 + ct*(v28 + ct*(v29 + v30*ct))) + v36*sa  &
   + sqrtsa*(v31 + ct*(v32 + ct*(v33 + ct*(v34 + v35*ct)))));
 
a1 = v37 + ct*(v38 + ct*(v39 + v40*ct)) + sa*(v41 + v42*ct)

a2 = v43 + ct*(v44 + v45*ct + v46*sa)

a3 = v47 + v48*ct

b0 = v01 + ct*(v02 + ct*(v03 + v04*ct))  &
         + sa*(v05 + ct*(v06 + v07*ct)  &
     + sqrtsa*(v08 + ct*(v09 + ct*(v10 + v11*ct))))
 
b1 = 0.5d0*(v12 + ct*(v13 + v14*ct) + sa*(v15 + v16*ct))

b2 = v17 + ct*(v18 + v19*ct) + v20*sa

b1sq = b1*b1 
sqrt_disc = sqrt(b1sq - b0*b2)

cn = a0 + (2*a3*b0*b1/b2 - a2*b0)/b2

cm = a1 + (4*a3*b1sq/b2 - a3*b0 - 2*a2*b1)/b2

ca = b1 - sqrt_disc
cb = b1 + sqrt_disc

part = (cn*b2 - cm*b1)/(b2*(cb - ca))

gsw_dynamic_enthalpy = db2pa*(p*(a2 - 2d0*a3*b1/b2 + 0.5d0*a3*p)/b2  &
                      + (cm/(2d0*b2))*log(1d0 + p*(2d0*b1 + b2*p)/b0)  &
                      + part*log(1d0 + (b2*p*(cb - ca))/(ca*(cb + b2*p))))

return
end function

!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
! freezing temperatures
!--------------------------------------------------------------------------

!==========================================================================
function gsw_ct_freezing(sa,p,saturation_fraction) 
!==========================================================================

! Calculates the Conservative Temperature at which of seawater freezes 
! from Absolute Salinity and pressure.
!
! sa     : Absolute Salinity                                 [g/kg]
! p      : sea pressure                                      [dbar]
! saturation_fraction : saturation fraction
!
! gsw_ct_freezing : Conservative Temperature freezing point  [deg C]

implicit none

integer, parameter :: r14 = selected_real_kind(14,30)

real (r14), parameter :: c0 = 0.017947064327968736d0, c1 = -6.076099099929818d0
real (r14), parameter :: c2 = 4.883198653547851d0, c3 = -11.88081601230542d0
real (r14), parameter :: c4 = 13.34658511480257d0, c5 = -8.722761043208607d0
real (r14), parameter :: c6 = 2.082038908808201d0, c7 = -7.389420998107497d0
real (r14), parameter :: c8 = -2.110913185058476d0, c9 = 0.2295491578006229d0
real (r14), parameter :: c10 = -0.9891538123307282d0, c11 = -0.08987150128406496d0
real (r14), parameter :: c12 = 0.3831132432071728d0, c13 = 1.054318231187074d0
real (r14), parameter :: c14 = 1.065556599652796d0, c15 = -0.7997496801694032d0
real (r14), parameter :: c16 = 0.3850133554097069d0, c17 = -2.078616693017569d0
real (r14), parameter :: c18 = 0.8756340772729538d0, c19 = -2.079022768390933d0
real (r14), parameter :: c20 = 1.596435439942262d0, c21 = 0.1338002171109174d0
real (r14), parameter :: c22 = 1.242891021876471d0

real (r14) :: sa, p, saturation_fraction, ct_freezing, gsw_ct_freezing
real (r14) :: sa_r, x, p_r, a, b

sa_r = sa*1d-2
x = sqrt(sa_r)
p_r = p*1d-4

ct_freezing = c0 &
 + sa_r*(c1 + x*(c2 + x*(c3 + x*(c4 + x*(c5 + c6*x))))) &
 + p_r*(c7 + p_r*(c8 + c9*p_r)) &
 + sa_r*p_r*(c10 + p_r*(c12 + p_r*(c15 + c21*sa_r)) + sa_r*(c13 + c17*p_r + c19*sa_r) &
 + x*(c11 + p_r*(c14 + c18*p_r)  + sa_r*(c16 + c20*p_r + c22*sa_r)))

! Adjust for the effects of dissolved air 
a = 0.014289763856964d0             ! Note that a = 0.502500117621/35.16504.
b = 0.057000649899720d0
ct_freezing = ct_freezing - saturation_fraction*(1d-3)*(2.4d0 - a*sa)*(1d0 + b*(1d0 - sa/35.16504d0))

if (p.gt.10000d0 .or. sa.gt.120d0 .or. (p+sa*71.428571428571402d0).gt.13571.42857142857d0) then
    ct_freezing = 9d15
end if

gsw_ct_freezing = ct_freezing

return
end function

!--------------------------------------------------------------------------

!==========================================================================
function gsw_t_freezing(sa,p,saturation_fraction) 
!==========================================================================

! Calculates the in-situ temperature at which of seawater freezes 
! from Absolute Salinity and pressure.
!
! sa     : Absolute Salinity                                 [g/kg]
! p      : sea pressure                                      [dbar]
! saturation_fraction : saturation fraction
!
! gsw_t_freezing : in-situ temperature freezing point  [deg C]

implicit none

integer, parameter :: r14 = selected_real_kind(14,30)

real (r14) :: sa, p, saturation_fraction, ct_freezing, gsw_ct_freezing
real (r14) :: gsw_t_from_ct, t_freezing, gsw_t_freezing

ct_freezing = gsw_CT_freezing(sa,p,saturation_fraction)
t_freezing = gsw_t_from_ct(sa,ct_freezing,p)

if (ct_freezing.gt.9d10) then
 t_freezing = 9d15
end if

gsw_t_freezing = t_freezing

return
end function

!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
! isobaric melting enthalpy and isobaric evaporation enthalpy
!--------------------------------------------------------------------------

!==========================================================================
function gsw_latentheat_melting(sa,p)  
!==========================================================================

! Calculates latent heat, or enthalpy, of melting.
!
! sa     : Absolute Salinity                               [g/kg]
! p      : sea pressure                                    [dbar]
! 
! gsw_latentheat_melting : latent heat of melting          [kg/m^3]

implicit none

integer, parameter :: r14 = selected_real_kind(14,30)

real (r14), parameter :: c0 =  3.334265169240710d5, c1 = -2.789444646733159d0
real (r14), parameter :: c2 = -1.822150156453350d4, c3 = -4.984585692734338d3
real (r14), parameter :: c4 = -7.371966528571920d1, c5 = -7.605802553358546d3
real (r14), parameter :: c6 =  1.195857305019339d3, c7 =  1.233720336206392d3
real (r14), parameter :: c8 =  2.294798676591890d2, c9 =  9.655751370889338d2
real (r14), parameter :: c10 = -5.792068522727968d2, c11 = -1.649446955902331d3
real (r14), parameter :: c12 = -1.029021448430547d3, c13 = -3.171558017172501d2
real (r14), parameter :: c14 = -1.751401389905041d2, c15 =  6.836527214265952d2
real (r14), parameter :: c16 =  1.078283734113611d3, c17 =  5.613896351265648d2
real (r14), parameter :: c18 =  6.968934948667265d2, c19 =  1.793032021946783d2
real (r14), parameter :: c20 =  8.692558481134256d1, c21 = -2.371103254714944d2
real (r14), parameter :: c22 = -5.775033277201674d2, c23 = -3.019749254648732d2
real (r14), parameter :: c24 = -6.420420579160927d2, c25 = -2.657570848596042d2
real (r14), parameter :: c26 = -1.646738151143109d1, c27 =  4.618228988300871d0

real (r14) :: sa, p, s_u, x, y, gsw_latentheat_melting

s_u = 40d0*(35.16504d0/35d0)
x = sqrt(sa/s_u)
y = p*1d-4

gsw_latentheat_melting = c0 + x*(c1 + c4*y + x*(c3   &
    + y*(c7 + c12*y) + x*(c6 + y*(c11 + y*(c17 + c24*y))  &
    + x*(c10  + y*(c16 + c23*y) + x*(c15 + c22*y + c21*x)))))  &
    + y*(c2 + y*(c5 + c8*x + y*(c9 + x*(c13 + c18*x)  &
    + y*(c14 + x*(c19 + c25*x) + y*(c20 + c26*x + c27*y)))))

return
end

!--------------------------------------------------------------------------

!==========================================================================
function gsw_latentheat_evap_ct(sa,ct) 
!==========================================================================

! Calculates latent heat, or enthalpy, of evaporation.
!
! sa     : Absolute Salinity                               [g/kg]
! ct     : Conservative Temperature                        [deg C]
! 
! gsw_latentheat_evaporation : latent heat of evaporation  [J/kg]

implicit none

integer, parameter :: r14 = selected_real_kind(14,30)

real (r14), parameter :: c0 =   2.499065844825125d6, c1 =  -1.544590633515099d-1
real (r14), parameter :: c2 =  -9.096800915831875d4, c3 =   1.665513670736000d2
real (r14), parameter :: c4 =   4.589984751248335d1, c5 =   1.894281502222415d1
real (r14), parameter :: c6 =   1.192559661490269d3, c7 =  -6.631757848479068d3
real (r14), parameter :: c8 =  -1.104989199195898d2, c9 =  -1.207006482532330d3
real (r14), parameter :: c10 = -3.148710097513822d3, c11 =  7.437431482069087d2
real (r14), parameter :: c12 =  2.519335841663499d3, c13 =  1.186568375570869d1
real (r14), parameter :: c14 =  5.731307337366114d2, c15 =  1.213387273240204d3
real (r14), parameter :: c16 =  1.062383995581363d3, c17 = -6.399956483223386d2
real (r14), parameter :: c18 = -1.541083032068263d3, c19 =  8.460780175632090d1
real (r14), parameter :: c20 = -3.233571307223379d2, c21 = -2.031538422351553d2
real (r14), parameter :: c22 =  4.351585544019463d1, c23 = -8.062279018001309d2
real (r14), parameter :: c24 =  7.510134932437941d2, c25 =  1.797443329095446d2
real (r14), parameter :: c26 = -2.389853928747630d1, c27 =  1.021046205356775d2

real (r14) :: sa, ct, s_u, x, y, gsw_latentheat_evap_ct

s_u = 40d0*(35.16504d0/35d0)
x = sqrt(sa/s_u)
y = ct/40

gsw_latentheat_evap_ct = c0 + x*(c1 + c4*y + x*(c3   &
    + y*(c7 + c12*y) + x*(c6 + y*(c11 + y*(c17 + c24*y)) &
    + x*(c10 + y*(c16 + c23*y) + x*(c15 + c22*y + c21*x)))))  &
    + y*(c2 + y*(c5 + c8*x + y*(c9 + x*(c13 + c18*x)  &
    + y*(c14 + x*(c19 + c25*x) + y*(c20 + c26*x + c27*y)))))

return
end

!--------------------------------------------------------------------------

!==========================================================================
function gsw_latentheat_evap_t(sa,t)  
!==========================================================================

! Calculates latent heat, or enthalpy, of evaporation.
!
! sa     : Absolute Salinity                               [g/kg]
! t      : in-situ temperature                             [deg C]
! 
! gsw_latentheat_evap_t : latent heat of evaporation       [J/kg]

implicit none

integer, parameter :: r14 = selected_real_kind(14,30)

real (r14) :: sa, t, ct, gsw_ct_from_pt, gsw_latentheat_evap_ct
real (r14) :: gsw_latentheat_evap_t

ct = gsw_ct_from_pt(sa,t)

gsw_latentheat_evap_t = gsw_latentheat_evap_ct(sa,ct)

return
end

!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
! basic thermodynamic properties in terms of in-situ t, based on the exact Gibbs function
!--------------------------------------------------------------------------

!==========================================================================
function gsw_rho_t_exact(sa,t,p) 
!==========================================================================

! Calculates in-situ density of seawater from Absolute Salinity and 
! in-situ temperature.
!
! sa     : Absolute Salinity                               [g/kg]
! t      : in-situ temperature                             [deg C]
! p      : sea pressure                                    [dbar]
! 
! gsw_rho_t_exact : in-situ density                        [kg/m^3]

implicit none

integer, parameter :: r14 = selected_real_kind(14,30)

integer :: n0, n1
real (r14) :: sa, t, p, gsw_rho_t_exact, gsw_gibbs

n0 = 0
n1 = 1

gsw_rho_t_exact = 1.d0/gsw_gibbs(n0,n0,n1,sa,t,p)

return
end

!--------------------------------------------------------------------------

!==========================================================================
function gsw_pot_rho_t_exact(sa,t,p,p_ref)  
!==========================================================================

! Calculates the potential density of seawater
!
! sa     : Absolute Salinity                               [g/kg]
! t      : in-situ temperature                             [deg C]
! p      : sea pressure                                    [dbar]
! p_ref  : reference sea pressure                          [dbar]
! 
! gsw_pot_rho_t_exact : potential density                  [kg/m^3]

implicit none

integer, parameter :: r14 = selected_real_kind(14,30)

real (r14) :: sa, t, p, p_ref, gsw_pot_rho_t_exact, pt
real (r14) :: gsw_pt_from_t, gsw_rho_t_exact

pt = gsw_pt_from_t(sa,t,p,p_ref)

gsw_pot_rho_t_exact = gsw_rho_t_exact(sa,pt,p_ref)

return
end

!--------------------------------------------------------------------------

!==========================================================================
function gsw_alpha_wrt_t_exact(sa,t,p) 
!==========================================================================

! Calculates thermal expansion coefficient of seawater with respect to 
! in-situ temperature
!
! sa     : Absolute Salinity                               [g/kg]
! t      : insitu temperature                              [deg C]
! p      : sea pressure                                    [dbar]
!
! gsw_alpha_wrt_t_exact : thermal expansion coefficient    [1/K]
!                         wrt (in-situ) temperature

implicit none

integer, parameter :: r14 = selected_real_kind(14,30)

integer :: n0, n1
real (r14) :: sa, t, p, gsw_gibbs, gsw_alpha_wrt_t_exact

n0 = 0
n1 = 1

gsw_alpha_wrt_t_exact = gsw_gibbs(n0,n1,n1,sa,t,p)/gsw_gibbs(n0,n0,n1,sa,t,p)

end

!--------------------------------------------------------------------------

!==========================================================================
function gsw_beta_const_t_exact(sa,t,p)  
!==========================================================================

! Calculates saline (haline) contraction coefficient of seawater at 
! constant in-situ temperature.
!
! sa     : Absolute Salinity                               [g/kg]
! t      : in-situ temperature                             [deg C]
! p      : sea pressure                                    [dbar]
! 
! gsw_beta_const_t_exact : haline contraction coefficient  [kg/g]

implicit none

integer, parameter :: r14 = selected_real_kind(14,30)

integer :: n0, n1
real (r14) :: sa, t, p, gsw_beta_const_t_exact, gsw_gibbs, uPS

n0 = 0
n1 = 1 

gsw_beta_const_t_exact = -gsw_gibbs(n1,n0,n1,sa,t,p)/gsw_gibbs(n0,n0,n1,sa,t,p)

return
end function

!--------------------------------------------------------------------------

!==========================================================================
function gsw_specvol_t_exact(sa,t,p)  
!==========================================================================

! Calulates the specific volume of seawater
!
! sa     : Absolute Salinity                               [g/kg]
! t      : in-situ temperature                             [deg C]
! p      : sea pressure                                    [dbar]
! 
! gsw_specvol_t_exact : specific volume                    [kg/m^3]

implicit none

integer, parameter :: r14 = selected_real_kind(14,30)

integer :: n0, n1
real (r14) :: sa, t, p, gsw_specvol_t_exact, gsw_gibbs

n0 = 0
n1 = 1

gsw_specvol_t_exact = gsw_gibbs(n0,n0,n1,sa,t,p)

return
end

!--------------------------------------------------------------------------

!==========================================================================
function gsw_sound_speed_t_exact(sa,t,p)  
!==========================================================================

! Calculates the sound speed of seawater
!
! sa     : Absolute Salinity                               [g/kg]
! t      : in-situ temperature                             [deg C]
! p      : sea pressure                                    [dbar]
! 
! gsw_sound_speed_t_exact : sound speed                    [m/s]

implicit none

integer, parameter :: r14 = selected_real_kind(14,30)

integer :: n0, n1, n2
real (r14) :: sa, t, p, g_tt, g_tp, gsw_sound_speed_t_exact, gsw_gibbs

n0 = 0
n1 = 1
n2 = 2

g_tt = gsw_gibbs(n0,n2,n0,sa,t,p)
g_tp = gsw_gibbs(n0,n1,n1,sa,t,p)

gsw_sound_speed_t_exact = gsw_gibbs(n0,n0,n1,sa,t,p) * &
               sqrt(g_tt/(g_tp*g_tp - g_tt*gsw_gibbs(n0,n0,n2,sa,t,p)))
                                    
return
end

!--------------------------------------------------------------------------

!==========================================================================
function gsw_kappa_t_exact(sa,t,p)  
!==========================================================================

! isentropic compressibility of seawater
!
! sa     : Absolute Salinity                               [g/kg]
! t      : in-situ temperature                             [deg C]
! p      : sea pressure                                    [dbar]
!
! gsw_kappa_t_exact : isentropic compressibility           [1/Pa]

implicit none

integer, parameter :: r14 = selected_real_kind(14,30)

integer :: n0, n1, n2
real (r14) :: sa, t, p, g_tt, g_tp, gsw_kappa_t_exact, gsw_gibbs

n0 = 0
n1 = 1
n2 = 2

g_tt = gsw_gibbs(n0,n2,n0,sa,t,p)
g_tp = gsw_gibbs(n0,n1,n1,sa,t,p)

gsw_kappa_t_exact = (g_tp*g_tp - g_tt*gsw_gibbs(n0,n0,n2,sa,t,p)) / &
                         (gsw_gibbs(n0,n0,n1,sa,t,p)*g_tt)

return
end

!--------------------------------------------------------------------------

!==========================================================================
function gsw_enthalpy_t_exact(sa,t,p) 
!==========================================================================

! Calculates the specific enthalpy of seawater
!
! sa     : Absolute Salinity                               [g/kg]
! t      : in-situ temperature                             [deg C]
! p      : sea pressure                                    [dbar]
! 
! gsw_enthalpy_t_exact : specific enthalpy                 [J/kg]

implicit none

integer, parameter :: r14 = selected_real_kind(14,30)

integer :: n0, n1
real (r14) :: sa, t, p, gsw_enthalpy_t_exact, gsw_gibbs

n0 = 0
n1 = 1

gsw_enthalpy_t_exact = gsw_gibbs(n0,n0,n0,sa,t,p) - (t+273.15d0)*gsw_gibbs(n0,n1,n0,sa,t,p)

return
end

!--------------------------------------------------------------------------

!==========================================================================
function gsw_entropy_t_exact(sa,t,p) 
!==========================================================================

! Calculates the specific entropy of seawater
!
! sa     : Absolute Salinity                               [g/kg]
! t      : in-situ temperature                             [deg C]
! p      : sea pressure                                    [dbar]
! 
! gsw_entropy_t_exact : specific entropy                   [J/(kg K)]

implicit none

integer, parameter :: r14 = selected_real_kind(14,30)

integer :: n0, n1
real (r14) :: sa, t, p, gsw_entropy_t_exact, gsw_gibbs

n0 = 0
n1 = 1

gsw_entropy_t_exact = -gsw_gibbs(n0,n1,n0,sa,t,p)

return
end

!--------------------------------------------------------------------------

!==========================================================================
function gsw_cp_t_exact(sa,t,p)    
!==========================================================================

! Calculates isobaric heat capacity of seawater
!
! sa     : Absolute Salinity                               [g/kg]
! t      : in-situ temperature                             [deg C]
! p      : sea pressure                                    [dbar]
! 
! gsw_cp_t_exact : heat capacity                           [J/(kg K)]

implicit none

integer, parameter :: r14 = selected_real_kind(14,30)

integer :: n0, n2
real (r14) :: sa, t, p, gsw_cp_t_exact, gsw_gibbs

n0 = 0
n2 = 2

gsw_cp_t_exact = -(t+273.15d0)*gsw_gibbs(n0,n2,n0,sa,t,p)

return
end function

!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
! Library functions of the GSW toolbox
!--------------------------------------------------------------------------

!==========================================================================
function gsw_gibbs(ns,nt,np,sa,t,p)
!==========================================================================

! seawater specific Gibbs free energy and derivatives up to order 2
!
! ns     : order of s derivative
! nt     : order of t derivative
! np     : order of p derivative
! sa     : Absolute Salinity                               [g/kg]
! t      : temperature                                     [deg C]
! p      : sea pressure                                    [dbar]
! 
! gsw_gibbs  : specific Gibbs energy or its derivative

implicit none
integer, parameter :: r14 = selected_real_kind(14,30)

integer :: ns, nt, np
real (r14) :: sa, t, p, gsw_gibbs
real (r14) :: sfac, x2, x, y, z, g03, g08

sfac = 0.0248826675584615d0              

x2 = sfac*sa
x = sqrt(x2)
y = t*0.025d0
z = p*1d-4

if(ns.eq.0 .and. nt.eq.0 .and. np.eq.0) then
          
  g03 = 101.342743139674d0 + z*(100015.695367145d0 + &
      z*(-2544.5765420363d0 + z*(284.517778446287d0 + &
      z*(-33.3146754253611d0 + (4.20263108803084d0 - 0.546428511471039d0*z)*z)))) + &
      y*(5.90578347909402d0 + z*(-270.983805184062d0 + &
      z*(776.153611613101d0 + z*(-196.51255088122d0 + (28.9796526294175d0 - 2.13290083518327d0*z)*z))) + &
      y*(-12357.785933039d0 + z*(1455.0364540468d0 + &
      z*(-756.558385769359d0 + z*(273.479662323528d0 + z*(-55.5604063817218d0 + 4.34420671917197d0*z)))) + &
      y*(736.741204151612d0 + z*(-672.50778314507d0 + &
      z*(499.360390819152d0 + z*(-239.545330654412d0 + (48.8012518593872d0 - 1.66307106208905d0*z)*z))) + &
      y*(-148.185936433658d0 + z*(397.968445406972d0 + &
      z*(-301.815380621876d0 + (152.196371733841d0 - 26.3748377232802d0*z)*z)) + &
      y*(58.0259125842571d0 + z*(-194.618310617595d0 + &
      z*(120.520654902025d0 + z*(-55.2723052340152d0 + 6.48190668077221d0*z))) + &
      y*(-18.9843846514172d0 + y*(3.05081646487967d0 - 9.63108119393062d0*z) + &
      z*(63.5113936641785d0 + z*(-22.2897317140459d0 + 8.17060541818112d0*z))))))))
          
  g08 = x2*(1416.27648484197d0 + z*(-3310.49154044839d0 + &
        z*(384.794152978599d0 + z*(-96.5324320107458d0 + (15.8408172766824d0 - 2.62480156590992d0*z)*z))) + &
        x*(-2432.14662381794d0 + x*(2025.80115603697d0 + &
        y*(543.835333000098d0 + y*(-68.5572509204491d0 + &
        y*(49.3667694856254d0 + y*(-17.1397577419788d0 + 2.49697009569508d0*y))) - 22.6683558512829d0*z) + &
        x*(-1091.66841042967d0 - 196.028306689776d0*y + &
        x*(374.60123787784d0 - 48.5891069025409d0*x + 36.7571622995805d0*y) + 36.0284195611086d0*z) + &
        z*(-54.7919133532887d0 + (-4.08193978912261d0 - 30.1755111971161d0*z)*z)) + &
        z*(199.459603073901d0 + z*(-52.2940909281335d0 + (68.0444942726459d0 - 3.41251932441282d0*z)*z)) + &
        y*(-493.407510141682d0 + z*(-175.292041186547d0 + (83.1923927801819d0 - 29.483064349429d0*z)*z) + &
        y*(-43.0664675978042d0 + z*(383.058066002476d0 + z*(-54.1917262517112d0 + 25.6398487389914d0*z)) + &
        y*(-10.0227370861875d0 - 460.319931801257d0*z + y*(0.875600661808945d0 + 234.565187611355d0*z))))) + &
        y*(168.072408311545d0 + z*(729.116529735046d0 + &
        z*(-343.956902961561d0 + z*(124.687671116248d0 + z*(-31.656964386073d0 + 7.04658803315449d0*z)))) + &
        y*(880.031352997204d0 + y*(-225.267649263401d0 + &
        y*(91.4260447751259d0 + y*(-21.6603240875311d0 + 2.13016970847183d0*y) + &
        z*(-297.728741987187d0 + (74.726141138756d0 - 36.4872919001588d0*z)*z)) + &
        z*(694.244814133268d0 + z*(-204.889641964903d0 + (113.561697840594d0 - 11.1282734326413d0*z)*z))) + &
        z*(-860.764303783977d0 + z*(337.409530269367d0 + &
        z*(-178.314556207638d0 + (44.2040358308d0 - 7.92001547211682d0*z)*z))))))
        
  if(sa.gt.0.d0) then
      g08 = g08 + x2*(5812.81456626732d0 + 851.226734946706d0*y)*log(x)
  endif

  gsw_gibbs = g03 + g08
  
elseif(ns.eq.1 .and. nt.eq.0 .and. np.eq.0) then
        
  g08 = 8645.36753595126d0 + z*(-6620.98308089678d0 + &
       z*(769.588305957198d0 + z*(-193.0648640214916d0 + (31.6816345533648d0 - 5.24960313181984d0*z)*z))) + &
       x*(-7296.43987145382d0 + x*(8103.20462414788d0 + &
       y*(2175.341332000392d0 + y*(-274.2290036817964d0 + &
       y*(197.4670779425016d0 + y*(-68.5590309679152d0 + 9.98788038278032d0*y))) - 90.6734234051316d0*z) + &
       x*(-5458.34205214835d0 - 980.14153344888d0*y + &
       x*(2247.60742726704d0 - 340.1237483177863d0*x + 220.542973797483d0*y) + 180.142097805543d0*z) + &
       z*(-219.1676534131548d0 + (-16.32775915649044d0 - 120.7020447884644d0*z)*z)) + &
       z*(598.378809221703d0 + z*(-156.8822727844005d0 + (204.1334828179377d0 - 10.23755797323846d0*z)*z)) + &
       y*(-1480.222530425046d0 + z*(-525.876123559641d0 + (249.57717834054571d0 - 88.449193048287d0*z)*z) + &
       y*(-129.1994027934126d0 + z*(1149.174198007428d0 + z*(-162.5751787551336d0 + 76.9195462169742d0*z)) + &
       y*(-30.0682112585625d0 - 1380.9597954037708d0*z + y*(2.626801985426835d0 + 703.695562834065d0*z))))) + &
       y*(1187.3715515697959d0 + z*(1458.233059470092d0 + &
       z*(-687.913805923122d0 + z*(249.375342232496d0 + z*(-63.313928772146d0 + 14.09317606630898d0*z)))) + &
       y*(1760.062705994408d0 + y*(-450.535298526802d0 + &
       y*(182.8520895502518d0 + y*(-43.3206481750622d0 + 4.26033941694366d0*y) + &
       z*(-595.457483974374d0 + (149.452282277512d0 - 72.9745838003176d0*z)*z)) + &
       z*(1388.489628266536d0 + z*(-409.779283929806d0 + (227.123395681188d0 - 22.2565468652826d0*z)*z))) + &
       z*(-1721.528607567954d0 + z*(674.819060538734d0 + &
       z*(-356.629112415276d0 + (88.4080716616d0 - 15.84003094423364d0*z)*z)))))
  
  if(sa.gt.0.d0) then
    g08 = g08 + (11625.62913253464d0 + 1702.453469893412d0*y)*log(x)
  else
    g08 = 0.d0
  endif
  
  gsw_gibbs = 0.5*sfac*g08

elseif(ns.eq.0 .and. nt.eq.1 .and. np.eq.0) then
               
  g03 = 5.90578347909402d0 + z*(-270.983805184062d0 + &
       z*(776.153611613101d0 + z*(-196.51255088122d0 + (28.9796526294175d0 - 2.13290083518327d0*z)*z))) + &
       y*(-24715.571866078d0 + z*(2910.0729080936d0 + &
       z*(-1513.116771538718d0 + z*(546.959324647056d0 + z*(-111.1208127634436d0 + 8.68841343834394d0*z)))) + &
       y*(2210.2236124548363d0 + z*(-2017.52334943521d0 + &
       z*(1498.081172457456d0 + z*(-718.6359919632359d0 + (146.4037555781616d0 - 4.9892131862671505d0*z)*z))) + &
       y*(-592.743745734632d0 + z*(1591.873781627888d0 + &
       z*(-1207.261522487504d0 + (608.785486935364d0 - 105.4993508931208d0*z)*z)) + &
       y*(290.12956292128547d0 + z*(-973.091553087975d0 + &
       z*(602.603274510125d0 + z*(-276.361526170076d0 + 32.40953340386105d0*z))) + &
       y*(-113.90630790850321d0 + y*(21.35571525415769d0 - 67.41756835751434d0*z) + &
       z*(381.06836198507096d0 + z*(-133.7383902842754d0 + 49.023632509086724d0*z)))))))
              
  g08 = x2*(168.072408311545d0 + z*(729.116529735046d0 + &
        z*(-343.956902961561d0 + z*(124.687671116248d0 + z*(-31.656964386073d0 + 7.04658803315449d0*z)))) + &
        x*(-493.407510141682d0 + x*(543.835333000098d0 + x*(-196.028306689776d0 + 36.7571622995805d0*x) + &
        y*(-137.1145018408982d0 + y*(148.10030845687618d0 + y*(-68.5590309679152d0 + 12.4848504784754d0*y))) - &
        22.6683558512829d0*z) + z*(-175.292041186547d0 + (83.1923927801819d0 - 29.483064349429d0*z)*z) + &
        y*(-86.1329351956084d0 + z*(766.116132004952d0 + z*(-108.3834525034224d0 + 51.2796974779828d0*z)) + &
        y*(-30.0682112585625d0 - 1380.9597954037708d0*z + y*(3.50240264723578d0 + 938.26075044542d0*z)))) + &
        y*(1760.062705994408d0 + y*(-675.802947790203d0 + &
        y*(365.7041791005036d0 + y*(-108.30162043765552d0 + 12.78101825083098d0*y) + &
        z*(-1190.914967948748d0 + (298.904564555024d0 - 145.9491676006352d0*z)*z)) + &
        z*(2082.7344423998043d0 + z*(-614.668925894709d0 + (340.685093521782d0 - 33.3848202979239d0*z)*z))) + &
        z*(-1721.528607567954d0 + z*(674.819060538734d0 + &
        z*(-356.629112415276d0 + (88.4080716616d0 - 15.84003094423364d0*z)*z)))))
      
  if(sa.gt.0.d0) then
    g08 = g08 + 851.226734946706d0*x2*log(x)
  end if
  
  gsw_gibbs = (g03 + g08)*0.025d0

elseif(ns.eq.0 .and. nt.eq.0 .and. np.eq.1) then
    
  g03 = 100015.695367145d0 + z*(-5089.1530840726d0 + &
        z*(853.5533353388611d0 + z*(-133.2587017014444d0 + (21.0131554401542d0 - 3.278571068826234d0*z)*z))) + &
        y*(-270.983805184062d0 + z*(1552.307223226202d0 + &
        z*(-589.53765264366d0 + (115.91861051767d0 - 10.664504175916349d0*z)*z)) + &
        y*(1455.0364540468d0 + z*(-1513.116771538718d0 + &
        z*(820.438986970584d0 + z*(-222.2416255268872d0 + 21.72103359585985d0*z))) + &
        y*(-672.50778314507d0 + z*(998.720781638304d0 + &
        z*(-718.6359919632359d0 + (195.2050074375488d0 - 8.31535531044525d0*z)*z)) + &
        y*(397.968445406972d0 + z*(-603.630761243752d0 + (456.589115201523d0 - 105.4993508931208d0*z)*z) + &
        y*(-194.618310617595d0 + y*(63.5113936641785d0 - 9.63108119393062d0*y + &
        z*(-44.5794634280918d0 + 24.511816254543362d0*z)) + &
        z*(241.04130980405d0 + z*(-165.8169157020456d0 + &
        25.92762672308884d0*z)))))))                                                           
  
  g08 = x2*(-3310.49154044839d0 + z*(769.588305957198d0 + &
        z*(-289.5972960322374d0 + (63.3632691067296d0 - 13.1240078295496d0*z)*z)) + &
        x*(199.459603073901d0 + x*(-54.7919133532887d0 + 36.0284195611086d0*x - 22.6683558512829d0*y + &
        (-8.16387957824522d0 - 90.52653359134831d0*z)*z) + &
        z*(-104.588181856267d0 + (204.1334828179377d0 - 13.65007729765128d0*z)*z) + &
        y*(-175.292041186547d0 + (166.3847855603638d0 - 88.449193048287d0*z)*z + &
        y*(383.058066002476d0 + y*(-460.319931801257d0 + 234.565187611355d0*y) + &
        z*(-108.3834525034224d0 + 76.9195462169742d0*z)))) + &
        y*(729.116529735046d0 + z*(-687.913805923122d0 + &
        z*(374.063013348744d0 + z*(-126.627857544292d0 + 35.23294016577245d0*z))) + &
        y*(-860.764303783977d0 + y*(694.244814133268d0 + &
        y*(-297.728741987187d0 + (149.452282277512d0 - 109.46187570047641d0*z)*z) + &
        z*(-409.779283929806d0 + (340.685093521782d0 - 44.5130937305652d0*z)*z)) + &
        z*(674.819060538734d0 + z*(-534.943668622914d0 + (176.8161433232d0 - 39.600077360584095d0*z)*z)))))
     
  gsw_gibbs = (g03 + g08)*1d-8

elseif(ns.eq.0 .and. nt.eq.2 .and. np.eq.0) then

  g03 = -24715.571866078d0 + z*(2910.0729080936d0 + z* &
       (-1513.116771538718d0 + z*(546.959324647056d0 + z*(-111.1208127634436d0 + 8.68841343834394d0*z)))) + &
       y*(4420.4472249096725d0 + z*(-4035.04669887042d0 + &
       z*(2996.162344914912d0 + z*(-1437.2719839264719d0 + (292.8075111563232d0 - 9.978426372534301d0*z)*z))) + &
       y*(-1778.231237203896d0 + z*(4775.621344883664d0 + &
       z*(-3621.784567462512d0 + (1826.356460806092d0 - 316.49805267936244d0*z)*z)) + &
       y*(1160.5182516851419d0 + z*(-3892.3662123519d0 + &
       z*(2410.4130980405d0 + z*(-1105.446104680304d0 + 129.6381336154442d0*z))) + &
       y*(-569.531539542516d0 + y*(128.13429152494615d0 - 404.50541014508605d0*z) + &
       z*(1905.341809925355d0 + z*(-668.691951421377d0 + 245.11816254543362d0*z))))))

  g08 = x2*(1760.062705994408d0 + x*(-86.1329351956084d0 + &
        x*(-137.1145018408982d0 + y*(296.20061691375236d0 + y*(-205.67709290374563d0 + 49.9394019139016d0*y))) + &
        z*(766.116132004952d0 + z*(-108.3834525034224d0 + 51.2796974779828d0*z)) + &
        y*(-60.136422517125d0 - 2761.9195908075417d0*z + y*(10.50720794170734d0 + 2814.78225133626d0*z))) + &
        y*(-1351.605895580406d0 + y*(1097.1125373015109d0 + y*(-433.20648175062206d0 + 63.905091254154904d0*y) + &
        z*(-3572.7449038462437d0 + (896.713693665072d0 - 437.84750280190565d0*z)*z)) + &
        z*(4165.4688847996085d0 + z*(-1229.337851789418d0 + (681.370187043564d0 - 66.7696405958478d0*z)*z))) + &
        z*(-1721.528607567954d0 + z*(674.819060538734d0 + &
        z*(-356.629112415276d0 + (88.4080716616d0 - 15.84003094423364d0*z)*z))))
     
  gsw_gibbs = (g03 + g08)*0.000625d0  

elseif(ns.eq.1 .and. nt.eq.0 .and. np.eq.1) then

  g08 = -6620.98308089678d0 + z*(1539.176611914396d0 + &
        z*(-579.1945920644748d0 + (126.7265382134592d0 - 26.2480156590992d0*z)*z)) + &
        x*(598.378809221703d0 + x*(-219.1676534131548d0 + 180.142097805543d0*x - 90.6734234051316d0*y + &
        (-32.65551831298088d0 - 362.10613436539325d0*z)*z) + &
        z*(-313.764545568801d0 + (612.4004484538132d0 - 40.95023189295384d0*z)*z) + &
        y*(-525.876123559641d0 + (499.15435668109143d0 - 265.347579144861d0*z)*z + &
        y*(1149.174198007428d0 + y*(-1380.9597954037708d0 + 703.695562834065d0*y) + &
        z*(-325.1503575102672d0 + 230.7586386509226d0*z)))) + &
        y*(1458.233059470092d0 + z*(-1375.827611846244d0 + &
        z*(748.126026697488d0 + z*(-253.255715088584d0 + 70.4658803315449d0*z))) + &
        y*(-1721.528607567954d0 + y*(1388.489628266536d0 + &
        y*(-595.457483974374d0 + (298.904564555024d0 - 218.92375140095282d0*z)*z) + &
        z*(-819.558567859612d0 + (681.370187043564d0 - 89.0261874611304d0*z)*z)) + &
        z*(1349.638121077468d0 + z*(-1069.887337245828d0 + (353.6322866464d0 - 79.20015472116819d0*z)*z))));    
                                                          
  g08 = g08
                                                                   
  gsw_gibbs = g08*sfac*0.5d-8
         
elseif(ns.eq.0 .and. nt.eq.1 .and. np.eq.1) then

  g03 = -270.983805184062d0 + z*(1552.307223226202d0 + z*(-589.53765264366d0 + (115.91861051767d0 - 10.664504175916349d0*z)*z)) + &
        y*(2910.0729080936d0 + z*(-3026.233543077436d0 + &
        z*(1640.877973941168d0 + z*(-444.4832510537744d0 + 43.4420671917197d0*z))) + &
        y*(-2017.52334943521d0 + z*(2996.162344914912d0 + &
        z*(-2155.907975889708d0 + (585.6150223126464d0 - 24.946065931335752d0*z)*z)) + &
        y*(1591.873781627888d0 + z*(-2414.523044975008d0 + (1826.356460806092d0 - 421.9974035724832d0*z)*z) + &
        y*(-973.091553087975d0 + z*(1205.20654902025d0 + z*(-829.084578510228d0 + 129.6381336154442d0*z)) + &
        y*(381.06836198507096d0 - 67.41756835751434d0*y + z*(-267.4767805685508d0 + 147.07089752726017d0*z))))))
    
  g08 = x2*(729.116529735046d0 + z*(-687.913805923122d0 + &
        z*(374.063013348744d0 + z*(-126.627857544292d0 + 35.23294016577245d0*z))) + &
        x*(-175.292041186547d0 - 22.6683558512829d0*x + (166.3847855603638d0 - 88.449193048287d0*z)*z + &
        y*(766.116132004952d0 + y*(-1380.9597954037708d0 + 938.26075044542d0*y) + &
        z*(-216.7669050068448d0 + 153.8390924339484d0*z))) + &
        y*(-1721.528607567954d0 + y*(2082.7344423998043d0 + &
        y*(-1190.914967948748d0 + (597.809129110048d0 - 437.84750280190565d0*z)*z) + &
        z*(-1229.337851789418d0 + (1022.055280565346d0 - 133.5392811916956d0*z)*z)) + &
        z*(1349.638121077468d0 + z*(-1069.887337245828d0 + (353.6322866464d0 - 79.20015472116819d0*z)*z))))
    
  gsw_gibbs = (g03 + g08)*2.5d-10

elseif(ns.eq.0 .and. nt.eq.0 .and. np.eq.2) then
           
  g03 = -5089.1530840726d0 + z*(1707.1066706777221d0 + &
      z*(-399.7761051043332d0 + (84.0526217606168d0 - 16.39285534413117d0*z)*z)) + &
     y*(1552.307223226202d0 + z*(-1179.07530528732d0 + (347.75583155301d0 - 42.658016703665396d0*z)*z) + &
      y*(-1513.116771538718d0 + z*(1640.877973941168d0 + z*(-666.7248765806615d0 + 86.8841343834394d0*z)) + &
        y*(998.720781638304d0 + z*(-1437.2719839264719d0 + (585.6150223126464d0 - 33.261421241781d0*z)*z) + &
         y*(-603.630761243752d0 + (913.178230403046d0 - 316.49805267936244d0*z)*z + &
           y*(241.04130980405d0 + y*(-44.5794634280918d0 + 49.023632509086724d0*z) + &
            z*(-331.6338314040912d0 + 77.78288016926652d0*z))))))
            
  g08 = x2*(769.588305957198d0 + z*(-579.1945920644748d0 + (190.08980732018878d0 - 52.4960313181984d0*z)*z) + &
      x*(-104.588181856267d0 + x*(-8.16387957824522d0 - 181.05306718269662d0*z) + &
       (408.2669656358754d0 - 40.95023189295384d0*z)*z + &
       y*(166.3847855603638d0 - 176.898386096574d0*z + y*(-108.3834525034224d0 + 153.8390924339484d0*z))) + &
      y*(-687.913805923122d0 + z*(748.126026697488d0 + z*(-379.883572632876d0 + 140.9317606630898d0*z)) + &
       y*(674.819060538734d0 + z*(-1069.887337245828d0 + (530.4484299696d0 - 158.40030944233638d0*z)*z) + &
         y*(-409.779283929806d0 + y*(149.452282277512d0 - 218.92375140095282d0*z) + &
          (681.370187043564d0 - 133.5392811916956d0*z)*z))))
    
  gsw_gibbs = (g03 + g08)*1d-16 

end if

return
end function

!--------------------------------------------------------------------------

!==========================================================================
function gsw_saar(p,long,lat)
!==========================================================================

! Calculates the Absolute Salinity Anomaly Ratio, SAAR.
!
! p      : sea pressure                                    [dbar]
! long   : longitude                                       [deg E]     
! lat    : latitude                                        [deg N]
!
! gsw_saar : Absolute Salinity Anomaly Ratio               [unitless]

implicit none

!integer, parameter :: int9 = selected_int_kind(9) 
integer, parameter :: r14 = selected_real_kind(14,30)

integer, parameter :: nx=91, ny=45, nz=45

integer :: indx0, indy0, indz0, i, j, icalled, k
integer :: nmean, flag_saar
integer, dimension(4) :: deli, delj

real (r14), dimension(4) :: saar, saar_old
real (r14), dimension(nz) :: p_ref
real (r14), dimension(ny) :: lats_ref
real (r14), dimension(nx) :: longs_ref
real (r14), dimension(ny,nx) :: ndepth_ref 
real (r14), dimension(nz,ny,nx) :: saar_ref
!real (r14), dimension(nz,ny,nx) :: delta_sa_ref
real (r14) :: p, lat, long, dlong, dlat
real (r14) :: gsw_saar, p0_original, lon0_in, sa_upper, sa_lower 
real (r14) :: r1, s1, t1, saar_mean, ndepth_max

data deli/0,1,1,0/, delj/0,0,1,1/

data icalled/0/, flag_saar/0/

save icalled, flag_saar, longs_ref, lats_ref, p_ref, ndepth_ref, saar_ref

gsw_saar = 9d15

if(lat .lt. -86d0 .or. lat .gt. 90d0) then
 gsw_saar = 9d15
 return
end if

if(long .lt. 0) then
 long = long + 360
end if

if(icalled.eq.0) then
  icalled = 1
   open(10,file='gsw_data_v3_0.dat',status='old',err=1)
   flag_saar = 1
   read(10,*) (longs_ref(i), i=1,nx)
   read(10,*) (lats_ref(i), i=1,ny)
   read(10,*) (p_ref(i), i=1,nz)
   read(10,*) ((ndepth_ref(j,i), j=1,ny), i=1,nx)
   read(10,*) (((saar_ref(k,j,i), k=1,nz), j=1,ny), i=1,nx)
   !read(10,*) (((delta_sa_ref(k,j,i), k=1,nz), j=1,ny), i=1,nx)
   close(10)
   go to 2
1  saar_ref = 9d15
   flag_saar = 0
2  continue
end if

if (flag_saar.eq.0) then
   write(*,*) "*** gsw_data_v3_0.dat is missing !!! ***"
   write(*,*) "Set the full path of gsw_data_v3_0.dat in gsw_saar"
end if

!Set gsw_saar = 9d15 and return if there is no data file present
if(flag_saar .eq. 0) then
 gsw_saar = 9d15
 return
endif

dlong = longs_ref(2)-longs_ref(1)
dlat = lats_ref(2)-lats_ref(1)

indx0 = floor(1 + (nx-1)*(long-longs_ref(1))/(longs_ref(nx)-longs_ref(1)))
if(indx0.eq.nx) then
   indx0 = nx-1
end if

indy0 = floor(1 + (ny-1)*(lat-lats_ref(1))/(lats_ref(ny)-lats_ref(1)))
if(indy0.eq.ny) then
   indy0 = ny-1
end if

ndepth_max = -1
do k = 1,4
   if(ndepth_ref(indy0+delj(k),indx0+deli(k)).gt.0.d0) then
      ndepth_max = max(ndepth_max,ndepth_ref(indy0+delj(k),indx0+deli(k)))
   end if
end do

if(ndepth_max.eq.-1.d0) then
  gsw_saar = 0d0 
   return
end if 

p0_original = p
if(p.gt.p_ref(int(ndepth_max))) then
 p = p_ref(int(ndepth_max))
end if
call indx(p_ref,nz,p,indz0)
    
r1 = (long-longs_ref(indx0))/(longs_ref(indx0+1)-longs_ref(indx0));
s1 = (lat-lats_ref(indy0))/(lats_ref(indy0+1)-lats_ref(indy0));
t1 = (p-p_ref(indz0))/(p_ref(indz0+1)-p_ref(indz0));

do k = 1,4
   saar(k) = saar_ref(indz0,indy0+delj(k),indx0+deli(k))
end do

if(260.d0.le.long.and.long.le.291.999d0.and.3.4d0.le.lat.and.lat.le.19.55d0) then
  saar_old = saar
  call gsw_add_barrier(saar_old,long,lat,longs_ref(indx0),lats_ref(indy0),dlong,dlat,saar)
else if(abs(sum(saar)) .ge. 1d10) then
  saar_old = saar
  call gsw_add_mean(saar_old,long,lat,saar)
end if

sa_upper = (1.d0-s1)*(saar(1) + r1*(saar(2)-saar(1))) + s1*(saar(4) + r1*(saar(3)-saar(4)))

do k = 1,4
   saar(k) = saar_ref(indz0+1,indy0+delj(k),indx0+deli(k))
end do

if(260.d0.le.long.and.long.le.291.999d0.and.3.4d0.le.lat.and.lat.le.19.55d0) then
   saar_old = saar
   call gsw_add_barrier(saar_old,long,lat,longs_ref(indx0),lats_ref(indy0),dlong,dlat,saar)
else if(abs(sum(saar)) .ge. 1d10) then 
   saar_old = saar
   call gsw_add_mean(saar_old,long,lat,saar)
end if

sa_lower = (1.d0-s1)*(saar(1) + r1*(saar(2)-saar(1))) + s1*(saar(4) + r1*(saar(3)-saar(4)))
if(abs(sa_lower) .ge. 1d10) then
  sa_lower = sa_upper
end if
gsw_saar = sa_upper + t1*(sa_lower-sa_upper)

if(abs(gsw_saar).ge.1d10) then
   gsw_saar = 9d15
endif

p = p0_original
  
return
end function

!--------------------------------------------------------------------------

!==========================================================================
function gsw_delta_sa_ref(p,long,lat)
!==========================================================================

! Calculates the Absolute Salinity Anomaly reference value, delta_SA_ref.
!
! p      : sea pressure                                    [dbar]
! long   : longiture                                       [deg E]     
! lat    : latitude                                        [deg N]
!
! gsw_delta_sa_ref : Absolute Salinity Anomaly reference value    [g/kg]

implicit none

!integer, parameter :: int9 = selected_int_kind(9) 
integer, parameter :: r14 = selected_real_kind(14,30)

integer, parameter :: nx=91, ny=45, nz=45

integer :: indx0, indy0, indz0, i, j, icalled2, k
integer :: nmean, flag_dsar
integer, dimension(4) :: deli, delj

real (r14), dimension(4) :: dsar, dsar_old
real (r14), dimension(nz) :: p_ref
real (r14), dimension(ny) :: lats_ref
real (r14), dimension(nx) :: longs_ref
real (r14), dimension(ny,nx) :: ndepth_ref 
real (r14), dimension(nz,ny,nx) :: saar_ref, delta_sa_ref
real (r14) :: p, lat, long, dlong, dlat
real (r14) :: gsw_delta_sa_ref, p0_original, lon0_in, sa_upper, sa_lower 
real (r14) :: r1, s1, t1, dsar_mean, ndepth_max

data deli/0,1,1,0/, delj/0,0,1,1/

data icalled2/0/, flag_dsar/0/

save icalled2, flag_dsar, longs_ref, lats_ref, p_ref, ndepth_ref, delta_sa_ref

gsw_delta_sa_ref = 9d15

if(lat .lt. -86d0 .or. lat .gt. 90d0) then
 gsw_delta_sa_ref = 9d15
 return
end if

if(long.lt.0) then
 long = long + 360
end if

if(icalled2.eq.0) then
   icalled2 = 1
   open(10,file='gsw_data_v3_0.dat',status='old',err=1)
   flag_dsar = 1
   read(10,*) (longs_ref(i), i=1,nx)
   read(10,*) (lats_ref(i), i=1,ny)
   read(10,*) (p_ref(i), i=1,nz)
   read(10,*) ((ndepth_ref(j,i), j=1,ny), i=1,nx)
   read(10,*) (((saar_ref(k,j,i), k=1,nz), j=1,ny), i=1,nx)
   read(10,*) (((delta_sa_ref(k,j,i), k=1,nz), j=1,ny), i=1,nx)
   close(10)
   go to 2
1  delta_sa_ref = 9d15
   flag_dsar = 0
2  continue
end if

if (flag_dsar.eq.0) then
   write(*,*) "*** gsw_data_v3_0.dat is missing !!! ***"
   write(*,*) "Set the full path of gsw_data_v3_0.dat in gsw_delta_sa_ref"
end if

!Set gsw_delta_sa_ref = 9d15 and return if there is no data set present
if(flag_dsar.eq.0) then
 gsw_delta_sa_ref = 9d15
 return
endif

dlong = longs_ref(2)-longs_ref(1)
dlat = lats_ref(2)-lats_ref(1)

indx0 = floor(1 + (nx-1)*(long-longs_ref(1))/(longs_ref(nx)-longs_ref(1)))
if(indx0.eq.nx) then
   indx0 = nx-1
end if

indy0 = floor(1 + (ny-1)*(lat-lats_ref(1))/(lats_ref(ny)-lats_ref(1)))
if(indy0.eq.ny) then
   indy0 = ny-1
end if

ndepth_max = -1
do k = 1,4
   if(ndepth_ref(indy0+delj(k),indx0+deli(k)).gt.0.d0) then
      ndepth_max = max(ndepth_max,ndepth_ref(indy0+delj(k),indx0+deli(k)))
   end if
end do

if(ndepth_max.eq.-1.d0) then
  gsw_delta_sa_ref = 0d0 
   return
end if 

p0_original = p
if(p.gt.p_ref(int(ndepth_max))) then
 p = p_ref(int(ndepth_max))
end if
call indx(p_ref,nz,p,indz0)
    
r1 = (long-longs_ref(indx0))/(longs_ref(indx0+1)-longs_ref(indx0));
s1 = (lat-lats_ref(indy0))/(lats_ref(indy0+1)-lats_ref(indy0));
t1 = (p-p_ref(indz0))/(p_ref(indz0+1)-p_ref(indz0));

do k = 1,4
   dsar(k) = delta_sa_ref(indz0,indy0+delj(k),indx0+deli(k))
end do

if(260.d0.le.long.and.long.le.291.999d0.and.3.4d0.le.lat.and.lat.le.19.55d0) then
   dsar_old = dsar
   call gsw_add_barrier(dsar_old,long,lat,longs_ref(indx0),lats_ref(indy0),dlong,dlat,dsar)
else if(abs(sum(dsar)).ge.1d10) then 
   dsar_old = dsar
   call gsw_add_mean(dsar_old,long,lat,dsar)
end if

sa_upper = (1.d0-s1)*(dsar(1) + r1*(dsar(2)-dsar(1))) + s1*(dsar(4) + r1*(dsar(3)-dsar(4)))

do k = 1,4
   dsar(k) = delta_sa_ref(indz0+1,indy0+delj(k),indx0+deli(k))
end do

if(260.d0.le.long.and.long.le.291.999d0.and.3.4d0.le.lat.and.lat.le.19.55d0) then
   dsar_old = dsar
   call gsw_add_barrier(dsar_old,long,lat,longs_ref(indx0),lats_ref(indy0),dlong,dlat,dsar)
else if(abs(sum(dsar)).ge.1d10) then 
   dsar_old = dsar
   call gsw_add_mean(dsar_old,long,lat,dsar)
end if

sa_lower = (1.d0-s1)*(dsar(1) + r1*(dsar(2)-dsar(1))) + s1*(dsar(4) + r1*(dsar(3)-dsar(4)))
if(abs(sa_lower).ge.1d10) then
  sa_lower = sa_upper
end if
gsw_delta_sa_ref = sa_upper + t1*(sa_lower-sa_upper)

if(abs(gsw_delta_sa_ref).ge.1d10) then
   gsw_delta_sa_ref = 9d15
endif

p = p0_original
  
return
end function

!--------------------------------------------------------------------------

!==========================================================================
subroutine gsw_add_barrier(input_data,long,lat,long_grid,lat_grid,dlong_grid,dlat_grid,output_data)
!==========================================================================

!  Adds a barrier through Central America (Panama) and then averages
!  over the appropriate side of the barrier
! 
!  data_in      :  data                                                     [unitless]
!  long         :  Longitudes of data in decimal degrees east               [ 0 ... +360 ]
!  lat          :  Latitudes of data in decimal degrees north               [ -90 ... +90 ]
!  longs_grid   :  Longitudes of regular grid in decimal degrees east       [ 0 ... +360 ]
!  lats_grid    :  Latitudes of regular grid in decimal degrees north       [ -90 ... +90 ]
!  dlongs_grid  :  Longitude difference of regular grid in decimal degrees  [ deg longitude ]
!  dlats_grid   :  Latitude difference of regular grid in decimal degrees   [ deg latitude ]
!
! gsw_add_barrier  : average of data depending on which side of the 
!                    Panama cannal it is on                                 [unitless]

implicit none

integer, parameter :: r14 = selected_real_kind(14,30)

integer, dimension(4) :: above_line
integer k, nmean, above_line0, kk
real (r14), dimension(4) :: input_data, output_data
real (r14), dimension(6) :: longs_pan, lats_pan
real (r14) :: long, lat, r, lats_line, long_grid, lat_grid
real (r14) :: dlong_grid, dlat_grid, data_mean

data longs_pan/260.0000d0, 272.5900d0, 276.5000d0, 278.6500d0, 280.7300d0, 292.000d0/ 
data  lats_pan/ 19.5500d0,  13.9700d0,   9.6000d0,   8.1000d0,   9.3300d0,   3.400d0/ 

call indx(longs_pan,6,long,k)                            !   the long/lat point
r = (long-longs_pan(k))/(longs_pan(k+1)-longs_pan(k))
lats_line = lats_pan(k) + r*(lats_pan(k+1)-lats_pan(k))

if(lats_line.le.lat) then
   above_line0 = 1
else
   above_line0 = 0
end if

call indx(longs_pan,6,long_grid,k)                                     !  the 1 and 4 long/lat points 
r = (long_grid-longs_pan(k))/(longs_pan(k+1)-longs_pan(k))
lats_line = lats_pan(k) + r*(lats_pan(k+1)-lats_pan(k))

if(lats_line.le.lat_grid) then
   above_line(1) = 1
else
   above_line(1) = 0
end if

if(lats_line.le.lat_grid+dlat_grid) then
   above_line(4) = 1
else
   above_line(4) = 0
end if

call indx(longs_pan,6,long_grid+dlong_grid,k)                              !  the 2 and 3 long/lat points 
r = (long_grid+dlong_grid-longs_pan(k))/(longs_pan(k+1)-longs_pan(k))
lats_line = lats_pan(k) + r*(lats_pan(k+1)-lats_pan(k))

if(lats_line.le.lat_grid) then
   above_line(2) = 1
else
   above_line(2) = 0
end if

if(lats_line.le.lat_grid+dlat_grid) then
   above_line(3) = 1
else
   above_line(3) = 0
end if

nmean = 0 
data_mean = 0.d0

do kk = 1,4
   if ((abs(input_data(kk)).le.100d0).and.above_line0.eq.above_line(kk)) then
      nmean = nmean+1
      data_mean = data_mean+input_data(kk)
   end if
end do

if(nmean .eq. 0)then
   data_mean = 0d0    !errorreturn
else
   data_mean = data_mean/nmean
endif

do kk = 1,4
   if((abs(input_data(kk)).ge.1d10).or.above_line0.ne.above_line(kk)) then
      output_data(kk) = data_mean
   else
      output_data(kk) = input_data(kk)
   end if
end do

return
end subroutine

!--------------------------------------------------------------------------

!==========================================================================
subroutine gsw_add_mean(data_in,long,lat,data_out)
!==========================================================================

! Replaces NaN's with non-nan mean of the 4 adjacent neighbours
!
! data_in   : data set of the 4 adjacent neighbours   
! long      : longitude
! lat       : latitude
!
! data_out : non-nan mean of the 4 adjacent neighbours     [unitless]

implicit none

integer, parameter :: int9 = selected_int_kind(9)
integer, parameter :: r14 = selected_real_kind(14,30)

integer :: k, nmean

real (r14), dimension(4) :: data_in, data_out 
real (r14) :: data_mean, long, lat

nmean = 0
data_mean = 0.d0

do k = 1,4
   if (abs(data_in(k)).le.100d0) then
      nmean = nmean+1
      data_mean = data_mean+data_in(k)
   end if
end do

if(nmean.eq.0)then
   data_mean = 0d0    !errorreturn
else
   data_mean = data_mean/nmean
endif

do k = 1,4
   if(abs(data_in(k)).ge.100d0) then
      data_out(k) = data_mean
   else
      data_out(k) = data_in(k)
   end if
end do

return
end subroutine

!--------------------------------------------------------------------------

!==========================================================================
function xinterp1(x,y,n,x0)
!==========================================================================

! Linearly interpolate a real array   
!
! x      : y array (Must be monotonic)               
! y      : y array     
! n      : length of X and Y arrays
! x0     : value to be interpolated
!
! xinterp1 : Linearly interpolated value

implicit none

integer, parameter :: r14 = selected_real_kind(14,30)

integer :: n, k

real (r14), dimension(n) :: x, y
real (r14) :: x0, r, xinterp1

call indx(x,n,x0,k)
r = (x0-x(k))/(x(k+1)-x(k))
xinterp1 = y(k) + r*(y(k+1)-y(k))

return
end function

!--------------------------------------------------------------------------

!==========================================================================
subroutine indx(x,n,z,k)
!==========================================================================

!  Finds the index of the value in a monotonically increasing array
!
!  x	 :  array of monotonically increasing values
!  n     :  length of the array
!  z     :  value to be indexed
!
!  K      : index K - if X(K) <= Z < X(K+1), or
!  N-1     		    - if Z = X(N)
!

integer, parameter :: r14 = selected_real_kind(14,30)

integer :: n, k, ku, kl, km

real (r14), dimension(n) :: x
real (r14) :: z

if(z.gt.x(1).and.z.lt.x(n)) then
   kl=1
   ku=n
   do while (ku-kl.gt.1)
      km=0.5d0*(ku+kl)
      if(z.gt.x(km)) then
         kl=km
      else
         ku=km
      endif
   end do
   k=kl
   if(z.eq.x(k+1)) then 
     k = k+1
   end if
elseif (z.le.x(1)) then
      k = 1
elseif (z.ge.x(n)) then
      k = n-1
else
      write(*,*) 'ERROR in subroutine indx : out of range'
      write(*,*) 'z = ', z, 'n = ', n, 'x = ',x
end if

return
end subroutine

!--------------------------------------------------------------------------

!==========================================================================
function gsw_fdelta(p,long,lat)
!==========================================================================

! Calculates fdelta. 
!
! p      : sea pressure                                    [dbar]
! long   : longitude                                       [deg E]     
! lat    : latitude                                        [deg N]
!
! gsw_fdelta : Absolute Salinty Anomaly                    [unitless]

implicit none

integer, parameter :: r14 = selected_real_kind(14,30)

real (r14) ::  long, lat, p, gsw_saar, saar, gsw_fdelta

saar = gsw_saar(p,long,lat)

gsw_fdelta = ((1d0 + 0.35d0)*saar)/(1d0 - 0.35d0*saar);

if (saar.gt.1d10) then
    gsw_fdelta = 9d15
end if

return
end function

!--------------------------------------------------------------------------

!==========================================================================
function gsw_sa_from_sp_baltic(sp,long,lat)
!==========================================================================

! For the Baltic Sea, calculates Absolute Salinity with a value
! computed analytically from Practical Salinity
!
! sp     : Practical Salinity                              [unitless]
! long   : longitude                                       [deg E]     
! lat    : latitude                                        [deg N]
! p      : sea pressure                                    [dbar]
!
! gsw_sa_from_sp_baltic : Absolute Salinity                [g/kg]                         

implicit none

integer, parameter :: r14 = selected_real_kind(14,30)

real (r14), dimension(2) :: xb_right, yb_right
real (r14), dimension(3) :: xb_left, yb_left
real (r14) :: sp, long, lat, gsw_sa_from_sp_baltic, xinterp1, xx_left, xx_right

data xb_left/12.6d0, 7.d0, 26.d0/, yb_left/50.d0, 59.d0, 69.d0/
data xb_right/45.d0, 26.d0/, yb_right/50.d0, 69.d0/

if(xb_left(2).lt.long .and. long.lt.xb_right(1) .and. yb_left(1).lt.lat .and. lat.lt.yb_left(3)) then
  
    xx_left = xinterp1(yb_left, xb_left, 3, lat)
    
    xx_right = xinterp1(yb_right, xb_right, 2, lat)
    
    if(xx_left.le.long .and. long.le.xx_right) then
        gsw_sa_from_sp_baltic =((35.16504d0 - 0.087d0)/35d0)*sp + 0.087d0
    else
        gsw_sa_from_sp_baltic = 9d15
    end if

else
    gsw_sa_from_sp_baltic = 9d15
end if

return
end

!--------------------------------------------------------------------------

!==========================================================================
function gsw_sp_from_sa_baltic(sa,long,lat)
!==========================================================================

! For the Baltic Sea, calculates Practical Salinity with a value
! computed analytically from Absolute Salinity
!
! sa     : Absolute Salinity                               [g/kg]
! long   : longitude                                       [deg E]     
! lat    : latitude                                        [deg N]
! p      : sea pressure                                    [dbar]
!
! gsw_sp_from_sa_baltic  : Practical Salinity              [unitless]

implicit none

integer, parameter :: r14 = selected_real_kind(14,30)

real (r14), dimension(2) :: xb_right, yb_right
real (r14), dimension(3) :: xb_left, yb_left
real (r14) :: sa, long, lat, gsw_sp_from_sa_baltic, xinterp1, xx_left, xx_right

data xb_left/12.6d0, 7.d0, 26.d0/, yb_left/50.d0, 59.d0, 69.d0/
data xb_right/45.d0, 26.d0/, yb_right/50.d0, 69.d0/

if(xb_left(2).lt.long .and. long.lt.xb_right(1) .and. yb_left(1).lt.lat .and. lat.lt.yb_left(3)) then
  
    xx_left = xinterp1(yb_left, xb_left, 3, lat)
    
    xx_right = xinterp1(yb_right, xb_right, 2, lat)
    
    if(xx_left.le.long .and. long.le.xx_right) then
        gsw_sp_from_sa_baltic = (35.d0/(35.16504d0 - 0.087d0))*(sa - 0.087d0)
    else
        gsw_sp_from_sa_baltic = 9d15
    end if
     
else
    gsw_sp_from_sa_baltic = 9d15
end if

return
end

!--------------------------------------------------------------------------

!==========================================================================
function gsw_entropy_part(sa,t,p)
!==========================================================================

! entropy minus the terms that are a function of only SA
!
! sa     : Absolute Salinity                               [g/kg]
! t      : in-situ temperature                             [deg C]
! p      : sea pressure                                    [dbar]
! 
! gsw_entropy_part : entropy part

implicit none

integer, parameter :: r14 = selected_real_kind(14,30)

real (r14) :: sa, t, p, sfac, x2, x, y, z, g03, g08, gsw_entropy_part

sfac = 0.0248826675584615d0

x2 = sfac*sa
x = sqrt(x2)
y = t*0.025d0
z = p*1d-4

g03 = z*(-270.983805184062d0 + &
    z*(776.153611613101d0 + z*(-196.51255088122d0 + (28.9796526294175d0 - 2.13290083518327d0*z)*z))) + &
    y*(-24715.571866078d0 + z*(2910.0729080936d0 + &
    z*(-1513.116771538718d0 + z*(546.959324647056d0 + z*(-111.1208127634436d0 + 8.68841343834394d0*z)))) + &
    y*(2210.2236124548363d0 + z*(-2017.52334943521d0 + &
    z*(1498.081172457456d0 + z*(-718.6359919632359d0 + (146.4037555781616d0 - 4.9892131862671505d0*z)*z))) + &
    y*(-592.743745734632d0 + z*(1591.873781627888d0 + &
    z*(-1207.261522487504d0 + (608.785486935364d0 - 105.4993508931208d0*z)*z)) + &
    y*(290.12956292128547d0 + z*(-973.091553087975d0 + &
    z*(602.603274510125d0 + z*(-276.361526170076d0 + 32.40953340386105d0*z))) + &
    y*(-113.90630790850321d0 + y*(21.35571525415769d0 - 67.41756835751434d0*z) + &
    z*(381.06836198507096d0 + z*(-133.7383902842754d0 + 49.023632509086724d0*z)))))))

g08 = x2*(z*(729.116529735046d0 + &
    z*(-343.956902961561d0 + z*(124.687671116248d0 + z*(-31.656964386073d0 + 7.04658803315449d0*z)))) + &
    x*( x*(y*(-137.1145018408982d0 + y*(148.10030845687618d0 + y*(-68.5590309679152d0 + 12.4848504784754d0*y))) - &
    22.6683558512829d0*z) + z*(-175.292041186547d0 + (83.1923927801819d0 - 29.483064349429d0*z)*z) + &
    y*(-86.1329351956084d0 + z*(766.116132004952d0 + z*(-108.3834525034224d0 + 51.2796974779828d0*z)) + &
    y*(-30.0682112585625d0 - 1380.9597954037708d0*z + y*(3.50240264723578d0 + 938.26075044542d0*z)))) + &
    y*(1760.062705994408d0 + y*(-675.802947790203d0 + &
    y*(365.7041791005036d0 + y*(-108.30162043765552d0 + 12.78101825083098d0*y) + &
    z*(-1190.914967948748d0 + (298.904564555024d0 - 145.9491676006352d0*z)*z)) + &
    z*(2082.7344423998043d0 + z*(-614.668925894709d0 + (340.685093521782d0 - 33.3848202979239d0*z)*z))) + &
    z*(-1721.528607567954d0 + z*(674.819060538734d0 + &
    z*(-356.629112415276d0 + (88.4080716616d0 - 15.84003094423364d0*z)*z)))))

gsw_entropy_part = -(g03 + g08)*0.025d0

return
end function

!--------------------------------------------------------------------------

!==========================================================================
function gsw_entropy_part_zerop(sa,pt0)
!==========================================================================

! entropy part evaluated at the sea surface
!
! sa     : Absolute Salinity                               [g/kg]
! pt0    : insitu temperature                              [deg C]
! 
! gsw_entropy_part_zerop : entropy part at the sea surface

implicit none

integer, parameter :: r14 = selected_real_kind(14,30)

real (r14) :: sa, pt0, sfac, x2, x, y, g03, g08, gsw_entropy_part_zerop

sfac = 0.0248826675584615d0

x2 = sfac*sa
x = sqrt(x2)
y = pt0*0.025d0

g03 = y*(-24715.571866078d0 + y*(2210.2236124548363d0 + &
    y*(-592.743745734632d0 + y*(290.12956292128547d0 + &
    y*(-113.90630790850321d0 + y*21.35571525415769d0)))))

g08 = x2*(x*(x*(y*(-137.1145018408982d0 + y*(148.10030845687618d0 + &
    y*(-68.5590309679152d0 + 12.4848504784754d0*y)))) + &
    y*(-86.1329351956084d0 + y*(-30.0682112585625d0 + y*3.50240264723578d0))) + &
    y*(1760.062705994408d0 + y*(-675.802947790203d0 + &
    y*(365.7041791005036d0 + y*(-108.30162043765552d0 + 12.78101825083098d0*y)))))

gsw_entropy_part_zerop = -(g03 + g08)*0.025d0

return
end function

!--------------------------------------------------------------------------

!==========================================================================
function gsw_gibbs_pt0_pt0(sa,pt0)
!==========================================================================

! gibbs_tt at (sa,pt,0)
!
! sa     : Absolute Salinity                            [g/kg]
! pt0    : potential temperature                        [deg C]
! 
! gsw_gibbs_pt0_pt0 : gibbs_tt at (sa,pt,0)         

implicit none

integer, parameter :: r14 = selected_real_kind(14,30)

real (r14) :: sa, pt0, sfac, x2, x, y, g03, g08, gsw_gibbs_pt0_pt0

sfac = 0.0248826675584615d0

x2 = sfac*sa
x = sqrt(x2)
y = pt0*0.025d0

g03 = -24715.571866078d0 + &
    y*(4420.4472249096725d0 + &
    y*(-1778.231237203896d0 + &
    y*(1160.5182516851419d0 + &
    y*(-569.531539542516d0 + y*128.13429152494615d0))))

g08 = x2*(1760.062705994408d0 + x*(-86.1329351956084d0 + &
    x*(-137.1145018408982d0 + y*(296.20061691375236d0 + &
    y*(-205.67709290374563d0 + 49.9394019139016d0*y))) + &
    y*(-60.136422517125d0 + y*10.50720794170734d0)) + &
    y*(-1351.605895580406d0 + y*(1097.1125373015109d0 + &
    y*(-433.20648175062206d0 + 63.905091254154904d0*y))))

gsw_gibbs_pt0_pt0 = (g03 + g08)*0.000625d0

return
end function

!--------------------------------------------------------------------------

!==========================================================================
function gsw_specvol_sso_0_p(p) 
!==========================================================================

!  This function calculates specifc volume at the Standard Ocean Salinty,
!  SSO, and at a Conservative Temperature of zero degrees C, as a function 
!  of pressure, p, in dbar, using a streamlined version of the 48-term CT
!  version of specific volume, that is, a streamlined version of the code
!  "gsw_specvol(SA,CT,p)".
!
! p      : sea pressure                                    [dbar]
! 
! gsw_specvol_sso_0_p : specvol(sso,0,p)

implicit none

integer, parameter :: r14 = selected_real_kind(14,30)

real (r14), parameter :: v01 =  9.998420897506056d2, v05 = -6.698001071123802d0
real (r14), parameter :: v08 = -3.988822378968490d-2, v12 = -2.233269627352527d-2
real (r14), parameter :: v15 = -1.806789763745328d-4, v17 = -3.087032500374211d-7
real (r14), parameter :: v20 =  1.550932729220080d-10, v21 =  1.0d0;
real (r14), parameter :: v26 = -7.521448093615448d-3, v31 = -3.303308871386421d-5
real (r14), parameter :: v36 =  5.419326551148740d-6, v37 = -2.742185394906099d-5
real (r14), parameter :: v41 = -1.105097577149576d-7, v43 = -1.119011592875110d-10
real (r14), parameter :: v47 = -1.200507748551599d-15

real (r14) :: sso, sqrtsso, p, gsw_specvol_sso_0_p

sso = 35.16504d0
sqrtsso = 5.930011804372737d0     ! sqrt(SSO) = 5.930011804372737

gsw_specvol_sso_0_p = (v21 + sso*(v26 + v36*sso + v31*sqrtsso)  &
             + p*(v37 + v41*sso + p*(v43 + v47*p )))/ &
             (v01 + sso*(v05 + v08*sqrtsso) &
             + p*(v12 + v15*sso + p*(v17 + v20*sso)))

return
end function

!--------------------------------------------------------------------------

!==========================================================================
function  gsw_hill_ratio_at_sp2(t)
!==========================================================================

!  Calculates the Hill ratio, which is the adjustment needed to apply for
!  Practical Salinities smaller than 2.  This ratio is defined at a 
!  Practical Salinity = 2 and in-situ temperature, t using PSS-78. The Hill
!  ratio is the ratio of 2 to the output of the Hill et al. (1986) formula
!  for Practical Salinity at the conductivity ratio, Rt, at which Practical
!  Salinity on the PSS-78 scale is exactly 2.

implicit none

integer, parameter :: r14 = selected_real_kind(14,30)

real (r14), parameter :: a0 =  0.0080d0, a1 = -0.1692d0, a2 = 25.3851d0
real (r14), parameter :: a3 = 14.0941d0, a4 = -7.0261d0, a5 =  2.7081d0
real (r14), parameter :: b0 =  0.0005d0, b1 = -0.0056d0, b2 = -0.0066d0
real (r14), parameter :: b3 = -0.0375d0, b4 =  0.0636d0, b5 = -0.0144d0
real (r14), parameter :: g0 = 2.641463563366498d-1, g1 = 2.007883247811176d-4
real (r14), parameter :: g2 = -4.107694432853053d-6, g3 = 8.401670882091225d-8
real (r14), parameter :: g4 = -1.711392021989210d-9, g5 = 3.374193893377380d-11
real (r14), parameter :: g6 = -5.923731174730784d-13, g7 = 8.057771569962299d-15
real (r14), parameter :: g8 = -7.054313817447962d-17, g9 = 2.859992717347235d-19
real (r14), parameter :: rk  =  0.0162d0, sp2 = 2d0

real (r14) :: t, t68, ft68, rtx0, dsp_drtx, sp_est, rtx, rtxm, x, part1, part2
real (r14) :: sqrty, sp_hill_raw_at_sp2, gsw_hill_ratio_at_sp2

t68 = t*1.00024d0
ft68 = (t68 - 15d0)/(1d0 + rk*(t68 - 15d0))

!--------------------------------------------------------------------------
! Find the initial estimates of Rtx (Rtx0) and of the derivative dSP_dRtx
! at SP = 2. 
!--------------------------------------------------------------------------
rtx0 = g0 + t68*(g1 + t68*(g2 + t68*(g3 + t68*(g4 + t68*(g5 &
         + t68*(g6 + t68*(g7 + t68*(g8 + t68*g9))))))))
     
dsp_drtx =  a1 + (2*a2 + (3*a3 + (4*a4 + 5*a5*rtx0)*rtx0)*rtx0)*rtx0 + &
    ft68*(b1 + (2*b2 + (3*b3 + (4*b4 + 5*b5*rtx0)*rtx0)*rtx0)*rtx0)    

!--------------------------------------------------------------------------
! Begin a single modified Newton-Raphson iteration to find Rt at SP = 2.
!--------------------------------------------------------------------------
sp_est = a0 + (a1 + (a2 + (a3 + (a4 + a5*rtx0)*rtx0)*rtx0)*rtx0)*rtx0 &
        + ft68*(b0 + (b1 + (b2+ (b3 + (b4 + b5*rtx0)*rtx0)*rtx0)*rtx0)*rtx0)
rtx = rtx0 - (sp_est - sp2)/dsp_drtx
rtxm = 0.5d0*(rtx + rtx0)
dsp_drtx =  a1 + (2*a2 + (3*a3 + (4*a4 + 5*a5*rtxm)*rtxm)*rtxm)*rtxm &
        + ft68*(b1 + (2*b2 + (3*b3 + (4*b4 + 5*b5*rtxm)*rtxm)*rtxm)*rtxm)
rtx = rtx0 - (sp_est - sp2)/dsp_drtx

! This is the end of one full iteration of the modified Newton-Raphson 
! iterative equation solver.  The error in Rtx at this point is equivalent 
! to an error in SP of 9e-16 psu.  
                                
x = 400d0*rtx*rtx
sqrty = 10d0*rtx
part1 = 1 + x*(1.5d0 + x) 
part2 = 1 + sqrty*(1 + sqrty*(1 + sqrty))
sp_hill_raw_at_sp2 = SP2 - a0/part1 - b0*ft68/part2

gsw_hill_ratio_at_sp2 = 2/sp_hill_raw_at_sp2

return
end function

!==========================================================================
