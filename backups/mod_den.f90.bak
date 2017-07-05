
! Description: dealing with density field
!
!      Author: OU Yuyuan <ouyuyuan@lasg.iap.ac.cn>
!     Created: 2015-10-03 07:41:56 BJT
! Last Change: 2016-05-12 19:43:33 BJT

module mod_den

  use mod_arrays, only: &
    gi2, git, &
    bnd, ts, ch

  use mod_con, only: deltat, deltas

  use mod_kind, only: wp

  use mod_param, only: ni, nj, nk, tc

  use mod_type, only: type_mat, type_gvar_m3d

  implicit none
  private

  public &
    den_rho, &
    den_alpha, &
    den_prho
    
contains !{{{1

subroutine den_alpha (alpha, ts, ch, mask)!{{{1
  ! dianose specific volume from (T, S) and 
  !   normalized bottom pressure ch

  type (type_gvar_m3d) :: ts
  real (kind=wp), dimension(ni,nj,nk) :: alpha
  real (kind=wp), dimension(ni,nj), intent(in) :: ch
  integer, dimension(ni,nj,nk), intent(in) :: mask

  integer :: i, j, k

  do k = 1, nk
  do j = 1, nj
  do i = 1, ni
    if ( mask(i,j,k) == 1 ) &
      alpha(i,j,k) = 1.0 /   &
        den_rho( ts%x(1)%v(i,j,k), &
                 ts%x(2)%v(i,j,k), &
                 ch(i,j) * git%pr(k) + bnd%pa%v(i,j) )
  end do
  end do
  end do

end subroutine den_alpha

pure function den_rho (t, s, p0) !{{{1
  ! Calculates the density of seawater using the
  !   standard equation of state recommended by unesco(1981).
  !	  Coefficients of rho0 is given by \citet{Millero1981}
  !	  Coefficients of k is according to \citet{Jackett1995}
  !
  ! t = potential temperature, Degrees centigrade
  ! s = salinity, practical salinity units
  ! p = pressure, Pa
  !
  ! output  dens: kg / m^3
  !
  ! references:
  !	
  !	Coefficients of K is according to Jackett and Mcdougail
  !	J. Atmos. & Ocean. Tech.    1995 Apr., P381-389
  !
  !	Coefficients of rho0 is given by Millero and A.Poisson
  !	Deep-Sea Res., 28A, 625-629
  !

  real (kind=wp), intent(in) :: &
    t  ,& ! potential temperature, degrees Celcius
    s  ,& ! salinity, practical salinity units
    p0    ! sea pressure (preclude surface atmaspheric pressure), in Pa

  real (kind=wp) :: den_rho

  real (kind=wp) :: p, p2, & 
    t2, t3, t4, t5, &
    s32, s2, &
    r0, r, a, b, k

  ! convert Pa to bar
  p = p0 * 1e-5 + 1.013
!  p = p0 * 1e-5

  p2  = p*p
  t2  = t*t
  t3  = t2*t
  t4  = t3*t
  t5  = t4*t
  s32 = s*sqrt(s)
  s2  = s*s

  r0 = 999.842594 + 6.793952e-2*t - 9.095290e-3*t2  &
       + 1.001685e-4*t3 - 1.120083e-6*t4 + 6.536336e-9*t5
  a  = 8.24493e-1 - 4.0899e-3*t + 7.6438e-5*t2 - 8.2467e-7*t3 + 5.3875e-9*t4
  b  = -5.72466e-3 + 1.0227e-4*t - 1.6546e-6*t2

  r  = r0 + a*s + b*s32 + 4.8314e-4 * s2

  k  =   1.965933e4  + 1.444304e2 *t - 1.706103   *t2 + 9.648704e-3*t3 - 4.190253e-5*t4  &
     + ( 5.284855e1  - 3.101089e-1*t + 6.283263e-3*t2 - 5.084188e-5*t3 ) * s &
     + ( 3.886640e-1 + 9.085835e-3*t - 4.619924e-4*t2 ) * s32  &
     + ( 3.186519    + 2.212276e-2*t - 2.984642e-4*t2 + 1.956415e-6*t3 ) * p &
     + ( 6.704388e-3 - 1.847318e-4*t + 2.059331e-7*t2 ) * s * p &
     +   1.480266e-4*s32*p  &
     + ( 2.102898e-4 - 1.202016e-5*t + 1.394680e-7*t2 ) * p2 &
     + (-2.040237e-6 + 6.128773e-8*t + 6.207323e-10*t2) * s * p2

  den_rho = r / ( 1.0 - p / k )
end function den_rho

subroutine den_prho( prho ) !{{{1
  ! calc. the partial derivative of density rho
  ! prho = {partial rho} / {partial t/s}
  ! potential density = prho%x(1)*(t-t0) + prho%x(2)*(s-s0) + den0
  !
  ! it will compare prho%x(1)*t+prho%x(2)*s rather than potential density 
  !   itself in subroutine int_convect
  !
  type (type_mat), dimension(ni,nj,nk) :: prho

  real (kind=wp) :: t, s, p, rhoa, rhob
  integer :: i, j, k

  do k = 1, nk
  do j = 1, nj
  do i = 1, ni
  if ( ts(tc)%x(1)%g%msk(i,j,k) > 0 ) then
    t = ts(tc)%x(1)%v(i,j,k)
    s = ts(tc)%x(2)%v(i,j,k)
    p = ch%tc(i,j) * gi2%pr(k)

    rhoa = den_rho( t + deltat, s, p )
    rhob = den_rho( t - deltat, s, p )
    prho(i,j,k)%x(1) = (rhoa - rhob) / (2*deltat)

    rhoa = den_rho( t, s + deltas, p )
    rhob = den_rho( t, s - deltas, p )
    prho(i,j,k)%x(2) = (rhoa - rhob) / (2*deltas)
  end if
  end do
  end do
  end do

end subroutine den_prho 

end module mod_den !{{{1
!-------------------------------------------------------{{{1
! vim:fdm=marker:fdl=0:
! vim:foldtext=getline(v\:foldstart).'...'.(v\:foldend-v\:foldstart):
