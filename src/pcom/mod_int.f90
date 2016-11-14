
! Description: prognostic/diagnostic of integration variables
!              subroutines called by the integrating process
!
!      Author: OU Yuyuan <ouyuyuan@lasg.iap.ac.cn>
!     Created: 2015-11-14 07:04:18 BJT
! Last Change: 2016-04-06 16:10:30 BJT

module mod_int

  use mod_den, only: den_rrho

  use mod_arrays, only: &
    hg1, hg2, hg3, hg4, &
    vg2, g1, g3, &
    gt, gu, hgt, hgu, vgt, vgw, &
    am, bnd, ts, &
    bphi, bgraphi, &
    lon, lat, &
    frc, fri, graphib, grapa, glo_lon, glo_lat, &
    rrho, rrhodp, pbt, fcor, &
    up, upb, wm, z, adv, prho, &
    arrays_init, arrays_cp_shape, arrays_free

  use mod_con, only: &
    a, g, rho0, rrho0, tice, gamma_b, gamma_t, gamma_s, &
    smag_c, pi

  use mod_kind, only: wp

  use mod_io, only: io_quick_output

  use mod_mympi, only: mympi_swpbnd, mympi_output, &
    mympi_quick_output

  use mod_op, only: op_ter, op_div, op_lap, op_gra, &
    op_vint, op_px1, op_px2, &
    op_vint_ns, op_vint_ew, op_int_tobot, &
    op_adv, op_fri, op_adv_ts, op_dif

  use mod_param, only: ni, nj, nk, nim, njm, nkp, &
    tp, tc, tpp, tctr, &
    nm, names, myid, mid

  use mod_type, only: type_mat, type_memo_r, type_accu_gm3d, &
    type_accu_gr3d, type_accu_gr2d, &
    type_stg, type_vstg, type_stg3d, type_bnd, &
    type_gvar_r3d, type_gvar_r2d, type_gvar_m3d, &
    type_gvar_m2d, type_bintg3, type_pbt

  implicit none
  private

  public &
    int_trop, &
    int_clin, &
    int_ts, &
    int_pgra, &
    int_readyc, &
    int_ssh, &
    int_bnd
    
contains !{{{1

subroutine int_pgra (rrhodp, bphi, bgraphi) !{{{1
  ! calc. pressure gradient Force due to free surface eta
  type (type_gvar_r3d) :: rrhodp 
  type (type_bintg3) :: bphi, bgraphi
  type (type_gvar_m3d) :: wkgm3d
  type (type_mat), dimension(ni,nj,nk) :: wkmat
  real (kind=wp),  dimension(ni,nj,nk) :: p3d, wk, wkb, wkc
  integer :: is

  ! indefinite integration from p to sea bottom
  call op_int_tobot( rrhodp%v, rrho%v, rrho%g )
  p3d = spread(spread(rrho%g%vg%p,1,nj),1,ni)

  call arrays_cp_shape( rrhodp, rrhodp%g%ew, rrhodp%g%ns, wkgm3d )

  call op_ter(wkb, rrhodp%v, rrhodp%g%hg, rrhodp%g%hg%ew)
  call op_ter(wkc, rrho%v, rrho%g%hg, rrho%g%hg%ew)
  wk  = wkb + p3d * wkc
  call op_vint_ns( wk, bphi%xn, bphi%xs )

  call op_ter(wkb, rrhodp%v, rrhodp%g%hg, rrhodp%g%hg%ns)
  call op_ter(wkc, rrho%v, rrho%g%hg, rrho%g%hg%ns)
  wk  = wkb + p3d * wkc
  call op_vint_ew( wk, bphi%ye, bphi%yw )

  call op_gra( wkgm3d, rrhodp )
  call op_vint_ns( wkgm3d%x(1)%v, bgraphi%xn, bgraphi%xs )
  call op_vint_ew( wkgm3d%x(2)%v, bgraphi%ye, bgraphi%yw )

  call arrays_free( wkgm3d ) ! to prevent memory leakage
end subroutine int_pgra

subroutine int_readyc (grapa, dub, adv, fri) !{{{1
  ! prepare for baroclinic integration
  type (type_gvar_m2d) :: grapa, dub 
  type (type_gvar_m3d) :: adv(3), fri

  type (type_mat), dimension(ni,nj,nk) :: wkm3d
  type (type_mat), dimension(ni,nj) :: wkm2d
  real (kind=wp), dimension(ni,nj) :: spbt ! square root of pbt
  integer :: k

  ! atmospheric gradient
  call op_gra( wkm2d, bnd%pa%v, bnd%pa%hg, bnd%pa%hg%ew, bnd%pa%hg%ns )
  call op_ter( grapa%x(1)%v, wkm2d%x(1), bnd%pa%hg%ew, grapa%x(1)%hg )
  call op_ter( grapa%x(2)%v, wkm2d%x(2), bnd%pa%hg%ns, grapa%x(2)%hg )

  ! the optional parameter 1.0 for interpolation is ressonable 
  !   and is neccessary, otherwise NaN will be created when divide spbt
  call op_ter( spbt, pbt%tc, pbt%hg, hgu, 1.0 )
  spbt = sqrt( spbt )
  where (gu%lev == 0) spbt = 1.0

  ! unaveraged velocity
  call upwelling( wm )

  ! advection of velocities, retain 3 time steps for time filter
  adv(tpp) = adv(tp)
  adv(tp)  = adv(tc)
  call op_adv( adv(tc)%x(1)%v, up(tc)%x(1)%v, spbt )
  call op_adv( adv(tc)%x(2)%v, up(tc)%x(2)%v, spbt )

  if ( tctr%i == 1 ) then
    adv(tpp) = adv(tc)
    adv(tp)  = adv(tc)
  else if ( tctr%i == 2 ) then
    adv(tpp) = adv(tp)
  end if

  call op_fri( wkm3d, spbt, up(tc), bnd%tau )
  fri%x(1)%v = wkm3d%x(1)
  fri%x(2)%v = wkm3d%x(2)

  ! vertical integration of tendency of baroclinic velocity
  dub%x(1)%v = 0.0
  dub%x(2)%v = 0.0
  do k = 1, nk
    wkm2d%x(1) = - adv(tc)%x(1)%v(:,:,k) - spbt*rrho0*grapa%x(1)%v + &
      fri%x(1)%v(:,:,k)
    wkm2d%x(2) = - adv(tc)%x(2)%v(:,:,k) - spbt*rrho0*grapa%x(2)%v + &
      fri%x(2)%v(:,:,k)
    dub%x(1)%v = dub%x(1)%v + wkm2d%x(1)*gu%vg%dp(k)*gu%msk(:,:,k)
    dub%x(2)%v = dub%x(2)%v + wkm2d%x(2)*gu%vg%dp(k)*gu%msk(:,:,k)
  end do
  dub%x(1)%v = dub%x(1)%v / gu%pb
  dub%x(2)%v = dub%x(2)%v / gu%pb

end subroutine int_readyc

subroutine int_trop (pbt, dub) !{{{1
  ! barotropic integrations for 2 baroclinic steps
  type (type_pbt) :: pbt
  type (type_gvar_m2d) :: dub  ! m/s^2, advection of upb

  type (type_mat), dimension(ni,nj) :: &
    supb, & ! scale barotropic velocity
    grapbt, & ! 1/m, bottom pressure gradient
    upbaccu, & ! m/s, accumulate upb in the whole barotropic cycle
    pupb, & ! m/s^2, tendency of upb, (p upb)/(p t)
    pgra, & ! N, pressure gradient force
    grapbtint, & ! N, bottom pressure gradient mul. geopotential int.
    pbtgraint, & ! N, bottom pressure mul. gradient of geopotential int.
    vs, & ! horizontal viscosity
    wkmat
  real (kind=wp), dimension(ni,nj) :: &
    spbt, & ! U-grid
    wk
  integer, dimension(ni,nj) :: seau ! sea mask of U-grid
  integer :: nt, nstep

  seau = gu%msk(:,:,1)
  upbaccu%x(1) = upb(tc)%x(1)%v
  upbaccu%x(2) = upb(tc)%x(2)%v

  pbt%bc  = pbt%tc
  pbt%bc2 = pbt%tc
  nstep   = nm%bc / nm%bt * 2

  do nt = 1, nstep
    spbt = 1.0
    call op_ter( wk, pbt%tp, pbt%hg, hgu )
    where ( gu%lev > 0 ) spbt = sqrt(wk)

    ! calculate pbt at predictor time step
    supb%x(1) = gu%pb*spbt * upb(tp)%x(1)%v
    supb%x(2) = gu%pb*spbt * upb(tp)%x(2)%v
    call op_div( wk, supb%x(1), supb%x(2), hgu, hgu, pbt%hg)
    wk = pbt%tp - wk * nm%bt * gamma_b / gt%pb
    where ( gt%lev > 0 ) pbt%tc = wk

    call op_ter( wk, pbt%tc, pbt%hg, hgu )
    where ( gu%lev > 0 ) spbt = sqrt(wk)

    call op_lap( wkmat%x(1), upb(tp)%x(1)%v, hgu, hgu )
    call op_lap( wkmat%x(2), upb(tp)%x(2)%v, hgu, hgu )
    vs%x(1) = seau*am%v(:,:,1)*wkmat%x(1)
    vs%x(2) = seau*am%v(:,:,1)*wkmat%x(2)
    if ( nt==1 ) then
      dub%x(1)%v = dub%x(1)%v - vs%x(1)
      dub%x(2)%v = dub%x(2)%v - vs%x(2)
    end if

    ! pressure gradient due to the change of free surface
    ! special interpolation from hg2 to hg3 for pressure gradient force
    call op_gra( grapbt, pbt%tc, pbt%hg, pbt%hg%ew, pbt%hg%ns )
    grapbtint(:,1:njm)%x(1) = 0.5 * ( &
       grapbt(:,1:njm)%x(1) * bphi%xn(:,1:njm) + &
       grapbt(:, 2:nj)%x(1) * bphi%xs(:, 2:nj) )
    grapbtint(1:nim,:)%x(2) = 0.5 * ( &
       grapbt(1:nim,:)%x(2) * bphi%ye(1:nim,:) + &
       grapbt(2:ni, :)%x(2) * bphi%yw(2:ni, :) )
    call mympi_swpbnd( grapbtint )

    call op_ter( wk, pbt%tc, pbt%hg, pbt%hg%ew )
    pbtgraint(:,1:njm)%x(1) = 0.5 * ( &
       wk(:,1:njm) * bgraphi%xn(:,1:njm) + &
       wk(:, 2:nj) * bgraphi%xs(:, 2:nj) )
    call op_ter( wk, pbt%tc, pbt%hg, pbt%hg%ns )
    pbtgraint(1:nim,:)%x(2) = 0.5 * ( &
       wk(1:nim,:) * bgraphi%ye(1:nim,:) + &
       wk(2:ni, :) * bgraphi%yw(2:ni, :) )
    call mympi_swpbnd( pbtgraint )

    pgra%x(1) = (grapbtint%x(1) + pbtgraint%x(1) + graphib%x(1)%v)*seau*spbt
    pgra%x(2) = (grapbtint%x(2) + pbtgraint%x(2) + graphib%x(2)%v)*seau*spbt

    ! advection + viscosiy + pressure gradient + coriolis
    ! forces in land set to zero
    pupb%x(1) = seau*( vs%x(1) + dub%x(1)%v - pgra%x(1) + 0.5*fcor%v*upb(tp)%x(2)%v )
    pupb%x(2) = seau*( vs%x(2) + dub%x(2)%v - pgra%x(2) - 0.5*fcor%v*upb(tp)%x(1)%v )

    ! Coriolis adjustment
    wkmat%x(1) = pupb%x(1) * nm%bt + upb(tp)%x(1)%v
    wkmat%x(2) = pupb%x(2) * nm%bt + upb(tp)%x(2)%v
    upb(tc)%x(1)%v = ( wkmat%x(1) + 0.5*fcor%v*nm%bt*wkmat%x(2) ) / &
                     ( 1 + (0.5*fcor%v*nm%bt)**2 )
    upb(tc)%x(2)%v = ( wkmat%x(2) - 0.5*fcor%v*nm%bt*wkmat%x(1) ) / &
                     ( 1 + (0.5*fcor%v*nm%bt)**2 )

    upb(tc)%x(1)%v = 0.5 * (upb(tc)%x(1)%v + upb(tp)%x(1)%v)
    upb(tc)%x(2)%v = 0.5 * (upb(tc)%x(2)%v + upb(tp)%x(2)%v)

    ! upb at next time step
    pupb%x(1) = seau*( vs%x(1) - pgra%x(1) + dub%x(1)%v + fcor%v*upb(tc)%x(2)%v )
    pupb%x(2) = seau*( vs%x(2) - pgra%x(2) + dub%x(2)%v - fcor%v*upb(tc)%x(1)%v )
    upb(tc)%x(1)%v = pupb%x(1) * nm%bt + upb(tp)%x(1)%v
    upb(tc)%x(2)%v = pupb%x(2) * nm%bt + upb(tp)%x(2)%v
    upb(tp) = upb(tc)

    ! calculate pbt at corrector time step
    supb%x(1) = gu%pb*spbt * upb(tc)%x(1)%v
    supb%x(2) = gu%pb*spbt * upb(tc)%x(2)%v
    call op_div( wk, supb%x(1), supb%x(2), hgu, hgu, pbt%hg)
    pbt%tc = pbt%tp - gt%msk(:,:,1)*wk*nm%bt/gt%pb
    pbt%tp = pbt%tc

    upbaccu%x(1) = upbaccu%x(1) + upb(tc)%x(1)%v
    upbaccu%x(2) = upbaccu%x(2) + upb(tc)%x(2)%v
    pbt%bc2 = pbt%bc2 + pbt%tc
    if ( nt <= nstep/2 ) pbt%bc = pbt%bc + pbt%tc
  end do

  upb(tc)%x(1)%v = upbaccu%x(1) / ( nstep + 1 )
  upb(tc)%x(2)%v = upbaccu%x(2) / ( nstep + 1 )
  upb(tp) = upb(tc)

  pbt%bc  = pbt%bc / ( nstep/2 + 1 )
  pbt%bc2 = pbt%bc2 / ( nstep + 1 )

  pbt%tc = pbt%bc2
  pbt%tp = pbt%tc
end subroutine int_trop

subroutine int_clin (up, acuv) !{{{1
  ! prediction of velocity of baroclinic mode
  type (type_gvar_m3d) :: up(2)
  type (type_accu_gm3d) :: acuv

  type (type_mat), dimension(ni,nj,nk) :: u ! m/s, horizontal currents
  type (type_mat), dimension(ni,nj) :: &
    pgra, & ! pressure gradient, (grid 2, grid 4)
    pup, & ! time tendency of up(:,:,k)
    wkmat
  real (kind=wp), dimension(ni,nj) :: &
    spbt, & ! , square root of pbt, grid 3
    wk
  real (kind=wp), dimension(ni,nj,nkp) :: w ! m/s, vertical velocity
  real (kind=wp) :: c_tc, c_tp, c_tpp  ! advection filter coefficients
  integer :: k

  c_tc  =  23/12
  c_tp  = -16/12
  c_tpp =   5/12

  spbt  = 1.0
  call op_ter( wk, pbt%bc, pbt%hg, gu%hg, 1.0 )
  where ( gu%lev > 0 ) spbt = sqrt(wk)

  do k = 1, nk
    ! pressure gradient
    wk = rrho%g%phib + pbt%bc * rrhodp%v(:,:,k)
    call op_gra( pgra, wk, hgt, hgt%ew, hgt%ns )

    wk = pbt%bc * vgt%p(k)
    call op_gra( wkmat, wk, hgt, hgt%ew, hgt%ns )

    call op_ter( wk, rrho%v(:,:,k), hgt, hgt%ew ) 
    pgra%x(1) = pgra%x(1) + wk*wkmat%x(1)
    call op_ter( wk, rrho%v(:,:,k), hgt, hgt%ns )
    pgra%x(2) = pgra%x(2) + wk*wkmat%x(2)

    call op_ter( wk, pgra%x(1), hgt%ew, hgu )
    pup%x(1) = - (wk + rrho0*grapa%x(1)%v)*spbt
    call op_ter( wk, pgra%x(2), hgt%ns, hgu )
    pup%x(2) = - (wk + rrho0*grapa%x(2)%v)*spbt

    ! advection
    pup%x(1) = pup%x(1) - ( c_tc*adv(tc)%x(1)%v(:,:,k) + &
                            c_tp*adv(tp)%x(1)%v(:,:,k) + &
                            c_tpp*adv(tpp)%x(1)%v(:,:,k) )
    pup%x(2) = pup%x(2) - ( c_tc*adv(tc)%x(2)%v(:,:,k) + &
                            c_tp*adv(tp)%x(2)%v(:,:,k) + &
                            c_tpp*adv(tpp)%x(2)%v(:,:,k) )

    ! friction
    pup%x(1) = pup%x(1) + fri%x(1)%v(:,:,k)
    pup%x(2) = pup%x(2) + fri%x(2)%v(:,:,k)

    ! Coriolis
    pup%x(1) = pup%x(1) + fcor%v * up(tc)%x(2)%v(:,:,k)
    pup%x(2) = pup%x(2) - fcor%v * up(tc)%x(1)%v(:,:,k)

    ! Coriolis adjustment
    wkmat = pup
    wk = 0.5 * fcor%v * nm%bc
    pup%x(1) = wkmat%x(1) + wk * wkmat%x(2)
    pup%x(2) = wkmat%x(2) - wk * wkmat%x(1)
    wk = 1 + wk**2
    pup%x(1) = pup%x(1)*gu%msk(:,:,k) / wk
    pup%x(2) = pup%x(2)*gu%msk(:,:,k) / wk

    ! prognose of up
    up(tc)%x(1)%v(:,:,k) = up(tp)%x(1)%v(:,:,k) + pup%x(1) * nm%bc
    up(tc)%x(2)%v(:,:,k) = up(tp)%x(2)%v(:,:,k) + pup%x(2) * nm%bc
  end do

 ! interaction between barotropic and baroclinic modes

  call op_vint( wkmat%x(1), up(tc)%x(1)%v, up(tc)%x(1)%g )
  call op_vint( wkmat%x(2), up(tc)%x(2)%v, up(tc)%x(2)%g )

  do k = 1, nk
    where ( gu%msk(:,:,k) > 0 )
      up(tc)%x(1)%v(:,:,k) = up(tc)%x(1)%v(:,:,k) - wkmat%x(1) + upb(tc)%x(1)%v
      up(tc)%x(2)%v(:,:,k) = up(tc)%x(2)%v(:,:,k) - wkmat%x(2) + upb(tc)%x(2)%v
    end where
  end do
  up(tp) = up(tc)

  ! diagnose unweighted 3d currents
  call op_ter( wk, pbt%bc2, pbt%hg, gu%hg, 1.0 )
  spbt = sqrt(wk) ! why not set land to 1.0, as previous does?

  do k = 1, nk
    u(:,:,k)%x(1) = up(tc)%x(1)%v(:,:,k) * spbt / pbt%bc2
    u(:,:,k)%x(2) = up(tc)%x(2)%v(:,:,k) * spbt / pbt%bc2
  end do

  ! accumulate ouput fields for time-average output
  acuv%var%x(1)%v = acuv%var%x(1)%v + u%x(1)
  acuv%var%x(2)%v = acuv%var%x(2)%v + u%x(2)
  acuv%n = acuv%n + 1

! calculate new am by Smagorinsky Scheme
  call smag( am%v, u%x(1), u%x(2) )

end subroutine int_clin

subroutine smag( am, u, v ) !{{{2
  ! lateral viscosity by Smagorinsky Scheme
  ! a = c * sqrt( dt^2 + ds^2 )
  real (kind=wp), dimension(ni,nj,nk) :: am
  ! unweighted horizontal currents
  real (kind=wp), dimension(ni,nj,nk), intent(in) :: u, v 

  real (kind=wp) :: c, dx, dy, dt, ds, am_max, wk
  integer :: i, j, k

  do k = 1, nk
  do j = 2, nj-1
  do i = 2, ni-1
    dx = a*( hgu%ew%rh1(i,j)  *hgu%ew%dx(i,j)%x(1) + &
             hgu%ew%rh1(i+1,j)*hgu%ew%dx(i+1,j)%x(1) )
    dy = a*( hgu%ns%dx(i,j)%x(2) + hgu%ns%dx(i,j+1)%x(2) )

    dt = 4*(u(i+1,j,k) - u(i-1,j,k))/dx - &
         4*(v(i,j+1,k) - v(i,j-1,k))/dy + &
         v(i,j,k) * hgu%tn(i,j) / a

    ds = 4*(v(i+1,j,k) - v(i-1,j,k))/dx + &
         4*(u(i,j+1,k) - u(i,j-1,k))/dy - &
         u(i,j,k) * hgu%tn(i,j) / a

    c = ( 0.5*smag_c/pi * (0.5*dx + 0.5*dy) )**2

    wk     = c * sqrt( dt*dt + ds*ds )
    am_max = ( dx*dx + dy*dy ) / (4*8*nm%bc)

    am(i,j,k) = min( wk, am_max )
  end do
  end do
  end do

  call mympi_swpbnd(am)

end subroutine smag

subroutine int_ts (ts, acts, acw) !{{{1
  ! prognose of temperature and salinity
  type (type_gvar_m3d) :: ts(2)
  type (type_accu_gm3d) :: acts
  type (type_accu_gr3d) :: acw

  type (type_mat), dimension(ni,nj,nk) :: &
    pts, & ! tendency of tracer
    wkm
  real (kind=wp), dimension(ni,nj) :: pbnd
  integer :: k

  call upwelling( wm )

  ! diagnostic the vertical velocity in z-coordinate
  do k = 1, nk
    acw%var%v(:,:,k) = acw%var%v(:,:,k) + &
      wm%v(:,:,k)*rrho%v(:,:,k) / (g*pbt%bc2)
  end do
  acw%n = acw%n + 1

  call op_adv_ts( wkm%x(1), ts(tc)%x(1)%v, wm%v )
  call op_adv_ts( wkm%x(2), ts(tc)%x(2)%v, wm%v )
  pts%x(1) = - wkm%x(1)
  pts%x(2) = - wkm%x(2)

  call op_dif( wkm%x(1), pbt%bc, ts(tp)%x(1)%v )
  call op_dif( wkm%x(2), pbt%bc, ts(tp)%x(2)%v )
  pts%x(1) = pts%x(1) + wkm%x(1)
  pts%x(2) = pts%x(2) + wkm%x(2)

  ! restored boundary condition
  pbnd = gamma_t*( bnd%ts%x(1)%v - ts(tc)%x(1)%v(:,:,1) ) &
         * ts(tc)%x(1)%g%msk(:,:,1) / ts(tc)%x(1)%g%vg%dp(1) 
  pts(:,:,1)%x(1) = pts(:,:,1)%x(1) + pbnd

  pbnd = gamma_s*( bnd%ts%x(2)%v - ts(tc)%x(2)%v(:,:,1) ) &
         * ts(tc)%x(2)%g%msk(:,:,1) * pbt%bc
  pts(:,:,1)%x(2) = pts(:,:,1)%x(2) + pbnd

  ! compute T&S at next time level
  ts(tc)%x(1)%v = ts(tp)%x(1)%v + pts%x(1) * nm%bc * &
                  ts(tc)%x(1)%g%msk / spread(pbt%bc, 3, nk)
  ts(tc)%x(2)%v = ts(tp)%x(2)%v + pts%x(2) * nm%bc * &
                  ts(tc)%x(2)%g%msk / spread(pbt%bc, 3, nk)

  ! set T > -1.5 temporarily since no seaice model
  where ( ts(tc)%x(1)%g%msk > 0 .and. ts(tc)%x(1)%v < tice ) &
    ts(tc)%x(1)%v = tice

  ! convective adjustment
  call convect( ts(tc) )

  ts(tp) = ts(tc)

  ! accumulate ouput fields for time-average output
  acts%var%x(1)%v = acts%var%x(1)%v + ts(tc)%x(1)%v
  acts%var%x(2)%v = acts%var%x(2)%v + ts(tc)%x(2)%v
  acts%n = acts%n + 1

  ! calc. density and surface height

end subroutine int_ts

subroutine convect (tra) !{{{2
  ! convective adjustment if unstable stratification occurs
  ! a full convective adjustment scheme, based on GFDL's version
  ! may be a new scheme can be adopt, like the one I think in my mind.
  type (type_gvar_m3d) :: tra ! tracers (temperature and salinity)
  real (kind=wp), dimension(nk) :: &
    rhoa, & ! density anomaly
    rhob ! same as roha, but reference to the level below
  integer :: i, j, kcnt, nlev

  do i = 1, ni
  do j = 1, nj
    nlev = tra%x(1)%g%lev(i,j)
    if ( nlev > 1 ) then 
      rhoa = 0.0
      rhob = 0.0
      rhoa(2:nlev)   = prho(i,j,2:nlev)%x(1)*tra%x(1)%v(i,j,2:nlev) + &
                       prho(i,j,2:nlev)%x(2)*tra%x(2)%v(i,j,2:nlev)
      rhob(1:nlev-1) = prho(i,j,2:nlev)%x(1)*tra%x(1)%v(i,j,1:nlev-1) + &
                       prho(i,j,2:nlev)%x(2)*tra%x(2)%v(i,j,1:nlev-1)
      kcnt = 1
      do while ( kcnt < nlev )
        ! adjust the upper most pair of unstable layers
        if ( rhob(kcnt) > rhoa(kcnt+1) ) then
          call adjust( tra, rhoa, rhob, kcnt, nlev, i, j )
        else
          kcnt = kcnt + 1
        end if
      end do
    end if
  end do
  end do

end subroutine convect

subroutine adjust( tra, rhoa, rhob, kcnt, nlev, i, j ) !{{{2
  ! mixing the (k, k+1) pair of water layers
  type (type_gvar_m3d) :: tra
  real (kind=wp), dimension(:) :: rhoa, rhob
  integer :: kcnt
  integer, intent(in) :: nlev, i, j

  integer :: na, nb

  ! mix the two unstable layers
  na = kcnt
  nb = kcnt + 1
  call mix( tra, rhoa, rhob, na, nb, i, j )

  ! mix the lower layers if needed
  do while ( nb < nlev )
    if ( rhob(nb) > rhoa(nb+1) ) then
      nb = nb + 1
      call mix( tra, rhoa, rhob, na, nb, i, j )
    else
      if ( (na <= 1) .or. (rhob(na-1) <= rhoa(na)) ) exit
      ! mix the upper layer before checking the below layer again
      na = na - 1
      call mix( tra, rhoa, rhob, na, nb, i, j )
    end if
  end do

  ! mix the upper layers if needed after checking all the lower layers
  do while ( na > 1 )
    if ( rhob(na-1) <= rhoa(na) ) exit
    na = na - 1
    call mix( tra, rhoa, rhob, na, nb, i, j )
  end do

  ! search downwards
  kcnt = nb

end subroutine adjust

subroutine mix (tra, rhoa, rhob, na, nb, i, j) !{{{2
  ! mix water layers through the na th layer to the nb th layer
  type (type_gvar_m3d) :: tra
  real (kind=wp), dimension(:) :: rhoa, rhob
  integer, intent(in) :: na, nb, i, j

  tra%x(1)%v(i,j,na:nb) = sum( tra%x(1)%v(i,j,na:nb)*vgt%dp(na:nb) ) / sum(vgt%dp(na:nb))
  tra%x(2)%v(i,j,na:nb) = sum( tra%x(2)%v(i,j,na:nb)*vgt%dp(na:nb) ) / sum(vgt%dp(na:nb))

  rhoa(na:nb)   = prho(i,j,na:nb)%x(1)  *tra%x(1)%v(i,j,na:nb)   + &
                  prho(i,j,na:nb)%x(2)  *tra%x(2)%v(i,j,na:nb)
  rhob(na:nb-1) = prho(i,j,na+1:nb)%x(1)*tra%x(1)%v(i,j,na:nb-1) + &
                  prho(i,j,na+1:nb)%x(2)*tra%x(2)%v(i,j,na:nb-1)

  if ( nb < nk ) & 
    rhob(nb) = prho(i,j,nb+1)%x(1)*tra%x(1)%v(i,j,nb) + &
               prho(i,j,nb+1)%x(2)*tra%x(2)%v(i,j,nb)
end subroutine mix

subroutine int_ssh (acssh, rrho) !{{{1
  ! diagnose sea surface height
  ! ssh = -h + int(1/(rho g)) dp
  type (type_accu_gr2d) :: acssh
  real (kind=wp), dimension(ni,nj,nk) :: rrho

  real (kind=wp), dimension(ni,nj) :: ssh
  integer :: i, j

  ! calc. density with different pbt and ts
  call den_rrho (rrho, ts(tp), pbt%bc, gt%msk)

  do j = 1, nj
  do i = 1, ni
    ssh(i,j) = sum(rrho(i,j,:) * vgt%dp*pbt%bc(i,j) * gt%msk(i,j,:))/g
  end do
  end do
  ssh = gt%phib / g + ssh

  acssh%var%v = acssh%var%v + ssh
  acssh%n = acssh%n + 1

end subroutine int_ssh

subroutine upwelling( wm ) !{{{1
  ! diagnostic vertical mass advection wm
  type (type_gvar_r3d) :: wm ! horizontally grid 1

  real (kind=wp), dimension(ni,nj,nk) :: wk3d, spbt3d
  real (kind=wp), dimension(ni,nj) :: &
    spbt, & ! square root of pbt, horizontally on grid 3
    wk2d
  integer :: k

  call op_ter( spbt, pbt%tp, pbt%hg, g3%hg )
  spbt3d = spread( sqrt(spbt), 3, nk )
  call op_div( wk3d, spbt3d*up(tc)%x(1)%v, spbt3d*up(tc)%x(2)%v, &
                up(tc)%x(1)%g%hg, up(tc)%x(2)%g%hg, wm%g%hg )
  call op_vint(wk2d, wk3d, g1)
  wk2d = - wk2d

  do k = 2, nk
    wm%v(:,:,k) = wm%v(:,:,k-1) + vg2%dp(k-1) * (wk2d + wk3d(:,:,k-1))
    wm%v(:,:,k) = wm%v(:,:,k) * g1%msk(:,:,k)
  end do
end subroutine upwelling

subroutine int_bnd (bnd) !{{{1
  ! interpolate monthly climatic atmospheric forcing fiels to daily 
  use mod_pro, only: pro_days_month

  type (type_bnd) :: bnd

  integer :: y, m, d, m1, m2, days
  real (kind=wp) :: f

  y = tctr%ct%y
  m = tctr%ct%m
  d = tctr%ct%d

  ! forcing is assuming at the 15 of each month
  if ( d < 15 ) then ! day 1-14
    m2 = m
    m1 = m - 1
    if ( m == 1 ) m1 = 12

    days = pro_days_month (y, m1)
    f = ( days - 15.0 + d )/( days - 15.0 + 15 )

  else ! day 15 - end of month
    m1 = m
    m2 = m + 1
    if ( m == 12 ) m2 = 1

    days = pro_days_month (y, m1)
    f = ( d - 15.0 )/( days - 15 + 15.0 )
  end if

  bnd%tau%x(1)%v = ( frc%tau%x(1)%v(:,:,m2) - frc%tau%x(1)%v(:,:,m1) )*f &
                   + frc%tau%x(1)%v(:,:,m1)
  bnd%tau%x(2)%v = ( frc%tau%x(2)%v(:,:,m2) - frc%tau%x(2)%v(:,:,m1) )*f &
                   + frc%tau%x(2)%v(:,:,m1)

  bnd%ts%x(1)%v  = ( frc%ts%x(1)%v(:,:,m2) - frc%ts%x(1)%v(:,:,m1) )*f &
                   + frc%ts%x(1)%v(:,:,m1)
  bnd%ts%x(2)%v  = ( frc%ts%x(2)%v(:,:,m2) - frc%ts%x(2)%v(:,:,m1) )*f &
                   + frc%ts%x(2)%v(:,:,m1)

  bnd%pa%v  = ( frc%pa%v(:,:,m2)  - frc%pa%v(:,:,m1)  )*f + frc%pa%v(:,:,m1)
  bnd%fw%v  = ( frc%fw%v(:,:,m2)  - frc%fw%v(:,:,m1)  )*f + frc%fw%v(:,:,m1)

end subroutine int_bnd


subroutine chk( ista ) !{{{1
  ! check state of allocate array 

  integer, intent(in) ::  ista

  if ( ista /= 0 ) then
    write(*,*) 'Allocate array failed. Stop'
    stop 2
  end if
end subroutine chk

end module mod_int !{{{1
!-------------------------------------------------------{{{1
! vim:fdm=marker:fdl=0:
! vim:foldtext=getline(v\:foldstart).'...'.(v\:foldend-v\:foldstart):
