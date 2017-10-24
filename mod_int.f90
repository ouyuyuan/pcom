
! Description: prognostic/diagnostic of integration variables
!              subroutines called by the integrating process
!
!      Author: OU Yuyuan <ouyuyuan@lasg.iap.ac.cn>
!     Created: 2015-11-14 07:04:18 BJT
! Last Change: 2017-10-24 10:24:24 BJT

module mod_int

  use mod_den, only: den_alpha, den_rho

  use mod_arrays, only: &
    g1j, g2j, g3j, g4j, &
    gi2, g12, g32, &
    gt, gu, gtj, guj, git, giw, &
    am, bnd, ts, &
    bphi, bgraphi, &
    badv, &
    frc, fri, graphih, grapa, &
    alpha, adp, ch, cor, &
    up, upb, z, adv, prho, &
    glo_lon, glo_lat

  use mod_con, only: &
    a, g, rho0, a0, tice, gamma_b, gamma_t, gamma_s, &
    smag_c, pi

  use mod_kind, only: wp, one

  use mod_mympi, only: mympi_swpbnd, mympi_quick_output

  use mod_op, only: op_ter, op_div, op_lap, op_gra, &
    op_vint, op_vint_ns, op_vint_ew, op_int_tobot, &
    op_adv, op_fri, op_adv_ts, op_dif

  use mod_param, only: ni, nj, nk, nim, njm, nkp, &
    tp, tc, tpp, nm

  use mod_type, only: type_mat, type_accu_gm3d, &
    type_accu_gr3d, type_accu_gr2d, &
    type_bnd, tctr, &
    type_gvar_r3d, type_gvar_r2d, type_gvar_m3d, &
    type_gvar_m2d, type_bintgu, type_ch, &
    type_days_month

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

subroutine int_pgra (adp, bphi, bgraphi) !{{{1
  ! calc. pressure gradient Force due to free surface eta
  type (type_gvar_r3d) :: adp 
  type (type_bintgu) :: bphi, bgraphi
  type (type_mat), dimension(ni,nj,nk) :: wkmat
  real (kind=wp),  dimension(ni,nj,nk) :: pr3d, wk, wkb, wkc

  ! indefinite integration from p to sea bottom
  call op_int_tobot( adp%v, alpha%v, alpha%g )
  pr3d = spread(spread(git%pr,1,nj),1,ni)

  call op_ter(wkb, adp%v, gtj, gtj%ew)
  call op_ter(wkc, alpha%v, gtj, gtj%ew)
  wk  = wkb + pr3d * wkc
  call op_vint_ns( wk, bphi%xn, bphi%xs)

  call op_ter(wkb, adp%v, gtj, gtj%ns)
  call op_ter(wkc, alpha%v, gtj, gtj%ns)
  wk  = wkb + pr3d * wkc
  call op_vint_ew( wk, bphi%ye, bphi%yw )

  call op_gra( wkmat, adp%v, gtj, gtj%ew, gtj%ns )
  call op_vint_ns( wkmat%x(1), bgraphi%xn, bgraphi%xs)
  call op_vint_ew( wkmat%x(2), bgraphi%ye, bgraphi%yw )

end subroutine int_pgra

subroutine int_readyc (grapa, badv, adv, fri, wm) !{{{1
  ! prepare for baroclinic integration
  type (type_gvar_m2d) :: grapa, badv 
  type (type_gvar_m3d) :: adv(3), fri
  type (type_gvar_r3d) :: wm

  type (type_mat), dimension(ni,nj) :: wkm2d
  real (kind=wp), dimension(ni,nj) :: sch ! square root of ch
  integer :: k

  ! atmospheric gradient
  call op_gra( wkm2d, bnd%pa%v, bnd%pa%hg, bnd%pa%hg%ew, bnd%pa%hg%ns )
  call op_ter( grapa%x(1)%v, wkm2d%x(1), bnd%pa%hg%ew, grapa%x(1)%hg )
  call op_ter( grapa%x(2)%v, wkm2d%x(2), bnd%pa%hg%ns, grapa%x(2)%hg )

  ! the optional parameter 1.0 for interpolation is ressonable 
  !   and is neccessary, otherwise NaN will be created when divide sch
  call op_ter( sch, ch%tc, ch%hg, guj, one )
  sch = sqrt( sch )
  where (gu%lev == 0) sch = 1.0

  ! unaveraged velocity
  call upwelling( wm )

  ! advection of velocities, retain 3 time steps for time filter
  adv(tpp) = adv(tp)
  adv(tp)  = adv(tc)
  call op_adv( adv(tc)%x(1)%v, up(tc)%x(1)%v, sch )
  call op_adv( adv(tc)%x(2)%v, up(tc)%x(2)%v, sch )

  if ( tctr%i == 1 ) then
    adv(tpp) = adv(tc)
    adv(tp)  = adv(tc)
  else if ( tctr%i == 2 ) then
    adv(tpp) = adv(tp)
  end if

  call op_fri( fri, sch, up(tc), bnd%tau )

  ! vertical integration of tendency of baroclinic velocity
  badv%x(1)%v = 0.0
  badv%x(2)%v = 0.0
  do k = 1, nk
    wkm2d%x(1) = - adv(tc)%x(1)%v(:,:,k) - sch*a0*grapa%x(1)%v + &
      fri%x(1)%v(:,:,k)
    wkm2d%x(2) = - adv(tc)%x(2)%v(:,:,k) - sch*a0*grapa%x(2)%v + &
      fri%x(2)%v(:,:,k)
    badv%x(1)%v = badv%x(1)%v + wkm2d%x(1)*gu%vg%dpr(k)*gu%msk(:,:,k)
    badv%x(2)%v = badv%x(2)%v + wkm2d%x(2)*gu%vg%dpr(k)*gu%msk(:,:,k)
  end do

  badv%x(1)%v = badv%x(1)%v / gu%prh
  badv%x(2)%v = badv%x(2)%v / gu%prh

end subroutine int_readyc

subroutine int_trop (ch, upb, acch) !{{{1
  ! barotropic integrations for 2 baroclinic steps
  type (type_ch) :: ch
  type (type_gvar_m2d) :: upb(2)
  type (type_accu_gr2d) :: acch

  type (type_mat), dimension(ni,nj) :: &
    supb, & ! scale barotropic velocity
    grach, & ! 1/m, bottom pressure gradient
    upbaccu, & ! m/s, accumulate upb in the whole barotropic cycle
    pupb, & ! m/s^2, tendency of upb, (p upb)/(p t)
    pgra, & ! N, pressure gradient force
    grachint, & ! N, bottom pressure gradient mul. geopotential int.
    chgraint, & ! N, bottom pressure mul. gradient of geopotential int.
    vs1, & ! horizontal viscosity when nt=1
    vs, & ! horizontal viscosity
    wkmat
  real (kind=wp), dimension(ni,nj) :: &
    sch, & ! U-grid
    wk
  integer, dimension(ni,nj) :: seau ! sea mask of U-grid
  integer :: nt, nstep

  seau = gu%msk(:,:,1)
  upbaccu%x(1) = upb(tc)%x(1)%v
  upbaccu%x(2) = upb(tc)%x(2)%v

  ch%bc  = ch%tc
  ch%bc2 = ch%tc
  nstep   = nm%bc / nm%bt * 2

  do nt = 1, nstep
    sch = 1.0
    call op_ter( wk, ch%tp, ch%hg, guj )
    where ( gu%lev > 0 ) sch = sqrt(wk)

    ! calculate ch at predictor time step
    supb%x(1) = gu%prh*sch * upb(tp)%x(1)%v
    supb%x(2) = gu%prh*sch * upb(tp)%x(2)%v
    call op_div( wk, supb%x(1), supb%x(2), guj, guj, ch%hg)
    wk = ch%tp - wk * nm%bt * gamma_b / gt%prh
    where ( gt%lev > 0 ) ch%tc = wk

    call op_ter( wk, ch%tc, ch%hg, guj )
    where ( gu%lev > 0 ) sch = sqrt(wk)

    call op_lap( wkmat%x(1), upb(tp)%x(1)%v, guj, guj )
    call op_lap( wkmat%x(2), upb(tp)%x(2)%v, guj, guj )
    vs%x(1) = seau*am%v(:,:,1)*wkmat%x(1)
    vs%x(2) = seau*am%v(:,:,1)*wkmat%x(2)
    ! vs is changed every nt step, it only minus the first nt
    ! so this two lines ia not "unnecessary"
    if ( nt==1 ) then
      vs1%x(1) = vs%x(1)
      vs1%x(2) = vs%x(2)
    end if

    ! pressure gradient due to the change of free surface
    ! special interpolation from g2j to g3j for pressure gradient force
    call op_gra( grach, ch%tc, ch%hg, ch%hg%ew, ch%hg%ns )
    grachint(:,1:njm)%x(1) = 0.5 * ( &
       grach(:,1:njm)%x(1) * bphi%xn(:,1:njm) + &
       grach(:, 2:nj)%x(1) * bphi%xs(:, 2:nj) )
    grachint(1:nim,:)%x(2) = 0.5 * ( &
       grach(1:nim,:)%x(2) * bphi%ye(1:nim,:) + &
       grach(2:ni, :)%x(2) * bphi%yw(2:ni, :) )
    call mympi_swpbnd( grachint )

    call op_ter( wk, ch%tc, ch%hg, ch%hg%ew )
    chgraint(:,1:njm)%x(1) = 0.5 * ( &
       wk(:,1:njm) * bgraphi%xn(:,1:njm) + &
       wk(:, 2:nj) * bgraphi%xs(:, 2:nj) )
    call op_ter( wk, ch%tc, ch%hg, ch%hg%ns )
    chgraint(1:nim,:)%x(2) = 0.5 * ( &
       wk(1:nim,:) * bgraphi%ye(1:nim,:) + &
       wk(2:ni, :) * bgraphi%yw(2:ni, :) )
    call mympi_swpbnd( chgraint )

    pgra%x(1) = (grachint%x(1) + chgraint%x(1) + &
      graphih%x(1)%v)*seau*sch
    pgra%x(2) = (grachint%x(2) + chgraint%x(2) + &
      graphih%x(2)%v)*seau*sch

    ! advection + viscosiy + pressure gradient + coriolis
    ! forces in land set to zero
    pupb%x(1) = seau*( vs%x(1)-vs1%x(1) + badv%x(1)%v - pgra%x(1) + &
                       0.5*cor%v*upb(tp)%x(2)%v )
    pupb%x(2) = seau*( vs%x(2)-vs1%x(2) + badv%x(2)%v - pgra%x(2) - &
                       0.5*cor%v*upb(tp)%x(1)%v )

    ! Coriolis adjustment
    wkmat%x(1) = pupb%x(1) * nm%bt + upb(tp)%x(1)%v
    wkmat%x(2) = pupb%x(2) * nm%bt + upb(tp)%x(2)%v
    upb(tc)%x(1)%v = ( wkmat%x(1) + 0.5*cor%v*nm%bt*wkmat%x(2) ) / &
                     ( 1 + (0.5*cor%v*nm%bt)**2 )
    upb(tc)%x(2)%v = ( wkmat%x(2) - 0.5*cor%v*nm%bt*wkmat%x(1) ) / &
                     ( 1 + (0.5*cor%v*nm%bt)**2 )

    upb(tc)%x(1)%v = 0.5 * (upb(tc)%x(1)%v + upb(tp)%x(1)%v)
    upb(tc)%x(2)%v = 0.5 * (upb(tc)%x(2)%v + upb(tp)%x(2)%v)

    ! upb at next time step
    pupb%x(1) = seau*( vs%x(1)-vs1%x(1) - pgra%x(1) + badv%x(1)%v + &
                       cor%v*upb(tc)%x(2)%v )
    pupb%x(2) = seau*( vs%x(2)-vs1%x(2) - pgra%x(2) + badv%x(2)%v - &
                       cor%v*upb(tc)%x(1)%v )
    upb(tc)%x(1)%v = pupb%x(1) * nm%bt + upb(tp)%x(1)%v
    upb(tc)%x(2)%v = pupb%x(2) * nm%bt + upb(tp)%x(2)%v
    upb(tp) = upb(tc)

    ! calculate ch at corrector time step
    supb%x(1) = gu%prh*sch * upb(tc)%x(1)%v
    supb%x(2) = gu%prh*sch * upb(tc)%x(2)%v
    call op_div( wk, supb%x(1), supb%x(2), guj, guj, ch%hg)
    ch%tc = ch%tp - gt%msk(:,:,1)*wk*nm%bt/gt%prh
    if ( nt <= nstep/2 ) ch%bc = ch%bc + ch%tc
    ch%tp = ch%tc

    upbaccu%x(1) = upbaccu%x(1) + upb(tc)%x(1)%v
    upbaccu%x(2) = upbaccu%x(2) + upb(tc)%x(2)%v
    ch%bc2 = ch%bc2 + ch%tc
  end do

  upb(tc)%x(1)%v = upbaccu%x(1) / ( nstep + 1 )
  upb(tc)%x(2)%v = upbaccu%x(2) / ( nstep + 1 )
  upb(tp) = upb(tc)

  ch%bc  = ch%bc / ( nstep/2 + 1 )
  ch%bc2 = ch%bc2 / ( nstep + 1 )

  ch%tc = ch%bc2
  ch%tp = ch%tc

  acch%var%v = acch%var%v + ch%tc
  acch%n = acch%n + 1

  ! sea bottom pressure (include pa) can be calculated by
  ! ( modification of 2.44 of Ou2016_phd,
  ! prh is independence of vertical coordinate
  ! and p, ch, pt all share the same horizontal grid )
  ! bottom pressure = ch%tc * (gt%prh + bnd%pa%v) + bnd%pa%v 

end subroutine int_trop

subroutine int_clin (up, acuv, am) !{{{1
  ! prediction of velocity of baroclinic mode
  type (type_gvar_m3d) :: up(2)
  type (type_accu_gm3d) :: acuv
  real (kind=wp), dimension(ni,nj,nk) :: am

  type (type_mat), dimension(ni,nj,nk) :: u ! m/s, horizontal currents
  type (type_mat), dimension(ni,nj) :: &
    pgra, & ! pressure gradient, (grid 2, grid 4)
    pup, & ! time tendency of up(:,:,k)
    wkmat
  real (kind=wp), dimension(ni,nj) :: &
    sch, & ! , square root of ch, grid 3
    wk
  real (kind=wp) :: c_tc, c_tp, c_tpp  ! advection filter coefficients
  integer :: k

  c_tc  =  23.0/12
  c_tp  = -16.0/12
  c_tpp =   5.0/12

  sch  = 1.0
  call op_ter( wk, ch%bc, ch%hg, gu%hg, one )
  where ( gu%lev > 0 ) sch = sqrt(wk)

  do k = 1, nk
    ! pressure gradient
    wk = gt%phih + ch%bc * adp%v(:,:,k)
    call op_gra( pgra, wk, gtj, gtj%ew, gtj%ns )

    wk = ch%bc
    call op_gra( wkmat, wk, gtj, gtj%ew, gtj%ns )

    call op_ter( wk, alpha%v(:,:,k), gtj, gtj%ew ) 
    pgra%x(1) = pgra%x(1) + wk*git%pr(k)*wkmat%x(1)
    call op_ter( wk, alpha%v(:,:,k), gtj, gtj%ns )
    pgra%x(2) = pgra%x(2) + wk*git%pr(k)*wkmat%x(2)

    call op_ter( wk, pgra%x(1), gtj%ew, guj )
    pup%x(1) = - (wk + a0*grapa%x(1)%v)*sch
    call op_ter( wk, pgra%x(2), gtj%ns, guj )
    pup%x(2) = - (wk + a0*grapa%x(2)%v)*sch

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
    pup%x(1) = pup%x(1) + cor%v * up(tc)%x(2)%v(:,:,k)
    pup%x(2) = pup%x(2) - cor%v * up(tc)%x(1)%v(:,:,k)

    ! Coriolis adjustment
    wkmat = pup
    wk = 0.5 * cor%v * nm%bc
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
  call op_ter( wk, ch%bc2, ch%hg, gu%hg, one )
  sch = sqrt(wk) ! why not set land to 1.0, as previous does?

  do k = 1, nk
    u(:,:,k)%x(1) = up(tc)%x(1)%v(:,:,k) * sch / ch%bc2
    u(:,:,k)%x(2) = up(tc)%x(2)%v(:,:,k) * sch / ch%bc2
  end do

  ! accumulate ouput fields for time-average output
  acuv%var%x(1)%v = acuv%var%x(1)%v + u%x(1)
  acuv%var%x(2)%v = acuv%var%x(2)%v + u%x(2)
  acuv%n = acuv%n + 1

! calculate new am by Smagorinsky Scheme
  call smag( am, u%x(1), u%x(2) )

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
    dx = a*( guj%ew%rh(i,j)  *guj%ew%dx(i,j)%x(1) + &
             guj%ew%rh(i+1,j)*guj%ew%dx(i+1,j)%x(1) )
    dy = a*( guj%ns%dx(i,j)%x(2) + guj%ns%dx(i,j+1)%x(2) )

    dt = 4*(u(i+1,j,k) - u(i-1,j,k))/dx - &
         4*(v(i,j+1,k) - v(i,j-1,k))/dy + &
         v(i,j,k) * guj%tn(i,j) / a

    ds = 4*(v(i+1,j,k) - v(i-1,j,k))/dx + &
         4*(u(i,j+1,k) - u(i,j-1,k))/dy - &
         u(i,j,k) * guj%tn(i,j) / a

    c = ( 0.5*smag_c/pi * (0.5*dx + 0.5*dy) )**2

    wk     = c * sqrt( dt*dt + ds*ds )
    am_max = ( dx*dx + dy*dy ) / (4*8*nm%bc)

    am(i,j,k) = min( wk, am_max )
  end do
  end do
  end do

  call mympi_swpbnd(am)

end subroutine smag

subroutine int_ts (ts, acts, acrho, acw, wm) !{{{1
  ! prognose of temperature and salinity
  type (type_gvar_m3d) :: ts(2)
  type (type_accu_gm3d) :: acts
  type (type_accu_gr3d) :: acrho, acw
  type (type_gvar_r3d) :: wm

  type (type_mat), dimension(ni,nj,nk) :: &
    pts, & ! tendency of tracer
    wkm
  real (kind=wp), dimension(ni,nj,nk) :: wk
  real (kind=wp), dimension(ni,nj) :: pbnd
  integer :: k

  call upwelling( wm )

  ! diagnostic the vertical velocity in z-coordinate
  do k = 1, nk
    acw%var%v(:,:,k) = acw%var%v(:,:,k) + &
      wm%v(:,:,k)*alpha%v(:,:,k) / (g*ch%bc2)
  end do
  acw%n = acw%n + 1

  call op_adv_ts( wkm%x(1), ts(tc)%x(1)%v, wm%v )
  call op_adv_ts( wkm%x(2), ts(tc)%x(2)%v, wm%v )
  pts%x(1) = - wkm%x(1)
  pts%x(2) = - wkm%x(2)

  call op_dif( wkm%x(1), ch%bc, ts(tp)%x(1)%v )
  call op_dif( wkm%x(2), ch%bc, ts(tp)%x(2)%v )
  pts%x(1) = pts%x(1) + wkm%x(1)
  pts%x(2) = pts%x(2) + wkm%x(2)

  ! restored boundary condition
  pbnd = gamma_t*( bnd%ts%x(1)%v - ts(tc)%x(1)%v(:,:,1) ) &
         * ts(tc)%x(1)%g%msk(:,:,1) / ts(tc)%x(1)%g%vg%dpr(1) 
  pts(:,:,1)%x(1) = pts(:,:,1)%x(1) + pbnd

  pbnd = gamma_s*( bnd%ts%x(2)%v - ts(tc)%x(2)%v(:,:,1) ) &
         * ts(tc)%x(2)%g%msk(:,:,1) * ch%bc
  pts(:,:,1)%x(2) = pts(:,:,1)%x(2) + pbnd

  ! compute T&S at next time level
  ts(tc)%x(1)%v = ts(tp)%x(1)%v + pts%x(1) * nm%bc * &
                  ts(tc)%x(1)%g%msk / spread(ch%bc, 3, nk)
  ts(tc)%x(2)%v = ts(tp)%x(2)%v + pts%x(2) * nm%bc * &
                  ts(tc)%x(2)%g%msk / spread(ch%bc, 3, nk)

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

  ! calc. density 
  call den_alpha (wk, ts(tp), ch%tc, gt%msk)
  
  acrho%var%v = acrho%var%v + 1.0/wk
  acrho%n = acrho%n + 1 

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

  tra%x(1)%v(i,j,na:nb) = sum( tra%x(1)%v(i,j,na:nb)*git%dpr(na:nb) ) / sum(git%dpr(na:nb))
  tra%x(2)%v(i,j,na:nb) = sum( tra%x(2)%v(i,j,na:nb)*git%dpr(na:nb) ) / sum(git%dpr(na:nb))

  rhoa(na:nb)   = prho(i,j,na:nb)%x(1)  *tra%x(1)%v(i,j,na:nb)   + &
                  prho(i,j,na:nb)%x(2)  *tra%x(2)%v(i,j,na:nb)
  rhob(na:nb-1) = prho(i,j,na+1:nb)%x(1)*tra%x(1)%v(i,j,na:nb-1) + &
                  prho(i,j,na+1:nb)%x(2)*tra%x(2)%v(i,j,na:nb-1)

  if ( nb < nk ) & 
    rhob(nb) = prho(i,j,nb+1)%x(1)*tra%x(1)%v(i,j,nb) + &
               prho(i,j,nb+1)%x(2)*tra%x(2)%v(i,j,nb)
end subroutine mix

subroutine int_ssh (acssh, alpha) !{{{1
  ! diagnose sea surface height
  ! ssh = -h + int(1/(rho g)) dpr
  type (type_accu_gr2d) :: acssh
  real (kind=wp), dimension(ni,nj,nk) :: alpha

  real (kind=wp), dimension(ni,nj) :: ssh
  integer :: i, j

  ! calc. density with different ch and ts
  call den_alpha (alpha, ts(tp), ch%bc, gt%msk)

  do j = 1, nj
  do i = 1, ni
    ssh(i,j) = sum(alpha(i,j,:) * git%dpr*ch%bc(i,j) * gt%msk(i,j,:))/g
  end do
  end do
  ssh = gt%phih / g + ssh

  acssh%var%v = acssh%var%v + ssh
  acssh%n = acssh%n + 1

end subroutine int_ssh

subroutine upwelling( wm ) !{{{1
  ! diagnostic vertical mass advection wm
  type (type_gvar_r3d) :: wm ! horizontal g1j

  real (kind=wp), dimension(ni,nj,nk) :: wk3d, sch3d
  real (kind=wp), dimension(ni,nj) :: &
    sch, & ! square root of ch, horizonta g3j
    wk2d
  integer :: k

  call op_ter( sch, ch%tp, ch%hg, g32%hg )
  sch3d = spread( sqrt(sch), 3, nk )
  call op_div( wk3d, sch3d*up(tc)%x(1)%v, sch3d*up(tc)%x(2)%v, &
                guj, guj, wm%g%hg )
  call op_vint(wk2d, wk3d, g12)
  wk2d = - wk2d

  do k = 2, nk
    wm%v(:,:,k) = wm%v(:,:,k-1) + gi2%dpr(k-1) * (wk2d + wk3d(:,:,k-1))
    wm%v(:,:,k) = wm%v(:,:,k) * g12%msk(:,:,k)
  end do
end subroutine upwelling

subroutine int_bnd (bnd) !{{{1
  ! interpolate monthly climatic atmospheric forcing fiels to daily 
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

    days = type_days_month (y, m1)
    f = ( days - 15.0 + d )/( days - 15.0 + 15 )

  else ! day 15 - end of month
    m1 = m
    m2 = m + 1
    if ( m == 12 ) m2 = 1

    days = type_days_month (y, m1)
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


end module mod_int !{{{1
!-------------------------------------------------------{{{1
! vim:fdm=marker:fdl=0:
! vim:foldtext=getline(v\:foldstart).'...'.(v\:foldend-v\:foldstart):
