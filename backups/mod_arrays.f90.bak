
! Description: global arrays
!
!      Author: OU Yuyuan <ouyuyuan@lasg.iap.ac.cn>
!     Created: 2015-09-26 15:40:39 BJT
! Last Change: 2017-12-02 10:19:52 BJT

module mod_arrays

  ! imported variables !{{{1
  use mod_kind, only: wp, one, zero

  use mod_param, only: &
    tc, tp, tpp, &
    ni, nj, nk, nkp, glo_ni, glo_nj, &
    myid, mid

  use mod_con, only: am_c, km_c

  use mod_type, only: &
    type_mat, &
    type_accu_gr2d, &
    type_bnd, type_frc, type_bintgu, &
    type_gvar_r4d, type_gvar_r2d, type_gvar_r3d, &
    type_gi, type_gj, type_gij, &
    type_eq_ts, type_eq_uv, type_eq_uvb, type_eq_w, &
    type_eq_ch

  implicit none
  public

  ! grid information !{{{1

  type (type_gi), target :: g1j, g2j, g3j, g4j
  type (type_gi), pointer :: gtj, guj
  type (type_gj), target :: gi1, gi2
  type (type_gj), pointer :: git, giw

  type (type_gij), target :: &
    g11, g21, g31, g41, &
    g12, g22, g32, g42, &
    sg1, sg2, sg3, sg4    ! pseudo 3d grid, 
                          ! only horizontally grid is meaningful
                          ! the third dimension is time, eg, frc
  type (type_gij), pointer :: gt, gu

  ! 'vector' variables !{{{1
  type (type_mat), allocatable, dimension(:,:,:) :: &
    prho ! (p rho)/(p t/s), partial derivative of density

  ! coordinates !{{{1

  ! rad, (lon, lat) position for g1j, g2j, g3j and g4j
  real (kind=wp), allocatable, dimension(:,:,:) :: hposx, hposy

  real (kind=wp), allocatable, dimension(:) :: &
    glo_lon, & ! degree, guj grid, global longitude
    lon,     & !              ... local ...
    glo_lat, & !     ... gtj grid,    ... latitude
    lat,     & !              ... local ...
    z      ! m, gi2 grid, global depths of each layer

  ! accumulated variables, for time-average output !{{{1
  real (kind=wp), allocatable, dimension(:,:,:), target :: &
    act, acs, & ! accumulated (T, S)
    acrho,    & ! kg/m^3, density of sea water
    acu, acv, & ! m/s, unweighted horizontal velocity
    acw   ! m/s, vertical velocity, bottom is assume to be zero, so it is nk layers, not nkp layers
  real (kind=wp), allocatable, dimension(:,:), target :: &
    acch     !  , normalized sea bottom pressure
  type (type_accu_gr2d) :: &
    acssh ! m, sea surface height
  
  ! containers !{{{1
  real (kind=wp), allocatable, dimension(:,:,:), target :: &
    ch  ! bottom pressure, normalized by initial values, the third dimension is for subcontainers

  ! grid variables !{{{1

  type (type_gvar_r4d), target :: &
    temp,  & ! potential temperature
    salt,  & ! salinity
    upres, & ! prognostic pressure weighted zonal velocity
    vpres, & ! prognostic pressure weighted longitinal velocity
    adu, adv  ! baroclinic advection terms (advection of u,v)

  type (type_gvar_r3d), target :: &
    alpha, & ! m^3/kg, specific volume
    adp, & ! m^3/kg * Pa, indefinite integration of alpha from p to prh
    am,   & ! m^2/s, horizontal momentum viscosity coefficient
    km,   & ! m^2/s, vertical momentum viscosity coefficient
    kh,   & ! m^2/s, vertical diffusion coefficient for tracers
    wm,   & ! N m/s, vertical mass advection
    frix, friy ! N, friction forces in momentum equation
  type (type_gvar_r2d), target :: &
    ubtc, vbtc, ubtp, vbtp, & ! current/previous timestep of barotropic velocities
    badvx, badvy, & ! tendency of barotropic velocity
    graphihx, & ! Pa/m, x-component of GRAdient of sea Bottom Pressure
    graphihy, & ! Pa/m, y-component of GRAdient of sea Bottom Pressure
    grapax, grapay, & ! Pa/m, gradient of overloading atmospheric pressure
    cor ! N, Coriolis force
        ! grapay = ( pay*10/rrho_0 ) in pcom 1.0
  type (type_frc) :: frc
  type (type_bnd) :: bnd

  ! equation variables !{{{1
  type (type_eq_ts) :: eqts
  type (type_eq_uv) :: equv
  type (type_eq_uvb) :: equvb
  type (type_eq_w) :: eqw
  type (type_eq_ch) :: eqch

  ! other !{{{1

  ! geopotential height and gradient of gepotential height integration
  type (type_bintgu) :: bphi, bgraphi
  ! bphi%ye = ( pbye*1.0e-4 / rdy ) in pcom 1.0
  ! bgraphi%ye = ( pbye*1.0e-2 * 2.0 ) in pcom 1.0

  ! coefficients in calculating frictional force
  ! cv2 = 2 tan(lat) / a = (2*cv2/rdxu*1.0e2) of version 1.0 
  real (kind=wp), allocatable, dimension(:,:) :: cv1, cv2

  ! interfaces !{{{1
  interface arrays_init
    module procedure init_gvar_r2d
    module procedure init_gvar_r3d
    module procedure init_gvar_r4d
    module procedure init_r3d
    module procedure init_r2d
  end interface arrays_init

  interface link_grid
    module procedure link_grid_gij
    module procedure link_grid_gi
    module procedure link_grid_ver_gij
    module procedure link_grid_ver_gj
  end interface

contains !{{{1

subroutine arrays_allocate () !{{{1
  ! allocate arrays define in this module, 
  !   but global arrays are allocate in mod_master

  integer :: is

  allocate(g1j%rh(ni,nj), stat=is); call chk(is)
  allocate(g1j%dx(ni,nj), stat=is); call chk(is)
  allocate(g1j%tn(ni,nj), stat=is); call chk(is)
  g1j%n = 1; g1j%i = 0; g1j%j = 0

  allocate(g2j%rh(ni,nj), stat=is); call chk(is)
  allocate(g2j%dx(ni,nj), stat=is); call chk(is)
  allocate(g2j%tn(ni,nj), stat=is); call chk(is)
  g2j%n = 2; g2j%i = 1; g2j%j = 0

  allocate(g3j%rh(ni,nj), stat=is); call chk(is)
  allocate(g3j%dx(ni,nj), stat=is); call chk(is)
  allocate(g3j%tn(ni,nj), stat=is); call chk(is)
  g3j%n = 3; g3j%i = 1; g3j%j = 1

  allocate(g4j%rh(ni,nj), stat=is); call chk(is)
  allocate(g4j%dx(ni,nj), stat=is); call chk(is)
  allocate(g4j%tn(ni,nj), stat=is); call chk(is)
  g4j%n = 4; g4j%i = 0; g4j%j = 1

  call link_grid( g1j, g2j, g3j, g4j )

  ! vertical stagger grids
  ! -- gi1 --
  ! -- gi2 --
  !    ...
  ! -- gi1 --
  ! -- gi2 --
  ! -- gi1 --
  allocate(gi1%z(nkp),  stat=is); call chk(is)
  allocate(gi1%pr(nkp),  stat=is); call chk(is)
  allocate(gi1%dz(nkp), stat=is); call chk(is)
  allocate(gi1%dpr(nkp), stat=is); call chk(is)
  gi1%n  = 1;   gi1%k  = 0
  gi1%z  = 0.0; gi1%pr  = 0.0
  gi1%dz = 0.0; gi1%dpr = 0.0

  allocate(gi2%z(nk),  stat=is); call chk(is)
  allocate(gi2%pr(nk),  stat=is); call chk(is)
  allocate(gi2%dz(nk), stat=is); call chk(is)
  allocate(gi2%dpr(nk), stat=is); call chk(is)
  gi2%n  = 2;   gi2%k  = 1
  gi2%z  = 0.0; gi2%pr  = 0.0
  gi2%dz = 0.0; gi2%dpr = 0.0

  ! 3d stagger grids
  allocate(g12%msk(ni,nj,nk), stat=is); call chk(is)
  allocate(g12%lev(ni,nj), stat=is); call chk(is)
  allocate(g12%phih(ni,nj), stat=is); call chk(is)
  allocate(g12%prh(ni,nj), stat=is); call chk(is)
  g12%msk = 0; g12%lev = 0
  g12%phih = 0.0; g12%prh = 0.0
  g12%hg => g1j; g12%vg => gi2

  g22%hg => g2j; g22%vg => gi2

  allocate(g32%msk(ni,nj,nk), stat=is); call chk(is)
  allocate(g32%lev(ni,nj), stat=is); call chk(is)
  allocate(g32%phih(ni,nj), stat=is); call chk(is)
  allocate(g32%prh(ni,nj), stat=is); call chk(is)
  g32%msk = 0; g32%lev = 0
  g32%phih = 0.0; g32%prh = 0.0
  g32%hg => g3j; g32%vg => gi2

  g42%hg => g4j; g42%vg => gi2

  call link_grid( g12, g22, g32, g42 )

  ! wg is only for linking useage
  g11%hg => g1j; g11%vg => gi1
  g21%hg => g2j; g21%vg => gi1
  g31%hg => g3j; g31%vg => gi1
  g41%hg => g4j; g41%vg => gi1

  call link_grid( g11, g21, g31, g41 )

  ! up/down linkage
  call link_grid( gi1, gi2 )

  call link_grid( g12, g11 )
  call link_grid( g22, g21 )
  call link_grid( g32, g31 )
  call link_grid( g42, g41 )

  ! pseudo 3d grids
  sg1%hg => g1j; sg1%vg => null()
  sg2%hg => g2j; sg2%vg => null()
  sg3%hg => g3j; sg3%vg => null()
  sg4%hg => g4j; sg4%vg => null()

  call link_grid( sg1, sg2, sg3, sg4 )

  if ( myid == mid ) then
    allocate(glo_lat(glo_nj), stat=is)
    call chk(is); glo_lat = 0.0
    allocate(glo_lon(glo_ni), stat=is)
    call chk(is); glo_lon = 0.0
  end if

  allocate( bphi%xn(ni,nj), stat=is ); call chk(is)
  allocate( bphi%xs(ni,nj), stat=is ); call chk(is)
  allocate( bphi%ye(ni,nj), stat=is ); call chk(is)
  allocate( bphi%yw(ni,nj), stat=is ); call chk(is)

  allocate( bgraphi%xn(ni,nj), stat=is ); call chk(is)
  allocate( bgraphi%xs(ni,nj), stat=is ); call chk(is)
  allocate( bgraphi%ye(ni,nj), stat=is ); call chk(is)
  allocate( bgraphi%yw(ni,nj), stat=is ); call chk(is)

  call arrays_init( act, ni,nj,nk, zero )
  call arrays_init( acs, ni,nj,nk, zero )

  call arrays_init( acu, ni,nj,nk, zero )
  call arrays_init( acv, ni,nj,nk, zero )
  call arrays_init( acw, ni,nj,nk, zero )

  call arrays_init( acrho, ni,nj,nk, zero )

  call arrays_init( acssh%var, ni,nj, zero, g1j )
  acssh%n = 0

  call arrays_init( acch, ni,nj, zero )

  call arrays_init(temp, ni,nj,nk,2, zero, g12)
  call arrays_init(salt, ni,nj,nk,2, zero, g12)

  call arrays_init(adu, ni,nj,nk,3, zero, g32)
  call arrays_init(adv, ni,nj,nk,3, zero, g32)

  call arrays_init(upres, ni,nj,nk,2, zero, g32)
  call arrays_init(vpres, ni,nj,nk,2, zero, g32)

  call arrays_init(frix, ni,nj,nk, zero, g32)
  call arrays_init(friy, ni,nj,nk, zero, g32)

  call arrays_init(frc%taux, ni,nj,12, zero, sg3)
  call arrays_init(frc%tauy, ni,nj,12, zero, sg3)
  call arrays_init(frc%t,  ni,nj,12, zero, sg1)
  call arrays_init(frc%s,  ni,nj,12, zero, sg1)
  call arrays_init(frc%pa,  ni,nj,12, zero, sg1)
  call arrays_init(frc%fw,  ni,nj,12, zero, sg1)

  call arrays_init(bnd%taux, ni,nj, zero, g1j)
  call arrays_init(bnd%tauy, ni,nj, zero, g1j)
  call arrays_init(bnd%t,    ni,nj, zero, g1j)
  call arrays_init(bnd%s,    ni,nj, zero, g1j)
  call arrays_init(bnd%pa,   ni,nj, zero, g1j)
  call arrays_init(bnd%fw,   ni,nj, zero, g1j)

  ! initiate as fresh water
  call arrays_init(alpha, ni,nj,nk, 1/1000.0*one, g12) 
  call arrays_init(adp, ni,nj,nk, zero, g12) 

  call arrays_init(wm, ni,nj,nkp, zero, g11)

  call arrays_init(am, ni,nj,nk,  am_c, g32)
  call arrays_init(km, ni,nj,nkp, km_c, g32%ud)
  call arrays_init(kh, ni,nj,nkp, zero,  g12%ud)

  call arrays_init(graphihx, ni,nj, zero, g3j)
  call arrays_init(graphihy, ni,nj, zero, g3j)

  call arrays_init(grapax, ni,nj, zero, g3j)
  call arrays_init(grapay, ni,nj, zero, g3j)

  call arrays_init(ch, ni,nj,4, one)

  call arrays_init(cor, ni,nj, zero, g3j)

  call arrays_init(ubtc, ni,nj, zero, g3j)
  call arrays_init(vbtc, ni,nj, zero, g3j)
  call arrays_init(ubtp, ni,nj, zero, g3j)
  call arrays_init(vbtp, ni,nj, zero, g3j)

  call arrays_init(badvx, ni,nj, zero, g3j)
  call arrays_init(badvy, ni,nj, zero, g3j)

  allocate(prho(ni,nj,nk), stat=is); call chk(is) 
  prho%x(1) = 0.0
  prho%x(2) = 0.0

  allocate(hposx(ni,nj,4), stat=is); call chk(is) 
  allocate(hposy(ni,nj,4), stat=is); call chk(is) 
  hposx = 0.0
  hposy = 0.0

  allocate(cv1(ni,nj), stat=is); call chk(is); cv1 = 0.0
  allocate(cv2(ni,nj), stat=is); call chk(is); cv2 = 0.0

  allocate(lat(nj), stat=is); call chk(is); lat = 0.0

  allocate(lon(ni), stat=is); call chk(is); lon = 0.0

  allocate(z  (nk), stat=is); call chk(is); z  = 0.0

  call init_eqts ()

  call init_equv ()

  call init_equvb ()

  call init_eqw ()

  call init_eqch ()
end subroutine arrays_allocate

subroutine init_eqts ()!{{{1
  ! initialize temperature eqution

  eqts%tc => temp%v(:,:,:,tc)
  eqts%sc => salt%v(:,:,:,tc)
  eqts%tp => temp%v(:,:,:,tp)
  eqts%sp => salt%v(:,:,:,tp)

  eqts%act => act
  eqts%acs => acs
  eqts%acr => acrho

  eqts%g => temp%g

  eqts%n = 0

end subroutine init_eqts

subroutine init_equv ()!{{{1
  ! initialize horizontal momentum eqution

  equv%uc => upres%v(:,:,:,tc)
  equv%vc => vpres%v(:,:,:,tc)

  equv%up => upres%v(:,:,:,tp)
  equv%vp => vpres%v(:,:,:,tp)

  equv%acu => acu
  equv%acv => acv

  equv%auc  => adu%v(:,:,:,tc)
  equv%aup  => adu%v(:,:,:,tp)
  equv%aupp => adu%v(:,:,:,tpp)

  equv%avc  => adv%v(:,:,:,tc)
  equv%avp  => adv%v(:,:,:,tp)
  equv%avpp => adv%v(:,:,:,tpp)

  equv%fx   => frix%v
  equv%fy   => friy%v

  equv%pax  => grapax%v
  equv%pay  => grapay%v

  equv%g => upres%g

  equv%n = 0

end subroutine init_equv

subroutine init_equvb ()!{{{1
  ! initialize barotropic momentum eqution

  equvb%uc  => ubtc%v
  equvb%vc  => vbtc%v
  equvb%up  => ubtp%v
  equvb%vp  => vbtp%v

  equvb%ut => badvx%v
  equvb%vt => badvy%v

  equvb%hg => g3j

end subroutine init_equvb

subroutine init_eqw ()!{{{1
  ! initialize dianostic equation for vertical velocity

  eqw%acw => acw

  eqw%g => g11

  eqw%n = 0

end subroutine init_eqw

subroutine init_eqch ()!{{{1
  ! initialize ch equation, i.e., the mass-conservation equation in PCOM

  eqch%chc  => ch(:,:,1)
  eqch%chp  => ch(:,:,2)
  eqch%mch  => ch(:,:,3)
  eqch%mch2 => ch(:,:,4)

  eqch%acch => acch

  eqch%hg => g1j

  eqch%n = 0

end subroutine init_eqch

subroutine init_gvar_r4d (var, d1, d2, d3, d4, ini, g)!{{{1
  ! initialize grid variables
  type (type_gvar_r4d) :: var
  integer, intent(in) :: d1, d2, d3, d4
  real (kind=wp) :: ini
  type (type_gij), target :: g

  integer :: is

  allocate(var%v(d1,d2,d3,d4), stat = is)
  call chk(is)

  var%v = ini
  var%g => g

end subroutine init_gvar_r4d

subroutine init_gvar_r3d (var, d1, d2, d3, ini, g)!{{{1
  ! initialize grid variables
  type (type_gvar_r3d) :: var
  integer, intent(in) :: d1, d2, d3
  real (kind=wp) :: ini
  type (type_gij), target :: g

  integer :: is

  allocate(var%v(d1,d2,d3), stat = is)
  call chk(is)

  var%v = ini
  var%g => g

end subroutine init_gvar_r3d

subroutine init_gvar_r2d (var, d1, d2, ini, hg)!{{{1
  ! initialize grid variables
  type (type_gvar_r2d) :: var
  integer, intent(in) :: d1, d2
  real (kind=wp) :: ini
  type (type_gi), target :: hg

  integer :: is

  allocate(var%v(d1,d2), stat = is)
  call chk(is)
  var%v = ini
  var%hg => hg

end subroutine init_gvar_r2d

subroutine init_r3d (var, d1, d2, d3, ini)!{{{1
  ! initialize 3d real variables
  real (kind=wp), dimension(:,:,:), allocatable :: var
  integer, intent(in) :: d1, d2, d3
  real (kind=wp) :: ini

  integer :: is

  allocate(var(d1,d2,d3), stat = is)
  call chk(is)
  var = ini

end subroutine init_r3d

subroutine init_r2d (var, d1, d2, ini)!{{{1
  ! initialize grid variables
  real (kind=wp), dimension(:,:), allocatable :: var
  integer, intent(in) :: d1, d2
  real (kind=wp) :: ini

  integer :: is

  allocate(var(d1,d2), stat = is)
  call chk(is)
  var = ini

end subroutine init_r2d

subroutine link_grid_gij (grid1, grid2, grid3, grid4) !{{{1
  ! horizontally link the 4 stagger grids
  ! staager grids arrangement
  !             west
  !           1-------4
  !           |       |
  !    south  |       |  north
  !           |       |
  !           2-------3
  !             east
  type (type_gij), target :: grid1, grid2, grid3, grid4

  grid1%ns => grid4
  grid1%ew => grid2
  grid1%di => grid3

  grid2%ns => grid3
  grid2%ew => grid1
  grid2%di => grid4

  grid3%ns => grid2
  grid3%ew => grid4
  grid3%di => grid1

  grid4%ns => grid1
  grid4%ew => grid3
  grid4%di => grid2

end subroutine link_grid_gij

subroutine link_grid_gi (grid1, grid2, grid3, grid4) !{{{1
  ! horizontally link the 4 stagger grids
  ! staager grids arrangement
  !             west
  !           1-------4
  !           |       |
  !    south  |       |  north
  !           |       |
  !           2-------3
  !             east
  type (type_gi), target :: grid1, grid2, grid3, grid4

  grid1%ns => grid4
  grid1%ew => grid2
  grid1%di => grid3

  grid2%ns => grid3
  grid2%ew => grid1
  grid2%di => grid4

  grid3%ns => grid2
  grid3%ew => grid4
  grid3%di => grid1

  grid4%ns => grid1
  grid4%ew => grid3
  grid4%di => grid2

end subroutine link_grid_gi

subroutine link_grid_ver_gij( grid1, grid2 ) !{{{1
  ! up/down linkage between grids
  type (type_gij), target:: grid1, grid2

  grid1%ud => grid2
  grid2%ud => grid1

end subroutine link_grid_ver_gij

subroutine link_grid_ver_gj( grid1, grid2 ) !{{{1
  ! up/down linkage between grids
  type (type_gj), target :: grid1, grid2

  grid1%ud => grid2
  grid2%ud => grid1

end subroutine link_grid_ver_gj

subroutine chk( ista ) !{{{1
  ! check state of allocate array 

  integer, intent(in) ::  ista

  if ( ista /= 0 ) then
    write(*,*) 'Allocate array failed. Stop'
    stop 2
  end if
end subroutine chk

end module mod_arrays!{{{1
!-------------------------------------------------------{{{1
! vim:fdm=marker:fdl=0:
! vim:foldtext=getline(v\:foldstart).'...'.(v\:foldend-v\:foldstart):
