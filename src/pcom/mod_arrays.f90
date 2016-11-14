
! Description: global arrays
!
!      Author: OU Yuyuan <ouyuyuan@lasg.iap.ac.cn>
!     Created: 2015-09-26 15:40:39 BJT
! Last Change: 2016-04-05 18:30:02 BJT

module mod_arrays

  ! imported variables !{{{1
  use mod_kind, only: sglp, wp

  use mod_param, only: &
    nm, tc, tp, tpp, &
    ni, nj, nk, nkp, glo_ni, glo_nj, &
    npro, &
    myid, mid

  use mod_type, only: &
    type_mat, type_memo_r, type_accu_gm3d, &
    type_accu_gr3d, type_accu_gr2d, &
    type_bnd, type_frc, type_bintg3, &
    type_gvar_m3d, type_gvar_m2d, type_gvar_r2d, type_gvar_r3d, &
    type_stg, type_vstg, type_stg3d, type_pbt

  implicit none
  public

  ! grid information !{{{1

  type (type_stg), target :: hg1, hg2, hg3, hg4
  type (type_stg), pointer :: hgt, hgu
  type (type_vstg), target :: vg1, vg2
  type (type_vstg), pointer :: vgt, vgw

  ! pseudo grid (ps_grd1, etc) are for 3d grid variable 
  !   whose last dimension is time
  type (type_stg3d), target :: &
    g1, g2, g3, g4, & ! vertical grid is g2
    wg1, wg2, wg3, wg4, & ! vertical grid is g1
    sg1, sg2, sg3, sg4    ! pseudo 3d grid, 
                          ! only horizontally grid is meaningful
                          ! the third dimension is time, eg, frc
  type (type_stg3d), pointer :: gt, gu

  ! 'vector' variables !{{{1
  type (type_mat), allocatable, dimension(:,:,:) :: &
    prho, & ! (p rho)/(p t/s), partial derivative of density
    hpos ! rad, (lon, lat) position for grid 1,2,3,4

  ! global coordinates !{{{1
  real (kind=wp), allocatable, dimension(:) :: &
    glo_lon, & ! degree, ih-grid, global longitude
    lon,     & !              ... local ...
    glo_lat, & !     ... jh-grid,    ... latitude
    lat,     & !              ... local ...
    z      ! m, k-grid, global depths of each layer


  ! accumulated variables, for time-average output !{{{1
  type (type_accu_gm3d) :: &
    acts, & ! accumulated (T, S)
    acuv ! ... unweighted horizontal velocity
  type (type_accu_gr3d) :: &
    acw   ! m/s, vertical velocity, bottom is assume to be zero, so it is nk layer, not nkp layer
  type (type_accu_gr2d) :: &
    acssh ! m, sea surface height

  ! grid variables !{{{1
  type (type_gvar_m3d) :: &
    adv(3), & ! baroclinic advection terms
    up(2), & ! prognostic pressure weighted 3d velocity
    ts(2), & ! potential temperature and salinity
    fri ! N, friction forces in momentum equation
  type (type_gvar_m2d) :: &
    graphib, & ! Pa/m, GRAdient of sea Bottom Pressure
    grapa, & ! Pa/m, gradient of overloading atmospheric pressure
    dub, & ! tendency of barotropic velocity
    upb(2) ! barotropic velocity
    ! grapa%x(2)%v = ( pay*10/rrho_0 ) in pcom 1.0
  type (type_gvar_r3d) :: &
    rrho, & ! m^3/kg, specific volume
    rrhodp, & ! m^3/kg * Pa, indefinite integration of rrho from p to pb
    am,   & ! , horizontal momentum viscosity coefficient
    km,   & ! , vertical momentum viscosity coefficient
    kh,   & ! , vertical diffusion coefficient for tracers
    wm      ! N m/s, vertical mass advection
  type (type_gvar_r2d) :: &
    fcor ! N, Coriolis force
  type (type_frc) :: frc
  type (type_bnd) :: bnd
  type (type_pbt) :: pbt ! bottom pressure, normalized by initial values

  ! other !{{{1

  ! geopotential height and gradient of gepotential height integration
  type (type_bintg3) :: bphi, bgraphi
  ! bphi%ye = ( pbye*1.0e-4 / rdy ) in pcom 1.0
  ! bgraphi%ye = ( pbye*1.0e-2 * 2.0 ) in pcom 1.0

  ! coefficients in calculating frictional force
  ! cv2 = 2 tan(lat) / a = (2*cv2/rdxu*1.0e2) of version 1.0 
  real (kind=wp), allocatable, dimension(:,:) :: cv1, cv2

  ! interfaces !{{{1
  interface arrays_cp_shape
    module procedure cp_shape_gr2d
    module procedure cp_shape_gr2d_gm2d
    module procedure cp_shape_gr3d
    module procedure cp_shape_gr3d_gm3d
    module procedure cp_shape_gr3d_gm3d_b
    module procedure cp_shape_gm2d
    module procedure cp_shape_gm3d
  end interface

  interface arrays_free
    module procedure free_gr2d
    module procedure free_gr3d
    module procedure free_gm2d
    module procedure free_gm3d
  end interface

  interface arrays_init
    module procedure init_gvar_r2d
    module procedure init_gvar_m2d
    module procedure init_gvar_r3d
    module procedure init_gvar_m3d
    module procedure init_r2d
  end interface arrays_init

  interface link_grid
    module procedure link_grid_stg3d
    module procedure link_grid_stg
    module procedure link_grid_ver_stg3d
    module procedure link_grid_ver_vstg
  end interface

contains !{{{1

subroutine arrays_allocate () !{{{1
  ! allocate arrays define in this module, 
  !   but global arrays are allocate in mod_master
  use mod_param, only: ni, nj, nk, missing_float
  use mod_con, only: am_c, km_c

  integer :: is

  allocate(hg1%rh1(ni,nj), stat=is); call chk(is)
  allocate(hg1%dx(ni,nj), stat=is); call chk(is)
  allocate(hg1%tn(ni,nj), stat=is); call chk(is)
  hg1%n = 1; hg1%i = 0; hg1%j = 0

  allocate(hg2%rh1(ni,nj), stat=is); call chk(is)
  allocate(hg2%dx(ni,nj), stat=is); call chk(is)
  allocate(hg2%tn(ni,nj), stat=is); call chk(is)
  hg2%n = 2; hg2%i = 1; hg2%j = 0

  allocate(hg3%rh1(ni,nj), stat=is); call chk(is)
  allocate(hg3%dx(ni,nj), stat=is); call chk(is)
  allocate(hg3%tn(ni,nj), stat=is); call chk(is)
  hg3%n = 3; hg3%i = 1; hg3%j = 1

  allocate(hg4%rh1(ni,nj), stat=is); call chk(is)
  allocate(hg4%dx(ni,nj), stat=is); call chk(is)
  allocate(hg4%tn(ni,nj), stat=is); call chk(is)
  hg4%n = 4; hg4%i = 0; hg4%j = 1

  call link_grid( hg1, hg2, hg3, hg4 )

  ! vertical stagger grids
  ! -- vg1 --
  ! -- vg2 --
  !    ...
  ! -- vg1 --
  ! -- vg2 --
  ! -- vg1 --
  allocate(vg1%z(nkp),  stat=is); call chk(is)
  allocate(vg1%p(nkp),  stat=is); call chk(is)
  allocate(vg1%dz(nkp), stat=is); call chk(is)
  allocate(vg1%dp(nkp), stat=is); call chk(is)
  vg1%n  = 1;   vg1%k  = 0
  vg1%z  = 0.0; vg1%p  = 0.0
  vg1%dz = 0.0; vg1%dp = 0.0

  allocate(vg2%z(nk),  stat=is); call chk(is)
  allocate(vg2%p(nk),  stat=is); call chk(is)
  allocate(vg2%dz(nk), stat=is); call chk(is)
  allocate(vg2%dp(nk), stat=is); call chk(is)
  vg2%n  = 2;   vg2%k  = 1
  vg2%z  = 0.0; vg2%p  = 0.0
  vg2%dz = 0.0; vg2%dp = 0.0

  ! 3d stagger grids
  allocate(g1%msk(ni,nj,nk), stat=is); call chk(is)
  allocate(g1%lev(ni,nj), stat=is); call chk(is)
  allocate(g1%phib(ni,nj), stat=is); call chk(is)
  allocate(g1%pb(ni,nj), stat=is); call chk(is)
  g1%msk = 0; g1%lev = 0
  g1%phib = 0.0; g1%pb = 0.0
  g1%hg => hg1; g1%vg => vg2

  g2%hg => hg2; g2%vg => vg2

  allocate(g3%msk(ni,nj,nk), stat=is); call chk(is)
  allocate(g3%lev(ni,nj), stat=is); call chk(is)
  allocate(g3%phib(ni,nj), stat=is); call chk(is)
  allocate(g3%pb(ni,nj), stat=is); call chk(is)
  g3%msk = 0; g3%lev = 0
  g3%phib = 0.0; g3%pb = 0.0
  g3%hg => hg3; g3%vg => vg2

  g4%hg => hg4; g4%vg => vg2

  call link_grid( g1, g2, g3, g4 )

  ! wg is only for linking useage
  wg1%hg => hg1; wg1%vg => vg1
  wg2%hg => hg2; wg2%vg => vg1
  wg3%hg => hg3; wg3%vg => vg1
  wg4%hg => hg4; wg4%vg => vg1

  call link_grid( wg1, wg2, wg3, wg4 )

  ! up/down linkage
  call link_grid( vg1, vg2 )

  call link_grid( g1, wg1 )
  call link_grid( g2, wg2 )
  call link_grid( g3, wg3 )
  call link_grid( g4, wg4 )

  ! pseudo 3d grids
  sg1%hg => hg1; sg1%vg => null()
  sg2%hg => hg2; sg2%vg => null()
  sg3%hg => hg3; sg3%vg => null()
  sg4%hg => hg4; sg4%vg => null()

  call link_grid( sg1, sg2, sg3, sg4 )

  if ( myid == 0 ) then
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

  call arrays_init( acts%var, ni,nj,nk, 0.0, g1 )
  acts%n = 0; acts%nrec = 0

  call arrays_init( acuv%var, ni,nj,nk, 0.0, g3 )
  acuv%n = 0; acuv%nrec = 0

  ! acw lies on g3 just for output, this is ad hoc
  call arrays_init( acw%var, ni,nj,nk, 0.0, g3 )
  acw%n = 0; acw%nrec = 0

  call arrays_init( acssh%var, ni,nj, 0.0, hg1 )
  acssh%n = 0; acssh%nrec = 0

  call arrays_init(ts(tc), ni,nj,nk, 0.0, g1)
  call arrays_init(ts(tp), ni,nj,nk, 0.0, g1)

  call arrays_init(adv(tc), ni,nj,nk, 0.0, g3)
  call arrays_init(adv(tp), ni,nj,nk, 0.0, g3)
  call arrays_init(adv(tpp), ni,nj,nk, 0.0, g3)

  call arrays_init(up(tc), ni,nj,nk, 0.0, g3)
  call arrays_init(up(tp), ni,nj,nk, 0.0, g3)

  call arrays_init(fri, ni,nj,nk, 0.0, g3)

  call arrays_init(frc%tau, ni,nj,12, 0.0, sg3)
  call arrays_init(frc%ts,  ni,nj,12, 0.0, sg1)
  call arrays_init(frc%pa,       ni,nj,12, 0.0, sg1)
  call arrays_init(frc%fw,       ni,nj,12, 0.0, sg1)

  call arrays_init(bnd%tau, ni,nj, 0.0, hg1)
  call arrays_init(bnd%ts,  ni,nj, 0.0, hg1)
  call arrays_init(bnd%pa,       ni,nj, 0.0, hg1)
  call arrays_init(bnd%fw,       ni,nj, 0.0, hg1)

  ! initiate as fresh water
  call arrays_init(rrho, ni,nj,nk, 1/1000.0, g1) 
  call arrays_init(rrhodp, ni,nj,nk, 0.0, g1) 

  call arrays_init(wm, ni,nj,nkp, 0.0, wg1)

  call arrays_init(am, ni,nj,nk,  am_c, g3)
  call arrays_init(km, ni,nj,nkp, km_c, g3%ud)
  call arrays_init(kh, ni,nj,nkp, 0.0,  g1%ud)

  call arrays_init(graphib, ni,nj, 0.0, hg3)

  call arrays_init(grapa, ni,nj, 0.0, hg3)

  pbt%hg => hg1
  call arrays_init(pbt%tc, ni,nj, 1.0)
  call arrays_init(pbt%tp, ni,nj, 1.0)
  call arrays_init(pbt%bc, ni,nj, 1.0)
  call arrays_init(pbt%bc2, ni,nj, 1.0)

  call arrays_init(fcor, ni,nj, 0.0, hg3)

  call arrays_init(upb(tc), ni,nj, 0.0, hg3)
  call arrays_init(upb(tp), ni,nj, 0.0, hg3)

  call arrays_init(dub, ni,nj, 0.0, hg3)

  allocate(prho(ni,nj,nk), stat=is); call chk(is) 
  prho%x(1) = 0.0
  prho%x(2) = 0.0

  allocate(hpos(ni,nj,4), stat=is); call chk(is) 
  hpos%x(1) = 0.0
  hpos%x(2) = 0.0

  allocate(cv1(ni,nj), stat=is); call chk(is); cv1 = 0.0
  allocate(cv2(ni,nj), stat=is); call chk(is); cv2 = 0.0

  allocate(lat(nj), stat=is); call chk(is); lat = 0.0

  allocate(lon(ni), stat=is); call chk(is); lon = 0.0

  allocate(z  (nk), stat=is); call chk(is); z  = 0.0
end subroutine arrays_allocate

subroutine cp_shape_gr3d (va, vb) !{{{1
  ! initialize vb as va
  type (type_gvar_r3d), intent(in) :: va
  type (type_gvar_r3d) :: vb

  integer :: d(3), is

  if ( .not.allocated(vb%v) ) then
    d = shape(va%v)
    allocate(vb%v(d(1),d(2),d(3)), stat=is); call chk(is)
  end if

  vb%v = 0.0
  vb%g => va%g

end subroutine cp_shape_gr3d

subroutine cp_shape_gr2d (va, vb) !{{{1
  ! initialize vb as va
  type (type_gvar_r2d), intent(in) :: va
  type (type_gvar_r2d) :: vb

  integer :: d(2), is

  if ( .not.allocated(vb%v) ) then
    d = shape(va%v)
    allocate(vb%v(d(1),d(2)), stat=is); call chk(is)
  end if
  vb%v = 0.0
  vb%hg => va%hg

end subroutine cp_shape_gr2d

subroutine cp_shape_gr2d_gm2d (va, vb) !{{{1
  ! initialize vb as va
  type (type_gvar_r2d), intent(in) :: va
  type (type_gvar_m2d) :: vb

  call cp_shape_gr2d( va, vb%x(1) )
  call cp_shape_gr2d( va, vb%x(2) )

end subroutine cp_shape_gr2d_gm2d

subroutine cp_shape_gm3d (va, vb) !{{{1
  ! initialize vb as va
  type (type_gvar_m3d), intent(in) :: va
  type (type_gvar_m3d) :: vb

  call cp_shape_gr3d( va%x(1), vb%x(1) )
  call cp_shape_gr3d( va%x(2), vb%x(2) )

end subroutine cp_shape_gm3d

subroutine cp_shape_gr3d_gm3d (va, ga, gb, vb) !{{{1
  ! initialize vb as va
  type (type_gvar_r3d), intent(in) :: va
  type (type_stg3d), target :: ga, gb
  type (type_gvar_m3d) :: vb

  call cp_shape_gr3d( va, vb%x(1) )
  vb%x(1)%g => ga

  call cp_shape_gr3d( va, vb%x(2) )
  vb%x(2)%g => gb

end subroutine cp_shape_gr3d_gm3d

subroutine cp_shape_gr3d_gm3d_b (va, vb) !{{{1
  ! initialize vb as va
  type (type_gvar_r3d), intent(in) :: va
  type (type_gvar_m3d) :: vb

  call cp_shape_gr3d( va, vb%x(1) )
  call cp_shape_gr3d( va, vb%x(2) )

end subroutine cp_shape_gr3d_gm3d_b

subroutine cp_shape_gm2d (va, vb) !{{{1
  ! initialize vb as va
  type (type_gvar_m2d), intent(in) :: va
  type (type_gvar_m2d) :: vb

  call cp_shape_gr2d( va%x(1), vb%x(1) )
  call cp_shape_gr2d( va%x(2), vb%x(2) )

end subroutine cp_shape_gm2d

subroutine free_gr2d (var) !{{{1
  ! free memory of var
  type (type_gvar_r2d) :: var

  if ( allocated(var%v) ) then
    var%hg => null()
    deallocate( var%v )
  end if

end subroutine free_gr2d

subroutine free_gr3d (var) !{{{1
  ! free memory of var
  type (type_gvar_r3d) :: var

  if ( allocated(var%v) ) then
    var%g => null()
    deallocate( var%v )
  end if

end subroutine free_gr3d

subroutine free_gm2d (var) !{{{1
  ! free memory of var
  type (type_gvar_m2d) :: var

  call free_gr2d( var%x(1) )
  call free_gr2d( var%x(2) )

end subroutine free_gm2d

subroutine free_gm3d (var) !{{{1
  ! free memory of var
  type (type_gvar_m3d) :: var

  call free_gr3d( var%x(1) )
  call free_gr3d( var%x(2) )

end subroutine free_gm3d

subroutine init_gvar_r3d (var, d1, d2, d3, ini, g)!{{{1
  ! initialize grid variables
  type (type_gvar_r3d) :: var
  integer, intent(in) :: d1, d2, d3
  real (kind=wp) :: ini
  type (type_stg3d), target :: g

  integer :: is

  allocate(var%v(d1,d2,d3), stat = is)
  call chk(is)

  var%v = ini
  var%g => g

end subroutine init_gvar_r3d

subroutine init_gvar_m3d (var, d1, d2, d3, ini, g) !{{{1
  ! initialize grid variable
  type (type_gvar_m3d) :: var
  integer, intent(in) :: d1, d2, d3
  real (kind=wp) :: ini
  type (type_stg3d), target :: g
  call init_gvar_r3d (var%x(1), d1, d2, d3, ini, g)
  call init_gvar_r3d (var%x(2), d1, d2, d3, ini, g)
end subroutine init_gvar_m3d

subroutine init_gvar_r2d (var, d1, d2, ini, hg)!{{{1
  ! initialize grid variables
  type (type_gvar_r2d) :: var
  integer, intent(in) :: d1, d2
  real (kind=wp) :: ini
  type (type_stg), target :: hg

  integer :: is

  allocate(var%v(d1,d2), stat = is)
  call chk(is)
  var%v = ini
  var%hg => hg

end subroutine init_gvar_r2d

subroutine init_gvar_m2d (var, d1, d2, ini, hg) !{{{1
  ! initialize grid variable
  type (type_gvar_m2d) :: var
  integer, intent(in) :: d1, d2
  real (kind=wp) :: ini
  type (type_stg), target :: hg

  call init_gvar_r2d( var%x(1), d1, d2, ini, hg )
  call init_gvar_r2d( var%x(2), d1, d2, ini, hg )

end subroutine init_gvar_m2d

subroutine link_grid_stg3d (grid1, grid2, grid3, grid4) !{{{1
  ! horizontally link the 4 stagger grids
  ! staager grids arrangement
  !             west
  !           1-------4
  !           |       |
  !    south  |       |  north
  !           |       |
  !           2-------3
  !             east
  type (type_stg3d), target :: grid1, grid2, grid3, grid4

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

end subroutine link_grid_stg3d

subroutine link_grid_stg (grid1, grid2, grid3, grid4) !{{{1
  ! horizontally link the 4 stagger grids
  ! staager grids arrangement
  !             west
  !           1-------4
  !           |       |
  !    south  |       |  north
  !           |       |
  !           2-------3
  !             east
  type (type_stg), target :: grid1, grid2, grid3, grid4

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

end subroutine link_grid_stg

subroutine link_grid_ver_stg3d( grid1, grid2 ) !{{{1
  ! up/down linkage between grids
  type (type_stg3d), target:: grid1, grid2

  grid1%ud => grid2
  grid2%ud => grid1

end subroutine link_grid_ver_stg3d

subroutine link_grid_ver_vstg( grid1, grid2 ) !{{{1
  ! up/down linkage between grids
  type (type_vstg), target :: grid1, grid2

  grid1%ud => grid2
  grid2%ud => grid1

end subroutine link_grid_ver_vstg

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
