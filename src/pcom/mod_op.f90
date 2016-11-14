
! Description: operators
!
!      Author: OU Yuyuan <ouyuyuan@lasg.iap.ac.cn>
!     Created: 2015-10-11 15:23:38 BJT
! Last Change: 2016-04-05 14:51:29 BJT

module mod_op

  ! imported variables !{{{1
  use mod_arrays, only: &
    g1, g2, g3, g4, hg1, hg2, hg3, hg4, &
    vg1, vg2, &
    gt, gu, hgt, hgu, vgt, vgw, &
    bnd, am, km, kh, up, wm, pbt, &
    cv1, cv2, &
    arrays_init, arrays_cp_shape, arrays_free, &
    glo_lon, glo_lat, z

  use mod_con, only: &
    a, g, rho0, cdbot, ah_c

  use mod_kind, only: wp

  use mod_mympi, only: mympi_swpbnd, mympi_quick_output

  use mod_param, only: &
    glo_nj, ni, nj, nk, nkp, &
    tc, tp, tctr, &
    nim, njm, myid, mid

  use mod_type, only: type_mat, &
    type_stg, type_vstg, type_stg3d, &
    type_gvar_m3d, type_gvar_m2d, type_gvar_r2d, type_gvar_r3d

  implicit none
  private

  ! exported variables !{{{1
  public  &
    op_int_tobot, &
    op_vint, &
    op_vint_ns, &
    op_vint_ew, &
    op_fri, &
    op_dif, &
    op_adv, &
    op_adv_ts, &
    op_lap, &
    op_div, &
    op_gra, &
    op_px1, &
    op_px2, &
    op_px3, &
    op_ter

  ! interfaces !{{{1

  interface op_div
    module procedure div_r3d
    module procedure div_r2d
  end interface

  interface op_gra
    module procedure gra_r2d
    module procedure gra_r2d_b
    module procedure gra_r3d
    module procedure gra_gr2d
    module procedure gra_gr3d
  end interface

  interface op_px1
    module procedure p1_r3d
    module procedure p1_r2d
  end interface

  interface op_px2
    module procedure p2_r3d
    module procedure p2_r2d
  end interface

  interface op_px3
    module procedure p3_r3d
    module procedure p3_r3d_ter
  end interface

  interface op_lap
    module procedure lap_r2d
    module procedure lap_gm2d
  end interface

  interface op_ter
    module procedure vter_r3d
    module procedure ter_r3d
    module procedure ter_r2d
    module procedure ter_gm2d
    module procedure ter_i3d
    module procedure ter_i2d
  end interface

  interface op_vint
    module procedure vint_r3d
  end interface

  interface op_int_tobot
    module procedure int_tobot
  end interface
contains !{{{1

subroutine fri_r3d (ans, spbt, v, vp, vc, tau) !{{{1
  ! horizontal frictional force
  ! ans is on the same grid as vp
  ! ans = 1/spbt*div(pbt*gra(v)) + cv1*vp + cv2*spbt*(p vc/ p x1) + vert.
  real (kind=wp), dimension(:,:,:) :: ans
  real (kind=wp), dimension(:,:), intent(in) :: spbt
  real (kind=wp), dimension(:,:,:), intent(in) :: &
    v, & ! unweighted velocity 
    vp, & ! pressure weighted velocity
    vc
  real (kind=wp), dimension(ni,nj), intent(in) :: tau

  type (type_mat), dimension(ni,nj,nk) :: wkmat
  real (kind=wp), dimension(ni,nj,nkp) :: wkp
  real (kind=wp), dimension(ni,nj,nk) :: wka, wkb, spbt3d, pbt3da, pbt3db
  real (kind=wp), dimension(ni,nj) :: wkr2d
  real (kind=wp) :: mag
  integer :: is, i, j, k

  spbt3d = spread(spbt, 3, nk)

  ! horizontal viscosity

  ! div(pbt*gra(v))
  call op_gra( wkmat, v, hgu, hgu%ew, hgu%ns )
  call op_ter( wkr2d, pbt%tc, pbt%hg, hgu%ew )
  pbt3da = spread( wkr2d, 3, nk )
  call op_ter( wkr2d, pbt%tc, pbt%hg, hgu%ns )
  pbt3db = spread( wkr2d, 3, nk )
  call div_r3d( wka, pbt3da*wkmat%x(1), pbt3db*wkmat%x(2), &
                hgu%ew, hgu%ns, hgu )
  ans = am%v/spbt3d * wka

  ! cv1*vp
  ans = ans + am%v*spread(cv1,3,nk)*vp

  ! cv2*spbt*(p vc/ p x1)
  ! oddly center difference of vc
  call px_r3d_center( wka, vc, hgu )
  ans = ans + am%v*spread(cv2,3,nk)*spbt3d*wka

  ! vertical viscosity

  ! oddly vertical difference respect to thickness (not pressure)
  call pz_r3d( wkp, vp, gu%msk )
  wkp = g*rho0*km%v * wkp / (spread(spbt,3,nk))**2

  ! surface boundary condition
  wkp(:,:,1) = g*tau(:,:) / spbt(:,:)

  do j = 1, nj
  do i = 1, ni
    k = gu%lev(i,j)
    if ( k > 0 ) then ! bottom drag
      mag = sqrt( v(i,j,k)**2 + vc(i,j,k)**2 )
      wkp(i,j,k+1) = v(i,j,k)*g*cdbot*mag/spbt(i,j)
    end if
  end do
  end do

  call p3_r3d( wka, wkp )

  ans = ans*gu%msk + wka

end subroutine fri_r3d

subroutine op_fri(ans, spbt, up, tau) !{{{1
  ! horizontal frictional force
  ! output is on the same grid as up
  type (type_mat), dimension(:,:,:) :: ans
  real (kind=wp), dimension(:,:), intent(in) :: spbt
  type (type_gvar_m3d), intent(in) :: up
  type (type_gvar_m2d), intent(in) :: tau

  real (kind=wp), dimension(ni,nj,nk) :: u, v, spbt3d

  spbt3d = spread( spbt, 3, nk )
  u = up%x(1)%v / spbt3d
  v = up%x(2)%v / spbt3d

  call fri_r3d( ans%x(1), spbt, u, up%x(1)%v, -v, tau%x(1)%v )
  call fri_r3d( ans%x(2), spbt, v, up%x(2)%v,  u, tau%x(2)%v )

end subroutine op_fri

subroutine op_dif (ans, pbt, var) !{{{1
  ! diffusion of tracers
  ! ans is on the same grid as var
  ! ans = div(pbt*gra(var)) + vert.
  real (kind=wp), dimension(:,:,:) :: ans
  real (kind=wp), dimension(:,:), intent(in) :: pbt
  real (kind=wp), dimension(:,:,:), intent(in) :: var ! tracer

  type (type_mat), dimension(ni,nj,nk) :: wkm
  real (kind=wp), dimension(ni,nj,nkp) :: wkp
  real (kind=wp), dimension(ni,nj,nk) :: wk, pbt3da, pbt3db
  real (kind=wp), dimension(ni,nj) :: wkr2d
  integer :: is, k

  ! horizontal diffusion

  ! div(pbt*gra(var))
  call p1_r3d_b( wkm%x(1), var, hgt, hgt%ew, gt%msk )
  wkm%x(1) = wkm%x(1) / ( a * spread(hgt%ew%rh1, 3, nk) )
  call p2_r3d_b( wkm%x(2), var, hgt, hgt%ns, gt%msk )
  wkm%x(2) = wkm%x(2) / a
  call op_ter( wkr2d, pbt, hgt, hgt%ew )
  pbt3da = spread( wkr2d, 3, nk )
  call op_ter( wkr2d, pbt, hgt, hgt%ns )
  pbt3db = spread( wkr2d, 3, nk )
  call div_r3d( wk, pbt3da*wkm%x(1), pbt3db*wkm%x(2), &
                hgt%ew, hgt%ns, hgt )
  ans = ah_c * wk

  ! vertical diffusion

  ! oddly vertical difference respect to thickness (not pressure)
  call pz_r3d( wkp, var, gt%msk )
  wkp = g*rho0 * kh%v * wkp

  call p3_r3d( wk, wkp )
  ans = ans + wk
end subroutine op_dif

subroutine op_adv (ans, var, spbt) !{{{1
  ! calc. advection of horizontal pressure weighted velocity var
  ! op_adv(var) = div(var*v) - 0.5*var*div(v), 
  !   v = (u, w) are unweighted 3d velocity
  real (kind=wp), dimension(:,:,:) :: ans
  real (kind=wp), dimension(:,:,:), intent(in) :: var
  real (kind=wp), dimension(:,:), intent(in) :: spbt

  type (type_mat), dimension(ni,nj,nk) :: wkm
  real (kind=wp),  dimension(ni,nj,nk) :: wka, wkb, cos3d, spbt3d, u, v
  real (kind=wp),  dimension(ni,nj,nkp):: wkap, wkbp, w

  ! unaveraged velocity, interpolate to proper grid
  spbt3d = spread( spbt, 3, nk )
  wka = up(tc)%x(1)%v / spbt3d
  call ter_r3d( u, wka, hgu, hgu%ew )

  cos3d = spread(hgu%rh1,3,nk)
  wka = up(tc)%x(2)%v / spbt3d
  call ter_r3d( v, wka*cos3d, hgu, hgu%ns )

  wkap = wm%v / spread(pbt%tc,3,nkp)
  call ter_r3d( w, wkap, hgt, hgu )

  ! var*u
  call ter_r3d( wka, var, hgu, hgu%ew )
  call ter_r3d( wkb, var, hgu, hgu%ns )
  wkm%x(1) = wka * u
  wkm%x(2) = wkb * v

  ! horizontal: div(var*u)
  call p1_r3d( wka, wkm%x(1), hgu%ew, hgu )
  call p2_r3d( wkb, wkm%x(2), hgu%ns, hgu )
  ans = (wka + wkb) / (a*cos3d)

  ! var*w
  call vter_r3d( wkap, var, vgt, vgw )
  wkap = wkap*w
  ! set to zero if any of the upper and lower layer of U-grid is land
  wkbp(:,:,1)   = 0
  wkbp(:,:,nkp) = 0
  wkbp(:,:,2:nkp-1) = wkap(:,:,2:nkp-1)*&
                      gu%msk(:,:,1:nk-1)*gu%msk(:,:,2:nk)

  call p3_r3d( wka, wkbp )
  ans = ans + wka

  ! div(v)
  call p1_r3d( wkm%x(1), u, hgu%ew, hgu )
  call p2_r3d( wkm%x(2), v, hgu%ns, hgu )
  call p3_r3d( wka, w )
  ans = ans - 0.5*var*( (wkm%x(1)+wkm%x(2))/(a*cos3d) + wka )

  ans = ans * gu%msk
end subroutine op_adv

subroutine op_adv_ts (ans, var, wm) !{{{1
  ! calc. 3d advection of tracer var
  ! op_adv_ts(var) = vm * grad( var )
  !   vm = (um, wm) are mass weighted 3d velocity
  real (kind=wp), dimension(ni,nj,nk) :: ans
  real (kind=wp), dimension(ni,nj,nk), intent(in) :: var
  real (kind=wp),  dimension(ni,nj,nkp), intent(in) :: wm

  type (type_mat), dimension(ni,nj,nk) :: wkm, um
  real (kind=wp),  dimension(ni,nj,nkp):: wkp
  real (kind=wp),  dimension(ni,nj,nk) :: wka, wkb, cos3d
  real (kind=wp),  dimension(ni,nj) :: spbt
  integer :: k

  call op_ter( spbt, pbt%tc, pbt%hg, hgu, 1.0 )
  spbt = sqrt(spbt)

  do k = 1, nk
    um(:,:,k)%x(1) = up(tc)%x(1)%v(:,:,k) * spbt
    um(:,:,k)%x(2) = up(tc)%x(2)%v(:,:,k) * spbt * hgu%rh1
  end do

  ! horizontal advection
  call p1_r3d( wkm%x(1), var, hgt, hgt%ew )
  call p2_r3d( wkm%x(2), var, hgt, hgt%ns )

  cos3d = spread(hgt%rh1,3,nk)

  call ter_r3d( wka, um%x(1), hgu, hgt%ew )
  call ter_r3d( wkb, wka * wkm%x(1), hgt%ew, hgt )
  ! not interpolate the cosine factor
  ans = wkb / (a*cos3d)

  call ter_r3d( wka, um%x(2), hgu, hgt%ns )
  call ter_r3d( wkb, wka * wkm%x(2), hgt%ns, hgt )
  ! not interpolate the cosine factor
  ans = ans + wkb / (a*cos3d)

  ! vertical advection
  call dx3_r3d_b( wkp, var )
  call vter_r3d( wka, wm * wkp, vgw, vgt ) 
  ! not interpolate the layer 'thickness'
  wka = wka / spread(spread(vgt%dp,1,nj), 1, ni)
  ans = ans + wka

end subroutine op_adv_ts

subroutine vter_r3d(ans, var, vga, vgb, dft) !{{{1
  ! interpolate vertically from grid vga to grid vgb
  real (kind=wp) :: ans(:,:,:)
  real (kind=wp), intent(in) :: var(:,:,:)
  type (type_vstg), target :: vga, vgb
  real (kind=wp), optional :: dft

  integer :: nda

  if ( vga%n == vgb%n ) &
    stop 'no need to interpolate in ans of mod_op'

  nda = size(vga%p)

  if ( size(var,3) /= nda ) &
    stop 'var and vga unmatch in ans of mod_op'

  if ( present(dft) ) then
    ans = dft
  else
    ans = 0.0
  end if

  ! vg2 to vg1
  if ( vga%n == 2 ) then
    ans(:,:,2:nk) = ( var(:,:,1:nk-1) + var(:,:,2:nk) ) * 0.5
  ! from vg1 to vg2
  else
    ans(:,:,1:nk) = ( var(:,:,1:nk) + var(:,:,2:nkp) ) * 0.5
  end if

  ! no need to swap boundary horizontally
end subroutine vter_r3d

subroutine ter_r3d(ans, var, hga, hgb, dft) !{{{1
  ! interpolate from grid hga to grid hgb
  real (kind=wp) :: ans(:,:,:)
  real (kind=wp), intent(in) :: var(:,:,:)
  type (type_stg), target :: hga, hgb
  real (kind=wp), optional :: dft

  if ( present(dft) ) then
    ans = dft
  else
    ans = 0.0
  end if

  if ( hga%n == hgb%n ) then
    stop 'no need to interpolate in ans of mod_op'
  else if ( associated(hga%ew, hgb) ) then
    ans(1+hga%i:nim+hga%i,:,:) = &
    0.5 * ( var(1:nim,:,:) + var(2:ni,:,:) )
  else if ( associated(hga%ns, hgb) ) then
    ans(:,1+hga%j:njm+hga%j,:) = &
      0.5 * ( var(:,1:njm,:) + var(:,2:nj,:) )
  else if ( associated(hga%di, hgb) ) then
    ans(1+hga%i:nim+hga%i,1+hga%j:njm+hga%j,:) = &
      ( var(1:nim,1:njm,:) + var(1:nim,2:nj,:) + &
        var(2:ni, 1:njm,:) + var(2:ni, 2:nj,:) ) * 0.25
  else
    print *, 'unhandled relative position in ter_r3d of mod_op'
    stop
  end if

  call mympi_swpbnd (ans)

end subroutine ter_r3d

subroutine ter_r2d(ans, var, hga, hgb, dft) !{{{1
  ! interpolate from grid hga to grid hgb
  real (kind=wp), dimension(ni,nj) :: ans
  real (kind=wp), intent(in) :: var(:,:)
  type (type_stg), target :: hga, hgb
  real (kind=wp), optional :: dft

  if ( present(dft) ) then
    ans = dft
  else
    ans = 0.0
  end if

  if ( hga%n == hgb%n ) then
    stop 'no need to interpolate in ans of mod_op'
  else if ( associated(hga%ew, hgb) ) then
    ans(1+hga%i:nim+hga%i,:) = &
    0.5 * ( var(1:nim,:) + var(2:ni,:) )
  else if ( associated(hga%ns, hgb) ) then
    ans(:,1+hga%j:njm+hga%j) = &
      0.5 * ( var(:,1:njm) + var(:,2:nj) )
  else if ( associated(hga%di, hgb) ) then
    ans(1+hga%i:nim+hga%i,1+hga%j:njm+hga%j) = &
      ( var(1:nim,1:njm) + var(1:nim,2:nj) + &
        var(2:ni, 1:njm) + var(2:ni, 2:nj) ) * 0.25
  else
    print *, 'unhandled relative position in ter_r2d of mod_op'
    stop
  end if

  call mympi_swpbnd (ans)

end subroutine ter_r2d

subroutine ter_gm2d (ans, var) !{{{1
  ! interpolate var to ans's grid
  type (type_gvar_m2d) :: ans
  type (type_gvar_m2d), intent(in) :: var

  call ter_r2d( ans%x(1)%v, var%x(1)%v, var%x(1)%hg, ans%x(1)%hg )
  call ter_r2d( ans%x(2)%v, var%x(2)%v, var%x(2)%hg, ans%x(2)%hg )

end subroutine ter_gm2d

subroutine ter_i3d(ans, var, hga, hgb, dft) !{{{1
  ! interpolate from grid hga to grid hgb
  integer :: ans(:,:,:)
  integer, intent(in) :: var(:,:,:)
  type (type_stg), target :: hga, hgb
  integer, optional :: dft

  if ( hga%n == hgb%n ) stop 'no need to interpolate in ans of mod_op'

  if ( present(dft) ) then
    ans = dft
  else
    ans = 0
  end if

  if      ( associated(hga%ew, hgb) ) then
    ans(1+hga%i:nim+hga%i,:,:) = &
    0.5 * ( var(1:nim,:,:) + var(2:ni,:,:) )
  else if ( associated(hga%ns, hgb) ) then
    ans(:,1+hga%j:njm+hga%j,:) = &
      0.5 * ( var(:,1:njm,:) + var(:,2:nj,:) )
  else if ( associated(hga%di, hgb) ) then
    ans(1+hga%i:nim+hga%i,1+hga%j:njm+hga%j,:) = &
      ( var(1:nim,1:njm,:) + var(1:nim,2:nj,:) + &
        var(2:ni, 1:njm,:) + var(2:ni, 2:nj,:) ) * 0.25
  else
    print *, 'unhandled relative position in ter_i3d of mod_op'
    stop
  end if

  call mympi_swpbnd (ans)

end subroutine ter_i3d

subroutine ter_i2d (ans, var, hga, hgb, dft) !{{{1
  ! interpolate from grid hga to grid hgb
  integer, dimension(ni,nj) :: ans
  integer, intent(in) :: var(:,:)
  type (type_stg), target :: hga, hgb
  integer, optional :: dft

  if ( hga%n == hgb%n ) stop 'no need to interpolate in ans of mod_op'

  if ( present(dft) ) then
    ans = dft
  else
    ans = 0
  end if

  if      ( associated(hga%ew, hgb) ) then
    ans(1+hga%i:nim+hga%i,:) = &
    0.5 * ( var(1:nim,:) + var(2:ni,:) )
  else if ( associated(hga%ns, hgb) ) then
    ans(:,1+hga%j:njm+hga%j) = &
      0.5 * ( var(:,1:njm) + var(:,2:nj) )
  else if ( associated(hga%di, hgb) ) then
    ans(1+hga%i:nim+hga%i,1+hga%j:njm+hga%j) = &
      ( var(1:nim,1:njm) + var(1:nim,2:nj) + &
        var(2:ni, 1:njm) + var(2:ni, 2:nj) ) * 0.25
  else
    print *, 'unhandled relative position in ter_i2d of mod_op'
  end if

  call mympi_swpbnd (ans)

end subroutine ter_i2d

subroutine op_ter_t2u( va, vb ) !{{{1
  ! interpolate from grid 1 to grid 3
  real (kind=wp), intent(in) :: va(:,:)
  real (kind=wp):: vb(:,:)

  integer :: i, j

  do i = 1, nim
  do j = 1, njm
    if ( g3%lev(i,j) > 0 ) vb(i,j) = &
      ( va(i,j)   + va(i,j+1) + &
        va(i+1,j) + va(i+1,j+1) ) * 0.25
  end do
  end do

  call mympi_swpbnd (vb)

end subroutine op_ter_t2u

subroutine dx1_r3d ( ans, var, ga, gb )  !{{{1
  ! var(x+dx) - var(x)
  ! var on grid ga, result on grid pb
  real (kind=wp) :: ans(:,:,:)
  real (kind=wp), intent(in) :: var(:,:,:)
  type (type_stg), target :: ga, gb
  real (kind=wp), allocatable, dimension(:,:,:) :: temp

  integer :: d3, is

  d3 = size(var, 3)
  allocate(temp(ni,nj,d3), stat = is); call chk(is)

  ! interpolate to the matching grid first
  if ( .not. associated(gb%ew, ga) ) then
    call ter_r3d(temp, var, ga, gb%ew )
  else
    temp = var
  end if

  ans = 0.0

  ans(1+ga%i:nim+ga%i,:,:) = &
    temp(2:ni,:,:) - temp(1:nim,:,:)

  call mympi_swpbnd( ans)

  deallocate(temp)

end subroutine dx1_r3d

subroutine dx2_r3d( ans, var, ga, gb )  !{{{1
  ! var(y+dy) - var(y)
  ! var on grid ga, result on grid pb
  real (kind=wp) :: ans(:,:,:)
  real (kind=wp), intent(in) :: var(:,:,:)
  type (type_stg), target :: ga, gb
  real (kind=wp), allocatable, dimension(:,:,:) :: temp

  integer :: d3, is

  d3 = size(var, 3)
  allocate(temp(ni,nj,d3), stat = is); call chk(is)

  ! interpolate to the matching grid first
  if ( associated(gb%ns, ga) ) then
    temp = var
  else
    call ter_r3d( temp, var, ga, gb%ns )
  end if

  ans = 0.0

  ans(:,1+ga%j:njm+ga%j,:) = temp(:,2:nj,:) - temp(:,1:njm,:)

  call mympi_swpbnd( ans)
    
  deallocate(temp)
end subroutine dx2_r3d

subroutine p1_r3d ( ans, var, ga, gb )  !{{{1
  ! (partial var) / (partial x1), not include metric effects
  ! var on grid ga, result on grid pb
  real (kind=wp) :: ans(:,:,:)
  real (kind=wp), intent(in) :: var(:,:,:)
  type (type_stg), target :: ga, gb
  real (kind=wp), allocatable, dimension(:,:,:) :: temp

  integer :: d3, is

  d3 = size(var, 3)
  allocate(temp(ni,nj,d3), stat = is); call chk(is)

  ! interpolate to the matching grid first
  if ( .not. associated(gb%ew, ga) ) then
    call ter_r3d(temp, var, ga, gb%ew )
  else
    temp = var
  end if

  ans = 0.0

  ans(1+ga%i:nim+ga%i,:,:) = &
    temp(2:ni,:,:) - temp(1:nim,:,:)

  ans = ans / spread(gb%dx%x(1),3,d3)

  call mympi_swpbnd(ans)

  deallocate(temp)

end subroutine p1_r3d

subroutine px_r3d ( ans, var, hga, hgb )  !{{{1
  ! physical derivative, (partial var) / (partial x)
  ! var on grid hga, result on grid pb
  real (kind=wp) :: ans(:,:,:)
  real (kind=wp), intent(in) :: var(:,:,:)
  type (type_stg), target :: hga, hgb

  integer :: d3, k, ia, ib

  d3 = size(var, 3)

  if ( associated(hgb%ew, hga) ) then
    ans = 0.0
    do k = 1, d3
      ia = 1 + hga%i
      ib = nim + hga%i
      ans(ia:ib,:,k) = (var(2:ni,:,k) - var(1:nim,:,k)) / &
                       (a * hgb%rh1(ia:ib,:) * hgb%dx(ia:ib,:)%x(1))
    end do
  else
    stop "grid unmatched in px_r3d of mod_op"
  end if

  call mympi_swpbnd(ans)

end subroutine px_r3d

subroutine px_r3d_center ( ans, var, hg )  !{{{1
  ! physical derivative, (partial var) / (partial x)
  ! var on grid hgu, result also on grid hgu
  ! center difference scheme
  real (kind=wp) :: ans(:,:,:)
  real (kind=wp), intent(in) :: var(:,:,:)
  type (type_stg), target :: hg

  real (kind=wp) :: dx
  integer :: d3, k, i, j

  d3 = size(var, 3)

  if ( associated(hgu, hg) ) then
    ans = 0.0
    do k = 1, d3
    do j = 1, nj
    do i = 2, ni - 1
      dx = a * hg%ew%rh1(i,j) * &
        ( hg%ew%dx(i,j)%x(1) + hg%ew%dx(i+1,j)%x(1) )
      ans(i,j,k) = ( var(i+1,j,k) - var(i-1,j,k) ) / dx
    end do
    end do
    end do
  else
    stop "grid is not proper in px_r3d_center"
  end if

  call mympi_swpbnd(ans)

end subroutine px_r3d_center

subroutine p1_r3d_b ( ans, var, ga, gb, mask )  !{{{1
  ! (partial var) / (partial x1), not include metric effects
  ! var on grid ga, result on grid pb
  ! set to zero if any of the two points is missing
  real (kind=wp) :: ans(:,:,:)
  real (kind=wp), dimension(:,:,:), intent(in) :: var
  integer, dimension(:,:,:), intent(in) :: mask
  type (type_stg), target :: ga, gb
  real (kind=wp), allocatable, dimension(:,:,:) :: temp

  integer :: d3, is

  d3 = size(var, 3)
  allocate(temp(ni,nj,d3), stat = is); call chk(is)

  ! interpolate to the matching grid first
  if ( .not. associated(gb%ew, ga) ) then
    call ter_r3d(temp, var, ga, gb%ew )
  else
    temp = var
  end if

  ans = 0.0

  ans(1+ga%i:nim+ga%i,:,:) = &
    ( temp(2:ni,:,:) - temp(1:nim,:,:) ) * &
    mask(2:ni,:,:)*mask(1:nim,:,:)

  ans = ans / spread(gb%dx%x(1),3,d3)

  call mympi_swpbnd( ans)

  deallocate(temp)

end subroutine p1_r3d_b

subroutine p1_r2d ( ans, var, ga, gb )  !{{{1
  ! (partial var) / (partial x1), not include metric effects
  ! var on grid ga, result on grid pb
  real (kind=wp), dimension(ni,nj) :: ans
  real (kind=wp), intent(in) :: var(:,:)
  type (type_stg), target :: ga, gb

  real (kind=wp), dimension(ni,nj) :: temp

  ! interpolate to the matching grid first
  if ( .not. associated(gb%ew, ga) ) then
    call ter_r2d( temp, var, ga, gb%ew )
  else
    temp = var
  end if

  ans = 0.0

  ans(1+ga%i:nim+ga%i,:) = temp(2:ni,:) - temp(1:nim,:)

  ans = ans / gb%dx%x(1)

  call mympi_swpbnd( ans)
end subroutine p1_r2d

subroutine p2_r3d( ans, var, ga, gb )  !{{{1
  ! (partial var) / (partial x2), not include metric effects
  ! var on grid ga, result on grid pb
  real (kind=wp) :: ans(:,:,:)
  real (kind=wp), dimension(:,:,:), intent(in) :: var
  type (type_stg), target :: ga, gb
  real (kind=wp), allocatable, dimension(:,:,:) :: temp

  integer :: d3, is

  d3 = size(var, 3)
  allocate(temp(ni,nj,d3), stat = is); call chk(is)

  ! interpolate to the matching grid first
  if ( associated(gb%ns, ga) ) then
    temp = var
  else
    call ter_r3d( temp, var, ga, gb%ns )
  end if

  ans = 0.0

  ans(:,1+ga%j:njm+ga%j,:) = &
    ( temp(:,2:nj,:) - temp(:,1:njm,:) )
    
  ans = ans / spread(gb%dx%x(2),3,d3)

  call mympi_swpbnd( ans)
    
  deallocate(temp)
end subroutine p2_r3d

subroutine p2_r3d_b( ans, var, ga, gb, mask )  !{{{1
  ! (partial var) / (partial x2), not include metric effects
  ! var on grid ga, result on grid pb
  ! set to zero if any of the tow points is missing
  real (kind=wp) :: ans(:,:,:)
  real (kind=wp), dimension(:,:,:), intent(in) :: var
  integer, dimension(:,:,:), intent(in) :: mask
  type (type_stg), target :: ga, gb
  real (kind=wp), allocatable, dimension(:,:,:) :: temp

  integer :: d3, is

  d3 = size(var, 3)
  allocate(temp(ni,nj,d3), stat = is); call chk(is)

  ! interpolate to the matching grid first
  if ( associated(gb%ns, ga) ) then
    temp = var
  else
    call ter_r3d( temp, var, ga, gb%ns )
  end if

  ans = 0.0

  ans(:,1+ga%j:njm+ga%j,:) = &
    ( temp(:,2:nj,:) - temp(:,1:njm,:) ) * &
    mask(:,2:nj,:)*mask(:,1:njm,:)
    
  ans = ans / spread(gb%dx%x(2),3,d3)

  call mympi_swpbnd( ans)
    
  deallocate(temp)
end subroutine p2_r3d_b

subroutine p2_r2d ( ans, var, ga, gb )  !{{{1
  ! (partial var) / (partial x2), not include metric effects
  ! var on grid ga, result on grid pb
  real (kind=wp), dimension(ni,nj) :: ans
  real (kind=wp), intent(in) :: var(:,:)
  type (type_stg), target :: ga, gb

  real (kind=wp), dimension(ni,nj) :: temp

  ! interpolate to the matching grid first
  if ( .not. associated(gb%ns, ga) ) then
    call ter_r2d( temp, var, ga, gb%ns )
  else
    temp = var
  end if

  ans = 0.0

  ans(:,1+ga%j:njm+ga%j) = temp(:,2:nj) - temp(:,1:njm)

  ans = ans / gb%dx%x(2)

  call mympi_swpbnd( ans)

end subroutine p2_r2d

subroutine p3_r3d ( ans, var )  !{{{1
  ! default (partial var) / (partial x3)
  ! upward is positive
  ! vertically, var on vg1, result on vg2
  real (kind=wp), dimension(ni,nj,nk) :: ans
  real (kind=wp), dimension(ni,nj,nkp) :: var

  integer :: d3

  d3 = size(var, 3)

  if (d3 == nkp) then
    ans = 0.0
    ans(:,:,1:nk) = var(:,:,1:nk) - var(:,:,2:nkp)
    ans = ans / spread(spread(vg2%dp,1,nj), 1, ni)
  else
    stop 'var should vertically on vg1 in p3_r3d in mod_op'
  end if

  ! no need to swap boundary horizontally
end subroutine p3_r3d

subroutine p3_r3d_b ( ans, var )  !{{{1
  ! default (partial var) / (partial x3)
  ! upward is positive
  ! vertically, var on vg2, result on vg1
  real (kind=wp), dimension(ni,nj,nkp) :: ans
  real (kind=wp), dimension(ni,nj,nk) :: var

  integer :: d3

  d3 = size(var, 3)

  if (d3 == nk) then
    ans = 0.0
    ans(:,:,2:nk) = var(:,:,1:nk-1) - var(:,:,2:nk)
    ans = ans / spread(spread(vg1%dp,1,nj), 1, ni)
  else
    stop 'var should vertically on vg2 in p3_r3d_b in mod_op'
  end if

  ! no need to swap boundary horizontally
end subroutine p3_r3d_b

subroutine p3_r3d_ter( ans, var, hga, hgb )  !{{{1
  ! (partial var) / (partial x3)
  ! upward is positive
  ! var horizontally on hga, vertically on vg1
  ! result horizontally in hgb, vertically on vg2
  real (kind=wp), dimension(ni,nj,nk) :: ans
  real (kind=wp), dimension(ni,nj,nkp) :: var
  type (type_stg), target :: hga, hgb

  real (kind=wp), dimension(ni,nj,nkp) :: temp
  integer :: d3

  d3 = size(var, 3)
  if (d3 /= nkp) &
    stop 'var should vertically on vg1 in ans in mod_op'

  call ter_r3d( temp, var, hga, hgb )

  ans = 0.0
  ans(:,:,1:nk) = temp(:,:,1:nk) - temp(:,:,2:nkp)
  ans = ans / spread(spread(vg2%dp,1,nj), 1, ni)

  ! no need to swap boundary horizontally
end subroutine p3_r3d_ter

subroutine pz_r3d ( ans, var, msk )  !{{{1
  ! (partial var) / (partial z), with respect to height, not pressure
  ! upward is positive
  ! vertically, var on vg2, result on vg1
  ! set to zero if any of the upper and lower layer of grid mask is land
  real (kind=wp), dimension(ni,nj,nkp) :: ans
  real (kind=wp), dimension(ni,nj,nk), intent(in) :: var
  integer, dimension(ni,nj,nk), intent(in) :: msk

  integer :: d3, i, j, k

  d3 = size(var, 3)

  if (d3 == nk) then
    ans = 0.0
    do k = 2, nk
    do j = 1, nj
    do i = 1, ni
      ans(i,j,k) = ( var(i,j,k-1) - var(i,j,k) ) * &
        msk(i,j,k-1)*msk(i,j,k) / vgw%dz(k)
    end do
    end do
    end do
  else
    stop 'var should vertically on vg2 in pz_r3d in mod_op'
  end if

  ! no need to swap boundary horizontally
end subroutine pz_r3d

subroutine dx3_r3d ( ans, var )  !{{{1
  ! default (partial var) / (partial x3)
  ! upward is positive
  ! vertically, var on vg1, result on vg2
  real (kind=wp), dimension(ni,nj,nk) :: ans
  real (kind=wp), dimension(ni,nj,nkp) :: var
  type (type_stg), target :: hga, hgb

  integer :: d3

  d3 = size(var, 3)

  if (d3 == nkp) then
    ans = 0.0
    ans(:,:,1:nk) = var(:,:,1:nk) - var(:,:,2:nkp)
  else
    stop 'var should vertically on vg1 in dx3_r3d in mod_op'
  end if

  ! no need to swap boundary horizontally
end subroutine dx3_r3d

subroutine dx3_r3d_b ( ans, var )  !{{{1
  ! default (partial var) / (partial x3)
  ! upward is positive
  ! vertically, var on vg2, result on vg1
  real (kind=wp), dimension(ni,nj,nkp) :: ans
  real (kind=wp), dimension(ni,nj,nk) :: var
  type (type_stg), target :: hga, hgb

  integer :: d3

  d3 = size(var, 3)

  if (d3 == nk) then
    ans = 0.0
    ans(:,:,2:nk) = var(:,:,1:nk-1) - var(:,:,2:nk)
  else
    stop 'var should vertically on vg2 in dx3_r3d_b in mod_op'
  end if

  ! no need to swap boundary horizontally
end subroutine dx3_r3d_b

subroutine div_r3d (ans, va, vb, ga, gb, gc) !{{{1
  ! horizontal divergence operator for 3d vector (va, vb)
  ! (va, vb) on grid (ga, gb), result on grid gc
  real (kind=wp), dimension(ni,nj,nk) :: ans
  real (kind=wp), dimension(ni,nj,nk), intent(in) :: va, vb
  type (type_stg), intent(in) :: ga, gb, gc

  real (kind=wp), dimension(ni,nj,nk) :: wka, wkb, wk

  wk = spread(gb%rh1,3,nk)
  call p1_r3d( wka, va, ga, gc )
  call p2_r3d( wkb, vb*wk, gb, gc)
  ans = wka + wkb

  wk = spread(gc%rh1,3,nk)
  ans = ans / (a * wk)

end subroutine div_r3d

subroutine div_r2d (ans, va, vb, ga, gb, gc) !{{{1
  ! horizontal divergence operator for 2d vector (va, vb)
  ! (va, vb) on grid (ga, gb), result on grid gc
  real (kind=wp), dimension(ni,nj) :: ans
  real (kind=wp), dimension(ni,nj), intent(in) :: va, vb
  type (type_stg), intent(in) :: ga, gb, gc

  real (kind=wp), dimension(ni,nj) :: wk

  call p1_r2d( wk, va, ga, gc )
  ans = wk
  call p2_r2d( wk, vb*gb%rh1, gb, gc )
  ans = ans + wk

  ans = ans / (a * gc%rh1)

end subroutine div_r2d

subroutine gra_r2d (ans, var, hga, hgb, hgc) !{{{1
  ! horizontal gradient operator
  ! var on grid hga, output on (hgb, hgc)
  type (type_mat) :: ans(ni,nj)
  real (kind=wp),  dimension(ni,nj), intent(in) :: var
  type (type_stg), target :: hga, hgb, hgc

  call p1_r2d( ans%x(1), var, hga, hgb )
  ans%x(1) = ans%x(1) / ( a * hgb%rh1 )
  call p2_r2d( ans%x(2), var, hga, hgc )
  ans%x(2) = ans%x(2) / a

end subroutine gra_r2d

subroutine gra_r2d_b (ans, var, hga, hgb, hgc) !{{{1
  ! horizontal gradient operator
  ! var on grid hga, output on (hgb, hgc)
  type (type_gvar_m2d) :: ans
  real (kind=wp),  dimension(ni,nj), intent(in) :: var
  type (type_stg), target :: hga, hgb, hgc

  ans%x(1)%hg => hgb; ans%x(2)%hg => hgc

  call p1_r2d( ans%x(1)%v, var, hga, hgb )
  ans%x(1)%v = ans%x(1)%v / ( a * hgb%rh1 )
  call p2_r2d( ans%x(2)%v, var, hga, hgc )
  ans%x(2)%v = ans%x(2)%v / a

end subroutine gra_r2d_b

subroutine gra_r3d (ans, var, ga, gb, gc) !{{{1
  ! horizontal gradient operator for 3d scalar
  ! var on grid ga, output on (gb, gc)
  type (type_mat) :: ans(:,:,:)
  real (kind=wp),  dimension(:,:,:), intent(in) :: var
  type (type_stg) :: ga, gb, gc

  integer :: k, d3, is

  call p1_r3d( ans%x(1), var, ga, gb)
  ans%x(1) = ans%x(1) / ( a * spread(gb%rh1,3,d3) )
  call p2_r3d( ans%x(2), var, ga, gc)
  ans%x(2) = ans%x(2) / a

end subroutine gra_r3d

subroutine gra_gr2d (ans, var) !{{{1
  ! gradient of var
  type (type_gvar_m2d) :: ans
  type (type_gvar_r2d), intent(in) :: var

  ans%x(1)%hg => var%hg%ew
  call p1_r2d( ans%x(1)%v, var%v, var%hg, ans%x(1)%hg )
  ans%x(1)%v = ans%x(1)%v / ( a * ans%x(1)%hg%rh1 )

  ans%x(2)%hg => var%hg%ns
  call p2_r2d( ans%x(2)%v, var%v, var%hg, ans%x(2)%hg )
  ans%x(2)%v = ans%x(2)%v / a

end subroutine gra_gr2d

subroutine gra_gr3d (ans, var) !{{{1
  ! gradient of var
  type (type_gvar_m3d) :: ans
  type (type_gvar_r3d), intent(in) :: var

  integer :: d3

  d3 = size(var%v, 3)

  ans%x(1)%g => var%g%ew
  call p1_r3d( ans%x(1)%v, var%v, var%g%hg, ans%x(1)%g%hg )
  ans%x(1)%v = ans%x(1)%v / ( a * spread(ans%x(1)%g%hg%rh1, 3, d3) )

  ans%x(2)%g => var%g%ns
  call p2_r3d( ans%x(2)%v, var%v, var%g%hg, ans%x(2)%g%hg )
  ans%x(2)%v = ans%x(2)%v / a

end subroutine gra_gr3d
subroutine lap_r2d (ans, var, ga, gb) !{{{1
  ! horizontal Laplacian operator
  ! var on grid ga, output on grid gb
  real (kind=wp), dimension(ni,nj) :: ans
  real (kind=wp), dimension(:,:), intent(in) :: var
  type (type_stg), intent(in) :: ga, gb

  real (kind=wp), dimension(ni,nj) :: wka, wkb
  
  call p1_r2d( wka, var, ga, gb%ew )
  wka = wka / gb%ew%rh1 
  call p1_r2d( wkb, wka, gb%ew, gb )
  ans = wkb

  call p2_r2d( wka, var, ga, gb%ns )
  wka = wka * gb%ns%rh1
  call p2_r2d( wkb, wka, gb%ns, gb )
  ans = ans + wkb

  ans = ans / ( a*a*gb%rh1 )

end subroutine lap_r2d

subroutine lap_gm2d (ans, var) !{{{1
  ! Laplacian operator
  type (type_gvar_m2d) :: ans
  type (type_gvar_m2d), intent(in) :: var

  ans%x(1)%hg => var%x(1)%hg
  call lap_r2d( ans%x(1)%v, var%x(1)%v, var%x(1)%hg, ans%x(1)%hg )

  ans%x(2)%hg => var%x(2)%hg
  call lap_r2d( ans%x(2)%v, var%x(2)%v, var%x(2)%hg, ans%x(2)%hg )
end subroutine lap_gm2d

subroutine int_tobot (ans, var, grd)!{{{1
  ! indefinite integration from current depth to the sea bottom
  ! assuming var is on the center of the layer
  real (kind=wp), dimension(ni,nj,nk) :: ans
  real (kind=wp), intent(in) :: var(:,:,:)
  type (type_stg3d), intent(in) :: grd

  real (kind=wp), dimension(ni,nj,nk) :: temp
  integer :: k

  temp = var * spread(spread(grd%vg%dp(1:nk),1,nj), 1, ni) * grd%msk

  do k = 1, nk-1
    ans(:,:,k) = temp(:,:,k)*0.5 + sum(temp(:,:,k+1:nk), 3)
  end do

  ans(:,:,nk) = temp(:,:,nk)*0.5
end subroutine int_tobot

subroutine op_vint_ns (var, ren, res)!{{{1
  ! vertical integration from sea surface to sea bottom
  ! var is horizontally on grid 2, vertically on vg2
  ! result on grid 3 (north and south)
  real (kind=wp), intent(in) :: var(:,:,:)
  real (kind=wp), dimension(:,:) :: ren, res

  real (kind=wp), dimension(ni,nj,nk) :: dp3d

  dp3d = spread( spread(g3%vg%dp,1,nj), 1, ni )
  ! divide by an integration facter pb
  ren = sum(var*dp3d*g3%msk, 3) / g3%pb
  ! shift mask and pb one grid northwards
  res = sum(var*dp3d*cshift(g3%msk,-1,2), 3) / cshift(g3%pb,-1,2)

  call mympi_swpbnd (ren, res) 

end subroutine op_vint_ns

subroutine op_vint_ew (var, ree, rew)!{{{1
  ! vertical integration from sea surface to sea bottom
  ! var is horizontally on grid 4, vertically on vg2
  ! result on grid 3 (east and west)
  real (kind=wp), intent(in) :: var(:,:,:)
  real (kind=wp), dimension(:,:) :: ree, rew

  real (kind=wp), dimension(ni,nj,nk) :: dp3d

  dp3d = spread( spread(g3%vg%dp,1,nj), 1, ni )
  ! divide by an integration facter pb
  ree = sum(var*dp3d*g3%msk, 3) / g3%pb
  ! shift mask and pb one grid northwards
  rew = sum(var*dp3d*cshift(g3%msk,-1,1), 3) / cshift(g3%pb,-1,1)

  call mympi_swpbnd (ree, rew) 

end subroutine op_vint_ew

subroutine vint_r3d (ans, var, grd)!{{{1
  ! vertical integration from sea surface to sea bottom
  ! result divided by a factor of pb
  real (kind=wp), dimension(ni,nj) :: ans
  real (kind=wp), intent(in) :: var(:,:,:)
  type (type_stg3d) :: grd

  real (kind=wp), dimension(ni,nj,nk) :: dp3d

  dp3d = spread( spread(grd%vg%dp(1:nk),1,nj), 1, ni )
  ans = sum(var * dp3d * grd%msk, 3) / grd%pb

end subroutine vint_r3d

subroutine chk( ista ) !{{{1
  ! check state of allocate array 

  integer, intent(in) ::  ista

  if ( ista /= 0 ) then
    write(*,*) 'Allocate array failed. Stop'
    stop 2
  end if
end subroutine chk

end module mod_op!{{{1
!-------------------------------------------------------{{{1
! vim:fdm=marker:fdl=0:
! vim:foldtext=getline(v\:foldstart).'...'.(v\:foldend-v\:foldstart):
