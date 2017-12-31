
! Description: operators
!
!      Author: OU Yuyuan <ouyuyuan@lasg.iap.ac.cn>
!     Created: 2015-10-11 15:23:38 BJT
! Last Change: 2017-12-31 09:47:48 BJT

module mod_op

  ! imported variables !{{{1
  use mod_debug, only: debug_var

  use mod_arrays, only: &
    g32, gi1, gi2, &
    gt, gu, gtj, guj, git, giw, &
    am, km, kh, wm, &
    cv1, cv2, &
    glo_lat, glo_lon, z, &
    equv, eqch

  use mod_con, only: &
    a, g, rho0, cdbot, ah_c

  use mod_kind, only: wp, one

  use mod_mympi, only: mympi_swpbnd

  use mod_param, only: &
    ni, nj, nk, nkp, &
    tc, tp, nim, njm

  use mod_type, only: type_mat, &
    tctr, type_gi, type_gj, type_gij, &
    type_gvar_r2d, type_gvar_r3d

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
    module procedure gra_r3d_b
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

  interface op_ter
    module procedure vter_r3d
    module procedure ter_r3d
    module procedure ter_r2d
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

subroutine fri_r3d ( ans, sch, v, vp, vc, tau ) !{{{1
  ! horizontal frictional force
  ! ans is on the same grid as vp
  ! ans = 1/sch*div(ch*gra(v)) + cv1*vp + cv2*sch*(p vc/ p x1) + vert.
  real (kind=wp), dimension(:,:,:) :: ans
  real (kind=wp), dimension(:,:), intent(in) :: sch
  real (kind=wp), dimension(:,:,:), intent(in) :: &
    v, & ! unweighted velocity 
    vp, & ! pressure weighted velocity
    vc
  real (kind=wp), dimension(ni,nj), intent(in) :: tau

  type (type_mat), dimension(ni,nj,nk) :: wkmat
  real (kind=wp), dimension(ni,nj,nkp) :: wkp
  real (kind=wp), dimension(ni,nj,nk) :: wka, sch3d, ch3da, ch3db
  real (kind=wp), dimension(ni,nj) :: wkr2d
  real (kind=wp) :: mag
  integer :: i, j, k

  sch3d = spread(sch, 3, nk)

  ! horizontal viscosity

  ! div(ch*gra(v))
  call op_gra( wkmat, v, guj, guj%ew, guj%ns )
  call op_ter( wkr2d, eqch%chc, eqch%hg, guj%ew )
  ch3da = spread( wkr2d, 3, nk )
  call op_ter( wkr2d, eqch%chc, eqch%hg, guj%ns )
  ch3db = spread( wkr2d, 3, nk )
  call div_r3d( wka, ch3da*wkmat%x(1), ch3db*wkmat%x(2), &
                guj%ew, guj%ns, guj )
  ans = am%v/sch3d * wka

  ! cv1*vp
  ans = ans + am%v*spread(cv1,3,nk)*vp

  ! cv2*sch*(p vc/ p x1)
  ! oddly center difference of vc
  call px_r3d_center( wka, vc, guj )
  ans = ans + am%v*spread(cv2,3,nk)*sch3d*wka

  ! vertical viscosity

  ! oddly vertical difference respect to thickness (not pressure)
  call pz_r3d( wkp, vp, gu%msk )
  wkp(:,:,1:nk) = g*rho0*km%v * wkp(:,:,1:nk) / (spread(sch,3,nk))**2

  ! surface boundary condition
  wkp(:,:,1) = g*tau(:,:) / sch(:,:)

  do j = 1, nj
  do i = 1, ni
    k = gu%lev(i,j)
    if ( k > 0 ) then ! bottom drag
      mag = sqrt( v(i,j,k)**2 + vc(i,j,k)**2 )
      wkp(i,j,k+1) = v(i,j,k)*g*cdbot*mag/sch(i,j)
    end if
  end do
  end do

  call p3_r3d( wka, wkp )

  ans = ans*gu%msk + wka

end subroutine fri_r3d

subroutine op_fri(fx, fy, sch, taux, tauy) !{{{1
  ! horizontal frictional force
  ! output is on the same grid as equv
  real (kind=wp), dimension(:,:,:), intent(in) :: fx, fy
  real (kind=wp), dimension(:,:), intent(in) :: sch, taux, tauy

  real (kind=wp), dimension(ni,nj,nk) :: u, v, sch3d

  sch3d = spread( sch, 3, nk )
  u = equv%uc / sch3d
  v = equv%vc / sch3d

  call fri_r3d( fx, sch, u, equv%uc, -v, taux )
  call fri_r3d( fy, sch, v, equv%vc,  u, tauy )

end subroutine op_fri

subroutine op_dif (ans, ch, var) !{{{1
  ! diffusion of tracers
  ! ans is on the same grid as var
  ! ans = div(ch*gra(var)) + vert.
  real (kind=wp), dimension(:,:,:) :: ans
  real (kind=wp), dimension(:,:), intent(in) :: ch
  real (kind=wp), dimension(:,:,:), intent(in) :: var ! tracer

  type (type_mat), dimension(ni,nj,nk) :: wkm
  real (kind=wp), dimension(ni,nj,nkp) :: wkp
  real (kind=wp), dimension(ni,nj,nk) :: wk, ch3da, ch3db
  real (kind=wp), dimension(ni,nj) :: wkr2d

  ! horizontal diffusion

  ! div(ch*gra(var))
  call p1_r3d_b( wkm%x(1), var, gtj, gtj%ew, gt%msk )
  wkm%x(1) = wkm%x(1) / ( a * spread(gtj%ew%rh, 3, nk) )
  call p2_r3d_b( wkm%x(2), var, gtj, gtj%ns, gt%msk )
  wkm%x(2) = wkm%x(2) / a
  call op_ter( wkr2d, ch, gtj, gtj%ew )
  ch3da = spread( wkr2d, 3, nk )
  call op_ter( wkr2d, ch, gtj, gtj%ns )
  ch3db = spread( wkr2d, 3, nk )
  call div_r3d( wk, ch3da*wkm%x(1), ch3db*wkm%x(2), &
                gtj%ew, gtj%ns, gtj )
  ans = ah_c * wk

  ! vertical diffusion

  ! oddly vertical difference respect to thickness (not pressure)
  call pz_r3d( wkp, var, gt%msk )
  wkp = g*rho0 * kh%v * wkp

  call p3_r3d( wk, wkp )
  ans = ans + wk
end subroutine op_dif

subroutine op_adv (ans, var, sch) !{{{1
  ! calc. advection of horizontal pressure weighted velocity var
  ! op_adv(var) = div(var*v) - 0.5*var*div(v), 
  !   v = (u, w) are unweighted 3d velocity
  real (kind=wp), dimension(:,:,:) :: ans
  real (kind=wp), dimension(:,:,:), intent(in) :: var
  real (kind=wp), dimension(:,:), intent(in) :: sch

  type (type_mat), dimension(ni,nj,nk) :: wkm
  real (kind=wp),  dimension(ni,nj,nk) :: wka, wkb, cos3d, sch3d, u, v
  real (kind=wp),  dimension(ni,nj,nkp):: wkap, wkbp, w

  ! unaveraged velocity, interpolate to proper grid
  sch3d = spread( sch, 3, nk )
  wka = equv%uc / sch3d
  call ter_r3d( u, wka, guj, guj%ew )

  cos3d = spread(guj%rh,3,nk)
  wka = equv%vc / sch3d
  call ter_r3d( v, wka*cos3d, guj, guj%ns )

  wkap = wm%v / spread(eqch%chc,3,nkp)
  call ter_r3d( w, wkap, gtj, guj )

  ! var*u
  call ter_r3d( wka, var, guj, guj%ew )
  call ter_r3d( wkb, var, guj, guj%ns )
  wkm%x(1) = wka * u
  wkm%x(2) = wkb * v

  ! horizontal: div(var*u)
  call p1_r3d( wka, wkm%x(1), guj%ew, guj )
  call p2_r3d( wkb, wkm%x(2), guj%ns, guj )
  ans = (wka + wkb) / (a*cos3d)

  ! var*w
  call vter_r3d( wkap, var, git, giw )
  wkap = wkap*w
  ! set to zero if any of the upper and lower layer of U-grid is land
  wkbp(:,:,1)   = 0
  wkbp(:,:,nkp) = 0
  wkbp(:,:,2:nkp-1) = wkap(:,:,2:nkp-1)*&
                      gu%msk(:,:,1:nk-1)*gu%msk(:,:,2:nk)

  call p3_r3d( wka, wkbp )
  ans = ans + wka

  ! div(v)
  call p1_r3d( wkm%x(1), u, guj%ew, guj )
  call p2_r3d( wkm%x(2), v, guj%ns, guj )
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
  real (kind=wp),  dimension(ni,nj) :: sch
  integer :: k

  call op_ter( sch, eqch%chc, eqch%hg, guj, one )
  sch = sqrt(sch)

  do k = 1, nk
    um(:,:,k)%x(1) = equv%uc(:,:,k) * sch
    um(:,:,k)%x(2) = equv%vc(:,:,k) * sch * guj%rh
  end do

  ! horizontal advection
  call p1_r3d( wkm%x(1), var, gtj, gtj%ew )
  call p2_r3d( wkm%x(2), var, gtj, gtj%ns )

  cos3d = spread(gtj%rh,3,nk)

  call ter_r3d( wka, um%x(1), guj, gtj%ew )
  call ter_r3d( wkb, wka * wkm%x(1), gtj%ew, gtj )
  ! not interpolate the cosine factor
  ans = wkb / (a*cos3d)

  call ter_r3d( wka, um%x(2), guj, gtj%ns )
  call ter_r3d( wkb, wka * wkm%x(2), gtj%ns, gtj )
  ! not interpolate the cosine factor
  ans = ans + wkb / (a*cos3d)

  ! vertical advection
  call dx3_r3d_b( wkp, var )
  call vter_r3d( wka, wm * wkp, giw, git ) 
  ! not interpolate the layer 'thickness'
  wka = wka / spread(spread(git%dpr,1,nj), 1, ni)
  ans = ans + wka

end subroutine op_adv_ts

subroutine vter_r3d(ans, var, vga, vgb, dft) !{{{1
  ! interpolate vertically from grid vga to grid vgb
  real (kind=wp) :: ans(:,:,:)
  real (kind=wp), intent(in) :: var(:,:,:)
  type (type_gj), target :: vga, vgb
  real (kind=wp), optional :: dft

  integer :: nda

  if ( vga%n == vgb%n ) &
    stop 'no need to interpolate in ans of mod_op'

  nda = size(vga%pr)

  if ( size(var,3) /= nda ) &
    stop 'var and vga unmatch in ans of mod_op'

  if ( present(dft) ) then
    ans = dft
  else
    ans = 0.0
  end if

  ! gi2 to gi1
  if ( vga%n == 2 ) then
    ans(:,:,2:nk) = ( var(:,:,1:nk-1) + var(:,:,2:nk) ) * 0.5
  ! from gi1 to gi2
  else
    ans(:,:,1:nk) = ( var(:,:,1:nk) + var(:,:,2:nkp) ) * 0.5
  end if

  ! no need to swap boundary horizontally
end subroutine vter_r3d

subroutine ter_r3d(ans, var, hga, hgb, dft) !{{{1
  ! interpolate from grid hga to grid hgb
  real (kind=wp) :: ans(:,:,:)
  real (kind=wp), intent(in) :: var(:,:,:)
  type (type_gi), target :: hga, hgb
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
  type (type_gi), target :: hga, hgb
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

subroutine ter_i3d(ans, var, hga, hgb, dft) !{{{1
  ! interpolate from grid hga to grid hgb
  integer :: ans(:,:,:)
  integer, intent(in) :: var(:,:,:)
  type (type_gi), target :: hga, hgb
  integer, optional :: dft

  if ( hga%n == hgb%n ) stop 'no need to interpolate in ans of mod_op'

  if ( present(dft) ) then
    ans = dft
  else
    ans = 0
  end if

  if      ( associated(hga%ew, hgb) ) then
    ans(1+hga%i:nim+hga%i,:,:) = &
     ( var(1:nim,:,:) + var(2:ni,:,:) ) / 2
  else if ( associated(hga%ns, hgb) ) then
    ans(:,1+hga%j:njm+hga%j,:) = &
       ( var(:,1:njm,:) + var(:,2:nj,:) ) / 2
  else if ( associated(hga%di, hgb) ) then
    ans(1+hga%i:nim+hga%i,1+hga%j:njm+hga%j,:) = &
      ( var(1:nim,1:njm,:) + var(1:nim,2:nj,:) + &
        var(2:ni, 1:njm,:) + var(2:ni, 2:nj,:) ) / 4
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
  type (type_gi), target :: hga, hgb
  integer, optional :: dft

  if ( hga%n == hgb%n ) stop 'no need to interpolate in ans of mod_op'

  if ( present(dft) ) then
    ans = dft
  else
    ans = 0
  end if

  if      ( associated(hga%ew, hgb) ) then
    ans(1+hga%i:nim+hga%i,:) = &
    ( var(1:nim,:) + var(2:ni,:) ) / 2
  else if ( associated(hga%ns, hgb) ) then
    ans(:,1+hga%j:njm+hga%j) = &
      ( var(:,1:njm) + var(:,2:nj) ) / 2
  else if ( associated(hga%di, hgb) ) then
    ans(1+hga%i:nim+hga%i,1+hga%j:njm+hga%j) = &
      ( var(1:nim,1:njm) + var(1:nim,2:nj) + &
        var(2:ni, 1:njm) + var(2:ni, 2:nj) ) / 4
  else
    print *, 'unhandled relative position in ter_i2d of mod_op'
  end if

  call mympi_swpbnd (ans)

end subroutine ter_i2d

subroutine p1_r3d ( ans, var, ga, gb )  !{{{1
  ! (partial var) / (partial x1), not include metric effects
  ! var on grid ga, result on grid gb
  real (kind=wp) :: ans(:,:,:)
  real (kind=wp), intent(in) :: var(:,:,:)
  type (type_gi), target :: ga, gb
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

subroutine px_r3d_center ( ans, var, hg )  !{{{1
  ! physical derivative, (partial var) / (partial x)
  ! var on grid guj, result also on grid guj
  ! center difference scheme
  real (kind=wp) :: ans(:,:,:)
  real (kind=wp), intent(in) :: var(:,:,:)
  type (type_gi), target :: hg

  real (kind=wp) :: dx
  integer :: d3, k, i, j

  d3 = size(var, 3)

  if ( associated(guj, hg) ) then
    ans = 0.0
    do k = 1, d3
    do j = 1, nj
    do i = 2, ni - 1
      dx = a * hg%ew%rh(i,j) * &
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
  ! var on grid ga, result on grid gb
  ! set to zero if any of the two points is missing
  real (kind=wp) :: ans(:,:,:)
  real (kind=wp), dimension(:,:,:), intent(in) :: var
  integer, dimension(:,:,:), intent(in) :: mask
  type (type_gi), target :: ga, gb
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
  ! var on grid ga, result on grid gb
  real (kind=wp), dimension(ni,nj) :: ans
  real (kind=wp), intent(in) :: var(:,:)
  type (type_gi), target :: ga, gb

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
  ! var on grid ga, result on grid gb
  real (kind=wp) :: ans(:,:,:)
  real (kind=wp), dimension(:,:,:), intent(in) :: var
  type (type_gi), target :: ga, gb
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
  ! var on grid ga, result on grid gb
  ! set to zero if any of the tow points is missing
  real (kind=wp) :: ans(:,:,:)
  real (kind=wp), dimension(:,:,:), intent(in) :: var
  integer, dimension(:,:,:), intent(in) :: mask
  type (type_gi), target :: ga, gb
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
  ! var on grid ga, result on grid gb
  real (kind=wp), dimension(ni,nj) :: ans
  real (kind=wp), intent(in) :: var(:,:)
  type (type_gi), target :: ga, gb

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

subroutine p3_r3d ( ans, var)  !{{{1
  ! default (partial var) / (partial x3)
  ! upward is positive
  ! vertically, var on gi1, result on gi2
  real (kind=wp), dimension(ni,nj,nk) :: ans
  real (kind=wp), dimension(ni,nj,nkp) :: var

  integer :: d3

  d3 = size(var, 3)

  if (d3 == nkp) then
    ans = 0.0
    ans(:,:,1:nk) = var(:,:,1:nk) - var(:,:,2:nkp)
    ans = ans / spread(spread(gi2%dpr,1,nj), 1, ni)
  else
    stop 'var should vertically on gi1 in p3_r3d in mod_op'
  end if

  ! no need to swap boundary horizontally
end subroutine p3_r3d

subroutine p3_r3d_ter( ans, var, hga, hgb )  !{{{1
  ! (partial var) / (partial x3)
  ! upward is positive
  ! var horizontally on hga, vertically on gi1
  ! result horizontally in hgb, vertically on gi2
  real (kind=wp), dimension(ni,nj,nk) :: ans
  real (kind=wp), dimension(ni,nj,nkp) :: var
  type (type_gi), target :: hga, hgb

  real (kind=wp), dimension(ni,nj,nkp) :: temp
  integer :: d3

  d3 = size(var, 3)
  if (d3 /= nkp) &
    stop 'var should vertically on gi1 in ans in mod_op'

  call ter_r3d( temp, var, hga, hgb )

  ans = 0.0
  ans(:,:,1:nk) = temp(:,:,1:nk) - temp(:,:,2:nkp)
  ans = ans / spread(spread(gi2%dpr,1,nj), 1, ni)

  ! no need to swap boundary horizontally
end subroutine p3_r3d_ter

subroutine pz_r3d ( ans, var, msk )  !{{{1
  ! (partial var) / (partial z), with respect to height, not pressure
  ! upward is positive
  ! vertically, var on gi2, result on gi1
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
        msk(i,j,k-1)*msk(i,j,k) / giw%dz(k)
    end do
    end do
    end do
  else
    stop 'var should vertically on gi2 in pz_r3d in mod_op'
  end if

  ! no need to swap boundary horizontally
end subroutine pz_r3d

subroutine dx3_r3d_b ( ans, var )  !{{{1
  ! default (partial var) / (partial x3)
  ! upward is positive
  ! vertically, var on gi2, result on gi1
  real (kind=wp), dimension(ni,nj,nkp) :: ans
  real (kind=wp), dimension(ni,nj,nk) :: var

  integer :: d3

  d3 = size(var, 3)

  if (d3 == nk) then
    ans = 0.0
    ans(:,:,2:nk) = var(:,:,1:nk-1) - var(:,:,2:nk)
  else
    stop 'var should vertically on gi2 in dx3_r3d_b in mod_op'
  end if

  ! no need to swap boundary horizontally
end subroutine dx3_r3d_b

subroutine div_r3d (ans, va, vb, ga, gb, gc) !{{{1
  ! horizontal divergence operator for 3d vector (va, vb)
  ! (va, vb) on grid (ga, gb), result on grid gc
  real (kind=wp), dimension(ni,nj,nk) :: ans
  real (kind=wp), dimension(ni,nj,nk), intent(in) :: va, vb
  type (type_gi), intent(in) :: ga, gb, gc

  real (kind=wp), dimension(ni,nj,nk) :: wka, wkb, wk

  wk = spread(gb%rh,3,nk)
  call p1_r3d( wka, va, ga, gc )
  call p2_r3d( wkb, vb*wk, gb, gc)
  ans = wka + wkb

  wk = spread(gc%rh,3,nk)
  ans = ans / (a * wk)

end subroutine div_r3d

subroutine div_r2d (ans, va, vb, ga, gb, gc) !{{{1
  ! horizontal divergence operator for 2d vector (va, vb)
  ! (va, vb) on grid (ga, gb), result on grid gc
  real (kind=wp), dimension(ni,nj) :: ans
  real (kind=wp), dimension(ni,nj), intent(in) :: va, vb
  type (type_gi), intent(in) :: ga, gb, gc

  real (kind=wp), dimension(ni,nj) :: wk

  call p1_r2d( wk, va, ga, gc )
  ans = wk
  call p2_r2d( wk, vb*gb%rh, gb, gc )
  ans = ans + wk

  ans = ans / (a * gc%rh)

end subroutine div_r2d

subroutine gra_r2d (ans, var, hga, hgb, hgc) !{{{1
  ! horizontal gradient operator
  ! var on grid hga, output on (hgb, hgc)
  type (type_mat) :: ans(ni,nj)
  real (kind=wp),  dimension(ni,nj), intent(in) :: var
  type (type_gi), target :: hga, hgb, hgc

  call p1_r2d( ans%x(1), var, hga, hgb )
  ans%x(1) = ans%x(1) / ( a * hgb%rh )
  call p2_r2d( ans%x(2), var, hga, hgc )
  ans%x(2) = ans%x(2) / a

end subroutine gra_r2d

subroutine gra_r2d_b (var, hg, grax, gray) !{{{1
  ! horizontal gradient operator
  ! var on grid hg, output on (hg%ew, hg%ns)
  real (kind=wp), dimension(ni,nj), intent(in) :: var
  type (type_gi), target :: hg
  real (kind=wp), dimension(ni,nj) :: grax, gray

  call p1_r2d( grax, var, hg, hg%ew )
  grax = grax / ( a * hg%ew%rh )
  call p2_r2d( gray, var, hg, hg%ns )
  gray = gray / a

end subroutine gra_r2d_b

subroutine gra_r3d_b (ansx, ansy, var, ga, gb, gc) !{{{1
  ! horizontal gradient operator for 3d scalar
  ! var on grid ga, output on (gb, gc)
  real (kind=wp),  dimension(:,:,:) :: ansx, ansy
  real (kind=wp),  dimension(:,:,:), intent(in) :: var
  type (type_gi) :: ga, gb, gc

  call p1_r3d( ansx, var, ga, gb)
  ansx = ansx / ( a * spread(gb%rh,3,size(var,3)) )
  call p2_r3d( ansy, var, ga, gc)
  ansy = ansy / a
end subroutine gra_r3d_b

subroutine gra_r3d (ans, var, ga, gb, gc) !{{{1
  ! horizontal gradient operator for 3d scalar
  ! var on grid ga, output on (gb, gc)
  type (type_mat) :: ans(:,:,:)
  real (kind=wp),  dimension(:,:,:), intent(in) :: var
  type (type_gi) :: ga, gb, gc

  call p1_r3d( ans%x(1), var, ga, gb)
  ans%x(1) = ans%x(1) / ( a * spread(gb%rh,3,size(var,3)) )
  call p2_r3d( ans%x(2), var, ga, gc)
  ans%x(2) = ans%x(2) / a
end subroutine gra_r3d

subroutine op_lap (ans, var, ga, gb) !{{{1
  ! horizontal Laplacian operator
  ! var on grid ga, output on grid gb
  real (kind=wp), dimension(ni,nj) :: ans
  real (kind=wp), dimension(:,:), intent(in) :: var
  type (type_gi), intent(in) :: ga, gb

  real (kind=wp), dimension(ni,nj) :: wka, wkb
  
  call p1_r2d( wka, var, ga, gb%ew )
  wka = wka / gb%ew%rh 
  call p1_r2d( wkb, wka, gb%ew, gb )
  ans = wkb

  call p2_r2d( wka, var, ga, gb%ns )
  wka = wka * gb%ns%rh
  call p2_r2d( wkb, wka, gb%ns, gb )
  ans = ans + wkb

  ans = ans / ( a*a*gb%rh )

end subroutine op_lap

subroutine int_tobot (ans, var, grd)!{{{1
  ! indefinite integration from current depth to the sea bottom
  ! assuming var is on the center of the layer
  real (kind=wp), dimension(ni,nj,nk) :: ans
  real (kind=wp), intent(in) :: var(:,:,:)
  type (type_gij), intent(in) :: grd

  real (kind=wp), dimension(ni,nj,nk) :: temp
  integer :: k

  temp = var * spread(spread(grd%vg%dpr(1:nk),1,nj), 1, ni) * grd%msk

  do k = 1, nk-1
    ans(:,:,k) = temp(:,:,k)*0.5 + sum(temp(:,:,k+1:nk), 3)
  end do

  ans(:,:,nk) = temp(:,:,nk)*0.5
end subroutine int_tobot

subroutine op_vint_ns (var, ren, res)!{{{1
  ! vertical integration from sea surface to sea bottom
  ! var is horizontally on g2j, vertically on gi2
  ! result on g3j (north and south)
  real (kind=wp), intent(in) :: var(:,:,:)
  real (kind=wp), dimension(:,:) :: ren, res

  real (kind=wp), dimension(ni,nj,nk) :: dpr3d

  dpr3d = spread( spread(gi2%dpr,1,nj), 1, ni )
  ! divide by an integration facter prh
  ren = sum(var*dpr3d*g32%msk, 3) / g32%prh
  ! shift mask and prh one grid northwards
  res = sum(var*dpr3d*cshift(g32%msk,-1,2), 3) / cshift(g32%prh,-1,2)

  call mympi_swpbnd (ren, res) 
end subroutine op_vint_ns

subroutine op_vint_ew (var, ree, rew)!{{{1
  ! vertical integration from sea surface to sea bottom
  ! var is horizontally on g4j, vertically on gi2
  ! result on g3j (east and west)
  real (kind=wp), intent(in) :: var(:,:,:)
  real (kind=wp), dimension(:,:) :: ree, rew

  real (kind=wp), dimension(ni,nj,nk) :: dpr3d

  dpr3d = spread( spread(g32%vg%dpr,1,nj), 1, ni )
  ! divide by an integration facter prh
  ree = sum(var*dpr3d*g32%msk, 3) / g32%prh
  ! shift mask and prh one grid northwards
  rew = sum(var*dpr3d*cshift(g32%msk,-1,1), 3) / cshift(g32%prh,-1,1)

  call mympi_swpbnd (ree, rew) 

end subroutine op_vint_ew

subroutine vint_r3d (ans, var, grd)!{{{1
  ! vertical integration from sea surface to sea bottom
  ! result divided by a factor of prh
  real (kind=wp), dimension(ni,nj) :: ans
  real (kind=wp), intent(in) :: var(:,:,:)
  type (type_gij) :: grd

  real (kind=wp), dimension(ni,nj,nk) :: dpr3d

  dpr3d = spread( spread(grd%vg%dpr(1:nk),1,nj), 1, ni )
  ans = sum(var * dpr3d * grd%msk, 3) / grd%prh

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
