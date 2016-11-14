
! Description: algebra rules of grid variables
!
!      Author: OU Yuyuan <ouyuyuan@lasg.iap.ac.cn>
!     Created: 2015-12-24 08:08:32 BJT
! Last Change: 2016-01-31 15:08:41 BJT

module mod_gvar

  ! imported variables !{{{1

  use mod_arrays, only: vg1, vg2, &
    pbt, cv1, cv2, &
    am, km, bnd, &
    glo_lon, glo_lat, z

  use mod_con, only: a, g, rho0, cdbot

  use mod_kind, only: wp

  use mod_mympi, only: mympi_quick_output

  use mod_param, only: &
    ni, nim, nj, njm, nk, nkp, tc, tp, &
    myid, mid

  use mod_op, only: &
    op_ter, op_div, op_int_tobot, &
    op_px1, op_px2, op_px3, op_lap

  use mod_type, only: &
    type_mat, &
    type_stg, &
    type_vstg, &
    type_stg3d, &
    type_gvar_r2d, &
    type_gvar_r3d, &
    type_gvar_m2d, &
    type_gvar_m3d, &
    type_add, &
    operator (/)

  implicit none
  private

!  public &
!    operator (+), &
!    operator (-), &
!    operator (*), &
!    operator (/), &
!    assignment (=)

  ! interfaces !{{{1

  interface gvar_shift_ew
    module procedure shift_ew_gr3d
    module procedure shift_ew_gr2d
  end interface

  interface gvar_shift_ns
    module procedure shift_ns_gr3d
    module procedure shift_ns_gr2d
  end interface

  interface gvar_lap
    module procedure lap_gr2d
    module procedure lap_gm2d
  end interface

  interface gvar_fri
    module procedure fri_gr3d
    module procedure fri_gm3d
  end interface

  interface gvar_adv
    module procedure adv_gr3d
    module procedure adv_gm3d
  end interface

  interface gvar_int_tobot
    module procedure int_tobot_gr3d
  end interface

  interface gvar_gra
    module procedure gra_gr3d
  end interface

  interface gvar_gra_old
    module procedure gra_r2d
    module procedure gra_gr2d
  end interface

  interface gvar_div
    module procedure div_gm2d_ter
    module procedure div_gm3d
    module procedure div_gm3d_ter
  end interface

  interface gvar_ter
    module procedure ter_gr2d
    module procedure ter_gr3d
    module procedure ter_gm2d
    module procedure ter_gm3d
  end interface

  interface operator(+)
    module procedure gr2d_plus_gr2d
    module procedure gm2d_plus_m2d
    module procedure gm2d_plus_gm2d
    module procedure gr3d_plus_gr3d
    module procedure gm3d_plus_gm3d
  end interface

  interface operator(-)
    module procedure r2d_minus_gr2d
    module procedure gr2d_minus_r2d
    module procedure gr3d_minus_r2d
    module procedure gr3d_minus_gr3d
    module procedure gm2d_minus_m2d
    module procedure gm3d_minus_m2d
  end interface

  interface operator(*)
    module procedure r_mul_gr2d
    module procedure r_mul_gr3d
    module procedure r_mul_gm2d
    module procedure r_mul_gm3d
    module procedure r2d_mul_gr2d
    module procedure r2d_mul_gr3d
    module procedure r3d_mul_gr3d
    module procedure r2d_mul_gm2d
    module procedure r2d_mul_gm3d
    module procedure gr2d_mul_gr2d
    module procedure gr2d_mul_gr3d
    module procedure gr3d_mul_gr3d
    module procedure gr2d_mul_gm3d
    module procedure gr3d_mul_gm3d
  end interface

  interface operator(/)
    module procedure gr2d_div_r
    module procedure gr2d_div_r2d
    module procedure gr3d_div_r
    module procedure gr3d_div_r2d
    module procedure gr3d_div_r3d
    module procedure gm2d_div_r
    module procedure gm2d_div_r2d
    module procedure gm3d_div_r2d
    module procedure gm3d_div_r3d
  end interface

  interface assignment(=)
    module procedure r2d_eq_gr2d
    module procedure r3d_eq_gr3d
    module procedure m2d_eq_gm2d
    module procedure m3d_eq_gm3d
    module procedure gr2d_eq_r
    module procedure gr2d_eq_r2d
    module procedure gr3d_eq_r
    module procedure gm2d_eq_r
    module procedure gm2d_eq_m2d
    module procedure gm3d_eq_r
  end interface

  interface gvar_assign
    module procedure m2d_assign_gm2d
    module procedure gm2d_assign_m2d
  end interface

  interface gvar_cp_shape
    module procedure cp_shape_gr2d
    module procedure cp_shape_gr2d_gm2d
    module procedure cp_shape_gr3d
    module procedure cp_shape_gr3d_gm3d
    module procedure cp_shape_gm2d
    module procedure cp_shape_gm3d
  end interface

  interface gvar_free
    module procedure free_gr2d
    module procedure free_gr3d
    module procedure free_gm2d
    module procedure free_gm3d
  end interface

contains !{{{1

function lap_gr2d (var) !{{{1
  ! default horizontal Laplacian operator
  type (type_gvar_r2d), intent(in) :: var
  type (type_gvar_r2d) :: lap_gr2d
  
  call cp_shape_gr2d( var, lap_gr2d )
  lap_gr2d%v = op_lap( var%v, var%hg, var%hg )

end function lap_gr2d

function lap_gm2d (var) !{{{1
  ! default horizontal Laplacian operator
  type (type_gvar_m2d), intent(in) :: var
  type (type_gvar_m2d) :: lap_gm2d
  
  lap_gm2d%x(1) = lap_gr2d( var%x(1) )
  lap_gm2d%x(2) = lap_gr2d( var%x(2) )

end function lap_gm2d

function fri_gr3d (spbt, va, vb, vc, tau) !{{{1
  ! horizontal frictional force
  ! output is on the same grid as vb
  real (kind=wp), dimension(:,:), intent(in) :: spbt
  type (type_gvar_r3d), intent(in) :: va, vb, vc
  real (kind=wp), dimension(ni,nj) :: tau
  type (type_gvar_r3d) :: fri_gr3d

  type (type_gvar_r3d) :: wk
  type (type_gvar_m3d) :: wkgm3d
  real (kind=wp) :: mag
  integer :: is, i, j, k

  call cp_shape_gr3d_gm3d( va, va%g%ew, va%g%ns, wkgm3d )
  call cp_shape_gr3d( va, wk )
  ! horizontal viscosity

  call gra_gr3d( va, wkgm3d )
  call px1_gr3d( shift_ew_gr3d(vc), wk )
  fri_gr3d = 1.0/spbt * div_gm3d( pbt(tc)*wkgm3d )
  call type_add( fri_gr3d, cv1*vb, cv2*spbt*wk )

  fri_gr3d = fri_gr3d * am

  ! vertical viscosity

  wk = g*rho0/spbt**2 * km * px3_gr3d_z( vb )

  ! surface boundary condition
  wk%v(:,:,1) = g*tau(:,:) / spbt(:,:)

  do i = 1, ni
  do j = 1, nj
    k = vb%g%lev(i,j)
    if ( k > 0 ) then ! bottom drag
      mag = sqrt( va%v(i,j,k)**2 + vc%v(i,j,k)**2 )
      wk%v(i,j,k+1) = g*cdbot*mag/spbt(i,j)
    end if
  end do
  end do

  call type_add( fri_gr3d, px3_gr3d(wk) )

  ! to prevent memory leakage
  call free_gr3d( wk ) 
  call free_gm3d( wkgm3d ) 

end function fri_gr3d

function fri_gm3d (spbt, va, vb) !{{{1
  ! horizontal frictional force
  ! output is on the same grid as vb
  real (kind=wp), dimension(:,:), intent(in) :: spbt
  type (type_gvar_m3d), intent(in) :: va, vb
  type (type_gvar_m3d) :: fri_gm3d

  fri_gm3d%x(1) = fri_gr3d( spbt, va%x(1), vb%x(1), &
    (-1.0)*va%x(2), bnd%tau%x(1)%v )
  fri_gm3d%x(2) = fri_gr3d( spbt, va%x(2), vb%x(2), &
    va%x(1), bnd%tau%x(2)%v )

end function fri_gm3d

function adv_gr3d (var, u, w) !{{{1
  ! calc. advection of horizontal pressure weighted velocity var
  ! adv(var) = div(var*v) - 0.5*var*div(v), 
  !   v = (u, w) are unweighted 3d velocity
  type (type_gvar_r3d), intent(in) :: var
  type (type_gvar_m3d) :: u
  type (type_gvar_r3d) :: w
  type (type_gvar_r3d) :: adv_gr3d

  adv_gr3d = &
    div_gm3d( var*shift_ewns_gm3d(u) ) + &
    px3_gr3d( shift_ud(var)*ter_gr3d(w, var%g%ud) ) - &
    0.5*var*( div_gm3d_ter(u, var%g) + px3_gr3d_ter(w, var%g) )

end function adv_gr3d

function adv_gm3d (var, u, w) !{{{1
  ! calc. advection of horizontal pressure weighted velocity var
  ! adv(var) = div(var*v) - 0.5*var*div(v), 
  !   v = (u, w) are unweighted 3d velocity
  type (type_gvar_m3d), intent(in) :: var
  type (type_gvar_m3d) :: u
  type (type_gvar_r3d) :: w
  type (type_gvar_m3d) :: adv_gm3d

  adv_gm3d%x(1) = adv_gr3d( var%x(1), u, w )
  adv_gm3d%x(2) = adv_gr3d( var%x(2), u, w )

end function adv_gm3d

function int_tobot_gr3d (var)!{{{1
  ! indefinite integration to sea bottom
  type (type_gvar_r3d) :: var
  type (type_gvar_r3d) :: int_tobot_gr3d

  call cp_shape_gr3d( var, int_tobot_gr3d )

  int_tobot_gr3d%v = op_int_tobot( var%v, var%g )

end function int_tobot_gr3d

function gra_r2d (var, hg) !{{{1
  ! horizontal gradient operator for 2d variable
  real (kind=wp) :: var(:,:)
  type (type_stg) :: hg
  type (type_gvar_m2d) :: gra_r2d

  integer :: d(2), is

  d = shape(var)
  allocate( gra_r2d%x(1)%v(d(1),d(2)), stat=is ); call chk(is)
  allocate( gra_r2d%x(2)%v(d(1),d(2)), stat=is ); call chk(is)
  gra_r2d%x(1)%hg => hg%ew
  gra_r2d%x(2)%hg => hg%ns

  gra_r2d%x(1)%v = op_px1(var, hg, gra_r2d%x(1)%hg) / &
                   ( a * gra_r2d%x(1)%hg%rh1 )
  gra_r2d%x(2)%v = op_px2(var, hg, gra_r2d%x(2)%hg) / a

end function gra_r2d

function gra_gr2d (var) !{{{1
  ! default horizontal gradient operator for 2d grid variable
  type (type_gvar_r2d) :: var
  type (type_gvar_m2d) :: gra_gr2d

  gra_gr2d%x(1) = px1_gr2d( var ) / ( a * var%hg%ew%rh1 )
  gra_gr2d%x(2) = px2_gr2d( var ) / a

end function gra_gr2d

subroutine gra_gr3d (va, vb) !{{{1
  ! default horizontal gradient operator
  type (type_gvar_r3d) :: va
  type (type_gvar_m3d) :: vb

  call px1_gr3d( va, vb%x(1) )
  call px2_gr3d( va, vb%x(2) )

  vb%x(1)%v = vb%x(1)%v / ( a * vb%x(1)%g%hg%rh1 )
  vb%x(2)%v = vb%x(2)%v / a

end subroutine gra_gr3d

function div_gm2d_ter (var, hg) !{{{1
  ! horizontal divergence operator for 2d grid variable
  ! result on the specific grid hg
  type (type_gvar_m2d), intent(in) :: var
  type (type_stg), target :: hg
  type (type_gvar_r2d) :: div_gm2d_ter

  call cp_shape_gr2d( var%x(1), div_gm2d_ter )
  div_gm2d_ter%hg => hg

  div_gm2d_ter%v = op_div( var%x(1)%v,  var%x(2)%v,  &
                           var%x(1)%hg, var%x(2)%hg, hg)

end function div_gm2d_ter

function div_gm3d (var) !{{{1
  ! default horizontal divergence operator for 3d grid variable
  type (type_gvar_m3d), intent(in) :: var
  type (type_gvar_r3d) :: div_gm3d

  if (.not.associated(var%x(1)%g%ew, var%x(2)%g%ns)) &
    stop 'grid of var in div_gm3d of mod_gvar is not suitable'

  call cp_shape_gr3d( var%x(1), div_gm3d )
  div_gm3d%g => var%x(1)%g%ew

  div_gm3d%v = op_div( var%x(1)%v,    var%x(2)%v, &
                           var%x(1)%g%hg, var%x(2)%g%hg, &
                           div_gm3d%g%hg )

end function div_gm3d

function div_gm3d_ter (var, g) !{{{1
  ! horizontal divergence operator for 3d grid variable
  ! result on the specific grid g
  type (type_gvar_m3d), intent(in) :: var
  type (type_stg3d), target :: g
  type (type_gvar_r3d) :: div_gm3d_ter

  if (.not. (associated(var%x(1)%g%vg, g%vg) .and. &
             associated(var%x(2)%g%vg, g%vg)) )&
    stop 'vertically grid is not proper in div_gm3d_ter of mod_gvar'

  call cp_shape_gr3d( var%x(1), div_gm3d_ter )
  div_gm3d_ter%g => g

  div_gm3d_ter%v = op_div( var%x(1)%v,    var%x(2)%v, &
                       var%x(1)%g%hg, var%x(2)%g%hg, g%hg )

end function div_gm3d_ter

subroutine px1_gr3d ( var, ans )  !{{{1
  ! (partial var) / (partial x1), not include metric effects
  type (type_gvar_r3d), intent(in) :: var
  type (type_gvar_r3d) :: ans

  ans%g => var%g%ew
  ans%v = op_px1( var%v, var%g%hg, ans%g%hg )

end subroutine px1_gr3d

function px1_gr2d ( var )  !{{{1
  ! (partial var) / (partial x1), not include metric effects
  type (type_gvar_r2d) :: var
  type (type_gvar_r2d) :: px1_gr2d

  call cp_shape_gr2d( var, px1_gr2d )
  px1_gr2d%hg => var%hg%ew

  px1_gr2d%v = op_px1( var%v, var%hg, px1_gr2d%hg )

end function px1_gr2d

subroutine px2_gr3d ( var, ans )  !{{{1
  ! (partial var) / (partial x1), not include metric effects
  type (type_gvar_r3d), intent(in) :: var
  type (type_gvar_r3d) :: ans

  ans%g => var%g%ns
  ans%v = op_px2( var%v, var%g%hg, ans%g%hg )

end subroutine px2_gr3d

function px2_gr2d ( var )  !{{{1
  ! (partial var) / (partial x1), not include metric effects
  type (type_gvar_r2d) :: var
  type (type_gvar_r2d) :: px2_gr2d

  call cp_shape_gr2d( var, px2_gr2d )
  px2_gr2d%hg => var%hg%ns

  px2_gr2d%v = op_px2( var%v, var%hg, px2_gr2d%hg )

end function px2_gr2d

function px3_gr3d ( var )  !{{{1
  ! default (partial var) / (partial x3)
  ! upward is positive
  type (type_gvar_r3d), intent(in) :: var
  type (type_gvar_r3d) :: px3_gr3d

  integer :: d1, d2, d3, is

  if ( associated(var%g%vg,vg1) ) then
    d1 = size(var%v, 1)
    d2 = size(var%v, 2)
    d3 = size(vg1%ud%p)
    allocate( px3_gr3d%v(d1,d2,d3), stat=is )
    call chk(is)
    px3_gr3d%g => var%g%ud
    px3_gr3d%v = op_px3( var%v )
  else
    stop 'var should be on vg1 in px3_gr3d in mod_gvar'
  end if

end function px3_gr3d

function px3_gr3d_z ( var )  !{{{1
  ! (partial var) / (partial z), z is thickness in meter
  ! this function is for the odd vertical momentum viscosity
  ! upward is positive
  type (type_gvar_r3d), intent(in) :: var
  type (type_gvar_r3d) :: px3_gr3d_z

  integer :: is

  if ( associated(var%g%vg,vg2) ) then
    allocate( px3_gr3d_z%v(ni,nj,nkp), stat=is ); call chk(is)
    px3_gr3d_z%g => var%g%ud
    px3_gr3d_z%v = 0.0
    ! the first and last layer is left to zero
    px3_gr3d_z%v(:,:,2:nk) = var%v(:,:,1:nk-1) - var%v(:,:,2:nk)
    px3_gr3d_z%v = px3_gr3d_z%v / &
      spread(spread(px3_gr3d_z%g%vg%dz,1,nj), 1, ni)
  else
    stop 'var should be on vg2 in px3_gr3d_z in mod_gvar'
  end if

end function px3_gr3d_z

function px3_gr3d_ter ( var, g )  !{{{1
  ! (partial var) / (partial x3)
  ! upward is positive
  ! result on grid g
  type (type_gvar_r3d), intent(in) :: var
  type (type_stg3d), target :: g
  type (type_gvar_r3d) :: px3_gr3d_ter

  integer :: d1, d2, d3, is

  if ( .not. (associated(var%g%vg,vg1) .and. &
              associated(g%vg,vg2)) ) &
    stop 'var should be on vg1, and g on vg2 in px3_gr3d_ter in od_gvar' 

  d1 = size(var%v, 1)
  d2 = size(var%v, 2)
  d3 = size(g%vg%p)
  allocate( px3_gr3d_ter%v(d1,d2,d3), stat=is )
  call chk(is)
  px3_gr3d_ter%g => g

  px3_gr3d_ter%v = op_px3( var%v, var%g%hg, g%hg )

end function px3_gr3d_ter

function shift_ud( var ) !{{{1
  ! interpolate var to its up-down grid
  type (type_gvar_r3d) :: var
  type (type_gvar_r3d) :: shift_ud

  call cp_shape_gr3d( var, shift_ud )
  shift_ud%g => var%g%ud

  shift_ud%v = op_ter(var%v, var%g%vg, shift_ud%g%vg) 

end function shift_ud

function shift_ew_gr2d( var ) !{{{1
  ! interpolate var to its east-west and north-south neight grid
  type (type_gvar_r2d) :: var
  type (type_gvar_r2d) :: shift_ew_gr2d

  call cp_shape_gr2d( var, shift_ew_gr2d )
  shift_ew_gr2d%hg => var%hg%ew

  shift_ew_gr2d%v = op_ter(var%v, var%hg, shift_ew_gr2d%hg)

end function shift_ew_gr2d

function shift_ew_gr3d( var ) !{{{1
  ! interpolate var to its east-west and north-south neight grid
  type (type_gvar_r3d) :: var
  type (type_gvar_r3d) :: shift_ew_gr3d

  call cp_shape_gr3d( var, shift_ew_gr3d )
  shift_ew_gr3d%g => var%g%ew

  shift_ew_gr3d%v = op_ter(var%v, var%g%hg, shift_ew_gr3d%g%hg)

end function shift_ew_gr3d

function shift_ns_gr2d( var ) !{{{1
  ! interpolate var to its east-west and north-south neight grid
  type (type_gvar_r2d) :: var
  type (type_gvar_r2d) :: shift_ns_gr2d

  call cp_shape_gr2d( var, shift_ns_gr2d )
  shift_ns_gr2d%hg => var%hg%ns

  shift_ns_gr2d%v = op_ter(var%v, var%hg, shift_ns_gr2d%hg)

end function shift_ns_gr2d

function shift_ns_gr3d( var ) !{{{1
  ! interpolate var to its east-west and north-south neight grid
  type (type_gvar_r3d) :: var
  type (type_gvar_r3d) :: shift_ns_gr3d

  call cp_shape_gr3d( var, shift_ns_gr3d )
  shift_ns_gr3d%g => var%g%ns

  shift_ns_gr3d%v = op_ter(var%v, var%g%hg, shift_ns_gr3d%g%hg)

end function shift_ns_gr3d

function shift_ewns_gm3d( var ) !{{{1
  ! interpolate var to its east-west and north-south neight grid
  type (type_gvar_m3d) :: var
  type (type_gvar_m3d) :: shift_ewns_gm3d

  shift_ewns_gm3d = ter_gm3d( var, var%x(1)%g%ew, var%x(2)%g%ns )

end function shift_ewns_gm3d

function ter_gr2d(var, hg, dft) !{{{1
  ! grid variable interpolation from va to vb
  type (type_gvar_r2d), intent(in) :: var
  type (type_stg), target :: hg
  real (kind=wp), optional :: dft
  type (type_gvar_r2d) :: ter_gr2d

  call cp_shape_gr2d( var, ter_gr2d )
  ter_gr2d%hg => hg

  if ( present(dft) ) then
    ter_gr2d%v = op_ter( var%v, var%hg, hg, dft )
  else
    ter_gr2d%v = op_ter( var%v, var%hg, hg )
  end if

end function ter_gr2d

function ter_gr3d(var, g) !{{{1
  ! horizontally interpolate va to grid g
  type (type_gvar_r3d), intent(in) :: var
  type (type_stg3d), target :: g
  type (type_gvar_r3d) :: ter_gr3d

  if ( associated(var%g%vg, g%vg) ) then
    call cp_shape_gr3d( var, ter_gr3d )
    ter_gr3d%g => g
    ter_gr3d%v = op_ter( var%v, var%g%hg, g%hg )
  else
    stop 'vg of g should be the same as var in ter_gr3d in mod_gvar'
  end if

end function ter_gr3d

function ter_gm2d(var, hga, hgb) !{{{1
  ! grid variable interpolation
  type (type_gvar_m2d), intent(in) :: var
  type (type_stg), target :: hga, hgb
  type (type_gvar_m2d) :: ter_gm2d

  ter_gm2d%x(1) = ter_gr2d( var%x(1), hga )
  ter_gm2d%x(2) = ter_gr2d( var%x(2), hgb )

end function ter_gm2d

function ter_gm3d(var, ga, gb) !{{{1
  ! grid variable interpolation
  type (type_gvar_m3d), intent(in) :: var
  type (type_stg3d) :: ga, gb
  type (type_gvar_m3d) :: ter_gm3d

  ter_gm3d%x(1) = ter_gr3d( var%x(1), ga )
  ter_gm3d%x(2) = ter_gr3d( var%x(2), gb )

end function ter_gm3d

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

function gr2d_plus_r2d (va, vb)!{{{1
  ! grid variable addition
  type (type_gvar_r2d), intent(in) :: va
  real (kind=wp), intent(in) :: vb(:,:)
  type (type_gvar_r2d) :: gr2d_plus_r2d

  call cp_shape_gr2d( va, gr2d_plus_r2d )
  gr2d_plus_r2d%v = va%v + vb

end function gr2d_plus_r2d

function gr2d_plus_gr2d (va, vb)!{{{1
  ! grid variable addition
  type (type_gvar_r2d), intent(in) :: va, vb
  type (type_gvar_r2d) :: gr2d_plus_gr2d

  if ( associated(va%hg, vb%hg) ) then
    call cp_shape_gr2d( vb, gr2d_plus_gr2d )
    gr2d_plus_gr2d%v = va%v + vb%v
  else
    stop 'va and vb should be on the same grid in gr2d_plus_gr2d in mod_gvar'
  end if

end function gr2d_plus_gr2d

function gm2d_plus_m2d (va, vb)!{{{1
  ! grid variable addition
  type (type_gvar_m2d), intent(in) :: va
  type (type_mat), intent(in) :: vb(:,:)
  type (type_gvar_m2d) :: gm2d_plus_m2d

  gm2d_plus_m2d%x(1) = gr2d_plus_r2d( va%x(1), vb%x(1) )
  gm2d_plus_m2d%x(2) = gr2d_plus_r2d( va%x(2), vb%x(2) )

end function gm2d_plus_m2d

function gm2d_plus_gm2d (va, vb)!{{{1
  ! grid variable addition
  type (type_gvar_m2d), intent(in) :: va, vb
  type (type_gvar_m2d) :: gm2d_plus_gm2d

  gm2d_plus_gm2d%x(1) = gr2d_plus_gr2d( va%x(1), vb%x(1) )
  gm2d_plus_gm2d%x(2) = gr2d_plus_gr2d( va%x(2), vb%x(2) )

end function gm2d_plus_gm2d

function gr3d_plus_gr3d (va, vb)!{{{1
  ! grid variable addition
  type (type_gvar_r3d), intent(in) :: va, vb
  type (type_gvar_r3d) :: gr3d_plus_gr3d

  if ( associated(va%g, vb%g) ) then
    call cp_shape_gr3d( vb, gr3d_plus_gr3d )
    gr3d_plus_gr3d%v = va%v + vb%v
  else
    stop 'va and vb should be on the same grid in gr3d_plus_gr3d in mod_gvar'
  end if

end function gr3d_plus_gr3d

function gm3d_plus_gm3d (va, vb)!{{{1
  ! grid variable addition
  type (type_gvar_m3d), intent(in) :: va, vb
  type (type_gvar_m3d) :: gm3d_plus_gm3d

  gm3d_plus_gm3d%x(1) = gr3d_plus_gr3d( va%x(1), vb%x(1) )
  gm3d_plus_gm3d%x(2) = gr3d_plus_gr3d( va%x(2), vb%x(2) )

end function gm3d_plus_gm3d

function r2d_minus_gr2d (va, vb)!{{{1
  ! grid variable addition
  real (kind=wp), intent(in) :: va(:,:)
  type (type_gvar_r2d), intent(in) :: vb
  type (type_gvar_r2d) :: r2d_minus_gr2d

  call cp_shape_gr2d( vb, r2d_minus_gr2d )
  r2d_minus_gr2d%v = va - vb%v

end function r2d_minus_gr2d

function gr3d_minus_gr3d (va, vb)!{{{1
  ! grid variable addition
  type (type_gvar_r3d), intent(in) :: va, vb
  type (type_gvar_r3d) :: gr3d_minus_gr3d

  if ( associated(va%g, vb%g) ) then
    call cp_shape_gr3d( vb, gr3d_minus_gr3d )
    gr3d_minus_gr3d%v = va%v - vb%v
  else
    stop 'va and vb should be on the same grid in gr3d_minus_gr3d in mod_gvar'
  end if

end function gr3d_minus_gr3d

function gr2d_minus_r2d (va, vb)!{{{1
  ! grid variable addition
  type (type_gvar_r2d), intent(in) :: va
  real (kind=wp), dimension(:,:), intent(in) :: vb
  type (type_gvar_r2d) :: gr2d_minus_r2d

  call cp_shape_gr2d( va, gr2d_minus_r2d )
  gr2d_minus_r2d%v = va%v - vb

end function gr2d_minus_r2d

function gr3d_minus_r2d (va, vb)!{{{1
  ! grid variable addition
  type (type_gvar_r3d), intent(in) :: va
  real (kind=wp), dimension(:,:), intent(in) :: vb
  type (type_gvar_r3d) :: gr3d_minus_r2d

  integer :: d3

  d3 = size(va%v,3)

  call cp_shape_gr3d( va, gr3d_minus_r2d )
  gr3d_minus_r2d%v = va%v - spread(vb,3,d3)

end function gr3d_minus_r2d

function gm2d_minus_m2d (va, vb)!{{{1
  ! grid variable addition
  type (type_gvar_m2d), intent(in) :: va
  type (type_mat), dimension(:,:), intent(in) :: vb
  type (type_gvar_m2d) :: gm2d_minus_m2d

  gm2d_minus_m2d%x(1) = va%x(1) - vb%x(1)
  gm2d_minus_m2d%x(2) = va%x(2) - vb%x(2)

end function gm2d_minus_m2d

function gm3d_minus_m2d (va, vb)!{{{1
  ! grid variable addition
  type (type_gvar_m3d), intent(in) :: va
  type (type_mat), dimension(:,:), intent(in) :: vb
  type (type_gvar_m3d) :: gm3d_minus_m2d

  gm3d_minus_m2d%x(1) = va%x(1) - vb%x(1)
  gm3d_minus_m2d%x(2) = va%x(2) - vb%x(2)

end function gm3d_minus_m2d

function r_mul_gr2d (va, vb)!{{{1
  ! multiplication
  ! result on the same grid as vb
  real (kind=wp), intent(in) :: va
  type (type_gvar_r2d), intent(in) :: vb
  type (type_gvar_r2d) :: r_mul_gr2d

  call cp_shape_gr2d( vb, r_mul_gr2d )
  r_mul_gr2d%v = va * vb%v

end function r_mul_gr2d

function r_mul_gr3d (va, vb)!{{{1
  ! multiplication
  ! result on the same grid as vb
  real (kind=wp), intent(in) :: va
  type (type_gvar_r3d), intent(in) :: vb
  type (type_gvar_r3d) :: r_mul_gr3d

  call cp_shape_gr3d( vb, r_mul_gr3d )
  r_mul_gr3d%v = va * vb%v

end function r_mul_gr3d

function r_mul_gm2d (va, vb)!{{{1
  ! multiplication
  ! result on the same grid as vb
  real (kind=wp), intent(in) :: va
  type (type_gvar_m2d), intent(in) :: vb
  type (type_gvar_m2d) :: r_mul_gm2d

  call cp_shape_gm2d( vb, r_mul_gm2d )
  r_mul_gm2d%x(1)%v = va * vb%x(1)%v
  r_mul_gm2d%x(2)%v = va * vb%x(2)%v

end function r_mul_gm2d

function r_mul_gm3d (va, vb)!{{{1
  ! multiplication
  ! result on the same grid as vb
  real (kind=wp), intent(in) :: va
  type (type_gvar_m3d), intent(in) :: vb
  type (type_gvar_m3d) :: r_mul_gm3d

  call cp_shape_gm3d( vb, r_mul_gm3d )
  r_mul_gm3d%x(1)%v = va * vb%x(1)%v
  r_mul_gm3d%x(2)%v = va * vb%x(2)%v

end function r_mul_gm3d

function r2d_mul_gr2d (va, vb)!{{{1
  ! multiplication
  ! result on the same grid as vb
  real (kind=wp),  dimension(:,:), intent(in) :: va
  type (type_gvar_r2d), intent(in) :: vb
  type (type_gvar_r2d) :: r2d_mul_gr2d

  call cp_shape_gr2d( vb, r2d_mul_gr2d )
  r2d_mul_gr2d%v = va * vb%v

end function r2d_mul_gr2d

function r2d_mul_gr3d (va, vb)!{{{1
  ! multiplication
  ! result on the same grid as vb
  real (kind=wp),  dimension(:,:), intent(in) :: va
  type (type_gvar_r3d), intent(in) :: vb
  type (type_gvar_r3d) :: r2d_mul_gr3d

  integer :: d3

  d3 = size(vb%v, 3)

  call cp_shape_gr3d( vb, r2d_mul_gr3d )
  r2d_mul_gr3d%v = spread(va,3,d3) * vb%v

end function r2d_mul_gr3d

function r3d_mul_gr3d (va, vb)!{{{1
  ! multiplication
  ! result on the same grid as vb
  real (kind=wp), dimension(:,:,:), intent(in) :: va
  type (type_gvar_r3d), intent(in) :: vb
  type (type_gvar_r3d) :: r3d_mul_gr3d

  call cp_shape_gr3d( vb, r3d_mul_gr3d )
  r3d_mul_gr3d%v = va * vb%v

end function r3d_mul_gr3d

function r2d_mul_gm2d (va, vb)!{{{1
  ! 2d real array multiply the same dimension matrix array
  ! result on the same grid as vb
  real (kind=wp),  dimension(:,:), intent(in) :: va
  type (type_gvar_m2d), intent(in) :: vb
  type (type_gvar_m2d) :: r2d_mul_gm2d

  r2d_mul_gm2d%x(1) = va * vb%x(1)
  r2d_mul_gm2d%x(2) = va * vb%x(2)

end function r2d_mul_gm2d

function r2d_mul_gm3d (va, vb)!{{{1
  ! 2d real array multiply the same dimension matrix array
  ! result on the same grid as vb
  real (kind=wp),  dimension(:,:), intent(in) :: va
  type (type_gvar_m3d), intent(in) :: vb
  type (type_gvar_m3d) :: r2d_mul_gm3d

  r2d_mul_gm3d%x(1) = va * vb%x(1)
  r2d_mul_gm3d%x(2) = va * vb%x(2)

end function r2d_mul_gm3d

function gr2d_mul_gr2d (va, vb)!{{{1
  ! grid variable multiplication
  ! result on the grid of vb
  type (type_gvar_r2d), intent(in) :: va, vb
  type (type_gvar_r2d) :: gr2d_mul_gr2d

  if ( associated(va%hg, vb%hg) ) then
    call cp_shape_gr2d(va, gr2d_mul_gr2d)
    gr2d_mul_gr2d%v = va%v * vb%v
  else
    stop 'va and vb should be on the same horizontal grid in gr2d_mul_gr2d in mod_gvar'
  end if

end function gr2d_mul_gr2d

function gr2d_mul_gr3d (va, vb)!{{{1
  ! grid variable multiplication
  ! result on the grid of vb
  type (type_gvar_r2d), intent(in) :: va
  type (type_gvar_r3d), intent(in) :: vb
  type (type_gvar_r3d) :: gr2d_mul_gr3d

  if ( associated(va%hg, vb%g%hg) ) then
    gr2d_mul_gr3d = r2d_mul_gr3d( va%v, vb )
  else
    stop 'va and vb should be on the same horizontal grid in gr2d_mul_gr3d in mod_gvar'
  end if

end function gr2d_mul_gr3d

function gr3d_mul_gr3d (va, vb)!{{{1
  ! grid variable multiplication
  type (type_gvar_r3d), intent(in) :: va, vb
  type (type_gvar_r3d) :: gr3d_mul_gr3d

  if ( associated(va%g, vb%g) ) then
    gr3d_mul_gr3d = r3d_mul_gr3d( va%v, vb )
  else
    stop 'va and vb should be on the same grid in gr3d_mul_gr3d in mod_gvar'
  end if

end function gr3d_mul_gr3d

function gr2d_mul_gm3d (va, vb)!{{{1
  ! interpolate before multiplication
  ! result is m3d, grid on vb's grid
  type (type_gvar_r2d), intent(in) :: va
  type (type_gvar_m3d), intent(in) :: vb
  type (type_gvar_m3d) :: gr2d_mul_gm3d

  gr2d_mul_gm3d%x(1) = &
    r2d_mul_gr3d( op_ter(va%v, va%hg, vb%x(1)%g%hg), vb%x(1))
  gr2d_mul_gm3d%x(2) = &
    r2d_mul_gr3d( op_ter(va%v, va%hg, vb%x(2)%g%hg), vb%x(2))

end function gr2d_mul_gm3d

function gr3d_mul_gm3d (va, vb)!{{{1
  ! interpolate before multiplication
  ! result is m3d, grid on vb's grid
  type (type_gvar_r3d), intent(in) :: va
  type (type_gvar_m3d), intent(in) :: vb
  type (type_gvar_m3d) :: gr3d_mul_gm3d

  gr3d_mul_gm3d%x(1) = &
    r3d_mul_gr3d( op_ter(va%v, va%g%hg, vb%x(1)%g%hg), vb%x(1))
  gr3d_mul_gm3d%x(2) = &
    r3d_mul_gr3d( op_ter(va%v, va%g%hg, vb%x(2)%g%hg), vb%x(2))

end function gr3d_mul_gm3d

function gr2d_div_r (va, vb)!{{{1
  ! reload divided operator
  type (type_gvar_r2d), intent(in) :: va
  real (kind=wp), intent(in) :: vb
  type (type_gvar_r2d) :: gr2d_div_r

  call cp_shape_gr2d( va, gr2d_div_r )

  gr2d_div_r%v = va%v / vb

end function gr2d_div_r

function gr2d_div_r2d (va, vb)!{{{1
  ! reload divided operator
  type (type_gvar_r2d), intent(in) :: va
  real (kind=wp), intent(in) :: vb(:,:)
  type (type_gvar_r2d) :: gr2d_div_r2d

  call cp_shape_gr2d( va, gr2d_div_r2d )

  gr2d_div_r2d%v = va%v / vb

end function gr2d_div_r2d

function gr3d_div_r (va, vb)!{{{1
  ! reload divided operator
  type (type_gvar_r3d), intent(in) :: va
  real (kind=wp), intent(in) :: vb
  type (type_gvar_r3d) :: gr3d_div_r

  call cp_shape_gr3d( va, gr3d_div_r )

  gr3d_div_r%v = va%v / vb

end function gr3d_div_r

function gr3d_div_r2d (va, vb)!{{{1
  ! reload divided operator
  type (type_gvar_r3d), intent(in) :: va
  real (kind=wp), intent(in) :: vb(:,:)
  type (type_gvar_r3d) :: gr3d_div_r2d

  integer :: d3

  d3 = size(va%v, 3)
  call cp_shape_gr3d( va, gr3d_div_r2d )

  gr3d_div_r2d%v = va%v / spread(vb,3,d3)

end function gr3d_div_r2d

function gr3d_div_r3d (va, vb)!{{{1
  ! reload divided operator
  type (type_gvar_r3d), intent(in) :: va
  real (kind=wp), intent(in) :: vb(:,:,:)
  type (type_gvar_r3d) :: gr3d_div_r3d

  call cp_shape_gr3d( va, gr3d_div_r3d )

  gr3d_div_r3d%v = va%v / vb

end function gr3d_div_r3d

function gm2d_div_r (va, vb)!{{{1
  ! reload divided operator
  type (type_gvar_m2d), intent(in) :: va
  real (kind=wp), intent(in) :: vb
  type (type_gvar_m2d) :: gm2d_div_r

  gm2d_div_r%x(1) = gr2d_div_r( va%x(1), vb )
  gm2d_div_r%x(2) = gr2d_div_r( va%x(2), vb )

end function gm2d_div_r

function gm2d_div_r2d (va, vb)!{{{1
  ! reload divided operator
  type (type_gvar_m2d), intent(in) :: va
  real (kind=wp), intent(in) :: vb(:,:)
  type (type_gvar_m2d) :: gm2d_div_r2d

  gm2d_div_r2d%x(1) = gr2d_div_r2d( va%x(1), vb )
  gm2d_div_r2d%x(2) = gr2d_div_r2d( va%x(2), vb )

end function gm2d_div_r2d

function gm3d_div_r2d (va, vb)!{{{1
  ! reload divided operator
  type (type_gvar_m3d), intent(in) :: va
  real (kind=wp), intent(in) :: vb(:,:)
  type (type_gvar_m3d) :: gm3d_div_r2d

  gm3d_div_r2d%x(1) = gr3d_div_r2d( va%x(1), vb )
  gm3d_div_r2d%x(2) = gr3d_div_r2d( va%x(2), vb )

end function gm3d_div_r2d

function gm3d_div_r3d (va, vb)!{{{1
  ! reload divided operator
  type (type_gvar_m3d), intent(in) :: va
  real (kind=wp), intent(in) :: vb(:,:,:)
  type (type_gvar_m3d) :: gm3d_div_r3d

  gm3d_div_r3d%x(1) = gr3d_div_r3d( va%x(1), vb )
  gm3d_div_r3d%x(2) = gr3d_div_r3d( va%x(2), vb )

end function gm3d_div_r3d

subroutine m2d_eq_gm2d (va, vb) !{{{1
  ! assignment operator
  type (type_mat), intent(inout) :: va(:,:)
  type (type_gvar_m2d), intent(in) :: vb

  va%x(1) = vb%x(1)%v
  va%x(2) = vb%x(2)%v
end subroutine m2d_eq_gm2d

subroutine m3d_eq_gm3d (va, vb) !{{{1
  ! assignment operator
  type (type_mat), dimension(:,:,:), intent(inout) :: va
  type (type_gvar_m3d), intent(in) :: vb

  va%x(1) = vb%x(1)%v
  va%x(2) = vb%x(2)%v
end subroutine m3d_eq_gm3d

subroutine r2d_eq_gr2d (va, vb) !{{{1
  ! assignment operator
  real (kind=wp), intent(inout) :: va(:,:)
  type (type_gvar_r2d), intent(in) :: vb

  va = vb%v
end subroutine r2d_eq_gr2d

subroutine r3d_eq_gr3d (va, vb) !{{{1
  ! assignment operator
  real (kind=wp), intent(inout) :: va(:,:,:)
  type (type_gvar_r3d), intent(in) :: vb

  va = vb%v
end subroutine r3d_eq_gr3d

subroutine gr2d_eq_r (va, vb) !{{{1
  ! assignment operator
  type (type_gvar_r2d), intent(inout) :: va
  real (kind=wp), intent(in) :: vb

  va%v = vb
end subroutine gr2d_eq_r

subroutine gr2d_eq_r2d (va, vb) !{{{1
  ! assignment operator
  type (type_gvar_r2d), intent(inout) :: va
  real (kind=wp), intent(in) :: vb(:,:)

  va%v = vb
end subroutine gr2d_eq_r2d

subroutine gr3d_eq_r (va, vb) !{{{1
  ! assignment operator
  type (type_gvar_r3d), intent(inout) :: va
  real (kind=wp), intent(in) :: vb

  va%v = vb
end subroutine gr3d_eq_r

subroutine gm2d_eq_r (va, vb) !{{{1
  ! assignment operator
  type (type_gvar_m2d), intent(inout) :: va
  real (kind=wp), intent(in) :: vb

  va%x(1) = vb
  va%x(2) = vb
end subroutine gm2d_eq_r

subroutine gm2d_eq_m2d (va, vb) !{{{1
  ! assignment operator
  type (type_gvar_m2d), intent(inout) :: va
  type (type_mat), dimension(:,:), intent(in) :: vb

  va%x(1)%v = vb%x(1)
  va%x(2)%v = vb%x(2)
end subroutine gm2d_eq_m2d

subroutine gm3d_eq_r (va, vb) !{{{1
  ! assignment operator
  type (type_gvar_m3d), intent(inout) :: va
  real (kind=wp), intent(in) :: vb

  va%x(1) = vb
  va%x(2) = vb
end subroutine gm3d_eq_r

subroutine m2d_assign_gm2d (va, vb) !{{{1
  ! assignment operator
  type (type_mat), intent(inout) :: va(:,:)
  type (type_gvar_m2d), intent(in) :: vb

  va%x(1) = vb%x(1)%v
  va%x(2) = vb%x(2)%v
end subroutine m2d_assign_gm2d

subroutine gm2d_assign_m2d (va, vb) !{{{1
  ! assignment operator
  type (type_gvar_m2d) :: va
  type (type_mat), dimension(:,:), intent(in) :: vb

  va%x(1)%v = vb%x(1)
  va%x(2)%v = vb%x(2)
end subroutine gm2d_assign_m2d

subroutine chk( ista ) !{{{1
  ! check state of allocate array 

  integer, intent(in) ::  ista

  if ( ista /= 0 ) then
    write(*,*) 'Allocate array failed. Stop'
    stop 2
  end if
end subroutine chk

end module mod_gvar!{{{1
!-------------------------------------------------------{{{1
! vim:fdm=marker:fdl=0:
! vim:foldtext=getline(v\:foldstart).'...'.(v\:foldend-v\:foldstart):
