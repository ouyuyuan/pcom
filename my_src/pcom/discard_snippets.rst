.. snippets have been discard but may be usefull in future

calc metric
===========

::
  subroutine calc_metric_old (mg1, mg2, mg3, mg4)
    ! calc. horizontal metric propertics
    ! (lon, lat) in I/O file is on grid 4

    type (type_met), dimension(ni,nj) :: mg1, mg2, mg3, mg4

    integer :: i, j
    type (type_mat) :: p
    real (kind=wp) :: dlon, dlat, clat

    ! inner area
    do j = 2, njm
    do i = 2, nim
      !westest/eastest longitude wrap up
      dlon = ( lon(i+1) - lon(i-1) )*0.5
      if ( dlon < 0 ) dlon = ( lon(i+1) + 360 - lon(i-1) )*0.5
      dlat = lat(j) - lat(j-1)
      p    = (/dlon, dlat/)
      mg1(i,j)%dx = p * torad

      dlon = lon(i+1) - lon(i)
      if ( dlon < 0 ) dlon = lon(i+1) + 360 - lon(i)
      dlat = lat(j) - lat(j-1)
      p    = (/dlon, dlat/)
      mg2(i,j)%dx = p * torad

      dlon = lon(i+1) - lon(i)
      if ( dlon < 0 ) dlon = lon(i+1) + 360 - lon(i)
      dlat = 0.5 * ( lat(j+1) - lat(j-1) )
      p    = (/dlon, dlat/)
      mg3(i,j)%dx = p * torad

      dlon = ( lon(i+1) - lon(i-1) )*0.5
      if ( dlon < 0 ) dlon = ( lon(i+1) + 360 - lon(i-1) )*0.5
      dlat = 0.5 * ( lat(j+1) - lat(j-1) )
      p    = (/dlon, dlat/)
      mg4(i,j)%dx = p * torad

      clat = 0.5 * ( lat(j) + lat(j-1) )
      p    = (/ cos( clat*torad ), 1.0/)
      mg1(i,j)%h = a * p

      mg2(i,j)%h = a * p

      clat = lat(j)
      p    = (/ cos( clat*torad ), 1.0/)
      mg3(i,j)%h = a * p

      mg4(i,j)%h = a * p
    end do
    end do
    call mympi_swpbnd (mg1, mg2, mg3, mg4)

    ! extropolate

    if (my%gs == 1)      call op_ext( mg1, mg2, mg3, mg4, 'southest' )

    if (my%gn == glo_nj) call op_ext( mg1, mg2, mg3, mg4, 'northest' )

    ! jacobian values

    do i = 1, ni
    do j = 1, nj
      mg1(i,j)%j  = product( mg1(i,j)%h%x )
      mg2(i,j)%j  = product( mg2(i,j)%h%x )
      mg3(i,j)%j  = product( mg3(i,j)%h%x )
      mg4(i,j)%j  = product( mg4(i,j)%h%x )
    end do
    end do
  end subroutine calc_metric_old

gradient force of atmospheric pressure
======================================
! fgpa: atmospheric pressure ::
  fgpa = rrho(:,:,1) * op_gra (bnd%pa)
  call mympi_swpbnd (fgpa)

debug sentences
===============
::
  if (myid==mid) print *, 'here a'
  if (myid==mid) print *, 'here b'
  if (myid==mid) print *, 'here c'
  if (myid==mid) print *, 'here d'
  if (myid==mid) print *, 'here e'

special treatment at northest/southest boundary
===============================================
::
  ! in the northest boundary, put actual value in the first row, 
  !   not missing value, this can save a lot of exception when 
  !   doing interpolation between T- and U-grids. And because the 
  !   the northest can be land, it will not produce problem if not
  !   include it in the general calculation. The similar applies to
  !   the southest row, but not the toppest and deepest layer. (since
  !   there is no 'upper land')
  if ( gp%i == 1 .or. gp%i == nm%py ) then
    gp % ni = gp % ni - 1
    gp % gn = gp%gs - (gp%ni - 1) + 1 ! only 1 row for boundary
  end if



distribute bottom pressure from some id to other ids
====================================================

::
  tag = 10
  do i = 2, ny-1
  do j = 2, nx-1
    k = levu(i, j)
    ! land
    if ( k == 0 ) then
      pbu(i,j) = 0.0
    ! domain contain sea bottom
    else if ( k >= gpw%gu .and. k <= gpw%gl ) then
      pbu(i,j) = pmn(k - gpw%gu + 1)
      call mympi_diffuse_z ( pbu(i,j), tag )
    ! domain above sea bottom
    else if ( k > gpw%gl ) then
      call mympi_recvsend_up ( pbu(i,j), tag )
    ! domain below sea bottom
    else if ( k < gpw%gu ) then
      call mympi_recvsend_down ( pbu(i,j), tag )
    end if
  end do
  end do

indefinite integral of vartical variable
========================================

::
    ! the top is for mpi-boundary
    pmn(2) = rmn(2) * g * grdu%dz(2)
    do k = 3, nz - 1
      pmn(k) = rmn(k) * g * grdu%dz(k) + pmn(k-1)
    end do
    call mympi_fix_integrate (pmn)

thickness increasement
======================

::
  ! glo_dz, solve the eq. ( dz(k-1) + dz(k) ) / 2 = z(k) - z(k-1)
  glo_dz(1) = 2*glo_z(1)
  do k = 2, glo_nz
    glo_dz(k) = 2 * ( glo_z(k) - glo_z(k-1) ) - glo_dz(k-1) 
  end do

half grid to integer grid interpolation
=======================================

::
  function op_xh2w (var)
    ! interpolate var from half grid to whole grid in the x-direction
    real (kind=wp), dimension(ny,nx), intent(in) :: var
    real (kind=wp), dimension(ny,nx) :: op_xh2w

    integer :: i, j

    op_xh2w = var
    forall ( i=2:ny-1, j=2:nx-1 ) &
        op_xh2w(i,j) = 0.5*( var(i,j-1) + var(i,j) )
    call mympi_swpbnd (op_xh2w)

    return
  end function op_xh2w

  function op_yh2w (var)
    ! interpolate var from half grid to whole grid in the y-direction
    real (kind=wp), dimension(ny,nx), intent(in) :: var
    real (kind=wp), dimension(ny,nx) :: op_yh2w

    integer :: i, j

    op_yh2w = var
    ! note i increase southwards
    forall ( i=2:ny-1, j=2:nx-1 ) &
        op_yh2w(i,j) = 0.5*( var(i+1,j) + var(i,j) )

    ! north south is not wrap up
    if ( gpw % gs == glo_ny ) op_yh2w(ny-1,:) = var(ny-1,:)

    call mympi_swpbnd (op_yh2w)

    return
  end function op_yh2w

subroutine output_grdvar () !{{{1
  ! output grid variables, mainly for code checing

  character (len=*), parameter :: fname = 'grdvar.nc'

  if ( myid == mid ) then
    call io_create_grdvar (fname)
    write(*, '(a)') 'created output file '//fname

    call io_write (fname, 'lon', glo_lon)
    call io_write (fname, 'lat', glo_lat)
    call io_write (fname, 'z',   z)

!     1d vertical variables
    call io_write (fname, 'vx1', vg1(1:nk)%z)
    call io_write (fname, 'vx2', vg2%z)
    call io_write (fname, 'vp1', vg1(1:nk)%p)
    call io_write (fname, 'vp2', vg2%p)
  end if

  call mpi_barrier (mpi_comm_world, err)

  call mympi_output (fname, 'tmsk', tmsk)
  call mympi_output (fname, 'umsk', umsk)

  call mympi_output (fname, 'itn', itn)
  call mympi_output (fname, 'iun', iun)

  call mympi_output (fname, 'phib', bot_old%phib)
  call mympi_output (fname, 'pbu', bot_old%pbu)
  call mympi_output (fname, 'pb', bot_old%pb)

  call mympi_output (fname, 'phibx', graphib_old%x(1))
  call mympi_output (fname, 'phiby', graphib_old%x(2))

  call mympi_output (fname, 'h1', met(:,:,1)%h%x(1))
  call mympi_output (fname, 'h2', met(:,:,2)%h%x(1))
  call mympi_output (fname, 'h3', met(:,:,3)%h%x(1))
  call mympi_output (fname, 'h4', met(:,:,4)%h%x(1))

  if ( myid == mid )write(*, '(a)') 'Done writing '//fname

end subroutine output_grdvar

subroutine inistat_old (ts_old, frc, bot_old, graphib_old) !{{{1
  ! prepare the initial state of the ocean
  type (type_mat), dimension(ni,nj,nk) :: ts_old
  type (type_frc), dimension(ni,nj,12) :: frc
  type (type_bot_old), dimension(ni,nj) :: bot_old
  type (type_mat), dimension(ni,nj) :: graphib_old

  real (kind=wp), allocatable, dimension(:,:,:) :: glo_pt, glo_sa
  type (type_frc), allocatable, dimension(:,:,:) :: glo_frc

  real (kind=wp) :: ra(ni,nj,nk)
  real (kind=wp), dimension(ni,nj) :: phib, pb
  real (kind=wp), dimension(nk) :: phimn

  integer :: is, m, i, j, k

  if ( myid == mid ) then
    ! 12 months forcing
    allocate(glo_frc(glo_ni, glo_nj, 12), stat=is); call chk(is)
    glo_frc%tau = 0.0; glo_frc%ts = 0.0
    glo_frc%pa = 0.0; glo_frc%fw = 0.0

    allocate(glo_pt(glo_ni, glo_nj, nk), stat=is); call chk(is)
    glo_pt = 0.0
    allocate(glo_sa(glo_ni, glo_nj, nk), stat=is); call chk(is)
    glo_sa = 0.0

    call io_read (nm%fi, 'pt',  glo_pt)
    call io_read (nm%fi, 'sa',  glo_sa)

    call io_read (nm%ff, 'taux', glo_frc%tau%x(1)) 
    call io_read (nm%ff, 'tauy', glo_frc%tau%x(2)) 

    call io_read (nm%ff, 'bct', glo_frc%ts%x(1)) 
    call io_read (nm%ff, 'bcs', glo_frc%ts%x(2)) 

    call io_read (nm%ff, 'pa', glo_frc%pa) 

    call io_read (nm%ff, 'fw', glo_frc%fw) 
  end if

  do m = 1, 12
    frc(:,:,m)%tau%x(1) = mympi_div (glo_frc(:,:,m)%tau%x(1))
    frc(:,:,m)%tau%x(2) = mympi_div (glo_frc(:,:,m)%tau%x(2))

    frc(:,:,m)%ts%x(1) = mympi_div (glo_frc(:,:,m)%ts%x(1))
    frc(:,:,m)%ts%x(2) = mympi_div (glo_frc(:,:,m)%ts%x(2))

    frc(:,:,m)%pa = mympi_div (glo_frc(:,:,m)%pa)
    frc(:,:,m)%fw  = mympi_div (glo_frc(:,:,m)%fw)

  end do

  ! forcing on land set to zero (no forcing)
  !!! potential bug, missing is not all at land of input file
!  where (spread(itn, 3, 12) == 0)
  where (frc%tau%x(1) == missing_float)
    frc%tau%x(1) = 0.0 
    frc%tau%x(2) = 0.0 
    frc%ts%x(1) = 0.0 
    frc%ts%x(2) = 0.0 
    frc%pa = 0.0 
    frc%fw = 0.0 
  end where

  ! initial (T, S)
  ts_old%x(1)  = mympi_div (glo_pt)
  ts_old%x(2)  = mympi_div (glo_sa)

  ! sea bottom pressure, geoptential height
  phimn(1) = - g * vg2(1)%dz
  do k = 2, nk
    phimn(k) = phimn(k-1)  - g * vg2(k)%dz
  end do

  do j = 1, nj
  do i = 1, ni
    k = itn(i, j)
    if ( k == 0 ) then
      bot_old(i,j)%phib = 0.0
    else
      bot_old(i,j)%phib = phimn(k)
      bot_old(i,j)%pb = vg1(k)%p
    end if
  end do
  end do

  graphib_old = op_gra_to3 (bot_old%phib)

  ra(:,:,:) = umsk * spread(spread(vg2%dp, 1, nj), 1, ni)
  bot_old(:,:)%pbu  = sum(ra(:,:,:), 3)

  ! in case of 1/pb, set to 0.1Pa in land
  where( itn == 0 ) bot_old%pb  = 0.1
  where( iun == 0 ) bot_old%pbu = 0.1

end subroutine inistat_old

function op_vint_r3d_old (var, topo, vg)!{{{1
  ! vertical integration from sea surface to sea bottom
  real (kind=wp), intent(in) :: var(:,:,:)
  integer, intent(in) :: topo(:,:)
  type (type_vg), intent(in) :: vg(:)

  real (kind=wp), dimension(ni,nj) :: op_vint_r3d_old

  integer :: i, j, n

  op_vint_r3d_old(:,:) = 0.0

  do j = 1, nj
  do i = 1, ni
    n = topo(i,j)
    if ( n > 0 ) then
      op_vint_r3d_old(i,j) = sum( var(i,j,1:n) * vg(1:n)%dp )
    end if
  end do
  end do

  ! scale factor resulting from integration of sigma coord.
  op_vint_r3d_old = op_vint_r3d_old / bot_old%pbu

end function op_vint_r3d_old

function ter_r2d_old(var, di) !{{{1
  ! interpolate from one grid to another
  ! di indicates directions
  real (kind=wp), intent(in) :: var(:,:)
  character (len=*), intent(in) :: di

  real (kind=wp), dimension(ni,nj) :: ter_r2d_old

  if      ( di == '1->2' ) then
    ter_r2d_old(1:nim,:) = 0.5 * ( var(1:nim,:) + var(2:ni,:) )
    ter_r2d_old(ni,:)    = 0.0

  else if ( di == '1->3' ) then
    ter_r2d_old(1:nim,1:njm) = &
      ( var(1:nim,1:njm) + var(1:nim,2:nj) + &
        var(2:ni, 1:njm) + var(2:ni, 2:nj) ) * 0.25
    ter_r2d_old(ni,:) = 0.0
    ter_r2d_old(:,nj) = 0.0

  else if ( di == '1->4' ) then
    ter_r2d_old(:,1:njm) = 0.5 * ( var(:,1:njm) + var(:,2:nj) )
    ter_r2d_old(:,nj)    = 0.0

  else if ( di == '2->3' ) then
    ter_r2d_old(:,1:njm) = 0.5 * ( var(:,1:njm) + var(:,2:nj) )
    ter_r2d_old(:,nj)    = 0.0

  else if ( di == '3->1' ) then
    ter_r2d_old(2:ni,2:nj) = &
      ( var(2:ni,2:nj)  + var(1:nim,2:nj) + &
        var(2:ni,1:njm) + var(1:nim,1:njm) ) * 0.25
    ter_r2d_old(1,:) = 0.0
    ter_r2d_old(:,1) = 0.0

  else if ( di == '3->2' ) then
    ter_r2d_old(:,2:nj) = 0.5 * ( var(:,2:nj) + var(:,1:njm) )
    ter_r2d_old(:,1)    = 0.0

  else if ( di == '3->4' ) then
    ter_r2d_old(2:ni,:) = 0.5 * ( var(2:ni,:) + var(1:nim,:) )
    ter_r2d_old(1,:)    = 0.0

  else if ( di == '4->3' ) then
    ter_r2d_old(1:nim,:) = 0.5 * ( var(1:nim,:) + var(2:ni,:) )
    ter_r2d_old(ni,:)    = 0.0

  else
    print *, 'unhandled interpolate indicator '//di//' in ter_r2d_old of mod_op'
  end if
end function ter_r2d_old

function op_px_r2d ( var )  !{{{1
  ! (partial var) / (partial x1), not include metric effects
  ! grid of var and direction of difference indicate by gd
  ! It is the caller's responsibility to swap bondaries
  type (type_gvar_r2d), intent(in) :: var
  type (type_gvar_r2d) :: op_px_r2d

  integer :: is, hg

  allocate(op_px_r2d%v(ni,nj), stat = is); call chk(is)

  op_px_r2d%v = 0.0

  op_px_r2d%vg = var%vg

  if      ( var%hg == 1 ) then
    op_px_r2d%hg = 2
  else if ( var%hg == 2 ) then
    op_px_r2d%hg = 1

  else if ( var%hg == 3 ) then
    op_px_r2d%hg = 4
  else if ( var%hg == 4 ) then
    op_px_r2d%hg = 3

  else
    stop 'unhandled value of var%hg in ter_r2d of mod_op'
  end if

  if ( var%hg == 1 .or. var%hg == 4 ) then
    op_px_r2d%v(1:nim,:) = var%v(2:ni,:) - var%v(1:nim,:)
  else
    op_px_r2d%v(2:ni,:)  = var%v(2:ni,:) - var%v(1:nim,:)
  end if

  op_px_r2d%v = op_px_r2d%v / met(:,:,op_px_r2d%hg)%dx%x(1)

end function op_px_r2d

function op_py_r2d ( var )  !{{{1
  ! (partial var) / (partial x2), not include metric effects
  ! grid of var and direction of difference indicate by gd
  ! It is the caller's responsibility to swap bondaries
  type (type_gvar_r2d), intent(in) :: var
  type (type_gvar_r2d) :: op_py_r2d

  integer :: is, hg

  allocate(op_py_r2d%v(ni,nj), stat = is); call chk(is)

  op_py_r2d%v = 0.0

  op_py_r2d%vg = var%vg

  if      ( var%hg == 1 ) then
    op_py_r2d%hg = 4
  else if ( var%hg == 4 ) then
    op_py_r2d%hg = 1

  else if ( var%hg == 2 ) then
    op_py_r2d%hg = 3
  else if ( var%hg == 3 ) then
    op_py_r2d%hg = 2

  else
    stop 'unhandled value of var%hg in ter_r2d of mod_op'
  end if

  if ( var%hg == 1 .or. var%hg == 2 ) then
    op_py_r2d%v(:,1:njm) = var%v(:,2:nj) - var%v(:,1:njm)
  else
    op_py_r2d%v(:,2:nj)  = var%v(:,2:nj) - var%v(:,1:njm)
  end if

  op_py_r2d%v = op_py_r2d%v / met(:,:,op_py_r2d%hg)%dx%x(2)

end function op_py_r2d

function op_dq_r3d ( var, gd )  !{{{1
  ! difference quotient
  ! grid of var and direction of difference indicate by gd
  real (kind=wp), dimension(:,:,:), intent(in) :: var
  character (len=*), intent(in) :: gd
  real (kind=wp), dimension(ni,nj,nk) :: op_dq_r3d

  integer :: k

  do k = 1, nk
    op_dq_r3d(:,:,k) = op_dq_r2d( var(:,:,k), gd )
  end do

end function op_dq_r3d

function op_dq_r2d ( var, gd )  !{{{1
  ! difference quotient
  ! grid of var and direction of difference indicate by gd
  real (kind=wp), dimension(:,:), intent(in) :: var
  character (len=*), intent(in) :: gd
  real (kind=wp), dimension(ni,nj) :: op_dq_r2d

  if      ( gd == 'g1x1' ) then
    ! shift eastward half-grid, result on grid 2 
    op_dq_r2d(1:nim,:) = &
    ( var(2:ni,:) - var(1:nim,:) ) / met(1:nim,:,2)%dx%x(1)
    op_dq_r2d(ni,:) = 0.0

  else if ( gd == 'g2x1' ) then
    ! shift westward half-grid, result on grid 1
    op_dq_r2d(2:ni,:) = &
    ( var(2:ni,:) - var(1:nim,:) ) / met(2:ni,:,1)%dx%x(1)
    op_dq_r2d(1,:) = 0.0

  else if ( gd == 'g3x1' ) then
    ! shift westward half-grid, result on grid 4
    op_dq_r2d(2:ni,:) = &
    ( var(2:ni,:) - var(1:nim,:) ) / met(2:ni,:,4)%dx%x(1)
    op_dq_r2d(1,:) = 0.0

  else if ( gd == 'g4x1' ) then
    ! shift eastward half-grid, result on grid 3
    ! note that the grid 'shift back' to grid 3, because it
    !   was shift forward from g3 to g4, so the left hand side
    !   is op_dq_r2d(1:nim,:), not op_dq_r2d(2:ni,:)
    op_dq_r2d(1:nim,:) = &
    ( var(2:ni,:) - var(1:nim,:) ) / met(1:nim,:,3)%dx%x(1)
    op_dq_r2d(ni,:) = 0.0

  else if ( gd == 'g1x2' ) then
    ! shift northward half-grid, result on grid 4
    op_dq_r2d(:,1:njm) = &
    ( var(:,2:nj) - var(:,1:njm) ) / met(:,1:njm,4)%dx%x(2)
    op_dq_r2d(:,nj) = 0.0

  else if ( gd == 'g2x2' ) then
    ! shift northward half-grid, result on grid 3
    op_dq_r2d(:,1:njm) = &
    ( var(:,2:nj) - var(:,1:njm) ) / met(:,1:njm,3)%dx%x(2)
    op_dq_r2d(:,nj) = 0.0

  else if ( gd == 'g3x2' ) then
    ! shift southward half-grid, result on grid 2
    op_dq_r2d(:,2:nj) = &
    ( var(:,2:nj) - var(:,1:njm) ) / met(:,2:nj,2)%dx%x(2)
    op_dq_r2d(:,1) = 0.0

  else if ( gd == 'g4x2' ) then
    ! shift southward half-grid, result on grid 1
    op_dq_r2d(:,2:nj) = &
    ( var(:,2:nj) - var(:,1:njm) ) / met(:,2:nj,1)%dx%x(2)
    op_dq_r2d(:,1) = 0.0

  else
    stop 'unhandled interpolate indicator in ter_r2d of mod_op'
  end if
end function op_dq_r2d

function gra_r2d_old (var, ga, gb, gc) !{{{1
  ! horizontal gradient operator for 2d scalar
  ! var on grid ga, output on gb, gc
  real (kind=wp),  dimension(ni,nj), intent(in) :: var
  type (type_stg), intent(in) :: ga, gb, gc

  type (type_mat), dimension(ni,nj) :: gra_r2d_old

  gra_r2d_old%x(1) = op_px1(var, ga, gb) / ( a * gb%rh1 )
  gra_r2d_old%x(2) = op_px2(var, ga, gc) / a

end function gra_r2d_old

function px3_gr3d_old ( var, hg, vg )  !{{{1
  ! (partial var) / (partial x3)
  ! upward is positive
  ! result horizontally in hg, vertically on vg (vg2)
  type (type_gvar_r3d), intent(in) :: var
  type (type_stg), target :: hg
  type (type_vstg), target :: vg
  real (kind=wp), dimension(ni,nj,nk) :: px3_gr3d_old

  real (kind=wp), dimension(ni,nj,nkp) :: temp
  integer :: d3

  if (.not.associated(var%g%vg, vg1)) &
    stop 'var should be on vg1 in px3_gr3d_old of mod_op'

  ! interpolate to the matching grid first
  if ( .not.associated(var%g%hg, hg) ) then
    temp = ter_r3d( var%v, var%g%hg, hg )
  else
    temp = var%v
  end if

  px3_gr3d_old = 0.0
  px3_gr3d_old(:,:,1:nk) = temp(:,:,1:nk) - temp(:,:,2:nkp)
  px3_gr3d_old = px3_gr3d_old / spread(spread(vg%dp,1,nj), 1, ni)

  ! no need to swap boundary horizontally
end function px3_gr3d_old

function op_adv_old (var) !{{{1
  ! calc. advection of horizontal pressure weighted velocity
  ! adv(var) = div(var v) - 0.5*var*div(v), v = (unu, unv, unw)
  ! var on the same grid as global up
  ! result g3
  type (type_gvar_r3d), intent(in) :: var
  real (kind=wp), dimension(ni,nj,nk) :: op_adv_old

  real (kind=wp), dimension(ni,nj,nk) :: unu, unv, u, v
  real (kind=wp), dimension(ni,nj,nkp) :: unw, w
  real (kind=wp), dimension(ni,nj) :: spbt
  integer :: is

  ! get unweighted velocity from global variables
  ! (they are on the same grid as weighted velocity)

  unw = wm%v / spread(pbt(tc)%v,3,nk)
  spbt = sqrt( op_ter(pbt(tc)%v, pbt(tc)%hg, g3%hg) )
  unu = up(tc)%x(1)%v / spread(spbt,3,nk) 
  unv = up(tc)%x(2)%v / spread(spbt,3,nk) 

  ! calc. -0.5*var*div(v)

  u = ter_r3d(unu, up(tc)%x(1)%g%hg, up(tc)%x(1)%g%hg%ew)
  v = ter_r3d(unv, up(tc)%x(2)%g%hg, up(tc)%x(2)%g%hg%ns)
  w = ter_r3d(unw, wm%g%hg, g3%hg)
  op_adv_old = - 0.5 * var%v * div_vec3d(u, v, w)

  ! calc. div(var v)

  u =  ter_r3d(var%v, up(tc)%x(1)%g%hg, up(tc)%x(1)%g%hg%ew) * u
  v =  ter_r3d(var%v, up(tc)%x(2)%g%hg, up(tc)%x(2)%g%hg%ns) * v
  w = vter_r3d(var%v, g3%vg, wm%g%vg)    * w
  op_adv_old = div_vec3d(u, v, w) + op_adv_old

end function op_adv_old

function div_vec3d (u, v, w) !{{{1
  ! 3d divergence of (u,v,w)
  ! care with the input parameters, this function does not 
  !   check the grids
  ! horizontally, (u,v,w) on (hg4, hg2, hg3), result on hg3
  ! vertically, (u,v,w) on (vg2, vg2, vg1), result on vg2
  real (kind=wp), dimension(ni,nj,nk) :: u, v
  real (kind=wp), dimension(ni,nj,nkp) :: w
  real (kind=wp), dimension(ni,nj,nk) :: div_vec3d

  div_vec3d = div_r3d(u, v, hg4, hg2, hg3) + &
              px3_r3d(w, hg3, hg3)

end function div_vec3d

subroutine ter_gm2d_old(va, vb) !{{{1
  ! grid variable interpolation from va to vb
  type (type_gvar_m2d), intent(in) :: va
  type (type_gvar_m2d) :: vb

  call ter_gr2d( va%x(1), vb%x(1) )
  call ter_gr2d( va%x(2), vb%x(2) )

end subroutine ter_gm2d_old

subroutine ter_gr2d_old(va, vb) !{{{1
  ! grid variable interpolation from va to vb
  type (type_gvar_r2d), intent(in) :: va
  type (type_gvar_r2d) :: vb

  vb%v = op_ter( va%v, va%hg, vb%hg )

end subroutine ter_gr2d_old

subroutine ter_gr3d(va, vb) !{{{1
  ! grid variable interpolation from va to vb
  type (type_gvar_r3d), intent(in) :: va
  type (type_gvar_r3d) :: vb

  if ( associated(va%g%vg, vb%g%vg) ) then 
    vb%v = op_ter( va%v, va%g%hg, vb%g%hg )
  else
    stop 'vg of va and vb in ter_gr3d of mod_gvar should be the same.'
  end if

end subroutine ter_gr3d

function gra_gr3d_old (var) !{{{1
  ! default horizontal gradient operator
  type (type_gvar_r3d) :: var
  type (type_gvar_m3d) :: gra_gr3d_old

  gra_gr3d_old%x(1) = px1_gr3d_old( var ) / ( a * var%g%hg%ew%rh1 )
  gra_gr3d_old%x(2) = px2_gr3d_old( var ) / a

end function gra_gr3d_old

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
  fri_gr3d = 1.0/spbt * div_gm3d( pbt(tc)*wkgm3d ) + &
             cv1*vb + cv2*spbt*wk
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

  fri_gr3d = fri_gr3d + px3_gr3d( wk )

  ! to prevent memory leakage
  call free_gr3d( wk ) 
  call free_gm3d( wkgm3d ) 

end function fri_gr3d

function px1_gr3d_old ( var )  !{{{1
  ! (partial var) / (partial x1), not include metric effects
  type (type_gvar_r3d) :: var
  type (type_gvar_r3d) :: px1_gr3d_old

  call cp_shape_gr3d( var, px1_gr3d_old )
  px1_gr3d_old%g => var%g%ew

  px1_gr3d_old%v = op_px1( var%v, var%g%hg, px1_gr3d_old%g%hg )

end function px1_gr3d_old

function px2_gr3d_old ( var )  !{{{1
  ! (partial var) / (partial x1), not include metric effects
  type (type_gvar_r3d) :: var
  type (type_gvar_r3d) :: px2_gr3d_old

  call cp_shape_gr3d( var, px2_gr3d_old )
  px2_gr3d_old%g => var%g%ns

  px2_gr3d_old%v = op_px2( var%v, var%g%hg, px2_gr3d_old%g%hg )

end function px2_gr3d_old

function r2d_plus_gr2d (va, vb)!{{{1
  ! grid variable addition
  real (kind=wp), intent(in) :: va(:,:)
  type (type_gvar_r2d), intent(in) :: vb
  type (type_gvar_r2d) :: r2d_plus_gr2d

  call cp_shape_gr2d( vb, r2d_plus_gr2d )
  r2d_plus_gr2d%v = va + vb%v

end function r2d_plus_gr2d

function m2d_plus_gm2d (va, vb)!{{{1
  ! grid variable addition
  type (type_mat), intent(in) :: va(:,:)
  type (type_gvar_m2d), intent(in) :: vb
  type (type_gvar_m2d) :: m2d_plus_gm2d

  m2d_plus_gm2d%x(1) = r2d_plus_gr2d( va%x(1), vb%x(1) )
  m2d_plus_gm2d%x(2) = r2d_plus_gr2d( va%x(2), vb%x(2) )

end function m2d_plus_gm2d

function vter_r3d_old(var, vga, vgb, dft) !{{{1
  ! interpolate vertically from grid vga to grid vgb
  real (kind=wp), intent(in) :: var(:,:,:)
  type (type_vstg), target :: vga, vgb
  real (kind=wp), optional :: dft
  real (kind=wp), allocatable:: vter_r3d_old(:,:,:)

  integer :: is, nda, ndb

  if ( vga%n == vgb%n ) &
    stop 'no need to interpolate in vter_r3d_old of mod_op'

  nda = size(vga%p)
  ndb = size(vgb%p)

  if ( size(var,3) /= nda ) &
    stop 'var and vga unmatch in vter_r3d_old of mod_op'

  allocate( vter_r3d_old(ni,nj,ndb), stat=is ); call chk(is)

  if ( present(dft) ) then
    vter_r3d_old = dft
  else
    vter_r3d_old = 0.0
  end if

  ! vg2 to vg1
  if ( vga%n == 2 ) then
    vter_r3d_old(:,:,2:nk) = ( var(:,:,1:nk-1) + var(:,:,2:nk) ) * 0.5
  ! from vg1 to vg2
  else
    vter_r3d_old(:,:,1:nk) = ( var(:,:,1:nk) + var(:,:,2:nkp) ) * 0.5
  end if

  ! no need to swap boundary horizontally

end function vter_r3d_old

function ter_r3d_old(var, hga, hgb, dft) !{{{1
  ! interpolate from grid hga to grid hgb
  real (kind=wp), intent(in) :: var(:,:,:)
  type (type_stg), target :: hga, hgb
  real (kind=wp), optional :: dft
  real (kind=wp), allocatable:: ter_r3d_old(:,:,:)

  integer :: k, d3, is

  d3 = size(var, 3)
  allocate( ter_r3d_old(ni,nj,d3), stat = is ); call chk(is)

  if ( present(dft) ) then
    ter_r3d_old = dft
  else
    ter_r3d_old = 0.0
  end if

  if ( hga%n == hgb%n ) then
    stop 'no need to interpolate in ter_r3d_old of mod_op'
  else if ( associated(hga%ew, hgb) ) then
    ter_r3d_old(1+hga%i:nim+hga%i,:,:) = &
    0.5 * ( var(1:nim,:,:) + var(2:ni,:,:) )
  else if ( associated(hga%ns, hgb) ) then
    ter_r3d_old(:,1+hga%j:njm+hga%j,:) = &
      0.5 * ( var(:,1:njm,:) + var(:,2:nj,:) )
  else if ( associated(hga%di, hgb) ) then
    ter_r3d_old(1+hga%i:nim+hga%i,1+hga%j:njm+hga%j,:) = &
      ( var(1:nim,1:njm,:) + var(1:nim,2:nj,:) + &
        var(2:ni, 1:njm,:) + var(2:ni, 2:nj,:) ) * 0.25
  else
    print *, 'unhandled relative position in ter_r3d_old of mod_op'
    stop
  end if

  call mympi_swpbnd (ter_r3d_old)

end function ter_r3d_old

function ter_r2d_old(var, hga, hgb, dft) !{{{1
  ! interpolate from grid hga to grid hgb
  real (kind=wp), intent(in) :: var(:,:)
  type (type_stg), target :: hga, hgb
  real (kind=wp), optional :: dft

  real (kind=wp), dimension(ni,nj) :: ter_r2d_old

  if ( present(dft) ) then
    ter_r2d_old = dft
  else
    ter_r2d_old = 0.0
  end if

  if ( hga%n == hgb%n ) then
    stop 'no need to interpolate in ter_r2d_old of mod_op'
  else if ( associated(hga%ew, hgb) ) then
    ter_r2d_old(1+hga%i:nim+hga%i,:) = &
    0.5 * ( var(1:nim,:) + var(2:ni,:) )
  else if ( associated(hga%ns, hgb) ) then
    ter_r2d_old(:,1+hga%j:njm+hga%j) = &
      0.5 * ( var(:,1:njm) + var(:,2:nj) )
  else if ( associated(hga%di, hgb) ) then
    ter_r2d_old(1+hga%i:nim+hga%i,1+hga%j:njm+hga%j) = &
      ( var(1:nim,1:njm) + var(1:nim,2:nj) + &
        var(2:ni, 1:njm) + var(2:ni, 2:nj) ) * 0.25
  else
    print *, 'unhandled relative position in ter_r2d_old of mod_op'
    stop
  end if

  call mympi_swpbnd (ter_r2d_old)

end function ter_r2d_old

function ter_i3d_old(var, hga, hgb, dft) !{{{1
  ! interpolate from grid hga to grid hgb
  integer, intent(in) :: var(:,:,:)
  type (type_stg), target :: hga, hgb
  integer, optional :: dft
  integer, allocatable:: ter_i3d_old(:,:,:)

  integer :: k, d3, is

  if ( hga%n == hgb%n ) stop 'no need to interpolate in ter_i3d_old of mod_op'

  d3 = size(var, 3)
  allocate( ter_i3d_old(ni,nj,d3), stat = is ); call chk(is)

  if ( present(dft) ) then
    ter_i3d_old = dft
  else
    ter_i3d_old = 0
  end if

  if      ( associated(hga%ew, hgb) ) then
    ter_i3d_old(1+hga%i:nim+hga%i,:,:) = &
    0.5 * ( var(1:nim,:,:) + var(2:ni,:,:) )
  else if ( associated(hga%ns, hgb) ) then
    ter_i3d_old(:,1+hga%j:njm+hga%j,:) = &
      0.5 * ( var(:,1:njm,:) + var(:,2:nj,:) )
  else if ( associated(hga%di, hgb) ) then
    ter_i3d_old(1+hga%i:nim+hga%i,1+hga%j:njm+hga%j,:) = &
      ( var(1:nim,1:njm,:) + var(1:nim,2:nj,:) + &
        var(2:ni, 1:njm,:) + var(2:ni, 2:nj,:) ) * 0.25
  else
    print *, 'unhandled relative position in ter_i3d_old of mod_op'
    stop
  end if

  call mympi_swpbnd (ter_i3d_old)

end function ter_i3d_old

function ter_i2d_old(var, hga, hgb, dft) !{{{1
  ! interpolate from grid hga to grid hgb
  integer, intent(in) :: var(:,:)
  type (type_stg), target :: hga, hgb
  integer, optional :: dft

  integer, dimension(ni,nj) :: ter_i2d_old

  if ( hga%n == hgb%n ) stop 'no need to interpolate in ter_i2d_old of mod_op'

  if ( present(dft) ) then
    ter_i2d_old = dft
  else
    ter_i2d_old = 0
  end if

  if      ( associated(hga%ew, hgb) ) then
    ter_i2d_old(1+hga%i:nim+hga%i,:) = &
    0.5 * ( var(1:nim,:) + var(2:ni,:) )
  else if ( associated(hga%ns, hgb) ) then
    ter_i2d_old(:,1+hga%j:njm+hga%j) = &
      0.5 * ( var(:,1:njm) + var(:,2:nj) )
  else if ( associated(hga%di, hgb) ) then
    ter_i2d_old(1+hga%i:nim+hga%i,1+hga%j:njm+hga%j) = &
      ( var(1:nim,1:njm) + var(1:nim,2:nj) + &
        var(2:ni, 1:njm) + var(2:ni, 2:nj) ) * 0.25
  else
    print *, 'unhandled relative position in ter_i2d_old of mod_op'
  end if

  call mympi_swpbnd (ter_i2d_old)

end function ter_i2d_old

function px1_r3d_old ( var, ga, gb )  !{{{1
  ! (partial var) / (partial x1), not include metric effects
  ! var on grid ga, result on grid pb
  real (kind=wp), intent(in) :: var(:,:,:)
  type (type_stg), target :: ga, gb
  real (kind=wp), allocatable, dimension(:,:,:) :: &
    px1_r3d_old, temp

  integer :: d3, is

  d3 = size(var, 3)

  allocate(px1_r3d_old(ni,nj,d3), stat = is); call chk(is)
  allocate(temp(ni,nj,d3), stat = is); call chk(is)

  ! interpolate to the matching grid first
  if ( .not. associated(gb%ew, ga) ) then
    call ter_r3d(temp, var, ga, gb%ew )
  else
    temp = var
  end if

  px1_r3d_old = 0.0

  px1_r3d_old(1+ga%i:nim+ga%i,:,:) = &
    temp(2:ni,:,:) - temp(1:nim,:,:)

  px1_r3d_old = px1_r3d_old / spread(gb%dx%x(1),3,d3)

  call mympi_swpbnd( px1_r3d_old)

end function px1_r3d_old

function px1_r2d_old ( var, ga, gb )  !{{{1
  ! (partial var) / (partial x1), not include metric effects
  ! var on grid ga, result on grid pb
  real (kind=wp), intent(in) :: var(:,:)
  type (type_stg), target :: ga, gb
  real (kind=wp), dimension(ni,nj) :: px1_r2d_old

  real (kind=wp), dimension(ni,nj) :: temp

  ! interpolate to the matching grid first
  if ( .not. associated(gb%ew, ga) ) then
    call ter_r2d( temp, var, ga, gb%ew )
  else
    temp = var
  end if

  px1_r2d_old = 0.0

  px1_r2d_old(1+ga%i:nim+ga%i,:) = temp(2:ni,:) - temp(1:nim,:)

  px1_r2d_old = px1_r2d_old / gb%dx%x(1)

  call mympi_swpbnd( px1_r2d_old)
end function px1_r2d_old

function px2_r3d_old ( var, ga, gb )  !{{{1
  ! (partial var) / (partial x2), not include metric effects
  ! var on grid ga, result on grid pb
  real (kind=wp), intent(in) :: var(:,:,:)
  type (type_stg), target :: ga, gb
  real (kind=wp), allocatable, dimension(:,:,:) :: &
    px2_r3d_old, temp

  integer :: d3, is

  d3 = size(var, 3)

  allocate(px2_r3d_old(ni,nj,d3), stat = is); call chk(is)
  allocate(temp(ni,nj,d3), stat = is); call chk(is)

  ! interpolate to the matching grid first
  if ( associated(gb%ns, ga) ) then
    temp = var
  else
    call ter_r3d( temp, var, ga, gb%ns )
  end if

  px2_r3d_old = 0.0

  px2_r3d_old(:,1+ga%j:njm+ga%j,:) = temp(:,2:nj,:) - temp(:,1:njm,:)
  px2_r3d_old = px2_r3d_old / spread(gb%dx%x(2),3,d3)

  call mympi_swpbnd( px2_r3d_old)
    
end function px2_r3d_old

function px2_r2d_old ( var, ga, gb )  !{{{1
  ! (partial var) / (partial x2), not include metric effects
  ! var on grid ga, result on grid pb
  real (kind=wp), intent(in) :: var(:,:)
  type (type_stg), target :: ga, gb
  real (kind=wp), dimension(ni,nj) :: px2_r2d_old

  real (kind=wp), dimension(ni,nj) :: temp

  ! interpolate to the matching grid first
  if ( .not. associated(gb%ns, ga) ) then
    call ter_r2d( temp, var, ga, gb%ns )
  else
    temp = var
  end if

  px2_r2d_old = 0.0

  px2_r2d_old(:,1+ga%j:njm+ga%j) = temp(:,2:nj) - temp(:,1:njm)

  px2_r2d_old = px2_r2d_old / gb%dx%x(2)

  call mympi_swpbnd( px2_r2d_old)

end function px2_r2d_old

function px3_r3d ( var )  !{{{1
  ! default (partial var) / (partial x3)
  ! upward is positive
  ! vertically, var on vg1, result on vg2
  real (kind=wp), dimension(ni,nj,nkp) :: var
  type (type_stg), target :: hga, hgb
  real (kind=wp), dimension(ni,nj,nk) :: px3_r3d

  integer :: d3

  d3 = size(var, 3)

  if (d3 == nkp) then
    px3_r3d = 0.0
    px3_r3d(:,:,1:nk) = var(:,:,1:nk) - var(:,:,2:nkp)
    px3_r3d = px3_r3d / spread(spread(vg2%dp,1,nj), 1, ni)
  else
    stop 'var should vertically on vg1 in px3_r3d in mod_op'
  end if

  ! no need to swap boundary horizontally
end function px3_r3d

function px3_r3d_ter ( var, hga, hgb )  !{{{1
  ! (partial var) / (partial x3)
  ! upward is positive
  ! var horizontally on hga, vertically on vg1
  ! result horizontally in hgb, vertically on vg2
  real (kind=wp), dimension(ni,nj,nkp) :: var
  type (type_stg), target :: hga, hgb
  real (kind=wp), dimension(ni,nj,nk) :: px3_r3d_ter

  real (kind=wp), dimension(ni,nj,nkp) :: temp
  integer :: d3

  d3 = size(var, 3)
  if (d3 /= nkp) &
    stop 'var should vertically on vg1 in px3_r3d_ter in mod_op'

  call ter_r3d( temp, var, hga, hgb )

  px3_r3d_ter = 0.0
  px3_r3d_ter(:,:,1:nk) = temp(:,:,1:nk) - temp(:,:,2:nkp)
  px3_r3d_ter = px3_r3d_ter / spread(spread(vg2%dp,1,nj), 1, ni)

  ! no need to swap boundary horizontally
end function px3_r3d_ter

function div_m3d (var, ga, gb, gc) !{{{1
  ! horizontal divergence operator for 3d vector
  ! var on grid (ga, gb), result on grid gc
  type (type_mat), dimension(:,:,:), intent(in) :: var
  type (type_stg), intent(in) :: ga, gb, gc
  real (kind=wp), dimension(ni,nj,nk) :: div_m3d

  div_m3d = div_r3d (var%x(1), var%x(2), ga, gb, gc)

end function div_m3d

function div_r3d (va, vb, ga, gb, gc) !{{{1
  ! horizontal divergence operator for 3d vector (va, vb)
  ! (va, vb) on grid (ga, gb), result on grid gc
  real (kind=wp), dimension(ni,nj,nk), intent(in) :: va, vb
  type (type_stg), intent(in) :: ga, gb, gc

  real (kind=wp), dimension(ni,nj,nk) :: div_r3d, wka, wkb, wk

  wk = spread(gb%rh1,3,nk)
  call px1_r3d( wka, va, ga, gc )
  call px2_r3d( wkb, vb*wk, gb, gc)
  div_r3d = wka + wkb

  wk = spread(gc%rh1,3,nk)
  div_r3d = div_r3d / (a * wk)

end function div_r3d

function div_m2d (var, ga, gb, gc) !{{{1
  ! horizontal divergence operator for 2d vector (va, vb)
  ! (va, vb) on grid (ga, gb), result on grid gc
  type (type_mat), dimension(:,:), intent(in) :: var
  type (type_stg), intent(in) :: ga, gb, gc
  real (kind=wp), dimension(ni,nj) :: div_m2d, wk

  call px1_r2d( wk, var%x(1), ga, gc )
  div_m2d = wk
  call px2_r2d( wk, var%x(2)*gb%rh1, gb, gc )
  div_m2d = div_m2d + wk

  div_m2d = div_m2d / (a * gc%rh1)

end function div_m2d

function div_r2d (va, vb, ga, gb, gc) !{{{1
  ! horizontal divergence operator for 2d vector (va, vb)
  ! (va, vb) on grid (ga, gb), result on grid gc
  real (kind=wp), dimension(ni,nj), intent(in) :: va, vb
  type (type_stg), intent(in) :: ga, gb, gc

  real (kind=wp), dimension(ni,nj) :: div_r2d, wk

  call px1_r2d( wk, va, ga, gc )
  div_r2d = wk
  call px2_r2d( wk, vb*gb%rh1, gb, gc )
  div_r2d = div_r2d + wk

  div_r2d = div_r2d / (a * gc%rh1)

end function div_r2d

function gra_r3d (var, ga, gb, gc) !{{{1
  ! horizontal gradient operator for 3d scalar
  ! var on grid ga, output on (gb, gc)
  real (kind=wp),  dimension(:,:,:), intent(in) :: var
  type (type_stg) :: ga, gb, gc
  type (type_mat), allocatable :: gra_r3d(:,:,:)

  integer :: k, d3, is

  d3 = size(var, 3)

  allocate( gra_r3d(ni,nj,d3), stat = is); call chk(is)

  call px1_r3d( gra_r3d%x(1), var, ga, gb)
  gra_r3d%x(1) = gra_r3d%x(1) / ( a * spread(gb%rh1,3,d3) )
  call px2_r3d( gra_r3d%x(2), var, ga, gc)
  gra_r3d%x(2) = gra_r3d%x(2) / a

end function gra_r3d

function gra_r2d (var, hga, hgb, hgc) !{{{1
  ! horizontal gradient operator
  ! var on grid hga, output on (hgb, hgc)
  real (kind=wp),  dimension(ni,nj), intent(in) :: var
  type (type_stg), target :: hga, hgb, hgc
  type (type_gvar_m2d) :: gra_r2d

  integer :: is

  allocate(gra_r2d%x(1)%v(ni,nj), stat=is); call chk(is)
  allocate(gra_r2d%x(2)%v(ni,nj), stat=is); call chk(is)
  gra_r2d%x(1)%hg => hgb; gra_r2d%x(2)%hg => hgc

  call px1_r2d( gra_r2d%x(1)%v, var, hga, hgb )
  gra_r2d%x(1)%v = gra_r2d%x(1)%v / ( a * hgb%rh1 )
  call px2_r2d( gra_r2d%x(2)%v, var, hga, hgc )
  gra_r2d%x(2)%v = gra_r2d%x(2)%v / a

end function gra_r2d

function lap_r2d (var, ga, gb) !{{{1
  ! horizontal Laplacian operator
  ! var on grid ga, output on grid gb
  real (kind=wp), dimension(:,:), intent(in) :: var
  type (type_stg), intent(in) :: ga, gb

  real (kind=wp), dimension(ni,nj) :: lap_r2d, wka, wkb
  
  call px1_r2d( wka, var, ga, gb%ew )
  wka = wka / gb%ew%rh1 
  call px1_r2d( wkb, wka, gb%ew, gb )
  lap_r2d = wkb

  call px2_r2d( wka, var, ga, gb%ns )
  wka = wka * gb%ns%rh1
  call px2_r2d( wkb, wka, gb%ns, gb )
  lap_r2d = lap_r2d + wkb

  lap_r2d = lap_r2d / ( a*a*gb%rh1 )

end function lap_r2d

function vint_r3d (var, grd)!{{{1
  ! vertical integration from sea surface to sea bottom
  ! result divided by a factor of pb
  real (kind=wp), intent(in) :: var(:,:,:)
  type (type_stg3d) :: grd

  real (kind=wp), dimension(ni,nj) :: vint_r3d
  real (kind=wp), dimension(ni,nj,nk) :: dp3d

  dp3d = spread( spread(grd%vg%dp(1:nk),1,nj), 1, ni )
  vint_r3d = sum(var * dp3d * grd%msk, 3) / grd%pb

end function vint_r3d

function vint_m3d (var, ga, gb)!{{{1
  ! vertical integration from the first to the bottom layer
  type (type_mat), intent(in) :: var(:,:,:)
  type (type_stg3d) :: ga, gb
  type (type_mat), dimension(ni,nj) :: vint_m3d

  vint_m3d%x(1) = vint_r3d( var%x(1), ga )
  vint_m3d%x(2) = vint_r3d( var%x(2), gb )

end function vint_m3d

function r3d_div_r2d (va, vb)!{{{1
  ! 2d real array multiply the same dimension matrix array
  real (kind=wp),  dimension(:,:,:), intent(in) :: va
  real (kind=wp),  dimension(:,:), intent(in) :: vb
  real (kind=wp), allocatable, dimension(:,:,:) :: r3d_div_r2d

  integer :: d(3), is

  d = shape(va)

  allocate( r3d_div_r2d(d(1),d(2),d(3)), stat=is )
  call chk(is)

  r3d_div_r2d = va / spread(vb,3,d(3))

end function r3d_div_r2d

function m2d_div_i (m2d, i)!{{{1
  ! reload divided operator
  type (type_mat), dimension(:,:), intent(in) :: m2d
  integer, intent(in) :: i

  type (type_mat), allocatable, dimension(:,:) :: m2d_div_i

  integer :: d1, d2, is

  d1 = size(m2d, 1)
  d2 = size(m2d, 2)

  allocate( m2d_div_i(d1, d2), stat=is ); call chk(is)

  m2d_div_i%x(1) = m2d%x(1) / i
  m2d_div_i%x(2) = m2d%x(2) / i

end function m2d_div_i

pure function m2d_div_r (m2d, r)!{{{1
  ! reload divided operator
  type (type_mat), dimension(:,:), intent(in) :: m2d
  real (kind=wp), intent(in) :: r

  type (type_mat), allocatable, dimension(:,:) :: m2d_div_r

  integer :: d1, d2

  d1 = size(m2d, 1)
  d2 = size(m2d, 2)

  allocate( m2d_div_r(d1, d2) )

  m2d_div_r%x(1) = m2d%x(1) / r
  m2d_div_r%x(2) = m2d%x(2) / r

end function m2d_div_r

pure function m2d_div_r2d (va, vb)!{{{1
  ! reload divided operator
  type (type_mat), dimension(:,:), intent(in) :: va
  real (kind=wp), intent(in) :: vb(:,:)

  type (type_mat), allocatable, dimension(:,:) :: m2d_div_r2d

  integer :: d1, d2

  d1 = size(va, 1)
  d2 = size(va, 2)

  allocate( m2d_div_r2d(d1, d2) )

  m2d_div_r2d%x(1) = va%x(1) / vb
  m2d_div_r2d%x(2) = va%x(2) / vb

end function m2d_div_r2d

pure function m3d_minus_m3d (va, vb)!{{{1
  ! 2d matrix minus 2d matrix
  type (type_mat), dimension(:,:,:), intent(in) :: va, vb

  type (type_mat), allocatable, dimension(:,:,:) :: m3d_minus_m3d

  integer :: d1, d2, d3

  d1 = size(va, 1)
  d2 = size(va, 2)
  d2 = size(va, 3)

  allocate( m3d_minus_m3d(d1, d2, d3) )

  m3d_minus_m3d%x(1) = va%x(1)- vb%x(1)
  m3d_minus_m3d%x(2) = va%x(2)- vb%x(2)

end function m3d_minus_m3d

pure function m2d_minus_m2d (va, vb)!{{{1
  ! 2d matrix minus 2d matrix
  type (type_mat), dimension(:,:), intent(in) :: va, vb

  type (type_mat), allocatable, dimension(:,:) :: m2d_minus_m2d

  integer :: d1, d2

  d1 = size(va, 1)
  d2 = size(va, 2)

  allocate( m2d_minus_m2d(d1, d2) )

  m2d_minus_m2d%x(1) = va%x(1)- vb%x(1)
  m2d_minus_m2d%x(2) = va%x(2)- vb%x(2)

end function m2d_minus_m2d

pure function m1d_minus_m1d (va, vb)!{{{1
  ! 2d matrix minus 2d matrix
  type (type_mat), dimension(:), intent(in) :: va, vb

  type (type_mat), allocatable, dimension(:) :: m1d_minus_m1d

  integer :: d1

  d1 = size(va)

  allocate( m1d_minus_m1d(d1) )

  m1d_minus_m1d%x(1) = va%x(1)- vb%x(1)
  m1d_minus_m1d%x(2) = va%x(2)- vb%x(2)

end function m1d_minus_m1d

pure function m_plus_m (va, vb)!{{{1
  ! 2d matrix minus 2d matrix
  type (type_mat), intent(in) :: va, vb

  type (type_mat) :: m_plus_m

  m_plus_m%x = va%x + vb%x

end function m_plus_m

pure function m2d_plus_m2d (va, vb)!{{{1
  ! 2d matrix minus 2d matrix
  type (type_mat), dimension(:,:), intent(in) :: va, vb

  type (type_mat), allocatable, dimension(:,:) :: m2d_plus_m2d

  integer :: d1, d2

  d1 = size(va, 1)
  d2 = size(va, 2)

  allocate( m2d_plus_m2d(d1, d2) )

  m2d_plus_m2d%x(1) = va%x(1) + vb%x(1)
  m2d_plus_m2d%x(2) = va%x(2) + vb%x(2)

end function m2d_plus_m2d

subroutine fri_m3d (ans, spbt, u, v, up, vp, taux, tauy) !{{{1
  ! horizontal frictional force
  ! output is on the same grid as vp
  type (type_mat), dimension(:,:,:) :: ans
  real (kind=wp), dimension(:,:,:), intent(in) :: u, v, up, vp
  real (kind=wp), dimension(:,:), intent(in) :: spbt, taux, tauy

  real (kind=wp), dimension(ni,nj,nk) :: wk

  wk = - v
  call fri_r3d( ans%x(1), spbt, u, up, wk, taux )
  call fri_r3d( ans%x(2), spbt, v, vp, u,  tauy )

end subroutine fri_m3d

subroutine fri_gr3d (ans, spbt, v, vp, vc, tau) !{{{1
  ! horizontal frictional force
  ! ans is on the same grid as vp
  ! ans = 1/spbt*div(pbt*gra(v)) + cv1*vp + cv2*spbt*(p vc/ p x1) + vert.
  type (type_gvar_r3d) :: ans
  real (kind=wp), dimension(:,:), intent(in) :: spbt
  type (type_gvar_r3d), intent(in) :: &
    v, & ! unweighted velocity 
    vp, & ! pressure weighted velocity
    vc
  real (kind=wp), dimension(ni,nj), intent(in) :: tau

  type (type_gvar_r3d) :: wkgr3d
  type (type_gvar_m3d) :: wkgm3d
  real (kind=wp), dimension(ni,nj,nk) :: wk, wkb, spbt3d, pbt3da, pbt3db
  real (kind=wp), dimension(ni,nj) :: wkr2d
  real (kind=wp) :: mag
  integer :: is, i, j, k

  call arrays_cp_shape( v, wkgm3d%x(1) )
  call arrays_cp_shape( v, wkgm3d%x(2) )
  ans%g => vp%g
  spbt3d = spread(spbt, 3, nk)

  ! horizontal viscosity

  ! div(pbt*gra(v))
  call op_gra( wkgm3d, v )
  call op_ter( wkr2d, pbt%tc, pbt%hg, wkgm3d%x(1)%g%hg )
  pbt3da = spread( wkr2d, 3, nk )
  call op_ter( wkr2d, pbt%tc, pbt%hg, wkgm3d%x(2)%g%hg )
  pbt3db = spread( wkr2d, 3, nk )
  call div_r3d( wk, pbt3da*wkgm3d%x(1)%v, pbt3db*wkgm3d%x(2)%v, &
                wkgm3d%x(1)%g%hg, wkgm3d%x(2)%g%hg, vp%g%hg )
  ans%v = am%v/spbt3d * wk

  ! cv1*vp
  ans%v = ans%v + am%v*spread(cv1,3,nk)*vp%v

  ! cv2*spbt*(p vc/ p x1)
  call ter_r3d( wk, vc%v, vc%g%hg, vc%g%hg%ew )
  call op_px1( wkb, wk, vc%g%hg%ew, vp%g%hg )
  ans%v = ans%v + am%v*spread(cv2,3,nk)*spbt3d*wkb

  ! vertical viscosity

  ! oddly vertical difference respect to thickness (not pressure)
  allocate( wkgr3d%v(ni,nj,nkp), stat = is ); call chk(is)
  wkgr3d%g => vp%g%ud
  wkgr3d%v(:,:,2:nk) = vp%v(:,:,1:nk-1) - vp%v(:,:,2:nk)
  wkgr3d%v = wkgr3d%v / spread(spread(wkgr3d%g%vg%dz,1,nj), 1, ni)
  wkgr3d%v = g*rho0/(spread(spbt,3,nk))**2 * km%v * wkgr3d%v

  ! surface boundary condition
  wkgr3d%v(:,:,1) = g*tau(:,:) / spbt(:,:)

  do i = 1, ni
  do j = 1, nj
    k = vp%g%lev(i,j)
    if ( k > 0 ) then ! bottom drag
      mag = sqrt( v%v(i,j,k)**2 + vc%v(i,j,k)**2 )
      wkgr3d%v(i,j,k+1) = g*cdbot*mag/spbt(i,j)
    end if
  end do
  end do

  call px3_r3d( wk, wkgr3d%v )
  ans%v = ans%v + wk

  ! to prevent memory leakage
  call arrays_free( wkgr3d ) 
  call arrays_free( wkgm3d ) 

end subroutine fri_gr3d

subroutine fri_gm3d (ans, spbt, v, vp) !{{{1
  ! horizontal frictional force
  ! output is on the same grid as vp
  real (kind=wp), dimension(:,:), intent(in) :: spbt
  type (type_gvar_m3d), intent(in) :: v, vp
  type (type_gvar_m3d) :: ans
  type (type_gvar_r3d) :: wk

  call arrays_cp_shape( v%x(2), wk )
  wk%v = - v%x(2)%v

  call fri_gr3d( ans%x(1), spbt, v%x(1), vp%x(1), &
        wk, bnd%tau%x(1)%v )
  call fri_gr3d( ans%x(2), spbt, v%x(2), vp%x(2), &
    v%x(1), bnd%tau%x(2)%v )

  call arrays_free( wk )
end subroutine fri_gm3d

subroutine adv_gm3d (ans, var, u, w) !{{{1
  ! calc. advection of horizontal pressure weighted velocity var
  ! adv_gm3d(var) = div(var*v) - 0.5*var*div(v), 
  !   v = (u, w) are unweighted 3d velocity
  type (type_gvar_m3d) :: ans
  type (type_gvar_m3d), intent(in) :: var, u
  type (type_gvar_r3d), intent(in) :: w

  call adv_gr3d(ans%x(1), var%x(1), u, w)
  call adv_gr3d(ans%x(2), var%x(2), u, w)
end subroutine adv_gm3d

subroutine adv_gr3d (ans, var, u, w) !{{{1
  ! calc. advection of horizontal pressure weighted velocity var
  ! adv_gr3d(var) = div(var*v) - 0.5*var*div(v), 
  !   v = (u, w) are unweighted 3d velocity
  type (type_gvar_r3d) :: ans
  type (type_gvar_r3d), intent(in) :: var, w
  type (type_gvar_m3d), intent(in) :: u

  type (type_mat), dimension(ni,nj,nk) :: wkmat
  real (kind=wp),  dimension(ni,nj,nkp):: wka, wkb
  real (kind=wp),  dimension(ni,nj,nk) :: wkc, wkd
  type (type_stg), pointer :: hg

  ans%g => var%g

  hg => var%g%hg

  ! var*u
  call ter_r3d( wkmat%x(1), var%v, hg, hg%ew )
  call ter_r3d( wkc,   u%x(1)%v, hg, hg%ew )
  wkmat%x(1) = wkmat%x(1) * wkc

  call ter_r3d( wkmat%x(2), var%v, hg, hg%ns )
  call ter_r3d( wkc,   u%x(2)%v, hg, hg%ns )
  wkmat%x(2) = wkmat%x(2) * wkc

  ! var*w
  call vter_r3d( wkb, var%v, var%g%vg, var%g%vg%ud )
  call ter_r3d( wka, w%v, w%g%hg, hg )
  wkb = wkb * wka

  ! div(v)
  call op_px3( wkc, w%v, w%g%hg, hg )
  call div_r3d( wkd, u%x(1)%v, u%x(2)%v, hg, hg, hg )
  wkc = wkc + wkd

  call px3_r3d( wkd, wkb )
  call div_r3d( ans%v, wkmat%x(1), wkmat%x(2), hg%ew, hg%ns, hg )
  ans%v = ans%v + wkd - 0.5*var%v*wkc
end subroutine adv_gr3d

subroutine op_adv_old (ans, var, u, v, w) !{{{1
  ! calc. advection of horizontal pressure weighted velocity var
  ! op_adv(var) = div(var*v) - 0.5*var*div(v), 
  !   v = (u, w) are unweighted 3d velocity
  real (kind=wp), dimension(:,:,:) :: ans, var, u, v, w

  type (type_mat), dimension(ni,nj,nk) :: wkm
  real (kind=wp),  dimension(ni,nj,nk) :: wka, wkb, cos3d
  real (kind=wp),  dimension(ni,nj,nkp):: wkap, wkbp

  ! var*u
  call ter_r3d( wka, var, hgu, hgu%ew )
  call ter_r3d( wkb,   u, hgu, hgu%ew )
  wkm%x(1) = wka * wkb

  cos3d = spread(hgu%rh1,3,nk)
  call ter_r3d( wka, var,     hgu, hgu%ns )
  call ter_r3d( wkb, v*cos3d, hgu, hgu%ns )
  wkm%x(2) = wka * wkb

  ! horizontal: div(var*u)
  call px1_r3d( wka, wkm%x(1), hgu%ew, hgu )
  call px2_r3d( wkb, wkm%x(2), hgu%ns, hgu )
  ans = (wka + wkb)*gu%msk / (a*cos3d)

  ! var*w
  call vter_r3d( wkap, var, vgt, vgw )
  call  ter_r3d( wkbp,   w, hgt, hgu )

  ! set to zero if any of the upper and lower layer of U-grid is land
  wkap = wkap*wkbp
  wkbp(:,:,1)   = 0
  wkbp(:,:,nkp) = 0
  wkbp(:,:,2:nkp-1) = wkap(:,:,2:nkp-1)*&
                      gu%msk(:,:,1:nk-1)*gu%msk(:,:,2:nk)

  call px3_r3d( wka, wkbp )
  ans = ans + wka

  ! div(v)
  call  op_px3( wka, w, hgt, hgu )
  call div_r3d_b( wkb, u, v, hgu, hgu, hgu, gu%msk )

  ans = ans - 0.5*var* (wka + wkb)
  ans = ans * gu%msk
end subroutine op_adv_old

subroutine p1_r3d_center ( ans, var, hg )  !{{{1
  ! (partial var) / (partial x1), not include metric effects
  ! var on grid hg, result also on grid hg (center scheme)
  real (kind=wp) :: ans(:,:,:)
  real (kind=wp), intent(in) :: var(:,:,:)
  type (type_stg), intent(in) :: hg

  integer :: i,j,k,d3

  d3 = size(var,3)

  ans = 0.0

  do i = 2, ni-1
  do j = 1, nj
  do k = 1, d3
    ans(i,j,k) = ( var(i+1,j,k) - var(i-1,j,k) ) / &
      ( hg%ew%dx(i,j)%x(1) + hg%ew%dx(i+1,j)%x(1) )
  end do
  end do
  end do

  call mympi_swpbnd(ans)

end subroutine p1_r3d_center

subroutine adv_gr3d_ts (ans, var, um, wm) !{{{1
  ! calc. 3d advection of tracer var
  ! adv_gr3d_ts(var) = vm * grad( var )
  !   vm = (um, wm) are mass weighted 3d velocity
  type (type_gvar_r3d) :: ans
  type (type_gvar_r3d), intent(in) :: var, wm
  type (type_gvar_m3d), intent(in) :: um

  type (type_gvar_m3d) :: wkgm
  type (type_mat), dimension(ni,nj,nk) :: wkmat
  real (kind=wp),  dimension(ni,nj,nkp):: wka, wkb
  real (kind=wp),  dimension(ni,nj,nk) :: wkc, wkd

  ans%g => var%g

  call arrays_cp_shape( var, wkgm )

  ! horizontal advection
!  call gra_gr3d( wkgm, var )
  wkgm%x(1)%g => var%g%ew
  call p1_r3d( wkgm%x(1)%v, var%v, var%g%hg, wkgm%x(1)%g%hg )
  wkgm%x(2)%g => var%g%ns
  call p2_r3d( wkgm%x(2)%v, var%v, var%g%hg, wkgm%x(2)%g%hg )

  call ter_r3d( wkc, um%x(1)%v, hgu, wkgm%x(1)%g%hg )
  call ter_r3d( wkd, wkc * wkgm%x(1)%v, wkgm%x(1)%g%hg, hgt )
  ! not interpolate the cosine factor
  wkd = wkd / (a * spread(hgt%rh1,3,nk))
  ans%v = wkd

  call ter_r3d( wkc, um%x(2)%v, hgu, wkgm%x(2)%g%hg )
  call ter_r3d( wkd, wkc * wkgm%x(2)%v, wkgm%x(2)%g%hg, hgt )
  ! not interpolate the cosine factor
  wkd = wkd / (a * spread(hgt%rh1,3,nk))
  ans%v = ans%v + wkd

  ! vertical advection
  call dx3_r3d_b( wka, var%v )
  call vter_r3d( wkc, wm%v * wka, wm%g%vg, ans%g%vg ) 
  ! not interpolate the layer 'thickness'
  wkc = wkc / spread(spread(ans%g%vg%dp,1,nj), 1, ni)
  ans%v = ans%v + wkc

  call arrays_free( wkgm )
end subroutine adv_gr3d_ts

subroutine adv_gm3d_ts (ans, var, um, wm) !{{{1
  ! calc. 3d advection of tracer var
  ! adv_gr3d_ts(var) = vm * grad( var )
  !   vm = (um, wm) are mass weighted 3d velocity
  type (type_gvar_m3d) :: ans
  type (type_gvar_m3d), intent(in) :: var, um
  type (type_gvar_r3d), intent(in) :: wm

  call adv_gr3d_ts(ans%x(1), var%x(1), um, wm)
  call adv_gr3d_ts(ans%x(2), var%x(2), um, wm)
end subroutine adv_gm3d_ts

subroutine dif_gr3d (ans, pbt, var) !{{{1
  ! diffusion of tracers
  ! ans is on the same grid as var
  ! ans = div(pbt*gra(var)) + vert.
  type (type_gvar_r3d) :: ans
  real (kind=wp), dimension(:,:), intent(in) :: pbt
  type (type_gvar_r3d), intent(in) :: var ! tracer

  type (type_gvar_r3d) :: wkgr3d
  type (type_gvar_m3d) :: wkgm3d
  real (kind=wp), dimension(ni,nj,nk) :: wk, pbt3da, pbt3db
  real (kind=wp), dimension(ni,nj) :: wkr2d
  integer :: is, k

  call arrays_cp_shape( var, wkgm3d%x(1) )
  call arrays_cp_shape( var, wkgm3d%x(2) )
  ans%g => var%g

  ! horizontal diffusion

  ! div(pbt*gra(var))
  call gra_gr3d_b( wkgm3d, var )
  call op_ter( wkr2d, pbt, hgt, wkgm3d%x(1)%g%hg )
  pbt3da = spread( wkr2d, 3, nk )
  call op_ter( wkr2d, pbt, hgt, wkgm3d%x(2)%g%hg )
  pbt3db = spread( wkr2d, 3, nk )
  call div_r3d( wk, pbt3da*wkgm3d%x(1)%v, pbt3db*wkgm3d%x(2)%v, &
                wkgm3d%x(1)%g%hg, wkgm3d%x(2)%g%hg, hgt )
  ans%v = ah_c * wk

  ! vertical diffusion

  ! oddly vertical difference respect to thickness (not pressure)
  allocate( wkgr3d%v(ni,nj,nkp), stat = is ); call chk(is)
  wkgr3d%g => var%g%ud
  wkgr3d%v(:,:,1) = 0.0
  wkgr3d%v(:,:,nkp) = 0.0
  wkgr3d%v(:,:,2:nk) = ( var%v(:,:,1:nk-1) - var%v(:,:,2:nk) ) * &
                       var%g%msk(:,:,1:nk-1) * var%g%msk(:,:,2:nk)
  wkgr3d%v = wkgr3d%v / spread(spread(wkgr3d%g%vg%dz,1,nj), 1, ni)
  wkgr3d%v = g*rho0 * kh%v * wkgr3d%v

  call p3_r3d( wk, wkgr3d%v )
  ans%v = ans%v + wk

  ! to prevent memory leakage
  call arrays_free( wkgr3d ) 
  call arrays_free( wkgm3d ) 

end subroutine dif_gr3d

subroutine dif_gm3d (ans, pbt, var) !{{{1
  ! horizontal diffusion for potential temperature and salinity
  ! output is on the same grid as var
  real (kind=wp), dimension(:,:), intent(in) :: pbt
  type (type_gvar_m3d), intent(in) :: var
  type (type_gvar_m3d) :: ans

  call dif_gr3d( ans%x(1), pbt, var%x(1))
  call dif_gr3d( ans%x(2), pbt, var%x(2))

end subroutine dif_gm3d

subroutine gra_gr3d_b (ans, var) !{{{1
  ! gradient of var
  ! set to zero if any of the two points is missing
  type (type_gvar_m3d) :: ans
  type (type_gvar_r3d), intent(in) :: var

  integer :: d3

  d3 = size(var%v, 3)

  ans%x(1)%g => var%g%ew
  call p1_r3d_b( ans%x(1)%v, var%v, var%g%hg, ans%x(1)%g%hg, var%g%msk )
  ans%x(1)%v = ans%x(1)%v / ( a * spread(ans%x(1)%g%hg%rh1, 3, d3) )

  ans%x(2)%g => var%g%ns
  call p2_r3d_b( ans%x(2)%v, var%v, var%g%hg, ans%x(2)%g%hg, var%g%msk )
  ans%x(2)%v = ans%x(2)%v / a

end subroutine gra_gr3d_b
subroutine div_r3d_b (ans, va, vb, ga, gb, gc, mask) !{{{1
  ! horizontal divergence operator for 3d vector (va, vb)
  ! (va, vb) on grid (ga, gb), result on grid gc
  ! mask out lands
  real (kind=wp), dimension(ni,nj,nk) :: ans
  real (kind=wp), dimension(ni,nj,nk), intent(in) :: va, vb
  type (type_stg), intent(in) :: ga, gb, gc
  integer, dimension(ni,nj,nk) :: mask

  real (kind=wp), dimension(ni,nj,nk) :: wka, wkb, wk

  wk = spread(gb%rh1,3,nk)
  call p1_r3d( wka, va, ga, gc )
  call p2_r3d( wkb, vb*wk, gb, gc)
!  ans = wka + wkb

  wk = spread(gc%rh1,3,nk)
  ans = ans*mask / (a * wk)

end subroutine div_r3d_b

subroutine test ( ) !{{{1

  call debug_output ()
!  if (myid == mid) &
!    call mympi_quick_output('temp.nc', 'var', km%v(:,:,2:nkp), glo_lon, glo_lat, z)
  if (myid==mid) stop 'finish debug in test of main'

end subroutine test

subroutine debug_output () !{{{1

  if (myid == mid) then
    call io_quick_output( 'check_output/vx1.nc', 'vx1', vg1%z, vg1%z )
    call io_quick_output( 'check_output/vx2.nc', 'vx2', vg2%z, vg2%z )
    call io_quick_output( 'check_output/vp1.nc', 'vp1', vg1%p, vg1%p )
    call io_quick_output( 'check_output/vp2.nc', 'vp2', vg2%p, vg2%p )
  end if

  call mympi_quick_output( 'check_output/tmask.nc', 'tmsk', &
    g1%msk, glo_lon, glo_lat, z )
  call mympi_quick_output( 'check_output/umask.nc', 'umsk', &
    g3%msk, glo_lon, glo_lat, z )
  call mympi_quick_output( 'check_output/itn.nc', 'itn', &
    g1%lev, glo_lon, glo_lat )
  call mympi_quick_output( 'check_output/iun.nc', 'iun', &
    g3%lev, glo_lon, glo_lat )

  call mympi_quick_output( 'check_output/phib.nc', 'phib', &
    g1%phib, glo_lon, glo_lat )
  call mympi_quick_output( 'check_output/pb.nc', 'pb', &
    g1%pb, glo_lon, glo_lat )
  call mympi_quick_output( 'check_output/pbu.nc', 'pbu', &
    g3%pb, glo_lon, glo_lat )

  call mympi_quick_output( 'check_output/phibx.nc', 'phibx', &
    graphib%x(1)%v, glo_lon, glo_lat )
  call mympi_quick_output( 'check_output/phiby.nc', 'phiby', &
    graphib%x(2)%v, glo_lon, glo_lat )

  call mympi_quick_output( 'check_output/h1g1.nc', 'h1g1', &
    hg1%rh1*a, glo_lon, glo_lat )
  call mympi_quick_output( 'check_output/h1g2.nc', 'h1g2', &
    hg2%rh1*a, glo_lon, glo_lat )
  call mympi_quick_output( 'check_output/h1g3.nc', 'h1g3', &
    hg3%rh1*a, glo_lon, glo_lat )
  call mympi_quick_output( 'check_output/h1g4.nc', 'h1g4', &
    hg4%rh1*a, glo_lon, glo_lat )

end subroutine debug_output

subroutine init_gvar_r3d_b (var, d1, d2, d3, ini, g)!{{{1
  ! initialize grid variables
  type (type_gvar_r3d) :: var
  integer, intent(in) :: d1, d2, d3
  real (kind=sglp) :: ini
  type (type_stg3d), target :: g

  integer :: is

  allocate(var%v(d1,d2,d3), stat = is)
  call chk(is)

  var%v = ini
  var%g => g

end subroutine init_gvar_r3d_b

subroutine init_gvar_m3d_b (var, d1, d2, d3, ini, g) !{{{1
  ! initialize grid variable
  type (type_gvar_m3d) :: var
  integer, intent(in) :: d1, d2, d3
  real (kind=sglp) :: ini
  type (type_stg3d), target :: g
  call init_gvar_r3d_b (var%x(1), d1, d2, d3, ini, g)
  call init_gvar_r3d_b (var%x(2), d1, d2, d3, ini, g)
end subroutine init_gvar_m3d_b

subroutine merge_out_r2d_rec (fname, varname, var, nrec) !{{{1
  ! merge 3d array from other domains to mid
  character (len=*), intent(in) :: fname, varname
  real (kind=wp), dimension(ni,nj), intent(in) :: var
  integer, intent(in) :: nrec

  real (kind=wp), allocatable, dimension(:,:) :: glo_var

  integer, parameter :: tag = 30
  type (type_my) :: d
  integer :: n, leng

  if (myid == mid) then

    allocate( glo_var(glo_ni, glo_nj), stat=is)
    call chk(is); glo_var = 0.0

    do n = 1, npro
      d = our(n)
      leng  = (d%ge-d%gw+1) * (d%gn-d%gs+1)
      if ( d%id == mid ) then
        glo_var(d%gw:d%ge, d%gs:d%gn) = &
          var(2:d%ni-1, 2:d%nj-1)
      else
        call mpi_recv (glo_var(d%gw:d%ge,d%gs:d%gn), &
          leng, mpi_real8, d%id, tag, mpi_comm_world, msta, err)
      end if
    end do

    call io_write (trim(fname), trim(varname), glo_var, nrec)
    deallocate(glo_var)

  else
    leng  = (my%ge-my%gw+1) * (my%gn-my%gs+1)
    call mpi_ssend (var(2:my%ni-1,2:my%nj-1), leng, &
      mpi_real8, mid, tag, mpi_comm_world, err)
  end if

end subroutine merge_out_r2d_rec

  type :: type_date
    character (len=19) :: d
  end type type_date

function time2date (time) !{{{1
  ! time to date
  type (type_time), intent(in)  :: time
  type (type_date) :: time2date

  write(time2date % d(1:4),   '(i4)') time % y
  write(time2date % d(6:7),   '(i2)') time % m
  write(time2date % d(9:10),  '(i2)') time % d
  write(time2date % d(12:13), '(i2)') time % h
  write(time2date % d(15:16), '(i2)') time % mi
  write(time2date % d(18:19), '(i2)') time % s

end function time2date
function date2time (date) !{{{1
  ! date to time, the string in date should be check if from user definition
  type (type_date), intent(in)  :: date
  type (type_time) :: date2time

  read(date % d(1:4),   '(i4)') date2time % y
  read(date % d(6:7),   '(i2)') date2time % m
  read(date % d(9:10),  '(i2)') date2time % d
  read(date % d(12:13), '(i2)') date2time % h
  read(date % d(15:16), '(i2)') date2time % mi
  read(date % d(18:19), '(i2)') date2time % s

  date2time % mm = pro_days_month (date2time % y, date2time % m)

end function date2time

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

!  interface arrays_cp_shape
!    module procedure cp_shape_gr2d
!    module procedure cp_shape_gr2d_gm2d
!    module procedure cp_shape_gr3d
!    module procedure cp_shape_gr3d_gm3d
!    module procedure cp_shape_gr3d_gm3d_b
!    module procedure cp_shape_gm2d
!    module procedure cp_shape_gm3d
!  end interface

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

  interface arrays_free
    module procedure free_gr2d
    module procedure free_gr3d
    module procedure free_gm2d
    module procedure free_gm3d
  end interface

subroutine io_create_grdvar (ncname) !{{{1
  ! for output grid variables, mainly for checking after 
  !   changing lots of code
  ! Note that this output file is mainly for code checking, 
  !   so the coordinates will not always corresspond to the 
  !   exact variables on which they underlying

  character (len=*), intent(in) :: ncname

  integer :: dimid1, dimid2, dimid3, i
  integer :: dimids(3)

  call check ( nf90_create (ncname, NF90_CLOBBER, ncid)  )

  !def dim. {{{2
  call check ( nf90_def_dim (ncid, 'lon', glo_ni, dimid1) )
  call check ( nf90_def_dim (ncid, 'lat', glo_nj, dimid2) )
  call check ( nf90_def_dim (ncid, 'z',   nk, dimid3) )

  !def global attr. {{{2
  call check ( nf90_put_att (ncid, NF90_GLOBAL, & 
    'created', "by subroutine io_create_grdvar in module mod_io") )

  ! def vars  !{{{2

  ! coordinates vars.
  call check ( nf90_def_var (ncid, "lon", nf90_float, &
    dimid1, varid) )
  call check ( nf90_put_att (ncid, varid, 'units', &
    'degree_east') )

  call check ( nf90_def_var (ncid, "lat", nf90_float, &
    dimid2, varid) )
  call check ( nf90_put_att (ncid, varid, 'units', &
    'degree_north') )

  call check ( nf90_def_var (ncid, "z", nf90_float, &
    dimid3, varid) )
  call check ( nf90_put_att (ncid, varid, 'units', &
    'm') )

  ! 3d vars.
  dimids = (/dimid1, dimid2, dimid3/)
  call check ( nf90_def_var (ncid, 'tmsk', nf90_int, &
    dimids, varid) )
  call check ( nf90_def_var (ncid, 'umsk', nf90_int, &
    dimids, varid) )

  ! 2d vars.
  call check ( nf90_def_var (ncid, 'itn', nf90_int, &
    (/dimid1, dimid2/), varid) )
  call check ( nf90_def_var (ncid, 'iun', nf90_int, &
    (/dimid1, dimid2/), varid) )
  call check ( nf90_def_var (ncid, 'phib', nf90_float, &
    (/dimid1, dimid2/), varid) )
  call check ( nf90_def_var (ncid, 'phibx', nf90_float, &
    (/dimid1, dimid2/), varid) )
  call check ( nf90_def_var (ncid, 'phiby', nf90_float, &
    (/dimid1, dimid2/), varid) )
  call check ( nf90_def_var (ncid, 'pbu', nf90_float, &
    (/dimid1, dimid2/), varid) )
  call check ( nf90_def_var (ncid, 'ph', nf90_float, &
    (/dimid1, dimid2/), varid) )
  call check ( nf90_def_var (ncid, 'h1', nf90_float, &
    (/dimid1, dimid2/), varid) )
  call check ( nf90_def_var (ncid, 'h2', nf90_float, &
    (/dimid1, dimid2/), varid) )
  call check ( nf90_def_var (ncid, 'h3', nf90_float, &
    (/dimid1, dimid2/), varid) )
  call check ( nf90_def_var (ncid, 'h4', nf90_float, &
    (/dimid1, dimid2/), varid) )

  ! 1d vars.
  call check ( nf90_def_var (ncid, 'vx1', nf90_float, &
    (/dimid3/), varid) )
  call check ( nf90_def_var (ncid, 'vx2', nf90_float, &
    (/dimid3/), varid) )
  call check ( nf90_def_var (ncid, 'vp1', nf90_float, &
    (/dimid3/), varid) )
  call check ( nf90_def_var (ncid, 'vp2', nf90_float, &
    (/dimid3/), varid) )

  !end def {{{2
  call check (nf90_enddef(ncid) )

  call check (nf90_close(ncid) )

end subroutine io_create_grdvar

