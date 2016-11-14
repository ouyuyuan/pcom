
! Description: 
!
!      Author: OU Yuyuan <ouyuyuan@lasg.iap.ac.cn>
!     Created: 2015-09-13 08:14:52 BJT
! Last Change: 2016-04-06 16:11:39 BJT

program main

  ! imported variables !{{{1
  !-------------------------------------------------------=-
  use mod_arrays, only: &
    adv, arrays_init, &
    cv1, cv2, &
    hg1, hg2, hg3, hg4, hgt, hgu, &
    sg1, sg2, sg3, sg4, &
    vg1, vg2, vgt, vgw, g1, g3, gt, gu, &
    acts, acuv, acw, acssh, am, arrays_allocate, &
    bnd, bphi, bgraphi, &
    dub, &
    frc, fcor, fri, &
    glo_lat, glo_lon, graphib, grapa, &
    hpos, &
    km, kh, &
    lon, lat, &
    rrho, rrhodp, &
    pbt, prho, &
    ts, &
    up, upb, &
    z, &
    arrays_cp_shape, arrays_free
    
  use mod_con, only: &
    rho0, g, am_c, km_c, gamma_b, a, torad, &
    omega

  use mod_den, only: den_rrho, den_rho, den_prho

  use mod_int, only: int_trop, int_bnd, &
    int_pgra, int_readyc, int_clin, int_ts, int_ssh

  use mod_io, only: io_create, io_create_grdvar, &
    io_get_dim_len, io_write, &
    io_quick_output, io_read

  use mod_kind, only: wp

  use mpi

  use mod_mympi, only: &
    mympi_div, mympi_divx, mympi_divy, &
    mympi_output, mympi_swpbnd, &
    mympi_merge, mympi_bcast, mympi_quick_output

  use mod_op, only: &
    op_lap, op_ter, op_gra

  use mod_param, only: &
    nm, my, missing_float, &
    glo_ni, glo_nj, &
    ni, nj, nk, nim, njm, nkp, &
    myid, npro, mid, names, ids, &
    param_set_my, tctr, &
    tc, tp

  use mod_pro, only: pro_print

  use mod_type, only: &
    type_accu_gm3d, type_accu_gr3d, type_accu_gr2d, &
    type_stg, type_vstg, type_stg3d, &
    type_str2time, &
    type_str2sec, type_mat, &
    type_memo_r, type_frc, &
    type_gvar_m2d, type_gvar_m3d, &
    type_gvar_r2d, type_gvar_r3d, &
    type_bintg3, &
    type_print, &
    operator (+), operator (<)

  ! local variables !{{{1
  implicit none

  integer, parameter :: fid_dia = 10

  integer :: err, leng, errorcode, i, j, k

  character (len=mpi_max_processor_name) :: hostname

  integer :: ia3(2,2,2)
  ! MPI environment !{{{1
  call mpi_init (err)
  if (err .ne. mpi_success) then 
    write(*,*) 'Error starting MPI program. Terminating. :-('
    call mpi_abort (mpi_comm_world, errorcode, err)
  end if

  call mpi_comm_rank (mpi_comm_world, myid, err)
  call mpi_comm_size (mpi_comm_world, npro, err)
  call mpi_get_processor_name (hostname, leng, err)
  write(*,'(a, i4, a, i4, a, a)') 'Number of tasks = ', npro, &
    ', My rank = ', myid, ', running on '//trim(hostname)//' :-)'

  ! preparation !{{{1
  if ( myid == mid ) open( fid_dia, file='diaginfo.txt' )

  ! initialize program environment
  call init ()

  ! setting viscosity and diffusity coefficients
  call set_cf( km%v, kh%v )

!   calc. variables depend on grid position
  call calc_hpos (hpos)

  ! setting stagger grids
  call set_hgrid (hg1, hg2, hg3, hg4, hgt, hgu) 
  call set_vgrid (vg1, vg2, vgt, vgw)
  call set_grid  (g1, g3, gt, gu)

  ! calc. grid variables
  call calc_cf( cv1, cv2 )

  ! Coriolis force
  fcor%v = 2.0*omega*sin(hpos(:,:,fcor%hg%n)%x(2))

  ! prepare initial state of the ocean
  call inistat(ts, frc, graphib)

  ! create output file
  if ( myid == mid ) then
    call io_create (nm%fo)
    write(fid_dia, '(a)') 'created output file '//nm%fo

    call io_write (trim(nm%fo), 'lon', glo_lon)
    call io_write (trim(nm%fo), 'lat', glo_lat)
    call io_write (trim(nm%fo), 'z',   z)
  end if

  ! integration cycle !{{{1

  ! time prepare !{{{2
  tctr%ct = type_str2time (nm%bd)
  tctr%pt = tctr%ct
  tctr%nt = ( type_str2sec(nm%ed) - type_str2sec(nm%bd) ) / nm%bc

  ! baroclinic integrate !{{{2

  tctr%t1 = mpi_wtime () ! time count start

  do i = 1, tctr%nt

    tctr%i = i
    ! linear interpolate monthly forcing to daily values of bnd
    !   at the initial of integration, or when the 
    !   integration enters a new day
    if ( (tctr%i == 1) .or. (tctr%ct%d /= tctr%pt%d) ) then 
      call int_bnd (bnd)
      call den_prho (prho)
    end if

    ! calc. the specific volume, reciprocal of density
    call den_rrho (rrho%v, ts(tc), pbt%tc, gt%msk)

    ! calc. pressure gradient forces
    call int_pgra (rrhodp, bphi, bgraphi)

    ! prepare for baroclinic integration
    call int_readyc (grapa, dub, adv, fri)

    ! integrate series time steps per baroclinic time step
    call int_trop (pbt, dub)

    ! prediction of baroclinic mode
    call int_clin (up, acuv)

    ! calc. sea surface height
    call int_ssh (acssh, rrho%v)

    ! prediction of temperature and salinity
    call int_ts (ts, acts, acw)

    tctr%pt = tctr%ct
    tctr%ct = tctr%ct + nm%bc

    call check_output (acts, acuv, acw, acssh)
  end do

  ! finish integration !{{{1
  if (myid == mid) then
    write (*, *) ''
    write (*, '(a)') 'Done model run.'
    write (*, '(a, i9, a)') &
      'Integrated for ', tctr%nt, ' baroclinic time steps.'
  end if

  call mpi_finalize (err)

contains  !{{{1

subroutine init () !{{{1
  ! initialize model environment
  use mod_param, only: param_set_nm

  integer :: is

  ! get namelist !{{{2

  ! although we can read in namelist in all processors, 
  ! but it seems the namelist file may be destroyed by 
  ! muliple-opening (happened once)
  if ( myid == mid ) call param_set_nm ()

  call mympi_bcast (nm%px)
  call mympi_bcast (nm%py)

  call mympi_bcast (nm%bd)
  call mympi_bcast (nm%ed)

  call mympi_bcast (nm%bt)
  call mympi_bcast (nm%bc)

  call mympi_bcast (nm%fi)
  call mympi_bcast (nm%ff)
  call mympi_bcast (nm%fo)

  call mympi_bcast (nm%per)

  ! determine dimensions !{{{2

  if ( myid == mid ) then
    call io_get_dim_len ( nm%fi, 'lon', glo_ni )
    call io_get_dim_len ( nm%fi, 'lat', glo_nj )
    call io_get_dim_len ( nm%fi, 'z',   nk )
  end if
  call mympi_bcast(glo_nj)
  call mympi_bcast(glo_ni)
  call mympi_bcast(nk)

  call param_set_my ()
  call mpi_barrier (mpi_comm_world, err)

  ni = my%ni
  nj = my%nj

  ! too many cpus are not allowed
  if ( ni <=2 ) stop 'cpu number in x direction exceeds grid points.'
  if ( nj <=2 ) stop 'cpu number in y direction exceeds grid points.'

  nim = ni - 1
  njm = nj - 1
  nkp = nk + 1

  if ( myid == mid ) call write_dim_info (fid_dia)

  ! allocate arrays !{{{2
  call arrays_allocate ()
  call mpi_barrier ( mpi_comm_world, err )

  ! get coordinate variables !{{{2

  if ( myid == mid ) then
    ! (glo_lon, glo_lat) is on grid 4
    call io_read (nm%fi, 'lon', glo_lon)
    call io_read (nm%fi, 'lat', glo_lat)
    call io_read (nm%fi,   'z', z)
  end if

  call mympi_bcast (z)

  ! (lon, lat) in on grid 4
  call mympi_divy (glo_lat, lat)
  if ( my%gs == 1 )     lat(1) = 2*lat(2) - lat(3)
  if ( my%gn == glo_nj) lat(nj) = 2*lat(njm) - lat(nj-2)

  call mympi_divx (glo_lon, lon)

end subroutine init

subroutine set_cf (km, kh) !{{{2
  ! setting the profile of vertical momentum viscosity coefficient
  ! reading the datasets of kh (vertical turbulent mixing coefficients), 
  !  \ref{Zhang2014}
  real (kind=wp), dimension(ni,nj,nkp) :: km, kh
  real (kind=wp) :: pf(12)

  real (kind=wp), allocatable, dimension(:,:,:) :: glo_kh
  integer :: k, is

  if (myid == mid) then
    allocate(glo_kh(glo_ni, glo_nj, nkp), stat=is); call chk(is)
    glo_kh = 0.0

    call io_read ( 'input/vmix.nc', 'vmix',  glo_kh(:,:,2:nkp) )
    where ( glo_kh > 9.0e+30 - 1.0 ) glo_kh = 0.1
    where ( glo_kh > 200.0 ) glo_kh = 200.0
    glo_kh = glo_kh * 1.0e-4 ! cm^2 / s to m^2/s
  end if

  km = km_c

  pf = (/ 0.1, 50.0, 50.0, 50.0, 40.0, &
         20.0, 10.0,  5.0,  2.0,  1.0, &
         0.5, 0.2/)
  pf = pf*1.0e-4 ! cm^2/s => m^2/s

  do k = 1, 12
    km(:,:,k) = pf(k)
  end do

  call mympi_div (glo_kh, kh)

end subroutine set_cf

subroutine calc_hpos (hpos) !{{{1
  ! calc. horizontal position on grid 1/2/3/4
  ! (lon, lat) in I/O file is on grid 4

  type (type_mat), dimension(ni,nj,4) :: hpos

  integer :: i, j
  type (type_mat) :: p
  real (kind=wp) :: loni, latj

  ! inner area !{{{2
  do j = 2, njm
  do i = 2, nim
    !westest/eastest longitude wrap up
    loni = lon(i)
    latj = (lat(j) + lat(j-1)) * 0.5
    p%x(1) = loni
    p%x(2) = latj
    hpos(i,j,1)%x(1) = p%x(1) * torad
    hpos(i,j,1)%x(2) = p%x(2) * torad

    loni = (lon(i) + lon(i+1)) * 0.5
    ! the 'eastest' point will be 360.0 degree, for keeping
    !   monotonous, we donot change it to 0 degree
    if ( lon(i+1) < lon(i) ) loni = (lon(i) + lon(i+1)+360) * 0.5
    latj = (lat(j) + lat(j-1)) * 0.5
    p%x(1) = loni
    p%x(2) = latj
    hpos(i,j,2)%x(1) = p%x(1) * torad
    hpos(i,j,2)%x(2) = p%x(2) * torad

    loni = (lon(i) + lon(i+1)) * 0.5
    if ( lon(i+1) < lon(i) ) loni = (lon(i) + lon(i+1)+360) * 0.5
    latj = lat(j)
    p%x(1) = loni
    p%x(2) = latj
    hpos(i,j,3)%x(1) = p%x(1) * torad
    hpos(i,j,3)%x(2) = p%x(2) * torad

    loni = lon(i)
    latj = lat(j)
    p%x(1) = loni
    p%x(2) = latj
    hpos(i,j,4)%x(1) = p%x(1) * torad
    hpos(i,j,4)%x(2) = p%x(2) * torad
  end do
  end do
  call mympi_swpbnd (hpos)

  ! extropolate  !{{{2

  if (my%gs == 1) then
    hpos(:,1,:)%x(1) = 2*hpos(:,2,:)%x(1) - hpos(:,3,:)%x(1)
    hpos(:,1,:)%x(2) = 2*hpos(:,2,:)%x(2) - hpos(:,3,:)%x(2)
  end if

  if (my%gn == glo_nj) then
    hpos(:,nj,:)%x(1) = 2*hpos(:,nj-1,:)%x(1) - hpos(:,nj-2,:)%x(1)
    hpos(:,nj,:)%x(2) = 2*hpos(:,nj-1,:)%x(2) - hpos(:,nj-2,:)%x(2)
  end if

end subroutine calc_hpos

subroutine set_hgrid (hg1, hg2, hg3, hg4, hgt, hgu) !{{{1
  ! setting horizontal stagger grid properties
  type (type_stg), target :: hg1, hg2, hg3, hg4
  type (type_stg), pointer :: hgt, hgu

  real (kind=wp) :: dx1, dx2
  integer :: is, i, j
  
  do j = 2, njm
  do i = 2, nim
    !westest/eastest longitude wrap up
    dx1 = hpos(i,j,2)%x(1) - hpos(i-1,j,2)%x(1)
    if ( dx1 < 0 ) dx1 = dx1 + 360.0*torad
    dx2 = hpos(i,j,4)%x(2) - hpos(i,j-1,4)%x(2)
    hg1%dx(i,j)%x(1) = dx1
    hg1%dx(i,j)%x(2) = dx2

    dx1 = hpos(i+1,j,1)%x(1) - hpos(i,j,1)%x(1)
    if ( dx1 < 0 ) dx1 = dx1 + 360.0*torad
    dx2 = hpos(i,j,3)%x(2) - hpos(i,j-1,3)%x(2)
    hg2%dx(i,j)%x(1) = dx1
    hg2%dx(i,j)%x(2) = dx2

    dx1 = hpos(i+1,j,4)%x(1) - hpos(i,j,4)%x(1)
    if ( dx1 < 0 ) dx1 = dx1 + 360.0*torad
    dx2 = hpos(i,j+1,2)%x(2) - hpos(i,j,2)%x(2)
    hg3%dx(i,j)%x(1) = dx1
    hg3%dx(i,j)%x(2) = dx2

    dx1 = hpos(i,j,3)%x(1) - hpos(i-1,j,3)%x(1)
    if ( dx1 < 0 ) dx1 = dx1 + 360.0*torad
    dx2 = hpos(i,j+1,1)%x(2) - hpos(i,j,1)%x(2)
    hg4%dx(i,j)%x(1) = dx1
    hg4%dx(i,j)%x(2) = dx2
  end do
  end do
  call mympi_swpbnd (hg1%dx)
  call mympi_swpbnd (hg2%dx)
  call mympi_swpbnd (hg3%dx)
  call mympi_swpbnd (hg4%dx)

  ! extropolate

  if (my%gs == 1) then
    hg1%dx(:,1)%x(1) = 2*hg1%dx(:,2)%x(1) - hg1%dx(:,3)%x(1)
    hg1%dx(:,1)%x(2) = 2*hg1%dx(:,2)%x(2) - hg1%dx(:,3)%x(2)

    hg2%dx(:,1)%x(1) = 2*hg2%dx(:,2)%x(1) - hg2%dx(:,3)%x(1)
    hg2%dx(:,1)%x(2) = 2*hg2%dx(:,2)%x(2) - hg2%dx(:,3)%x(2)

    hg3%dx(:,1)%x(1) = 2*hg3%dx(:,2)%x(1) - hg3%dx(:,3)%x(1)
    hg3%dx(:,1)%x(2) = 2*hg3%dx(:,2)%x(2) - hg3%dx(:,3)%x(2)

    hg4%dx(:,1)%x(1) = 2*hg4%dx(:,2)%x(1) - hg4%dx(:,3)%x(1)
    hg4%dx(:,1)%x(2) = 2*hg4%dx(:,2)%x(2) - hg4%dx(:,3)%x(2)
  end if

  if (my%gn == glo_nj) then
    hg1%dx(:,nj)%x(1) = 2*hg1%dx(:,nj-1)%x(1) - hg1%dx(:,nj-2)%x(1)
    hg1%dx(:,nj)%x(2) = 2*hg1%dx(:,nj-1)%x(2) - hg1%dx(:,nj-2)%x(2)

    hg2%dx(:,nj)%x(1) = 2*hg2%dx(:,nj-1)%x(1) - hg2%dx(:,nj-2)%x(1)
    hg2%dx(:,nj)%x(2) = 2*hg2%dx(:,nj-1)%x(2) - hg2%dx(:,nj-2)%x(2)

    hg3%dx(:,nj)%x(1) = 2*hg3%dx(:,nj-1)%x(1) - hg3%dx(:,nj-2)%x(1)
    hg3%dx(:,nj)%x(2) = 2*hg3%dx(:,nj-1)%x(2) - hg3%dx(:,nj-2)%x(2)

    hg4%dx(:,nj)%x(1) = 2*hg4%dx(:,nj-1)%x(1) - hg4%dx(:,nj-2)%x(1)
    hg4%dx(:,nj)%x(2) = 2*hg4%dx(:,nj-1)%x(2) - hg4%dx(:,nj-2)%x(2)
  end if

  ! Lame coefficients and triangle function

  do j = 1, nj
  do i = 1, ni
    hg1%rh1(i,j) = cos( hpos(i,j,1)%x(2) )
    hg2%rh1(i,j) = cos( hpos(i,j,2)%x(2) )
    hg3%rh1(i,j) = cos( hpos(i,j,3)%x(2) )
    hg4%rh1(i,j) = cos( hpos(i,j,4)%x(2) )

    hg1%tn(i,j) = tan( hpos(i,j,1)%x(2) )
    hg2%tn(i,j) = tan( hpos(i,j,2)%x(2) )
    hg3%tn(i,j) = tan( hpos(i,j,3)%x(2) )
    hg4%tn(i,j) = tan( hpos(i,j,4)%x(2) )
  end do
  end do

  hgt => hg1
  hgu => hg3
end subroutine set_hgrid

subroutine set_vgrid (vg1, vg2, vgt, vgw) !{{{1
  ! set vertical stagger grid properties
  type (type_vstg), target :: vg1, vg2
  type (type_vstg), pointer :: vgt, vgw

  real (kind=wp), dimension(nk) :: dz, rmn, pmn, tmn, smn
  real (kind=wp) :: pmnbot
  integer :: i, j, k, is

  ! calc. pressure and height !{{{2
  if ( myid == mid ) then
    call io_read (nm%fi, 'tmn', tmn)
    call io_read (nm%fi, 'smn', smn)
  end if
  call mympi_bcast (tmn)
  call mympi_bcast (smn)

  ! dz, solve the eq. ( dz(k-1) + dz(k) ) / 2 = z(k) - z(k-1)
  dz(1) = 2*z(1)
  do k = 2, nk
    dz(k) = 2 * ( z(k) - z(k-1) ) - dz(k-1) 
  end do

  pmnbot = 1.0
  pmn(:) = 0.0
  rmn(:) = rho0

  do while ( abs(pmn(nk) - pmnbot) .gt. 1e-6 )
    pmnbot = pmn(nk)

    rmn(1) = den_rho( tmn(1), smn(1), pmn(1)*0.5 )
    do k = 2, nk
      rmn(k) = den_rho( tmn(k), smn(k), (pmn(k-1) + pmn(k))*0.5 )
    end do

    pmn(1) = rmn(1) * g * dz(1)
    do k = 2, nk
      pmn(k) = rmn(k) * g * dz(k) + pmn(k-1)
    end do
  end do

  ! set vertical grid info. at the initial state !{{{2
  ! z
  vg1%z(1) = 0.0 ! sea surface
  do k = 2, nkp
    vg1%z(k) = vg1%z(k-1) + dz(k-1)
  end do

  vg2%z(:) = z(:)

  ! p
  vg1%p(1) = 0.0
  vg1%p(2:nkp) = pmn(:)

  vg2%p(1) = pmn(1)*0.5
  do k = 2, nk
    vg2%p(k) = ( pmn(k) + pmn(k-1) )*0.5
  end do

  ! dz
  vg1%dz = 0.0
  vg1%dz(1) = 0.5*dz(1)
  do k = 2, nk
    vg1%dz(k) = z(k) - z(k-1)
  end do
  vg1%dz(nkp) = vg1%dz(nk)

  vg2%dz(:) = dz(:)

  ! dp
  vg1%dp(1) = 0.5*pmn(1)
  do k = 2, nk
    vg1%dp(k) = vg2%p(k) - vg2%p(k-1)
  end do
  vg1%dp(nkp) = vg1%dp(nk)

  vg2%dp(1) = pmn(1)
  do k = 2, nk
    vg2%dp(k) = pmn(k) - pmn(k-1)
  end do

  vgw => vg1
  vgt => vg2

end subroutine set_vgrid

subroutine set_grid (g1, g3, gt, gu) !{{{1
  ! 3d model stagger grid
  type (type_stg3d),target :: g1, g3
  type (type_stg3d), pointer :: gt, gu

  integer, allocatable, dimension(:,:) :: glo_itn
  real (kind=wp), dimension(nk) :: phimn
  integer :: i, j, k, is

  if ( myid == mid ) then
    allocate(glo_itn(glo_ni, glo_nj), stat=is); call chk(is)
    glo_itn = 0
    call io_read (nm%fi, 'itn', glo_itn)
  end if
  call mympi_div (glo_itn, g1%lev)

  ! set to land first
  g1%msk(:,:,:) = 0 
  g3%msk(:,:,:) = 0 

  do k = 1, nk
  do j = 2, njm
  do i = 2, nim
    if ( g1%lev(i,j) >= k ) g1%msk(i,j,k) = 1
  end do
  end do
  end do
  call mympi_swpbnd ( g1%msk )

  ! a water U-grid should be surounded by T-grids
  do k = 1, nk
  do j = 1, njm
  do i = 1, nim
    g3%msk(i,j,k) = g1%msk(i,j,k) * g1%msk(i+1,j,k)  &
                  * g1%msk(i,j+1,k) * g1%msk(i+1,j+1,k)
  end do
  end do
  end do
  call mympi_swpbnd ( g3%msk )

  g3%lev(:,:) = sum(g3%msk(:,:,:), 3)

  ! sea bottom pressure, geoptential height
  phimn(1) = - g * vg2%dz(1)
  do k = 2, nk
    phimn(k) = phimn(k-1)  - g * vg2%dz(k)
  end do

  do j = 1, nj
  do i = 1, ni
    k = g1%lev(i, j)
    if ( k == 0 ) then
      g1%phib(i,j) = 0.0
    else
      g1%phib(i,j) = phimn(k)
      g1%pb(i,j) = vg1%p(k+1)
    end if

    k = g3%lev(i, j)
    if ( k == 0 ) then
      g3%phib(i,j) = 0.0
    else
      g3%phib(i,j) = phimn(k)
      g3%pb(i,j) = vg1%p(k+1)
    end if
  end do
  end do

  ! in case of 1/pb, set to 0.1Pa in land
  where( g1%lev == 0 ) g1%pb = 0.1
  where( g3%lev == 0 ) g3%pb = 0.1

  gt=>g1
  gu=>g3
end subroutine set_grid

subroutine calc_cf (cv1, cv2) !{{{1
  ! calc. coefficents of frictional force at hg3 grid
  real (kind=wp), dimension(ni,nj) :: cv1, cv2
  real (kind=wp), dimension(ni,nj) :: temp

  ! pcom 1.0 use this inderct method instead of directly use tan function
  ! these two ways differ a little bit
!  temp = tan( hpos(:,:,3)%x(2) )
  temp = sqrt(1-cos(hpos(:,:,3)%x(2))**2) / cos(hpos(:,:,3)%x(2))
  cv1 = (1 - temp*temp) / a**2
  cv2 = 2*temp / a

end subroutine calc_cf

subroutine inistat (ts, frc, graphib) !{{{1
  ! prepare the initial state of the ocean
  type (type_gvar_m3d) :: ts(2)
  type (type_frc) :: frc
  type (type_gvar_m2d) :: graphib, wk

  type (type_stg), pointer :: hg1, hg2
  real (kind=wp), allocatable, dimension(:,:,:) :: glo_pt, glo_sa
  type (type_frc) :: glo_frc
  integer :: is, m, i, j, k

  if ( myid == mid ) then
    ! 12 months forcing
    call arrays_init(glo_frc%tau%x(1), glo_ni,glo_nj,12, 0.0, frc%tau%x(1)%g)
    call arrays_init(glo_frc%tau%x(2), glo_ni,glo_nj,12, 0.0, frc%tau%x(2)%g)
    call arrays_init(glo_frc%ts%x(1),  glo_ni,glo_nj,12, 0.0, frc%ts%x(1)%g)
    call arrays_init(glo_frc%ts%x(2),  glo_ni,glo_nj,12, 0.0, frc%ts%x(2)%g)
    call arrays_init(glo_frc%pa,       glo_ni,glo_nj,12, 0.0, frc%pa%g)
    call arrays_init(glo_frc%fw,       glo_ni,glo_nj,12, 0.0, frc%fw%g)

    allocate(glo_pt(glo_ni, glo_nj, nk), stat=is); call chk(is)
    glo_pt = 0.0
    allocate(glo_sa(glo_ni, glo_nj, nk), stat=is); call chk(is)
    glo_sa = 0.0

    call io_read (nm%fi, 'pt',  glo_pt)
    call io_read (nm%fi, 'sa',  glo_sa)

    call io_read (nm%ff, 'taux', glo_frc%tau%x(1)%v) 
    call io_read (nm%ff, 'tauy', glo_frc%tau%x(2)%v) 

    call io_read (nm%ff, 'bct', glo_frc%ts%x(1)%v) 
    call io_read (nm%ff, 'bcs', glo_frc%ts%x(2)%v) 

    call io_read (nm%ff, 'pa', glo_frc%pa%v) 

    call io_read (nm%ff, 'fw', glo_frc%fw%v) 
  end if

  call mympi_div (glo_frc%tau, frc%tau)
  call mympi_div (glo_frc%ts, frc%ts)
  call mympi_div (glo_frc%pa, frc%pa)
  call mympi_div (glo_frc%fw, frc%fw)

  ! forcing on land set to zero (no forcing)
  !!! potential bug, missing is not all at land of input file
!  where (spread(g1%lev, 3, 12) == 0)
  where (frc%tau%x(1)%v == missing_float)
    frc%tau%x(1)%v = 0.0 
    frc%tau%x(2)%v = 0.0 
    frc%ts%x(1)%v = 0.0 
    frc%ts%x(2)%v = 0.0 
    frc%pa%v = 0.0 
    frc%fw%v = 0.0 
  end where

  ! initial (T, S)
  call mympi_div (glo_pt, ts(tc)%x(1))
  call mympi_div (glo_sa, ts(tc)%x(2))
  ts(tp) = ts(tc)

  !! why not average before differentiate?
  call arrays_cp_shape( graphib, wk )
  call op_gra( wk, g1%phib, g1%hg, g1%hg%ew, g1%hg%ns )
  call op_ter( graphib%x(1)%v, wk%x(1)%v, wk%x(1)%hg, graphib%x(1)%hg )
  call op_ter( graphib%x(2)%v, wk%x(2)%v, wk%x(2)%hg, graphib%x(2)%hg )

  call arrays_free( wk )
end subroutine inistat

subroutine write_dim_info (fid) !{{{1

  integer, intent(in) :: fid

    write( fid, '(a, i4, a, i2, a, i3, a)' )   &
      'glo_ni = ', glo_ni, ' with ', nm%px,          &
      ' processors in the x-direction, and ni approximate ', &
      ni-2, ' for the sub-domain.'

    write( fid, '(a, i4, a, i2, a, i3, a)' )   &
      'glo_nj = ', glo_nj, ' with ', nm%py,          &
      ' processors in the y-direction, and nj approximate ', &
      nj-2, ' for the sub-domain.'

    write( fid, '(a, i4)' ) 'nk = ', nk

    write( fid, * ) ''

end subroutine write_dim_info

subroutine check_output (acts, acuv, acw, acssh) !{{{1
  type (type_accu_gm3d) :: acts, acuv
  type (type_accu_gr3d) :: acw
  type (type_accu_gr2d) :: acssh

  ! the integration enters a new day, print time infos.
  if ( tctr%ct%d /= tctr%pt%d ) then
    tctr%t2 = mpi_wtime ()

    if ( myid==mid ) &
    write(*, '(a, i0.4,a,i0.2,a,i0.2, a, f5.2, a, i3)') &
      'integrate for ', tctr%pt%y, '-', tctr%pt%m, &
      '-', tctr%pt%d, ' use ', tctr%t2 - tctr%t1, &
      ' seconds on processor ', myid

    tctr%t1 = tctr%t2
  end if

  ! output if time is proper
  if ( ((tctr%ct%d/=tctr%pt%d).and.(nm%per.eq.'day'))   .or. &
       ((tctr%ct%m/=tctr%pt%m).and.(nm%per.eq.'month')) .or. &
       ((tctr%ct%y/=tctr%pt%y).and.(nm%per.eq.'year')) ) then

    call mympi_output (nm%fo, names%pt, names%sa, acts)

    call mympi_output (nm%fo, names%u, names%v, acuv, gu%msk)

    call mympi_output (nm%fo, names%w, acw, gu%msk)

    ! calc. fluctuation of ssh, minus global mean
    call mympi_output (nm%fo, names%ssh, acssh, &
      acssh%var%hg%rh1*gt%msk(:,:,1), gt%msk(:,:,1))

  end if

end subroutine check_output

subroutine chk( ista ) !{{{1
  ! check state of allocate array 

  integer, intent(in) ::  ista

  if ( ista /= 0 ) then
    write(*,*) 'Allocate array failed. Stop'
    stop 2
  end if
end subroutine chk

end program main !{{{1
!-------------------------------------------------------{{{1
! vim:fdm=marker:fdl=0:
! vim:foldtext=getline(v\:foldstart).'...'.(v\:foldend-v\:foldstart):
