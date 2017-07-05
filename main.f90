
! Description: 
!
!      Author: OU Yuyuan <ouyuyuan@lasg.iap.ac.cn>
!     Created: 2015-09-13 08:14:52 BJT
! Last Change: 2017-07-05 21:43:17 BJT

program main

  ! imported variables !{{{1
  !-------------------------------------------------------=-
  use mod_arrays, only: &
    adv, arrays_init, &
    cv1, cv2, &
    g1j, g2j, g3j, g4j, gtj, guj, &
    gi1, gi2, git, giw, g12, g32, gt, gu, &
    acts, acuv, acw, acssh, am, arrays_allocate, &
    bnd, bphi, bgraphi, &
    badv, &
    frc, cor, fri, &
    glo_lat, glo_lon, graphih, grapa, &
    hpos, &
    km, kh, &
    lon, lat, &
    alpha, adp, &
    ch, prho, &
    ts, &
    up, upb, &
    wm, &
    z
    
  use mod_con, only: &
    rho0, g, km_c, a, torad, omega

  use mod_den, only: den_alpha, den_rho, den_prho

  use mod_int, only: int_trop, int_bnd, &
    int_pgra, int_readyc, int_clin, int_ts, int_ssh

  use mod_io, only: io_create, &
    io_get_dim_len, io_write, io_read, &
    io_quick_output

  use mod_kind, only: wp, zero, lint

  use mpi

  use mod_mympi, only: &
    mympi_div, mympi_divx, mympi_divy, &
    mympi_output, mympi_swpbnd, mympi_bcast, &
    mympi_quick_output

  use mod_op, only: op_ter, op_gra

  use mod_param, only: &
    nm, my, missing_float, &
    glo_ni, glo_nj, &
    ni, nj, nk, nim, njm, nkp, &
    myid, npro, mid, names, &
    param_set_my, param_set_nm, &
    print_my, &
    tc, tp

  use mod_type, only: &
    type_accu_gm3d, type_accu_gr3d, type_accu_gr2d, &
    type_gi, type_gj, type_gij, &
    type_str2time, &
    type_str2sec, type_mat, &
    type_frc, tctr, type_tctr, &
    type_gvar_m2d, type_gvar_m3d, &
    type_gvar_r2d, type_gvar_r3d, &
    type_bintgu, &
    type_check_date, &
    operator (+), operator (<)

  ! local variables !{{{1
  implicit none

  integer, parameter :: fid_dia = 10

  integer :: err, leng, errorcode

  integer (kind=lint) :: i

  character (len=mpi_max_processor_name) :: hostname

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
  call set_hgrid (g1j, g2j, g3j, g4j, gtj, guj) 
  call set_vgrid (gi1, gi2, git, giw)
  call set_grid  (g12, g32, gt, gu)

  ! calc. grid variables
  call calc_cf( cv1, cv2 )

  ! Coriolis force
  cor%v = 2.0*omega*sin(hpos(:,:,cor%hg%n)%x(2))

  ! prepare initial state of the ocean
  call inistat(ts, frc, graphih)

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
    call den_alpha (alpha%v, ts(tc), ch%tc, gt%msk)

    ! calc. pressure gradient forces
    call int_pgra (adp, bphi, bgraphi)

    ! prepare for baroclinic integration
    call int_readyc (grapa, badv, adv, fri, wm)
!    call mympi_quick_output('output/test.nc', 'var', &
!      badv%x(1)%v, glo_lon, glo_lat) !DEBUG

    ! integrate series time steps per baroclinic time step
    call int_trop (ch, upb)

    ! prediction of baroclinic mode
    call int_clin (up, acuv, am%v)

    ! calc. sea surface height
    call int_ssh (acssh, alpha%v)

    ! prediction of temperature and salinity
    call int_ts (ts, acts, acw, wm)

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

  ! get namelist !{{{2

  ! although we can read in namelist in all processors, 
  ! but it seems the namelist file may be destroyed by 
  ! muliple-opening (happened once)
  if ( myid == mid ) then
    call param_set_nm (nm)
    call type_check_date(nm%bd)
    call type_check_date(nm%ed)
  end if

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

  call param_set_my (my)
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

subroutine set_hgrid (g1j, g2j, g3j, g4j, gtj, guj) !{{{1
  ! setting horizontal stagger grid properties
  type (type_gi), target :: g1j, g2j, g3j, g4j
  type (type_gi), pointer :: gtj, guj

  real (kind=wp) :: dx1, dx2
  integer :: i, j
  
  do j = 2, njm
  do i = 2, nim
    !westest/eastest longitude wrap up
    dx1 = hpos(i,j,2)%x(1) - hpos(i-1,j,2)%x(1)
    if ( dx1 < 0 ) dx1 = dx1 + 360.0*torad
    dx2 = hpos(i,j,4)%x(2) - hpos(i,j-1,4)%x(2)
    g1j%dx(i,j)%x(1) = dx1
    g1j%dx(i,j)%x(2) = dx2

    dx1 = hpos(i+1,j,1)%x(1) - hpos(i,j,1)%x(1)
    if ( dx1 < 0 ) dx1 = dx1 + 360.0*torad
    dx2 = hpos(i,j,3)%x(2) - hpos(i,j-1,3)%x(2)
    g2j%dx(i,j)%x(1) = dx1
    g2j%dx(i,j)%x(2) = dx2

    dx1 = hpos(i+1,j,4)%x(1) - hpos(i,j,4)%x(1)
    if ( dx1 < 0 ) dx1 = dx1 + 360.0*torad
    dx2 = hpos(i,j+1,2)%x(2) - hpos(i,j,2)%x(2)
    g3j%dx(i,j)%x(1) = dx1
    g3j%dx(i,j)%x(2) = dx2

    dx1 = hpos(i,j,3)%x(1) - hpos(i-1,j,3)%x(1)
    if ( dx1 < 0 ) dx1 = dx1 + 360.0*torad
    dx2 = hpos(i,j+1,1)%x(2) - hpos(i,j,1)%x(2)
    g4j%dx(i,j)%x(1) = dx1
    g4j%dx(i,j)%x(2) = dx2
  end do
  end do
  call mympi_swpbnd (g1j%dx)
  call mympi_swpbnd (g2j%dx)
  call mympi_swpbnd (g3j%dx)
  call mympi_swpbnd (g4j%dx)

  ! extropolate

  if (my%gs == 1) then
    g1j%dx(:,1)%x(1) = 2*g1j%dx(:,2)%x(1) - g1j%dx(:,3)%x(1)
    g1j%dx(:,1)%x(2) = 2*g1j%dx(:,2)%x(2) - g1j%dx(:,3)%x(2)

    g2j%dx(:,1)%x(1) = 2*g2j%dx(:,2)%x(1) - g2j%dx(:,3)%x(1)
    g2j%dx(:,1)%x(2) = 2*g2j%dx(:,2)%x(2) - g2j%dx(:,3)%x(2)

    g3j%dx(:,1)%x(1) = 2*g3j%dx(:,2)%x(1) - g3j%dx(:,3)%x(1)
    g3j%dx(:,1)%x(2) = 2*g3j%dx(:,2)%x(2) - g3j%dx(:,3)%x(2)

    g4j%dx(:,1)%x(1) = 2*g4j%dx(:,2)%x(1) - g4j%dx(:,3)%x(1)
    g4j%dx(:,1)%x(2) = 2*g4j%dx(:,2)%x(2) - g4j%dx(:,3)%x(2)
  end if

  if (my%gn == glo_nj) then
    g1j%dx(:,nj)%x(1) = 2*g1j%dx(:,nj-1)%x(1) - g1j%dx(:,nj-2)%x(1)
    g1j%dx(:,nj)%x(2) = 2*g1j%dx(:,nj-1)%x(2) - g1j%dx(:,nj-2)%x(2)

    g2j%dx(:,nj)%x(1) = 2*g2j%dx(:,nj-1)%x(1) - g2j%dx(:,nj-2)%x(1)
    g2j%dx(:,nj)%x(2) = 2*g2j%dx(:,nj-1)%x(2) - g2j%dx(:,nj-2)%x(2)

    g3j%dx(:,nj)%x(1) = 2*g3j%dx(:,nj-1)%x(1) - g3j%dx(:,nj-2)%x(1)
    g3j%dx(:,nj)%x(2) = 2*g3j%dx(:,nj-1)%x(2) - g3j%dx(:,nj-2)%x(2)

    g4j%dx(:,nj)%x(1) = 2*g4j%dx(:,nj-1)%x(1) - g4j%dx(:,nj-2)%x(1)
    g4j%dx(:,nj)%x(2) = 2*g4j%dx(:,nj-1)%x(2) - g4j%dx(:,nj-2)%x(2)
  end if

  ! Lame coefficients and triangle function

  do j = 1, nj
  do i = 1, ni
    g1j%rh(i,j) = cos( hpos(i,j,1)%x(2) )
    g2j%rh(i,j) = cos( hpos(i,j,2)%x(2) )
    g3j%rh(i,j) = cos( hpos(i,j,3)%x(2) )
    g4j%rh(i,j) = cos( hpos(i,j,4)%x(2) )

    g1j%tn(i,j) = tan( hpos(i,j,1)%x(2) )
    g2j%tn(i,j) = tan( hpos(i,j,2)%x(2) )
    g3j%tn(i,j) = tan( hpos(i,j,3)%x(2) )
    g4j%tn(i,j) = tan( hpos(i,j,4)%x(2) )
  end do
  end do

  gtj => g1j
  guj => g3j
end subroutine set_hgrid

subroutine set_vgrid (gi1, gi2, git, giw) !{{{1
  ! set vertical stagger grid properties
  type (type_gj), target :: gi1, gi2
  type (type_gj), pointer :: git, giw

  real (kind=wp), dimension(nk) :: dz, rmn, pmn, tmn, smn
  real (kind=wp) :: pmnbot
  integer :: k

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
  gi1%z(1) = 0.0 ! sea surface
  do k = 2, nkp
    gi1%z(k) = gi1%z(k-1) + dz(k-1)
  end do

  gi2%z(:) = z(:)

  ! pr
  gi1%pr(1) = 0.0
  gi1%pr(2:nkp) = pmn(:)

  gi2%pr(1) = pmn(1)*0.5
  do k = 2, nk
    gi2%pr(k) = ( pmn(k) + pmn(k-1) )*0.5
  end do

  ! dz
  gi1%dz = 0.0
  gi1%dz(1) = 0.5*dz(1)
  do k = 2, nk
    gi1%dz(k) = z(k) - z(k-1)
  end do
  gi1%dz(nkp) = gi1%dz(nk)

  gi2%dz(:) = dz(:)

  ! dpr
  gi1%dpr(1) = 0.5*pmn(1)
  do k = 2, nk
    gi1%dpr(k) = gi2%pr(k) - gi2%pr(k-1)
  end do
  gi1%dpr(nkp) = gi1%dpr(nk)

  gi2%dpr(1) = pmn(1)
  do k = 2, nk
    gi2%dpr(k) = pmn(k) - pmn(k-1)
  end do

  giw => gi1
  git => gi2

end subroutine set_vgrid

subroutine set_grid (g12, g32, gt, gu) !{{{1
  ! 3d model stagger grid
  type (type_gij),target :: g12, g32
  type (type_gij), pointer :: gt, gu

  integer, allocatable, dimension(:,:) :: glo_itn
  real (kind=wp), dimension(nk) :: phimn
  integer :: i, j, k, is

  if ( myid == mid ) then
    allocate(glo_itn(glo_ni, glo_nj), stat=is); call chk(is)
    glo_itn = 0
    call io_read (nm%fi, 'itn', glo_itn)
  end if
  call mympi_div (glo_itn, g12%lev)

  ! set to land first
  g12%msk(:,:,:) = 0 
  g32%msk(:,:,:) = 0 

  do k = 1, nk
  do j = 2, njm
  do i = 2, nim
    if ( g12%lev(i,j) >= k ) g12%msk(i,j,k) = 1
  end do
  end do
  end do
  call mympi_swpbnd ( g12%msk )

  ! a water U-grid should be surounded by T-grids
  do k = 1, nk
  do j = 1, njm
  do i = 1, nim
    g32%msk(i,j,k) = g12%msk(i,j,k) * g12%msk(i+1,j,k)  &
                  * g12%msk(i,j+1,k) * g12%msk(i+1,j+1,k)
  end do
  end do
  end do
  call mympi_swpbnd ( g32%msk )

  g32%lev(:,:) = sum(g32%msk(:,:,:), 3)

  ! sea bottom pressure, geoptential height
  phimn(1) = - g * gi2%dz(1)
  do k = 2, nk
    phimn(k) = phimn(k-1)  - g * gi2%dz(k)
  end do

  do j = 1, nj
  do i = 1, ni
    k = g12%lev(i, j)
    if ( k == 0 ) then
      g12%phih(i,j) = 0.0
    else
      g12%phih(i,j) = phimn(k)
      g12%prh(i,j) = gi1%pr(k+1)
    end if

    k = g32%lev(i, j)
    if ( k == 0 ) then
      g32%phih(i,j) = 0.0
    else
      g32%phih(i,j) = phimn(k)
      g32%prh(i,j) = gi1%pr(k+1)
    end if
  end do
  end do

  ! in case of 1/prh, set to 0.1Pa in land
  where( g12%lev == 0 ) g12%prh = 0.1
  where( g32%lev == 0 ) g32%prh = 0.1

  gt=>g12
  gu=>g32
end subroutine set_grid

subroutine calc_cf (cv1, cv2) !{{{1
  ! calc. coefficents of frictional force at g3j grid
  real (kind=wp), dimension(ni,nj) :: cv1, cv2
  real (kind=wp), dimension(ni,nj) :: temp

  ! pcom 1.0 use this inderct method instead of directly use tan function
  ! these two ways differ a little bit
!  temp = tan( hpos(:,:,3)%x(2) )
  temp = sqrt(1-cos(hpos(:,:,3)%x(2))**2) / cos(hpos(:,:,3)%x(2))
  cv1 = (1 - temp*temp) / a**2
  cv2 = 2*temp / a

end subroutine calc_cf

subroutine inistat (ts, frc, graphih) !{{{1
  ! prepare the initial state of the ocean
  type (type_gvar_m3d) :: ts(2)
  type (type_frc) :: frc
  type (type_gvar_m2d) :: graphih

  type (type_mat), dimension(ni,nj) :: wk
  real (kind=wp), allocatable, dimension(:,:,:) :: glo_pt, glo_sa
  type (type_frc) :: glo_frc
  integer :: is

  if ( myid == mid ) then
    ! 12 months forcing
    call arrays_init(glo_frc%tau%x(1), glo_ni,glo_nj,12, zero, frc%tau%x(1)%g)
    call arrays_init(glo_frc%tau%x(2), glo_ni,glo_nj,12, zero, frc%tau%x(2)%g)
    call arrays_init(glo_frc%ts%x(1),  glo_ni,glo_nj,12, zero, frc%ts%x(1)%g)
    call arrays_init(glo_frc%ts%x(2),  glo_ni,glo_nj,12, zero, frc%ts%x(2)%g)
    call arrays_init(glo_frc%pa,       glo_ni,glo_nj,12, zero, frc%pa%g)
    call arrays_init(glo_frc%fw,       glo_ni,glo_nj,12, zero, frc%fw%g)

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
!  where (spread(g12%lev, 3, 12) == 0)
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
  call op_gra( wk, gt%phih, gtj, gtj%ew, gtj%ns)
  call op_ter( graphih%x(1)%v, wk%x(1), gtj%ew, graphih%x(1)%hg )
  call op_ter( graphih%x(2)%v, wk%x(2), gtj%ns, graphih%x(2)%hg )
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

  ! the print elapsed time infos if neccessary
  if ( nm%per.eq.'hour' ) then 
    if ( tctr%ct%h /= tctr%pt%h ) call print_time_per_hour (tctr) 
  else 
    if ( tctr%ct%d /= tctr%pt%d ) call print_time_per_day (tctr)
  end if

  ! output if time is proper
  if ( ((tctr%ct%h/=tctr%pt%h).and.(nm%per.eq.'hour'))  .or. &
       ((tctr%ct%d/=tctr%pt%d).and.(nm%per.eq.'day'))   .or. &
       ((tctr%ct%m/=tctr%pt%m).and.(nm%per.eq.'month')) .or. &
       ((tctr%ct%y/=tctr%pt%y).and.(nm%per.eq.'year')) ) then

    call mympi_output (nm%fo, names%pt, names%sa, acts)

    call mympi_output (nm%fo, names%u, names%v, acuv, gu%msk)

    call mympi_output (nm%fo, names%w, acw, gu%msk)

    ! calc. fluctuation of ssh, minus global mean
    call mympi_output (nm%fo, names%ssh, acssh, &
      acssh%var%hg%rh*gt%msk(:,:,1), gt%msk(:,:,1))

  end if

end subroutine check_output

subroutine print_time_per_hour (tctr)
    type (type_tctr):: tctr

    tctr%t2 = mpi_wtime ()

    if ( myid==mid ) & 
      write(*, '(a, i0.4,a,i0.2,a,i0.2,a,i0.2, a, f8.2, a, i3)') &
      'integrate for ', tctr%pt%y, '-', tctr%pt%m, &
      '-', tctr%pt%d, ' ', tctr%ct%h, ':00:00 use ', tctr%t2 - tctr%t1, &
      ' seconds on processor ', myid

    tctr%t1 = tctr%t2
end subroutine 

subroutine print_time_per_day (tctr)
    type (type_tctr):: tctr

    tctr%t2 = mpi_wtime ()

    if ( myid==mid ) & 
      write(*, '(a, i0.4,a,i0.2,a,i0.2, a, f8.2, a, i3)') &
      'integrate for ', tctr%pt%y, '-', tctr%pt%m, &
      '-', tctr%pt%d, ' use ', tctr%t2 - tctr%t1, &
      ' seconds on processor ', myid

    tctr%t1 = tctr%t2
end subroutine 

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