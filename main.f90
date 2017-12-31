
! Description: 
!
!      Author: OU Yuyuan <ouyuyuan@lasg.iap.ac.cn>
!     Created: 2015-09-13 08:14:52 BJT
! Last Change: 2017-12-31 16:29:09 BJT

program main

  ! imported variables !{{{1
  !-------------------------------------------------------=-
  use mod_arrays, only: &
    arrays_init, &
    cv1, cv2, &
    g1j, g2j, g3j, g4j, gtj, guj, &
    gi1, gi2, git, giw, g12, g32, gt, gu, &
    acssh, & 
    am, arrays_allocate, &
    bnd, bphi, bgraphi, &
    frc, cor, &
    glo_lat, glo_lon, graphihx, graphihy, &
    hposx, hposy, &
    km, kh, &
    lon, lat, &
    alpha, adp, &
    eqts, equv, equvb, eqw, eqch, &
    prho, &
    wm, &
    z
    
  use mod_con, only: &
    rho0, g, km_c, a, torad, omega

  use mod_den, only: den_alpha, den_rho, den_prho

  use mod_int, only: int_trop, int_bnd, &
    int_pgra, int_readyc, int_clin, int_ts, &
    int_ssh

  use mod_io, only: io_get_dim_len, io_read, io_create_rst

  use mod_kind, only: wp, zero, lint, one

  use mpi

  use mod_mympi, only: &
    mympi_div, mympi_divx, mympi_divy, &
    mympi_output, mympi_swpbnd, mympi_bcast, &
    mympi_input

  use mod_op, only: op_ter, op_gra

  use mod_param, only: &
    nm, my, missing_float, &
    glo_ni, glo_nj, &
    ni, nj, nk, nim, njm, nkp, vars_info, &
    myid, npro, mid, &
    param_set_my, param_set_io, &
    print_my, &
    rst_info, &
    tc, tp

  use mod_type, only: &
    type_accu_gr2d, &
    type_gi, type_gj, type_gij, &
    type_str2time, type_time2str, type_str2sec, &
    type_frc, tctr, type_tctr, type_rst_info, &
    type_gvar_r2d, type_gvar_r3d, &
    type_eq_ts, type_eq_uv, type_eq_ch, &
    operator (+), operator (<)

  use mod_debug, only: debug_var

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
  call calc_hpos ()

  ! setting stagger grids
  call set_hgrid (g1j, g2j, g3j, g4j, gtj, guj) 
  call set_vgrid (gi1, gi2, git, giw)
  call set_grid  (g12, g32, gt, gu)

  ! calc. grid variables
  call calc_cf( cv1, cv2 )

  ! Coriolis force
  cor%v = 2.0*omega*sin(hposy(:,:,cor%hg%n))

  ! prepare initial state of the ocean
  call inistat(eqts, frc, graphihx, graphihy)

  ! time prepare !{{{2
  tctr%ct = type_str2time (nm%bd)
  tctr%pt = tctr%ct
  tctr%nt = ( type_str2sec(nm%ed) - type_str2sec(nm%bd) ) / nm%bc
  tctr%ist= 1

  if (nm%rst == 1) then ! restart from a moving ocean !{{{2
    ! read in variables from restart file
    call mympi_input ( rst_info )

    ! modify tctr
    tctr%ct = type_str2time (rst_info%cdate)
    tctr%pt = type_str2time (rst_info%pdate)
    tctr%ist = ( type_str2sec(rst_info%cdate) - type_str2sec(nm%bd) ) / nm%bc + 1

    if (myid == mid) print *, "I'm going for a restart-run from "//&
      trim(rst_info%fname)//" at "//trim(rst_info%cdate)//"."
  end if

  ! output sea bottom pressure of reference state !{{{2
  if ( nm%rst /= 1 ) call mympi_output (vars_info%prh, gt%prh, gt%msk(:,:,1))

  ! integration cycle !{{{1

  tctr%t1 = mpi_wtime () ! time count start

  ! baroclinic integrate !{{{2

  do i = tctr%ist, tctr%nt
    tctr%i = i

    ! linear interpolate monthly forcing to daily values of bnd
    !   at the initial of integration, or when the 
    !   integration enters a new day
    if ( ( i == 1 ) .or. (tctr%ct%d /= tctr%pt%d) ) then 
      call int_bnd ( bnd )
      call den_prho ( prho )
    end if

    ! calc. the specific volume, reciprocal of density
    call den_alpha ( alpha%v, eqts%tc, eqts%sc, eqch%chc, gt%msk )

    ! calc. pressure gradient forces
    call int_pgra (adp, bphi, bgraphi)

    ! prepare for baroclinic integration
    call int_readyc (equv, equvb, wm)

    ! integrate series time steps per baroclinic time step
    call int_trop (equvb, eqch)

    ! prediction of baroclinic mode
    call int_clin (equv, am%v)

    ! calc. sea surface height
    call int_ssh (acssh, alpha%v)

    ! prediction of temperature and salinity
    call int_ts (eqts, eqw, wm)

    tctr%pt = tctr%ct
    tctr%ct = tctr%ct + nm%bc

    ! output at proper time
    if ( ((tctr%ct%h/=tctr%pt%h).and.(nm%out_per.eq.'hour'))  .or. &
         ((tctr%ct%d/=tctr%pt%d).and.(nm%out_per.eq.'day'))   .or. &
         ((tctr%ct%m/=tctr%pt%m).and.(nm%out_per.eq.'month')) .or. &
         ((tctr%ct%y/=tctr%pt%y).and.(nm%out_per.eq.'year')) ) then
      call output ()
    end if
  end do

  ! finish integration !{{{1
  if (myid == mid) then
    write (*, *) ''
    write (*, '(a)') 'Done model run.'
    write (*, '(a, i9, a)') &
      'Integrated for total ', tctr%nt-tctr%ist+1, ' baroclinic time step(s).'
  end if

  call mpi_finalize (err)

contains  !{{{1

subroutine init () !{{{1
  ! initialize model environment
  
  ! get namelist
  call mympi_input (nm)

  ! set nc_variable infomations
  call param_set_io ()

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

subroutine set_cf (km, kh) !{{{1
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

subroutine calc_hpos () !{{{1
  ! calc. horizontal position (hposx, hposy) on grid 1/2/3/4
  ! (lon, lat) in I/O file is on grid 4

  real (kind=wp) :: px, py, loni, latj
  integer :: i, j

  ! inner area !{{{2
  do j = 2, njm
  do i = 2, nim
    !westest/eastest longitude wrap up
    loni = lon(i)
    latj = (lat(j) + lat(j-1)) * 0.5
    px = loni
    py = latj
    hposx(i,j,1) = px * torad
    hposy(i,j,1) = py * torad

    loni = (lon(i) + lon(i+1)) * 0.5
    ! the 'eastest' point will be 360.0 degree, for keeping
    !   monotonous, we donot change it to 0 degree
    if ( lon(i+1) < lon(i) ) loni = (lon(i) + lon(i+1)+360) * 0.5
    latj = (lat(j) + lat(j-1)) * 0.5
    px = loni
    py = latj
    hposx(i,j,2) = px * torad
    hposy(i,j,2) = py * torad

    loni = (lon(i) + lon(i+1)) * 0.5
    if ( lon(i+1) < lon(i) ) loni = (lon(i) + lon(i+1)+360) * 0.5
    latj = lat(j)
    px = loni
    py = latj
    hposx(i,j,3) = px * torad
    hposy(i,j,3) = py * torad

    loni = lon(i)
    latj = lat(j)
    px = loni
    py = latj
    hposx(i,j,4) = px * torad
    hposy(i,j,4) = py * torad
  end do
  end do
  call mympi_swpbnd (hposx)
  call mympi_swpbnd (hposy)

  ! extropolate  !{{{2

  if (my%gs == 1) then
    hposx(:,1,:) = 2*hposx(:,2,:) - hposx(:,3,:)
    hposy(:,1,:) = 2*hposy(:,2,:) - hposy(:,3,:)
  end if

  if (my%gn == glo_nj) then
    hposx(:,nj,:) = 2*hposx(:,nj-1,:) - hposx(:,nj-2,:)
    hposy(:,nj,:) = 2*hposy(:,nj-1,:) - hposy(:,nj-2,:)
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
    dx1 = hposx(i,j,2) - hposx(i-1,j,2)
    if ( dx1 < 0 ) dx1 = dx1 + 360.0*torad
    dx2 = hposy(i,j,4) - hposy(i,j-1,4)
    g1j%dx(i,j)%x(1) = dx1
    g1j%dx(i,j)%x(2) = dx2

    dx1 = hposx(i+1,j,1) - hposx(i,j,1)
    if ( dx1 < 0 ) dx1 = dx1 + 360.0*torad
    dx2 = hposy(i,j,3) - hposy(i,j-1,3)
    g2j%dx(i,j)%x(1) = dx1
    g2j%dx(i,j)%x(2) = dx2

    dx1 = hposx(i+1,j,4) - hposx(i,j,4)
    if ( dx1 < 0 ) dx1 = dx1 + 360.0*torad
    dx2 = hposy(i,j+1,2) - hposy(i,j,2)
    g3j%dx(i,j)%x(1) = dx1
    g3j%dx(i,j)%x(2) = dx2

    dx1 = hposx(i,j,3) - hposx(i-1,j,3)
    if ( dx1 < 0 ) dx1 = dx1 + 360.0*torad
    dx2 = hposy(i,j+1,1) - hposy(i,j,1)
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
    g1j%rh(i,j) = cos( hposy(i,j,1) )
    g2j%rh(i,j) = cos( hposy(i,j,2) )
    g3j%rh(i,j) = cos( hposy(i,j,3) )
    g4j%rh(i,j) = cos( hposy(i,j,4) )

    g1j%tn(i,j) = tan( hposy(i,j,1) )
    g2j%tn(i,j) = tan( hposy(i,j,2) )
    g3j%tn(i,j) = tan( hposy(i,j,3) )
    g4j%tn(i,j) = tan( hposy(i,j,4) )
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
!  temp = tan( hposy(:,:,3) )
  temp = sqrt(1-cos(hposy(:,:,3))**2) / cos(hposy(:,:,3))
  cv1 = (1 - temp*temp) / a**2
  cv2 = 2*temp / a

end subroutine calc_cf

subroutine inistat (eqts, frc, graphihx, graphihy) !{{{1
  ! prepare the initial state of the ocean
  type (type_eq_ts) :: eqts
  type (type_frc) :: frc
  type (type_gvar_r2d) :: graphihx, graphihy

  real (kind=wp), dimension(ni,nj) :: wkx, wky

  call mympi_input (frc)

  ! forcing on land set to zero (no forcing)
  !!! potential bug, missing is not all at land of input file
!  where (spread(g12%lev, 3, 12) == 0)
  where (frc%taux%v == missing_float)
    frc%taux%v = 0.0 
    frc%tauy%v = 0.0 
    frc%t%v = 0.0 
    frc%s%v = 0.0 
    frc%pa%v = 0.0 
    frc%fw%v = 0.0 
  end where

  call mympi_input (eqts)

  !! why not average before differentiate?
  call op_gra( gt%phih, gtj, wkx, wky)
  call op_ter( graphihx%v, wkx, gtj%ew, graphihx%hg )
  call op_ter( graphihy%v, wky, gtj%ns, graphihy%hg )

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

subroutine output () !{{{1
  ! time averaged model output

  integer, save :: nrec = 0

  ! print elapsed time infos if neccessary
  if ( nm%out_per.eq.'hour' ) then 
    if ( tctr%ct%h /= tctr%pt%h ) call print_time_per_hour (tctr) 
  else 
    if ( tctr%ct%d /= tctr%pt%d ) call print_time_per_day (tctr)
  end if

  call mympi_output (equv, eqw)

  call mympi_output (eqts)

  ! output ch
  eqch%acch = eqch%acch / eqch%n
  where (gt%msk(:,:,1) == 0)
    eqch%acch = missing_float
  end where
  call mympi_output (vars_info%ch, eqch%acch)
  eqch%acch  = 0.0 ! reset accumulated value
  eqch%n     = 0 ! reset counter

  ! calc. fluctuation of ssh, minus global mean
  call mympi_output (vars_info%ssh, &
    acssh, acssh%var%hg%rh*gt%msk(:,:,1), gt%msk(:,:,1))

  nrec = nrec + 1

  ! output restart file
  if ( mod(nrec, nm%rst_per) == 0 ) call mympi_output (rst_info)

end subroutine output

subroutine print_time_per_hour (tctr) !{{{1
    type (type_tctr):: tctr

    tctr%t2 = mpi_wtime ()

    if ( myid==mid ) & 
      write(*, '(a, i0.4,a,i0.2,a,i0.2,a,i0.2, a, f8.2, a, i3)') &
      'integrated for ', tctr%pt%y, '-', tctr%pt%m, &
      '-', tctr%pt%d, ' ', tctr%ct%h, ':00:00,  used ', tctr%t2 - tctr%t1, &
      ' seconds on processor ', myid

    tctr%t1 = tctr%t2
end subroutine print_time_per_hour

subroutine print_time_per_day (tctr) !{{{1
    type (type_tctr):: tctr

    tctr%t2 = mpi_wtime ()

    if ( myid==mid ) & 
      write(*, '(a, i0.4,a,i0.2,a,i0.2, a, f8.2, a, i3)') &
      'integrated for ', tctr%pt%y, '-', tctr%pt%m, &
      '-', tctr%pt%d, ', used ', tctr%t2 - tctr%t1, &
      ' seconds on processor ', myid

    tctr%t1 = tctr%t2
end subroutine print_time_per_day

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
