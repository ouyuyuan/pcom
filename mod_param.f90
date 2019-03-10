
! Description: global parameters
!
!      Author: OU Yuyuan <ouyuyuan@lasg.iap.ac.cn>
!     Created: 2015-02-26 08:20:12 BJT
! Last Change: 2019-03-10 08:55:47 BJT

module mod_param

  use mpi
  use mod_kind, only: sglp, wp
  use mod_type, only: type_vars_info, type_rst_info, type_var_info

  implicit none
  public

  ! exported variables !{{{1

  ! mpi envioronment
  type :: type_my
    integer :: id ! rank of the group
    integer :: i, j ! coordinate in the topology of processors
    integer :: ni, nj ! number of grid points in east-west/north-south directions
    integer :: s, n, w, e ! south/north/west/east/upper/lower neighbor ids
    integer :: gn, gs, gw, ge ! indices in the global array
  end type type_my
  type (type_my) :: my
  type (type_my), allocatable :: our(:)

  ! output variable informations
  type (type_vars_info) :: vars_info

  ! restart file
  type (type_rst_info) :: rst_info

  ! namelist parameters
  ! if you change this type definition, remember also change the following:
  !   1. read_namelist () in main.f90
  !   2. file namelist
  type :: type_nm
    sequence
    ! number of Processors in horizontal directions
    integer :: py, px

    ! Begin/End Date for integration, in the form of yyyy-mm-dd hh:mm:ss
    character (len=80) :: bd, ed
    ! BaroTropic/BaroClinic time step
    integer :: bt, bc

    ! output one restart file after output every 'rst_per' records
    integer :: rst_per, rst

    ! Filename of Initial/Forcing file
    ! directory name of output files
    character (len=80) :: fi, ff, ff_ws, ff_pa, od
    ! output PER month or PER year 
    character (len=80) :: out_per
  end type type_nm
  type (type_nm) :: nm

  integer, parameter :: missing_int = -999
  real (kind=sglp), parameter :: missing_float = -9.99999e30
  real (kind=wp), parameter :: missing_double = -9.99999e30

  ! parallel
  integer, parameter :: mid = 0 ! Master ID

  integer :: &
    myid,  & ! rank in the mpi_comm_world
    npro 

  integer, allocatable, dimension(:,:) :: ids

  ! dimension
  integer :: glo_ni, glo_nj, ni, nj, nk, &
    nim, njm, nkp ! minus/plus 1, frequently used, for simplity

  ! other
  integer, parameter :: tc = 1 ! current time value
  integer, parameter :: tp = 2 ! previous time value
  integer, parameter :: tpp = 3 ! pre-previous time value

  ! local variables !{{{1
  integer :: is, err

contains !{{{1

subroutine param_set_io ( ) !{{{1
  ! set output infomations

  call set_var_info (vars_info%pt, 'pt', &
    'potential temperature', 'degrees Celcius', 'mean')
  call set_var_info (vars_info%sa, 'sa', &
    'salinity', 'psu', 'mean')
  call set_var_info (vars_info%u, 'u', &
    'zonal velocity', 'm/s', 'mean')
  call set_var_info (vars_info%v, 'v', &
    'meridional velocity', 'm/s', 'mean')
  call set_var_info (vars_info%w, 'w', &
    'vertical velocity', 'm/s', 'mean')
  call set_var_info (vars_info%rho, 'rho', &
    'sea water density', 'kg/m^3', 'mean')
  call set_var_info (vars_info%ssh, 'ssh', &
    'sea surface height', 'm', 'mean')
  call set_var_info (vars_info%ch, 'ch', &
    'normalized sea bottom pressure', '', 'mean')
  call set_var_info (vars_info%prh, 'prh', &
    'reference state: sea bottom pressure (T-grid)', 'Pa', 'mean')

  ! set restart infomations
  if (nm%rst == 1) then 
    rst_info%fname = trim(nm%od)//'restart.nc'
  end if

  call set_var_info (rst_info%chc, 'chc', &
    vars_info%ch%longname, vars_info%ch%units, 'restart')
  call set_var_info (rst_info%chp, 'chp', &
    vars_info%ch%longname, vars_info%ch%units, 'restart')

  call set_var_info (rst_info%bnd_taux , 'bnd_taux' , ' ' , ' ' , 'restart')
  call set_var_info (rst_info%bnd_tauy , 'bnd_tauy' , ' ' , ' ' , 'restart')
  call set_var_info (rst_info%bnd_t    , 'bnd_t'    , ' ' , ' ' , 'restart')
  call set_var_info (rst_info%bnd_s    , 'bnd_s'    , ' ' , ' ' , 'restart')
  call set_var_info (rst_info%bnd_pa   , 'bnd_pa'   , ' ' , ' ' , 'restart')
  call set_var_info (rst_info%bnd_fw   , 'bnd_fw'   , ' ' , ' ' , 'restart')

  call set_var_info (rst_info%buc, 'buc', ' ', ' ', 'restart')
  call set_var_info (rst_info%bvc, 'bvc', ' ', ' ', 'restart')
  call set_var_info (rst_info%bup, 'bup', ' ', ' ', 'restart')
  call set_var_info (rst_info%bvp, 'bvp', ' ', ' ', 'restart')

  call set_var_info (rst_info%prpt, 'prpt', ' ', ' ', 'restart')
  call set_var_info (rst_info%prps, 'prps', ' ', ' ', 'restart')

  call set_var_info (rst_info%tc, 'tc', &
    vars_info%pt%longname, vars_info%pt%units, 'restart')
  call set_var_info (rst_info%sc, 'sc', &
    vars_info%sa%longname, vars_info%sa%units, 'restart')
  call set_var_info (rst_info%tp, 'tp', &
    vars_info%pt%longname, vars_info%pt%units, 'restart')
  call set_var_info (rst_info%sp, 'sp', &
    vars_info%sa%longname, vars_info%sa%units, 'restart')

  call set_var_info (rst_info%uc, 'uc', &
    vars_info%u%longname, vars_info%u%units, 'restart')
  call set_var_info (rst_info%vc, 'vc', &
    vars_info%v%longname, vars_info%v%units, 'restart')
  call set_var_info (rst_info%up, 'up', &
    vars_info%u%longname, vars_info%u%units, 'restart')
  call set_var_info (rst_info%vp, 'vp', &
    vars_info%v%longname, vars_info%v%units, 'restart')
  
  call set_var_info (rst_info%auc, 'auc', &
    'advection of velocity', 'm/s^2', 'restart')
  call set_var_info (rst_info%avc, 'avc', &
    'advection of velocity', 'm/s^2', 'restart')
  call set_var_info (rst_info%aup, 'aup', &
    'advection of velocity', 'm/s^2', 'restart')
  call set_var_info (rst_info%avp, 'avp', &
    'advection of velocity', 'm/s^2', 'restart')
  call set_var_info (rst_info%aupp, 'aupp', &
    'advection of velocity', 'm/s^2', 'restart')
  call set_var_info (rst_info%avpp, 'avpp', &
    'advection of velocity', 'm/s^2', 'restart')

  call set_var_info (rst_info%am, 'am', &
    'horizontal mixing coeficient', '', 'restart')
end subroutine param_set_io

subroutine set_var_info ( var_info, varname, longname, units, vartype ) !{{{1
  ! set variable infomation
  type (type_var_info) :: var_info
  character (len=*), intent(in) :: varname, longname, units, vartype

  var_info%name     = varname
  var_info%longname = longname
  var_info%units    = units
  var_info%vartype  = vartype

end subroutine set_var_info

subroutine print_my (one) !{{{1
  ! print one define in mod_param
  type (type_my) :: one
  
  write(*,'(a,i3, 6(a,i2), 4(a, i3))') &
    'id=',  one%id, &
    ', i=',     one%i, ', j=',     one%j, &
    ', west=',  one%w, ', east=',  one%e, &
    ', south=', one%s, ', north=', one%n, &
    ', gw=', one%gw, ', ge=', one%ge, &
    ', gs=', one%gs, ', gn=', one%gn 
  
end subroutine print_my

subroutine param_set_nm ( nm ) !{{{1
  ! set the parameters in the namelist structure
  type (type_nm) :: nm

  integer, parameter :: fid_nam = 20

  ! namelist
  integer :: npy, npx

  character (len=80) :: bdate, edate
  integer :: dtbt, dtbc, restart_per, restart_run

  character (len=80) :: fname_ini, fname_frc, fname_frc_ws, fname_frc_pa, out_dir, out_per

  ! it can't use the syntax % for namelist /.../ statesment
  !    so we redefine namelist parameters

  namelist /mpi_ctrl/ npy, npx
  namelist /time_ctrl/ bdate, edate, dtbt, dtbc
  namelist /io_ctrl/ fname_ini, fname_frc, fname_frc_ws, fname_frc_pa, out_dir, out_per, &
           restart_per, restart_run

  open(fid_nam, file='namelist')
  read(fid_nam, mpi_ctrl)
  read(fid_nam, time_ctrl)
  read(fid_nam, io_ctrl)
  close(fid_nam)

  ! check integrated time step

  if ( mod(24*60*60, dtbt) /= 0 ) stop &
    'dtbt should be divided by seconds in a day'

  if ( mod(24*60*60, dtbc) /= 0 ) stop &
    'dtbc should be divided by seconds in a day'

  if ( mod(dtbc, dtbt) /=0 ) stop &
    'dtbc should be divived by dtbt'

  nm%py = npy
  nm%px = npx

  nm%bd = bdate
  nm%ed = edate
  nm%bt = dtbt
  nm%bc = dtbc

  nm%fi      = fname_ini
  nm%ff      = fname_frc
  nm%ff_ws   = fname_frc_ws
  nm%ff_pa   = fname_frc_pa
  nm%od      = out_dir
  nm%out_per = out_per
  nm%rst_per = restart_per
  nm%rst     = restart_run

end subroutine param_set_nm

subroutine param_set_my ( my )  !{{{1
  ! set mpi groups
  type (type_my) :: my

  integer :: i

  call set_ids( ids )

  ! set my
  my%id = myid
  call get_neighbor (my, ids)

  ! set our for mid
  if ( myid == mid) then
    allocate(our(npro), stat=is); call chk(is)
    do i = 1, npro
      our(i)%id = i - 1
      call get_neighbor (our(i), ids)
    end do
  end if
end subroutine param_set_my

subroutine get_neighbor (one, ids)  !{{{1
  ! get neighbor ids and environments in the global array

  integer, intent(in) :: ids(:,:)
  type (type_my) :: one

  integer :: ind(2), n, m, d1, d2, i, j

  d1 = size(ids, 1)
  d2 = size(ids, 2)

  ! get postion in the ids topology !{{{2
  ind(:) = maxloc(ids(:,:), mask = ids(:,:) <= one%id)

  one%i = ind(1)
  one%j = ind(2)

  i = one%i
  j = one%j

  ! north/south !{{{2
  ! get neighbors !{{{3
  if ( one%j /= 1 )  one%s = ids(i, j-1)
  if ( one%j == 1 )  one%s = mpi_proc_null

  if ( one%j /= d2 ) one%n = ids(i, j+1)
  if ( one%j == d2 ) one%n = mpi_proc_null

  ! get positions in the global array !{{{3
  ! the first and last element in each dimension is for boundary
  n = glo_nj / nm%py
  m = mod(glo_nj, nm%py)

  ! put the remainder m grids on the first m processors
  if ( one%j <= m ) then
    one%nj = n + 2 + 1
    one%gn = one%j * (n+1)
  else
    one%nj = n + 2
    one%gn = one%j * n + m
  end if
  one%gs = one%gn - (one%nj - 2) + 1

  ! east-west !{{{2
  ! east-west is wrap up
  ! get neighbors !{{{3
  if ( one%i /= 1 )  one%w = ids(i-1, j)
  if ( one%i == 1 )  one%w = ids(d1,  j)

  if ( one%i /= d1 ) one%e = ids(i+1, j)
  if ( one%i == d1 ) one%e = ids(1,   j)

  ! get positions in the global array !{{{3
  n = glo_ni / nm%px
  m = mod(glo_ni, nm%px)

  if ( one%i <= m ) then
    one%ni = n + 2 + 1
    one%ge = one%i * (n+1)
  else
    one%ni = n + 2
    one%ge = one%i * n + m
  end if
  one%gw = one%ge - (one%ni - 2) + 1

end subroutine get_neighbor

subroutine set_ids ( ids ) !{{{1
  ! if nm%px = 4, nm%py = 2, then ids is:
  !   (NOTE that east-west direction is wrap up, 
  !     and the map is 'a reverse matrix')
  !
  !    0 4
  !    1 5
  !    2 6
  !    3 7
  !
  integer, allocatable, dimension(:,:) :: ids

  integer, allocatable, dimension(:) :: ia
  integer :: i

  allocate(ids(nm%px, nm%py), stat=is);
  call chk(is); ids = 0

  allocate(ia(npro), stat=is); call chk(is); ia = 0

  forall ( i=1:npro ) ia(i) = i - 1
  ids = reshape(ia, shape(ids))

end subroutine set_ids

subroutine chk( ista ) !{{{1
  ! check state of allocate array 

  integer, intent(in) ::  ista

  if ( ista /= 0 ) then
    write(*,*) 'Allocate array failed. Stop'
    stop 2
  end if
end subroutine chk

end module mod_param !{{{1
!-------------------------------------------------------{{{1
! vim:fdm=marker:fdl=0:
! vim:foldtext=getline(v\:foldstart).'...'.(v\:foldend-v\:foldstart):
