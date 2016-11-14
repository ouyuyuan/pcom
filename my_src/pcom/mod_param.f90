
! Description: global parameters
!
!      Author: OU Yuyuan <ouyuyuan@lasg.iap.ac.cn>
!     Created: 2015-02-26 08:20:12 BJT
! Last Change: 2016-04-05 19:42:16 BJT

module mod_param

  use mpi
  use mod_kind, only: sglp, wp
  use mod_pro, only: pro_print
  use mod_type, only: type_time, type_mat

  implicit none
  public

  ! exported variables !{{{1

  ! mpi group, the whole world
  type :: type_my
    integer :: id ! communicator/group/rank of the group
    integer :: i, j ! coordinate in the topology of processors
    integer :: ni, nj ! number of grid points in 3d directions
    integer :: s, n, w, e ! south/north/west/east/upper/lower neighbor ids
    integer :: gn, gs, gw, ge ! indices in the global array
  end type type_my
  type (type_my) :: my
  type (type_my), allocatable :: our(:)

  ! output variable names
  type :: type_names
    character (len=2) :: pt = 'pt'
    character (len=2) :: sa = 'sa'
    character (len=1) :: u = 'u'
    character (len=1) :: v = 'v'
    character (len=1) :: w = 'w'
    character (len=len('ssh')) :: ssh = 'ssh'
  end type type_names
  type (type_names) :: names 

  ! units of variables
  type :: type_units
    character( len=len('degrees Celcius') ) :: pt = 'degrees Celcius'
    character( len=len('psu') ) :: sa = 'psu'
    character( len=len('m/s') ) :: u = 'm/s'
    character( len=len('m/s') ) :: v = 'm/s'
    character( len=len('m/s') ) :: w = 'm/s'
    character( len=len('m') ) :: ssh = 'm'
  end type type_units
  type (type_units) :: units

  ! namelist parameters
  ! if you change this type definition, remember also change the following:
  !   1. read_namelist () in main.f90
  !   2. file namelist
  type :: type_nm
    sequence
    ! number of Processors in 3d directions
    integer :: py, px

    ! Begin/End Date for integration, in the form of yyyy-mm-dd hh:mm:ss
    character (len=80) :: bd, ed
    ! BaroTropic/BaroClinic time step
    integer :: bt, bc 

    ! Filename of Initial/Forcing/Output file
    character (len=80) :: fi, ff, fo
    ! output PER month or PER year 
    character (len=80) :: per
  end type type_nm
  type (type_nm) :: nm

  integer, parameter :: missing_int = -999
  real (kind=sglp), parameter :: missing_float = -9.99999e30

  ! parallel
  integer, parameter :: mid = 0 ! Master ID

  integer :: &
    myid,  & ! rank in the mpi_comm_world
    npro 

  integer, allocatable, dimension(:,:) :: ids

  ! dimension
  integer :: glo_ni, glo_nj, ni, nj, nk, &
    nim, njm, nkp ! minus/plus 1, frequently used, for simplity

  ! integration time control !{{{1
  type :: type_tctr
    type (type_time) :: pt, ct ! previous/current baroclinic step
    real (kind=wp)   :: t1, t2 ! mpi time barrier
    integer*8        :: nt     ! total baroclinic time steps
    integer*8        :: i      ! current baroclinic time steps
  end type
  type (type_tctr) :: tctr

  ! other
  integer, parameter :: tc = 1 ! current time value
  integer, parameter :: tp = 2 ! previous time value
  integer, parameter :: tpp = 3 ! pre-previous time value

  ! local variables !{{{1
  integer :: is, err

contains !{{{1

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

subroutine param_set_nm () !{{{1
  ! set the parameters in the namelist structure

  integer, parameter :: fid_nam = 20

  ! namelist
  integer :: npy, npx

  character (len=80) :: bdate, edate
  integer :: dtbt, dtbc

  character (len=80) :: fname_ini, fname_frc, fname_out, out_per

  ! it can't use the syntax % for namelist /.../ statesment
  !    so we redefine namelist parameters

  namelist /mpi_ctrl/ npy, npx
  namelist /time_ctrl/ bdate, edate, dtbt, dtbc
  namelist /io_ctrl/ fname_ini, fname_frc, fname_out, out_per

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

  nm%fi = fname_ini
  nm%ff = fname_frc
  nm%fo = fname_out
  nm%per= out_per

end subroutine param_set_nm

subroutine param_set_my ()  !{{{1
  ! set mpi groups
  integer :: i

  call set_ids ()

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

subroutine set_ids () !{{{1
  ! if nm%px = 4, nm%py = 2, then ids is:
  !   (NOTE that east-west direction is wrap up, 
  !     and the map is 'a reverse matrix')
  !
  !    0 4
  !    1 5
  !    2 6
  !    3 7
  !
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
