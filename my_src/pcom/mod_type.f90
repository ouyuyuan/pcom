
! Description: data-type definition
!
!      Author: OU Yuyuan <ouyuyuan@lasg.iap.ac.cn>
!     Created: 2015-02-26 08:20:12 BJT
! Last Change: 2016-04-05 18:27:10 BJT

module mod_type

  use mod_kind, only: sglp, wp
  use mod_pro, only: pro_days_month, pro_days_year

  implicit none
  public

  ! u, v are in half grid, T, S in whole grid

  ! memo !{{{1
  type :: type_memo_r
    real (kind=wp) :: p, c ! previous/current values
  end type type_memo_r

  ! matrix structure !{{{1
  ! for horizontal (u,v) and tracer (theta,s)
  type :: type_mat
    real (kind=wp) :: x(2)
  end type type_mat

  ! time structure !{{{1
  ! yyyy-mm-dd hh:mm:ss
  type :: type_date
    character (len=19) :: d
  end type type_date

  type :: type_time
    integer :: y, m, d, h, mi, s
    integer :: mm ! maxima day of the current month
  end type type_time

  ! stagger grid  !{{{1
  type :: type_stg
    integer :: n ! grid number
    integer :: i, j ! x/y position of (i,j) grid
    ! rescaled Lame coefficients, h1 / a
    !   rh2 and h3 is 1 for current version of pcom
    real (kind=wp), allocatable :: rh1(:,:)
    real (kind=wp), allocatable :: tn(:,:) ! tan( lat )
    ! coordinate increaments
    type (type_mat), allocatable :: dx(:,:)
    ! north-south / east-west / diagonal neighbors
    type (type_stg), pointer :: ns, ew, di
  end type type_stg

  type :: type_vstg
    integer :: n
    integer :: k
    ! m, grid thickness
    real (kind=wp), allocatable :: dz(:)
    ! pa, vertical coordinate increments
    real (kind=wp), allocatable :: dp(:)
    ! m, geometry height
    real (kind=wp), allocatable :: z(:)
    ! Pa, pressure, vertical coordinate
    real (kind=wp), allocatable :: p(:)
    type (type_vstg), pointer :: ud ! up/down neighbor grid
  end type type_vstg

  type :: type_stg3d
    integer, allocatable :: msk(:,:,:)
    integer, allocatable :: lev(:,:)
    ! initial geopotential height at sea bottom
    real (kind=wp), allocatable :: phib(:,:)
    ! initial pressure at sea bottom (exclude atmopheric pressure)
    real (kind=wp), allocatable :: pb(:,:)
    type (type_stg3d), pointer :: ns, ew, di, ud
    type (type_stg),  pointer :: hg
    type (type_vstg), pointer :: vg
  end type type_stg3d

  ! grid variables !{{{1

  type type_gvar_r2d
    real (kind=wp), allocatable :: v(:,:) ! data values
    type (type_stg), pointer :: hg
  end type type_gvar_r2d

  type type_gvar_r3d
    real (kind=wp), allocatable :: v(:,:,:) ! data values
    type (type_stg3d), pointer :: g
  end type type_gvar_r3d

  type type_gvar_m2d
    type (type_gvar_r2d) :: x(2)
  end type type_gvar_m2d

  type type_gvar_m3d
    type (type_gvar_r3d) :: x(2)
  end type type_gvar_m3d

  ! compound variables !{{{1

  ! forcing
  type :: type_frc
    type (type_gvar_m3d) :: tau ! wind stress (taux, tauy)
    type (type_gvar_m3d) :: ts  ! climatic mean surface (temp., salinity)
    type (type_gvar_r3d) :: pa  ! sea level atmospheric pressure
    type (type_gvar_r3d) :: fw  ! fresh water flux (evaporation - precp.)
  end type type_frc

  ! surface boundary
  type :: type_bnd
    type (type_gvar_m2d) :: tau
    type (type_gvar_m2d) :: ts
    type (type_gvar_r2d) :: pa
    type (type_gvar_r2d) :: fw
  end type type_bnd

  ! barotropic integration at g3
  type :: type_bintg3 
    real (kind=wp), allocatable, dimension(:,:) :: xn, xs ! north/south of g2 at g3
    real (kind=wp), allocatable, dimension(:,:) :: ye, yw ! east / west of g4 at g3
  end type type_bintg3

  ! normalized sea bottom pressure when using stagger time scheme
  type :: type_pbt
    real (kind=wp), allocatable :: tp(:,:) ! previous time step values
    real (kind=wp), allocatable :: tc(:,:) ! current  time step values
    ! mean values in 1 baroclinic step ( = pbt_st(:,:,2) of v1.0 )
    real (kind=wp), allocatable :: bc(:,:) 
    ! mean values in 2 baroclinic step ( = pbt_st(:,:,3) of v1.0 )
    real (kind=wp), allocatable :: bc2(:,:)
    type (type_stg), pointer :: hg
  end type type_pbt

  ! accumulated variables !{{{1
  ! for time-average output
  type :: type_accu_gm3d
    type (type_gvar_m3d) :: var
    integer :: n, nrec
  end type type_accu_gm3d

  type :: type_accu_gr3d
    type (type_gvar_r3d) :: var
    integer :: n, nrec
  end type type_accu_gr3d

  type :: type_accu_gr2d
    type (type_gvar_r2d) :: var
    integer :: n, nrec
  end type type_accu_gr2d

  !interface !{{{1

  interface operator(+)
    module procedure time_plus_integer
  end interface

  interface operator(<)
    module procedure time_lt_string
  end interface

  interface type_check_date
    module procedure check_date_string
  end interface

  interface type_print
    module procedure print_type_time
  end interface

contains !{{{1

subroutine print_type_time(var) !{{{1
  ! print self defined type of type_time
  
  type (type_time) :: var

  ! print as the form of yyyy-mm-dd hh:mm:ss
  write(*,'(i0.4,a,i0.2,a,i0.2,x,i0.2,a,i0.2,a,i0.2)') &
    var % y, '-', var % m,  '-', var % d, &
    var % h, ':', var % mi, ':', var % s
  
end subroutine print_type_time

function time_lt_string (t, s) !{{{1
  ! user should check s first to see whether it is a proper date string

  type (type_time), intent(in) :: t
  character (len=*), intent(in) :: s
  logical :: time_lt_string

  character (len=14) :: s1, s2

  write( s1(1:4),  '(i4)') t % y
  write( s1(5:6),  '(i2)') t % m
  write( s1(7:8),  '(i2)') t % d
  write( s1(9:10), '(i2)') t % h
  write( s1(11:12),'(i2)') t % mi
  write( s1(13:14),'(i2)') t % s

  s2(1:4)   = s(1:4)
  s2(5:6)   = s(6:7)
  s2(7:8)   = s(9:10)
  s2(9:10)  = s(12:13)
  s2(11:12) = s(15:16)
  s2(13:44) = s(18:19)

  time_lt_string = llt(s1, s2)

end function time_lt_string

function type_time2sec ( t ) !{{{1
  ! how many seconds of the current time from 0001-01-01 00:00:00
  ! Note: this algorithm donot consider any history 'mistakes' for calendar
  !   the result of this routine have been compare to NCL 6.1.0's build-in
  !   function (see ou_string2time in ~/archive/ncl.ncl). But NCL thinks 
  !   there are 2*24*60*60 seconds for 0100-03-01 00:00:00 since 
  !   0100-02-28 00:00:00, obviously it treat 100 as a leap year (but also the
  !   built-in function isleapyear does not treat 100 as a leap year), so I
  !   think NCL has bugs in determine seconds from a date.
  !                  OU Niansen  2015-09-24

  type (type_time) :: t
  real (kind=wp) :: type_time2sec

  integer :: i

  type_time2sec = 0

  do i = 1, t%y - 1
    type_time2sec = type_time2sec + pro_days_year (i) * 24*60*60 
  end do

  do i = 1, t%m - 1
    type_time2sec = type_time2sec + pro_days_month (t%y, i) * 24*60*60 
  end do

  ! day start at 1
  type_time2sec = type_time2sec + (t%d - 1) * 24*60*60 

  ! hour/minute/second are start at 0 
  type_time2sec = type_time2sec + t%h * 60*60 + t%mi * 60 + t%s

end function type_time2sec

function type_str2sec ( str ) !{{{1
  ! how many seconds from string like "0001-01-01 00:00:00"
  character (len=*), intent(in) :: str
  integer*8 :: type_str2sec

  type (type_time) :: t

  t = type_str2time(str)
  type_str2sec = type_time2sec ( t )

end function type_str2sec

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

function type_str2time (str) !{{{1
  ! type_time plus an integer (in seconds)
  character (len=*), intent(in) :: str
  type (type_time) :: type_str2time

  integer :: leng, i, n
  character (len=80) :: str_num

  leng = len_trim(str)

  ! select digits character 

  ! will not successfully trim blanks without this line
  str_num = repeat('', len(str_num))

  n = 0
  do i = 1, leng
    if ( lge(str(i:i),'0') .and. lle(str(i:i),'9') ) then
      n = n + 1
      str_num(n:n) = str(i:i)
    end if
  end do

  if ( len_trim(str_num) /= 14 ) stop 'str in type_str2time should be in the form of yyyy-mm-dd hh:mm:ss'

  read(str_num(1:4), '(i4)') type_str2time % y
  if ( type_str2time % y <= 0 ) &
    stop 'year should be greater than 0 in the input string in function type_str2time'

  read(str_num(5:6), '(i2)') type_str2time % m
  type_str2time % mm = pro_days_month (type_str2time%y, type_str2time%m)
  if ( type_str2time % m <= 0 .or. type_str2time % m > type_str2time % mm ) then
    write (*, '(a,i2,a)') &
      'month should be between 1-', type_str2time % mm, &
      ' in the input string for the specified year in function type_str2time '
    stop
  end if

  read(str_num(7:8), '(i2)') type_str2time % d
  if ( type_str2time % d <= 0 .or. type_str2time % d > 31 ) &
    stop 'day should be between 1-31 in the input string in function type_str2time '

  read(str_num(9:10), '(i2)') type_str2time % h
  if ( type_str2time % h < 0 .or. type_str2time % h > 23 ) &
    stop 'hour should be between 0-23 in the input string in function type_str2time '

  read(str_num(11:12), '(i2)') type_str2time % mi
  if ( type_str2time % mi < 0 .or. type_str2time % mi > 59 ) &
    stop 'miniute should be between 0-59 in the input string in function type_str2time '

  read(str_num(13:14), '(i2)') type_str2time % s
  if ( type_str2time % s < 0 .or. type_str2time % s > 59 ) &
    stop 'second should be between 0-59 in the input string in function type_str2time '

end function type_str2time

function time_plus_integer (t, dt) !{{{1
  ! type_time plus an integer (in seconds)
  type (type_time), intent(in) :: t
  integer, intent(in) :: dt

  type (type_time) :: time_plus_integer

  if (dt > 24*60*60) stop 'dt should less than a day in function time_plus_integer'

  time_plus_integer = t
  time_plus_integer % s = time_plus_integer % s + dt

  if (time_plus_integer % s >= 60) then
    time_plus_integer % mi = time_plus_integer % mi + time_plus_integer % s / 60
    time_plus_integer % s  = mod(time_plus_integer % s, 60)
  end if

  if (time_plus_integer % mi >= 60) then
    time_plus_integer % h  = time_plus_integer % h + time_plus_integer % mi / 60
    time_plus_integer % mi = mod(time_plus_integer % mi, 60)
  end if

  if (time_plus_integer % h >= 24) then
    time_plus_integer % d = time_plus_integer % d + time_plus_integer % h / 24
    time_plus_integer % h = mod(time_plus_integer % h, 24)
  end if

  if (time_plus_integer % d > time_plus_integer % mm) then
    ! dt is no more than one day
    time_plus_integer % m = time_plus_integer % m + 1
    time_plus_integer % d = 1
  end if

  if (time_plus_integer % m > 12) then
    time_plus_integer % y = time_plus_integer % y + 1
    time_plus_integer % m = 1
  end if

  time_plus_integer % mm = pro_days_month (time_plus_integer%y, time_plus_integer%m)

  return

end function time_plus_integer

subroutine check_date_string (str) !{{{1
  ! check a string for whether it represents a proper date
  character (len=*), intent(in) :: str

  integer :: leng, i, n, &
    y, m, d, h, mi, s, mm
  character (len=80) :: str_num

  leng = len_trim(str)

  ! select digits character 

  ! will not successfully trim blanks without this line
  str_num = repeat('', len(str_num))

  n = 0
  do i = 1, leng
    if ( lge(str(i:i),'0') .and. lle(str(i:i),'9') ) then
      n = n + 1
      str_num(n:n) = str(i:i)
    end if
  end do

  if ( len_trim(str_num) /= 14 ) stop 'str in check_date_string should be in the form of yyyy-mm-dd hh:mm:ss'

  read(str_num(1:4), '(i4)') y
  if ( y <= 0 ) &
    stop 'year should be greater than 0 in the input string in subroutine check_date_string'

  read(str_num(5:6), '(i2)') m
  mm = pro_days_month (y, m)
  if ( m <= 0 .or. m > mm ) then
    write (*, '(a,i2,a)') &
      'month should be between 1-', mm, &
      ' in the input string for the specified year in subroutine check_date_string '
    stop
  end if

  read(str_num(7:8), '(i2)') d
  if ( d <= 0 .or. d > 31 ) &
    stop 'day should be between 1-31 in the input string in subroutine check_date_string '

  read(str_num(9:10), '(i2)') h
  if ( h < 0 .or. h > 23 ) &
    stop 'hour should be between 0-23 in the input string in subroutine check_date_string '

  read(str_num(11:12), '(i2)') mi
  if ( mi < 0 .or. mi > 59 ) &
    stop 'miniute should be between 0-59 in the input string in subroutine check_date_string '

  read(str_num(13:14), '(i2)') s
  if ( s < 0 .or. s > 59 ) &
    stop 'second should be between 0-59 in the input string in subroutine check_date_string '

end subroutine check_date_string

subroutine chk( ista ) !{{{1
  ! check state of allocate array 

  integer, intent(in) ::  ista

  if ( ista /= 0 ) then
    write(*,*) 'Allocate array failed. Stop'
    stop 2
  end if
end subroutine chk

end module mod_type !{{{1
!-------------------------------------------------------{{{1
! vim:fdm=marker:fdl=0:
! vim:foldtext=getline(v\:foldstart).'...'.(v\:foldend-v\:foldstart):
