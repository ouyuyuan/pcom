
! Description: debug module 
!
!      Author: OU Yuyuan <ouyuyuan@lasg.iap.ac.cn>
!     Created: 2017-12-05 09:21:11 BJT
! Last Change: 2017-12-29 10:50:41 BJT

module mod_debug

  ! imported variables !{{{1
  use mod_arrays, only: eqch, eqts, equv
  
  use mod_kind, only: wp

  use mod_param, only: &
    ni, nj, nk, myid, mid, nm, rst_info

  use mod_type, only: &
    tctr, type_time2str, &
    type_eq_uv, type_eq_uvb, &
    type_var_info

  use mod_mympi, only: mympi_input, mympi_output

  implicit none
  private

  ! exported variables !{{{1
  public  debug_var

  interface debug_var
    module procedure debug
    module procedure debug_equv
    module procedure debug_equvb
    module procedure debug_r3d
    module procedure debug_r2d
  end interface

contains !{{{1

subroutine debug ( opt ) !{{{1
  ! debug variable according to option "opt"
  character (len=*), intent(in) :: opt

  if (opt == "equal to restart?") then
    call check_equal_to_rst ( )
  else
    print *, "unreconized option in debug_var"
    stop
  end if

end subroutine debug

subroutine debug_equv ( opt, var, varname ) !{{{1
  ! debug variable according to option "opt"
  character (len=*), intent(in) :: opt
  type (type_eq_uv), intent(in) :: var
  character (len=*), intent(in) :: varname

  if (opt == "nan?") then
    call check_nan_equv ( var, varname )
  else
    print *, "unreconized option in debug_var"
    stop
  end if
end subroutine debug_equv

subroutine debug_equvb ( opt, var, varname ) !{{{1
  ! debug variable according to option "opt"
  character (len=*), intent(in) :: opt
  type (type_eq_uvb), intent(in) :: var
  character (len=*), intent(in) :: varname

  if (opt == "nan?") then
    call check_nan_equvb ( var, varname )
  else
    print *, "unreconized option in debug_var"
    stop
  end if
end subroutine debug_equvb

subroutine debug_r3d ( opt, var, varname ) !{{{1
  ! debug variable according to option "opt"
  character (len=*), intent(in) :: opt
  real (kind=wp), dimension(:,:,:), intent(in) :: var
  character (len=*), intent(in) :: varname

  type (type_var_info) :: var_info

  var_info%vartype  = 'debug'
  var_info%units    = ''
  var_info%longname = ''
  var_info%name     = trim(varname)//'_norst'
  if (nm%rst == 1) var_info%name = trim(varname)//'_rst'

  if (opt == "nan?") then
    call check_nan_r3d ( var, varname )
  else if (opt == "zero?") then
    call check_zero_r3d ( var, varname )
  else if (opt == "txt") then
    if (myid == mid) call output_txt_r3d ( var, varname )
  else if (opt == "binary") then
    if ( myid == mid ) call output_binary_r3d ( var, varname )
  else if (opt == "equal to no restart?") then
    call mympi_output ( var_info, var )
    if ( nm%rst == 1 ) call check_equal_to_norst_r3d ( var, varname )
  else
    print *, "unreconized option in debug_var"
    stop
  end if

end subroutine debug_r3d
subroutine debug_r2d ( opt, var, varname ) !{{{1
  ! debug variable according to option "opt"
  character (len=*), intent(in) :: opt
  real (kind=wp), dimension(:,:), intent(in) :: var
  character (len=*), intent(in) :: varname

  type (type_var_info) :: var_info

  var_info%vartype  = 'debug'
  var_info%units    = ''
  var_info%longname = ''
  var_info%name     = trim(varname)//'_norst'
  if (nm%rst == 1) var_info%name = trim(varname)//'_rst'

  if (opt == "nan?") then
    call check_nan_r2d ( var, varname )
  else if (opt == "zero?") then
    call check_zero_r2d ( var, varname )
  else if (opt == "txt") then
    if (myid == mid) call output_txt_r2d (var, varname)
  else if (opt == "binary") then
    if (myid == mid) call output_binary_r2d (var, varname)
  else if (opt == "equal to restart?") then
    call check_equal_to_rst_r2d ( var, varname )
  else if (opt == "equal to no restart?") then
    call mympi_output ( var_info, var )
    if ( nm%rst == 1 ) call check_equal_to_norst_r2d ( var, varname )
  end if

end subroutine debug_r2d
subroutine check_equal_to_rst ( ) !{{{1
  ! check whether the variables equal to the restart file
  call check_equal_to_rst_r2d ( eqch%chc, "chc" )
  call check_equal_to_rst_r2d ( eqch%chp, "chp" )

  call check_equal_to_rst_r3d ( eqts%tc, "tc" )
  call check_equal_to_rst_r3d ( eqts%sc, "sc" )
  call check_equal_to_rst_r3d ( eqts%tp, "tp" )
  call check_equal_to_rst_r3d ( eqts%sp, "sp" )

  call check_equal_to_rst_r3d ( equv%uc, "uc" )
  call check_equal_to_rst_r3d ( equv%vc, "vc" )
  call check_equal_to_rst_r3d ( equv%up, "up" )
  call check_equal_to_rst_r3d ( equv%vp, "vp" )

  call check_equal_to_rst_r3d ( equv%auc, "auc" )
  call check_equal_to_rst_r3d ( equv%avc, "avc" )
  call check_equal_to_rst_r3d ( equv%aup, "aup" )
  call check_equal_to_rst_r3d ( equv%avp, "avp" )
  call check_equal_to_rst_r3d ( equv%aupp, "aupp" )
  call check_equal_to_rst_r3d ( equv%avpp, "avpp" )
end subroutine check_equal_to_rst

subroutine check_equal_to_norst_r3d ( var, varname ) !{{{1
  ! check whether inequal to norst happens
  real (kind=wp), dimension(:,:,:), intent(in) :: var
  character ( len = * ) :: varname

  real (kind=wp), allocatable, dimension(:,:,:) :: norst_var
  integer, allocatable, dimension(:,:,:) :: inequal
  integer :: d1, d2, d3, is, idx
  character (len = len('0000-00-00 00:00:00')) :: str
  character ( len = 80 ) :: fname

  fname = trim(nm%od)//"debug/"//trim(varname)//'_norst.nc'

  str = trim(type_time2str(tctr%ct))
  idx = len('0000-00-00') + 1
  str (idx:idx) = '_'
  idx = len('0000-00-00_00:00')
  fname = trim(nm%od)//"debug/"//trim(varname)//'_norst'//'_'//str(1:idx)//'.nc'

  d1 = size ( var, 1 )
  d2 = size ( var, 2 )
  d3 = size ( var, 3 )

  allocate ( norst_var(d1,d2,d3), stat = is ); call chk(is)
  allocate ( inequal(d1,d2,d3), stat = is ); call chk(is)
  norst_var = var
  inequal = 0

  call mympi_input ( fname, trim(varname)//'_norst', norst_var )

  where ( var /= norst_var )
    inequal = 1
  end where 

  call print_inequal_to_norestart (sum(inequal), varname)

  deallocate (inequal)
  deallocate (norst_var)

end subroutine check_equal_to_norst_r3d

subroutine check_equal_to_norst_r2d ( var, varname ) !{{{1
  ! check whether inequal to norst happens
  real (kind=wp), dimension(:,:), intent(in) :: var
  character ( len = * ) :: varname

  real (kind=wp), allocatable, dimension(:,:) :: norst_var
  integer, allocatable, dimension(:,:) :: inequal
  integer :: d1, d2, is, idx
  character (len = len('0000-00-00 00:00:00')) :: str
  character ( len = 80 ) :: fname

  fname = trim(nm%od)//"debug/"//trim(varname)//'_norst.nc'

  str = trim(type_time2str(tctr%ct))
  idx = len('0000-00-00') + 1
  str (idx:idx) = '_'
  idx = len('0000-00-00_00:00')
  fname = trim(nm%od)//"debug/"//trim(varname)//'_norst'//'_'//str(1:idx)//'.nc'

  d1 = size ( var, 1 )
  d2 = size ( var, 2 )

  allocate ( norst_var(d1,d2), stat = is ); call chk(is)
  allocate ( inequal(d1,d2), stat = is ); call chk(is)
  norst_var = var
  inequal = 0

  call mympi_input ( fname, trim(varname)//'_norst', norst_var )

  where ( var /= norst_var )
    inequal = 1
  end where 

  call print_inequal_to_norestart (sum(inequal), varname)

  deallocate (inequal)
  deallocate (norst_var)

end subroutine check_equal_to_norst_r2d

subroutine check_equal_to_rst_r3d (var, varname) !{{{1
  ! check whether equal_to_rst happens
  real (kind=wp), dimension(:,:,:), intent(in) :: var
  character (len=*), intent(in) :: varname

  real (kind=wp), allocatable, dimension(:,:,:) :: rst_var
  integer, allocatable, dimension(:,:,:) :: inequal
  integer :: d1, d2, d3, is

  d1 = size(var,1)
  d2 = size(var,2)
  d3 = size(var,3)

  allocate ( rst_var(d1,d2,d3), stat = is ); call chk(is)
  allocate ( inequal(d1,d2,d3), stat = is ); call chk(is)
  rst_var = var
  inequal = 0

  call mympi_input (rst_info%fname, varname, rst_var)

!  where (abs(var-rst_var) .gt. 1.0e-12)
  where ( var /= rst_var )
    inequal = 1
  end where 

  call print_inequal_to_restart (sum(inequal), varname)

  deallocate (inequal)
  deallocate (rst_var)

end subroutine check_equal_to_rst_r3d

subroutine check_equal_to_rst_r2d (var, varname) !{{{1
  ! check whether equal_to_rst happens
  real (kind=wp), dimension(:,:), intent(in) :: var
  character (len=*), intent(in) :: varname

  real (kind=wp), allocatable, dimension(:,:) :: rst_var
  integer, allocatable, dimension(:,:) :: inequal
  integer :: d1, d2, is

  d1 = size(var,1)
  d2 = size(var,2)

  allocate ( rst_var(d1,d2), stat = is ); call chk(is)
  allocate ( inequal(d1,d2), stat = is ); call chk(is)
  rst_var = var
  inequal = 0

  call mympi_input (rst_info%fname, varname, rst_var)

  where ( var /= rst_var )
    inequal = 1
  end where 

  call print_inequal_to_restart (sum(inequal), varname)

  deallocate (inequal)
  deallocate (rst_var)

end subroutine check_equal_to_rst_r2d

subroutine check_zero_r3d (var, varname) !{{{1
  ! check whether zero happens
  real (kind=wp), dimension(:,:,:), intent(in) :: var
  character (len=*), intent(in) :: varname

  integer, allocatable, dimension(:,:,:) :: zeros
  integer :: d1, d2, d3, is

  d1 = size(var,1)
  d2 = size(var,2)
  d3 = size(var,3)

  allocate ( zeros(d1,d2,d3), stat = is ); call chk(is)
  zeros = 0

  where (abs(var) .le. 1.0e-12)
    zeros = 1
  end where 

  call print_zeros (sum(zeros), varname)

  deallocate (zeros)

end subroutine check_zero_r3d

subroutine check_zero_r2d (var, varname) !{{{1
  ! check whether zero happens
  real (kind=wp), dimension(:,:), intent(in) :: var
  character (len=*), intent(in) :: varname

  integer, allocatable, dimension(:,:) :: zeros
  integer :: d1, d2, is

  d1 = size(var,1)
  d2 = size(var,2)

  allocate ( zeros(d1,d2), stat = is ); call chk(is)
  zeros = 0

  where (abs(var) .le. 1.0e-12)
    zeros = 1
  end where 

  call print_zeros (sum(zeros), varname)

  deallocate (zeros)

end subroutine check_zero_r2d

subroutine check_nan_equv (var, varname) !{{{1
  ! check whether NaN happens
  type (type_eq_uv) :: var
  character (len=*), intent(in) :: varname

  call check_nan_r3d (var%uc, "uc of "//varname)
  call check_nan_r3d (var%vc, "vc of "//varname)

  call check_nan_r3d (var%auc, "auc of "//varname)
  call check_nan_r3d (var%avc, "avc of "//varname)

  call check_nan_r3d (var%fx, "fx of "//varname)
  call check_nan_r3d (var%fy, "fy of "//varname)

  call check_nan_r2d (var%pax, "pax of "//varname)
  call check_nan_r2d (var%pay, "pay of "//varname)

end subroutine check_nan_equv

subroutine check_nan_equvb (var, varname) !{{{1
  ! check whether NaN happens
  type (type_eq_uvb) :: var
  character (len=*), intent(in) :: varname

  call check_nan_r2d (var%uc, "uc of "//varname)
  call check_nan_r2d (var%vc, "vc of "//varname)

  call check_nan_r2d (var%ut, "ut of "//varname)
  call check_nan_r2d (var%vt, "vt of "//varname)

end subroutine check_nan_equvb

subroutine check_nan_r3d (var, varname) !{{{1
  ! check whether NaN happens
  real (kind=wp), dimension(:,:,:), intent(in) :: var
  character (len=*), intent(in) :: varname

  integer, dimension(ni,nj,nk) :: nan

  nan = 0

  where (var .ne. var)
    nan = 1
  end where 

  call print_nans (sum(nan), varname)

end subroutine check_nan_r3d

subroutine check_nan_r2d (var, varname) !{{{1
  ! check whether NaN happens
  real (kind=wp), dimension(:,:), intent(in) :: var
  character (len=*), intent(in) :: varname

  integer, dimension(ni,nj) :: nan

  nan = 0

  where (var .ne. var)
    nan = 1
  end where 

  call print_nans (sum(nan), varname)

end subroutine check_nan_r2d

subroutine output_txt_r3d (var, varname) !{{{1
! output data to txt format for debug
  real (kind=wp), dimension(:,:,:), intent(in) :: var
  character (len=*), intent(in) :: varname

  integer, parameter :: out_unit = 10
  character (len=80) :: fname

  if (nm%rst == 1) then
    fname = trim(nm%od)//"debug/"//trim(varname)//'_rst.txt'
  else
    fname = trim(nm%od)//"debug/"//trim(varname)//'_norst.txt'
  end if

  open (unit=out_unit, file=trim(fname), action="write")
  write(out_unit, "(10f10.3)") var
  close(out_unit)

end subroutine output_txt_r3d

subroutine output_txt_r2d (var, varname) !{{{1
! output data to txt format for debug
  real (kind=wp), dimension(:,:), intent(in) :: var
  character (len=*), intent(in) :: varname

  integer, parameter :: out_unit = 10
  character (len=80) :: fname

  if (nm%rst == 1) then
    fname = trim(nm%od)//"debug/"//trim(varname)//'_rst.txt'
  else
    fname = trim(nm%od)//"debug/"//trim(varname)//'_norst.txt'
  end if

  open (unit=out_unit, file=trim(fname), action="write")
  write(out_unit, "(10f10.3)") var
  close(out_unit)

end subroutine output_txt_r2d

subroutine output_binary_r3d (var, varname) !{{{1
! output data to binary format for debug
  real (kind=wp), dimension(:,:,:), intent(in) :: var
  character (len=*), intent(in) :: varname

  integer, parameter :: out_unit = 10
  character (len=80) :: fname

  if (nm%rst == 1) then
    fname = trim(nm%od)//"debug/"//trim(varname)//'_rst.bin'
  else
    fname = trim(nm%od)//"debug/"//trim(varname)//'_norst.bin'
  end if

  open ( unit=out_unit, file=trim(fname), form="unformatted" )
  write ( out_unit ) var
  close ( out_unit )

end subroutine output_binary_r3d

subroutine output_binary_r2d (var, varname) !{{{1
! output data to binary format for debug
  real (kind=wp), dimension(:,:), intent(in) :: var
  character (len=*), intent(in) :: varname

  integer, parameter :: out_unit = 10
  character (len=80) :: fname

  if (nm%rst == 1) then
    fname = trim(nm%od)//"debug/"//trim(varname)//'_rst.bin'
  else
    fname = trim(nm%od)//"debug/"//trim(varname)//'_norst.bin'
  end if

  open ( unit=out_unit, file=trim(fname), access="direct" )
  write ( out_unit ) var
  close ( out_unit )

end subroutine output_binary_r2d

subroutine print_inequal_to_restart (inequalsum, varname) !{{{1
  ! check the numbers of non-defined values
  integer, intent(in) :: inequalsum
  character (len=*), intent(in) :: varname

  if (inequalsum .gt. 0) then
    write(*,'(a,i8,a,i2,a,i2)') "DEBUG: different to restart file: "//varname//&
      ", inequal grids:", inequalsum, ", id:", myid, &
      ", "//trim(type_time2str(tctr%ct))//", i:",tctr%i
    stop "Find inequals when debug."
  else
    write(*,'(a,i2,a,i2)') "DEBUG: Looks OK: "//varname//&
    ", id:", myid, ", "//trim(type_time2str(tctr%ct))//", i:",tctr%i
  end if

end subroutine print_inequal_to_restart

subroutine print_inequal_to_norestart (inequalsum, varname) !{{{1
  ! check the numbers of non-defined values
  integer, intent(in) :: inequalsum
  character (len=*), intent(in) :: varname

  if (inequalsum .gt. 0) then
    write(*,'(a,i8,a,i2,a,i2)') "DEBUG: different to no restart file: "//varname//&
      ", inequal grids:", inequalsum, ", id:", myid, &
      ", "//trim(type_time2str(tctr%ct))//", i:",tctr%i
    stop "Find inequals when debug."
  else
    write(*,'(a,i2,a,i2)') "DEBUG: Looks OK: "//varname//&
    ", id:", myid, ", "//trim(type_time2str(tctr%ct))//", i:",tctr%i
  end if

end subroutine print_inequal_to_norestart

subroutine print_zeros (zerosum, varname) !{{{1
  ! check the numbers of non-defined values
  integer, intent(in) :: zerosum
  character (len=*), intent(in) :: varname

  if (zerosum .gt. 0) then
    write(*,'(a,i6,a,i2,a,i2)') "DEBUG: divided by zero: "//varname//&
      ", zero grids:", zerosum, ", id:", myid, &
      ", "//trim(type_time2str(tctr%ct))//", i:",tctr%i
    stop "Find divided by ZERO when debug."
  else
    write(*,'(a,i2,a,i2)') "DEBUG: Looks OK: "//varname//&
    ", id:", myid, ", "//trim(type_time2str(tctr%ct))//", i:",tctr%i
  end if

end subroutine print_zeros

subroutine print_nans (nansum, varname) !{{{1
  ! check the numbers of non-defined values
  integer, intent(in) :: nansum
  character (len=*), intent(in) :: varname

  if (nansum .gt. 0) then
    write(*,'(a,i6,a,i2,a,i2)') "DEBUG: Something wrong: "//varname//&
      ", undefined grids:", nansum, ", id:", myid, &
      ", "//trim(type_time2str(tctr%ct))//", i:",tctr%i
    stop "Find NaN value when debug."
  else
    write(*,'(a,i2,a,i2)') "DEBUG: Looks OK: "//varname//&
    ", id:", myid, ", "//trim(type_time2str(tctr%ct))//", i:",tctr%i
  end if

end subroutine print_nans

subroutine chk( ista ) !{{{1
  ! check state of allocate array 

  integer, intent(in) ::  ista

  if ( ista /= 0 ) then
    write(*,*) 'Allocate array failed. Stop'
    stop 2
  end if
end subroutine chk

end module mod_debug!{{{1
!-------------------------------------------------------{{{1
! vim:fdm=marker:fdl=0:
! vim:foldtext=getline(v\:foldstart).'...'.(v\:foldend-v\:foldstart):
