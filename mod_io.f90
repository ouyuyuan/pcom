
! Description: basic NetCDF input/output interface
!
!      Author: OU Yuyuan <ouyuyuan@lasg.iap.ac.cn>
!     Created: 2015-03-06 10:38:13 BJT
! Last Change: 2017-12-31 16:10:58 BJT

module mod_io !{{{1 
!-------------------------------------------------------{{{1
!imported variables-------------------------------------{{{1
  use netcdf
  use mod_kind, only: wp

  use mod_type, only: tctr, type_var_info, type_time2str, &
    type_rst_info
  use mod_param, only: glo_ni, glo_nj, nk, &
    missing_float, missing_double, &
    nm, vars_info, rst_info

  use mod_arrays, only: glo_lon, glo_lat, z

  implicit none
  private

!public subroutines-------------------------------------{{{1

  public &
    io_get_dim_len,   &
    io_write,         &
    io_create_rst,    &
    io_read

!internal variables-------------------------------------{{{1

  integer :: ncid, varid, ndim1, ndim2, ndim3, &
             id_lon, id_lat, id_z, id_time

!interfaces---------------------------------------------{{{1
!-----------------------------------------------------------
  interface io_write
    module procedure write_r3d
    module procedure write_r2d
  end interface

  interface io_read
    module procedure read_r3d
    module procedure read_r2d
    module procedure read_r1d
    module procedure read_i3d
    module procedure read_i2d
    module procedure read_scalar
    module procedure read_glo_att
  end interface

contains !{{{1
!-------------------------------------------------------{{{1

subroutine io_create_rst ( rst_info ) !{{{1
  ! create file for restart run
  type (type_rst_info) :: rst_info

  integer :: idx

  rst_info%cdate = type_time2str(tctr%ct)
  rst_info%pdate = type_time2str(tctr%pt)

  rst_info%fname = trim(nm%od)//'restart_'//trim(rst_info%cdate)//'.nc'
  idx = len(trim(nm%od)//'restart_0000-00-00') + 1
  rst_info%fname(idx:idx) = '_'

  call create_nc (rst_info%fname, ncid)

  call def_dim_3d (id_lon, id_lat, id_z, id_time)

  call put_glo_att ('created', "by subroutine io_create_rst in module mod_io")
  call put_glo_att ('cdate', rst_info%cdate)
  call put_glo_att ('pdate', rst_info%pdate)

  ! def vars
  call def_var_r2d (rst_info%chc)
  call def_var_r2d (rst_info%chp)

  call def_var_r2d (rst_info%bnd_taux)
  call def_var_r2d (rst_info%bnd_tauy)
  call def_var_r2d (rst_info%bnd_t)
  call def_var_r2d (rst_info%bnd_s)
  call def_var_r2d (rst_info%bnd_pa)
  call def_var_r2d (rst_info%bnd_fw)

  call def_var_r2d (rst_info%buc)
  call def_var_r2d (rst_info%bvc)
  call def_var_r2d (rst_info%bup)
  call def_var_r2d (rst_info%bvp)

  call def_var_r3d (rst_info%tc)
  call def_var_r3d (rst_info%sc)
  call def_var_r3d (rst_info%tp)
  call def_var_r3d (rst_info%sp)

  call def_var_r3d (rst_info%prpt)
  call def_var_r3d (rst_info%prps)

  call def_var_r3d (rst_info%uc)
  call def_var_r3d (rst_info%vc)
  call def_var_r3d (rst_info%up)
  call def_var_r3d (rst_info%vp)

  call def_var_r3d (rst_info%auc)
  call def_var_r3d (rst_info%avc)
  call def_var_r3d (rst_info%aup)
  call def_var_r3d (rst_info%avp)
  call def_var_r3d (rst_info%aupp)
  call def_var_r3d (rst_info%avpp)

  call def_var_r3d ( rst_info%am )

  call check (nf90_enddef(ncid) )

  call put_dim_3d ( )

  call check (nf90_close(ncid) )

end subroutine io_create_rst

subroutine write_r3d(var_info, var) !{{{1

  type (type_var_info), intent(in) :: var_info
  real (kind=wp), intent(in) :: var(:,:,:)

  character (len = 80) :: ncname

  ! determine nc filename
  call get_ncname (var_info, ncname)

  ! create nc file for every output
  if (trim(var_info%vartype) .ne. 'restart') then
    call create_nc ( ncname, ncid )
    call def_dim_3d ( id_lon, id_lat, id_z, id_time )
    call put_glo_att ( 'created', "by subroutine write_r3d in module mod_io")
    call def_var_r3d ( var_info )
    call check (nf90_enddef(ncid) )
    call put_dim_3d ( )
  else
    call check (nf90_open(trim(ncname), nf90_write, ncid)  )
  end if

  call put_var_r3d (var_info%name, var)

  call increase_time ( )

  call check (nf90_close(ncid) )

end subroutine write_r3d

subroutine write_r2d(var_info, var) !{{{1

  type (type_var_info), intent(in) :: var_info
  real (kind=wp), intent(in) :: var(:,:)

  character (len = 80) :: ncname

  ! determine nc filename
  call get_ncname (var_info, ncname)

  ! create nc file for every output
  if (trim(var_info%vartype) .ne. 'restart') then
    call create_nc ( ncname, ncid )
    call def_dim_2d ( id_lon, id_lat, id_time )
    call put_glo_att ( 'created', "by subroutine write_r2d in module mod_io")
    call def_var_r2d ( var_info )
    call check (nf90_enddef(ncid) )
    call put_dim_2d ( )
  else
    call check (nf90_open(trim(ncname), nf90_write, ncid)  )
  end if

  call put_var_r2d (var_info%name, var)

  call increase_time ( )

  call check (nf90_close(ncid) )

end subroutine write_r2d

subroutine read_r3d(ncname, varname, var) !{{{1
  character (len=*), intent(in) :: ncname, varname
  real (kind=wp) :: var(:,:,:)

  call check (nf90_open(ncname, NF90_NOWRITE, ncid) )

  call check (nf90_inq_varid(ncid, varname, varid) )

  call check (nf90_get_var(ncid, varid, var) )

  call check (nf90_close(ncid) )

  write(*,'(a)') 'got '//trim(varname)// ' from '//trim(ncname)

end subroutine read_r3d 

subroutine read_r2d(ncname, varname, var) !{{{1
  character (len=*), intent(in) :: ncname, varname
  real (kind=wp) :: var(:,:)

  call check (nf90_open(ncname, NF90_NOWRITE, ncid) )

  call check (nf90_inq_varid(ncid, varname, varid) )

  call check (nf90_get_var(ncid, varid, var) )

  call check (nf90_close(ncid) )

  write(*,'(a)') 'got '//trim(varname)// ' from '//trim(ncname)

end subroutine read_r2d 

subroutine read_r1d(ncname, varname, var) !{{{1
  character (len=*), intent(in) :: ncname, varname
  real (kind=wp) :: var(:)

  call check (nf90_open(ncname, NF90_NOWRITE, ncid) )

  call check (nf90_inq_varid(ncid, varname, varid) )

  call check (nf90_get_var(ncid, varid, var) )

  call check (nf90_close(ncid) )

  write(*,'(a)') 'got '//trim(varname)// ' from '//trim(ncname)
end subroutine read_r1d 

subroutine read_i3d(ncname, varname, var) !{{{1
  character (len=*), intent(in) :: ncname, varname
  integer :: var(:,:,:)

  call check (nf90_open(ncname, NF90_NOWRITE, ncid) )

  call check (nf90_inq_varid(ncid, varname, varid) )

  call check (nf90_get_var(ncid, varid, var) )

  call check (nf90_close(ncid) )

  write(*,'(a)') 'got '//trim(varname)// ' from '//trim(ncname)

end subroutine read_i3d 

subroutine read_i2d(ncname, varname, var) !{{{1
  character (len=*), intent(in) :: ncname, varname
  integer :: var(:,:)

  call check (nf90_open(ncname, NF90_NOWRITE, ncid) )

  call check (nf90_inq_varid(ncid, varname, varid) )

  call check (nf90_get_var(ncid, varid, var) )

  call check (nf90_close(ncid) )

  write(*,'(a)') 'got '//trim(varname)// ' from '//trim(ncname)

end subroutine read_i2d 

subroutine read_scalar(ncname, varname, var) !{{{1
  character (len=*), intent(in) :: ncname, varname
  real (kind=wp) :: var

  call check (nf90_open(ncname, NF90_NOWRITE, ncid) )

  call check (nf90_inq_varid(ncid, varname, varid) )

  call check (nf90_get_var(ncid, varid, var) )

  call check (nf90_close(ncid) )
end subroutine read_scalar

subroutine read_glo_att(ncname, att, var) !{{{1
  ! read a global attribute
  character (len=*), intent(in) :: ncname, att
  character (len=*) :: var

  call check (nf90_open(ncname, NF90_NOWRITE, ncid) )

  call check (nf90_get_att(ncid, NF90_GLOBAL, att, var) )

  call check (nf90_close(ncid) )
end subroutine read_glo_att

subroutine io_get_dim_len(ncname, dimname, var) !{{{1
! get the length of a dimension with 'dimname' from a 
! NetCDF file with 'ncname'
  character (len = *), intent(in) :: ncname, dimname
  integer, intent(inout) :: var
  integer :: dimid

  call check ( nf90_open(ncname, NF90_NOWRITE, ncid) )

  call check ( nf90_inq_dimid(ncid, dimname, dimid) )

  call check ( nf90_inquire_dimension(ncid, dimid, len = var) )

  call check ( nf90_close(ncid) )

end subroutine io_get_dim_len

subroutine get_ncname (var_info, ncname) !{{{1
! create nc file name per output time
  type (type_var_info), intent (in) :: var_info
  character (len=*) :: ncname

  character (len = len('0000-00-00 00:00:00')) :: str
  integer :: idx

  str = trim(type_time2str(tctr%ct))
  idx = len('0000-00-00') + 1
  str (idx:idx) = '_'
  idx = len('0000-00-00_00:00')

  if ( trim ( var_info%vartype ) .eq. 'restart') then
    ncname = rst_info%fname
  else if ( trim ( var_info%vartype ) .eq. 'debug') then
    ncname = trim(nm%od)//"debug/"//trim(var_info%name)//'_'//str(1:idx)//'.nc'
  else
    ncname = trim(nm%od)//trim(var_info%name)//'_'//str(1:idx)//'.nc'
  end if

end subroutine get_ncname

subroutine def_dim_2d (id_lon, id_lat, id_time) !{{{1
! dimension definition
  integer :: id_lon, id_lat, id_time

  call check ( nf90_def_dim (ncid, 'lon', glo_ni, id_lon) )
  call check ( nf90_def_dim (ncid, 'lat', glo_nj, id_lat) )
  call check ( nf90_def_dim (ncid, 'time', NF90_UNLIMITED, id_time) )

  call check ( nf90_def_var (ncid, "lon", nf90_float, &
    id_lon, varid) )
  call check ( nf90_put_att (ncid, varid, 'long_name', &
    'longitude') )
  call check ( nf90_put_att (ncid, varid, 'units', &
    'degree_east') )

  call check ( nf90_def_var (ncid, "lat", nf90_float, &
    id_lat, varid) )
  call check ( nf90_put_att (ncid, varid, 'long_name', &
    'latitude') )
  call check ( nf90_put_att (ncid, varid, 'units', &
    'degree_north') )

  call check ( nf90_def_var (ncid, "time", nf90_float, &
    id_time, varid) )
  call check ( nf90_put_att (ncid, varid, 'units', &
    'hours since 0001-01-01 00:00:00') )

end subroutine def_dim_2d

subroutine def_dim_3d (id_lon, id_lat, id_z, id_time) !{{{1
! dimension definition
  integer :: id_lon, id_lat, id_z, id_time

  call def_dim_2d (id_lon, id_lat, id_time)

  call check ( nf90_def_dim (ncid, 'z',   nk, id_z) )

  call check ( nf90_def_var (ncid, "z", nf90_float, &
    id_z, varid) )
  call check ( nf90_put_att (ncid, varid, 'long_name', &
    'depth') )
  call check ( nf90_put_att (ncid, varid, 'units', &
    'm') )

end subroutine def_dim_3d

subroutine def_var_r2d (var_info) !{{{1
  ! define a 2d variable
  type (type_var_info), intent(in) :: var_info

  integer :: dimids(3)

  dimids = (/id_lon, id_lat, id_time/)

  if ( var_info%vartype == 'mean' ) then
    call check ( nf90_def_var (ncid, var_info%name, &
      nf90_float, dimids, varid) )
    call check ( nf90_put_att (ncid, varid, '_FillValue', &
      missing_float) )
  else
    call check ( nf90_def_var (ncid, var_info%name, &
      nf90_double, dimids, varid) )
    call check ( nf90_put_att (ncid, varid, '_FillValue', &
      missing_double) )
  end if

  call check ( nf90_put_att (ncid, varid, 'long_name', &
    var_info%longname) )
  call check ( nf90_put_att (ncid, varid, 'units', &
    var_info%units) )

end subroutine def_var_r2d

subroutine def_var_r3d (var_info) !{{{1
  ! define a 2d variable
  type (type_var_info), intent(in) :: var_info

  integer :: dimids(4)

  dimids = (/id_lon, id_lat, id_z, id_time/)

  if ( var_info%vartype == 'mean' ) then
    call check ( nf90_def_var (ncid, var_info%name, &
      nf90_float, dimids, varid) )
    call check ( nf90_put_att (ncid, varid, '_FillValue', &
      missing_float) )
  else
    call check ( nf90_def_var (ncid, var_info%name, &
      nf90_double, dimids, varid) )
    call check ( nf90_put_att (ncid, varid, '_FillValue', &
      missing_double) )
  end if

  call check ( nf90_put_att (ncid, varid, 'long_name', &
    var_info%longname) )
  call check ( nf90_put_att (ncid, varid, 'units', &
    var_info%units) )

end subroutine def_var_r3d

subroutine put_var_r3d (varname, var) !{{{1
  ! define a 2d variable
  character (len=*), intent(in) :: varname
  real (kind=wp), intent(in) :: var(:,:,:)

  integer :: stt(4), cnt(4)

  ndim1 = size(var, 1)
  ndim2 = size(var, 2)
  ndim3 = size(var, 3)

  stt  = (/1, 1, 1, 1/)
  cnt  = (/ndim1, ndim2, ndim3, 1/)

  call check (nf90_inq_varid(ncid, trim(varname), varid) )
  call check (nf90_put_var(ncid, varid, var, start=stt, count=cnt) )

end subroutine put_var_r3d

subroutine put_var_r2d (varname, var) !{{{1
  ! define a 2d variable
  character (len=*), intent(in) :: varname
  real (kind=wp), intent(in) :: var(:,:)

  integer :: stt(3), cnt(3)

  ndim1 = size(var, 1)
  ndim2 = size(var, 2)

  stt  = (/1, 1, 1/)
  cnt  = (/ndim1, ndim2, 1/)

  call check (nf90_inq_varid(ncid, trim(varname), varid) )
  call check (nf90_put_var(ncid, varid, var, start=stt, count=cnt) )

end subroutine put_var_r2d

subroutine increase_time () !{{{1
  ! increase time value
  real (kind=wp) :: hour

  hour = real(tctr%i) * real(nm%bc) / 3600.0
  call check (nf90_inq_varid(ncid, 'time', varid) )
  call check (nf90_put_var(ncid, varid, hour, start=(/1/)) )

end subroutine increase_time

subroutine put_glo_att (attname, str) !{{{1
! write global attribute of nc file
  character (len=*), intent(in) :: attname, str

  call check ( nf90_put_att (ncid, NF90_GLOBAL, attname, trim(str) ) )

end subroutine put_glo_att

subroutine create_nc (ncname, ncid) !{{{1
  ! create nc file
  character (len=*), intent(in) :: ncname
  integer :: ncid

  call check ( nf90_create (trim(ncname), NF90_CLOBBER, ncid)  )

end subroutine create_nc

subroutine put_dim_2d ( ) !{{{1
  ! write dimension coordinate

  call check (nf90_inq_varid(ncid, 'lon', varid) )
  call check (nf90_put_var(ncid, varid, glo_lon) )

  call check (nf90_inq_varid(ncid, 'lat', varid) )
  call check (nf90_put_var(ncid, varid, glo_lat) )

end subroutine put_dim_2d

subroutine put_dim_3d ( ) !{{{1
  ! write dimension coordinates

  call put_dim_2d ( )

  call check (nf90_inq_varid(ncid, 'z', varid) )
  call check (nf90_put_var(ncid, varid, z) )

end subroutine put_dim_3d

subroutine check(status) !{{{1
!-----------------------------------------------------------
! check netcdf call
!-----------------------------------------------------------
  integer, intent (in) :: status

  if(status /= nf90_noerr) then 
    print *, trim(nf90_strerror(status))
    stop
  end if
end subroutine check  

end module mod_io !{{{1
!-------------------------------------------------------{{{1
! vim:fdm=marker:fdl=0:
! vim:foldtext=getline(v\:foldstart).'...'.(v\:foldend-v\:foldstart):
