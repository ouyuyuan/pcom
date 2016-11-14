      !--------------------------------------------------------
      subroutine netcdf_read_cor(ncname,lat,lon,z,nx,ny,nz)
      !--------------------------------------------------------
      use netcdf
      implicit none
      integer nx,ny,nz
      character (len = *) :: ncname
      character (len = *), parameter :: lonname = "lon"
      character (len = *), parameter :: latname = "lat"
      character (len = *), parameter :: zname = "z"
      integer ncid
      integer varid,lon_varid, lat_varid, z_varid
      real lat(ny), lon(nx), z(nz)

      call check( nf90_open(ncname, nf90_nowrite, ncid) )

      call check( nf90_inq_varid(ncid, latname, lat_varid) )
      call check( nf90_inq_varid(ncid, lonname, lon_varid) )
      call check( nf90_inq_varid(ncid, zname,   z_varid  ) )

      call check( nf90_get_var(ncid, lat_varid, lat) )
      call check( nf90_get_var(ncid, lon_varid, lon) )
      call check( nf90_get_var(ncid, z_varid,   z  ) )

      call check( nf90_close(ncid) )
      print *,"*** SUCCESS Get coordinate information from ", ncname, "!"

      return
      end subroutine netcdf_read_cor

      !--------------------------------------------------------
      subroutine netcdf_read_var1d(ncname,varname,var,nx,missvalue)
      !--------------------------------------------------------
      use netcdf
      implicit none
      integer nx
      character (len = *), parameter :: missname = "missing_value"
      character (len = *) :: varname
      character (len = *) :: ncname
      integer ncid
      integer varid
      real var(nx),missvalue

      call check( nf90_open(ncname, nf90_nowrite, ncid) )

      call check( nf90_inq_varid(ncid, varname, varid) )
      call check( nf90_get_att(ncid, varid, missname, missvalue))
      call check( nf90_get_var(ncid, varid, var) )
      call check( nf90_close(ncid) )
      print *,"*** SUCCESS Get 1dvar ",varname," from file ", ncname, "!"

      return
      end subroutine netcdf_read_var1d

      !--------------------------------------------------------
      subroutine netcdf_read_var2d(ncname,varname,var,nx,ny,missvalue)
      !--------------------------------------------------------
      use netcdf
      implicit none
      integer nx,ny
      character (len = *), parameter :: missname = "missing_value"
      character (len = *) :: varname
      character (len = *) :: ncname
      integer ncid
      integer varid
      real var(nx,ny),missvalue

      call check( nf90_open(ncname, nf90_nowrite, ncid) )

      call check( nf90_inq_varid(ncid, varname, varid) )
      call check( nf90_get_att(ncid, varid, missname, missvalue))
      call check( nf90_get_var(ncid, varid, var) )
      call check( nf90_close(ncid) )
      print *,"*** SUCCESS Get 2dvar ",varname," from file ", ncname, "!"

      return
      end subroutine netcdf_read_var2d

      !--------------------------------------------------------
      subroutine netcdf_re_var2d_s(ncname,varname,var,nx,ny,missvalue,ts)
      !--------------------------------------------------------
      use netcdf
      implicit none
      integer nx,ny,ts
      character (len = *), parameter :: missname = "missing_value"
      character (len = *) :: varname
      character (len = *) :: ncname
      integer ncid,varid
      integer start(3),count(3)
      real var(nx,ny),missvalue

      count = (/nx,ny,1/)
      start = (/1,1,1/)
      start(3) = ts
      
      call check( nf90_open(ncname, nf90_nowrite, ncid) )

      call check( nf90_inq_varid(ncid, varname, varid) )
      call check( nf90_get_att(ncid, varid, missname, missvalue))
      call check( nf90_get_var(ncid, varid, var, start, count) )
      call check( nf90_close(ncid) )
      print *,"*** SUCCESS Get 2dvar slice ",varname," from file ", ncname, " month ",ts

      return
      end subroutine netcdf_re_var2d_s

      !--------------------------------------------------------
      subroutine netcdf_read_var(ncname,varname,var,nx,ny,nz,missvalue)
      !--------------------------------------------------------
      use netcdf
      implicit none
      integer nx,ny,nz
      character (len = *), parameter :: missname = "missing_value"
      character (len = *) :: varname
      character (len = *) :: ncname
      integer ncid
      integer varid
      real var(nx,ny,nz),missvalue

      call check( nf90_open(ncname, nf90_nowrite, ncid) )

      call check( nf90_inq_varid(ncid, varname, varid) )
      call check( nf90_get_att(ncid, varid, missname, missvalue))
      call check( nf90_get_var(ncid, varid, var) )
      call check( nf90_close(ncid) )
      print *,"*** SUCCESS Get 3dvar ",varname," from file ", ncname, "!"

      return
      end subroutine netcdf_read_var

      !--------------------------------------------------------
      subroutine netcdf_read_itn(var,nx,ny)
      !--------------------------------------------------------
      use netcdf
      implicit none
      integer nx,ny
      character (len = *), parameter :: varname = "itn"
      character (len = *), parameter :: ncname = "input/pcom_ini.nc"
      integer ncid
      integer varid
      integer var(nx,ny)

      call check( nf90_open(ncname, nf90_nowrite, ncid) )

      call check( nf90_inq_varid(ncid, varname, varid) )
      call check( nf90_get_var(ncid, varid, var) )
      call check( nf90_close(ncid) )
      print *,"*** SUCCESS Get var ",varname," from file ", ncname, "!"

      return
      end subroutine netcdf_read_itn


      !--------------------------------------------------------
      subroutine netcdf_write_cor(ncname,lat,lon,z,nx,ny,nz,month,t_units)
      !--------------------------------------------------------
      use netcdf
      implicit none
      character (len = *), parameter :: latname = "lat"
      character (len = *), parameter :: lonname = "lon"
      character (len = *), parameter :: zname = "z"
      character (len = *), parameter :: tname = "t"
      character (len = *), parameter :: units = "units"
      character (len = *), parameter :: lat_units = "degrees_north"
      character (len = *), parameter :: lon_units = "degrees_east"
      character (len = *), parameter :: z_units = "m"
      character (len = *), parameter :: long_name = "long_name"
      character (len = *), parameter :: lat_longname = "latitude"
      character (len = *), parameter :: lon_longname = "longitude"
      character (len = *), parameter :: z_longname = "sea depth"
      character (len = *), parameter :: t_longname = "Time"
      character (len = *) :: ncname
      character (len = *) :: t_units
      integer :: ncid,month
      integer :: nx,ny,nz
      real :: lat(ny), lon(nx), z(nz)
      integer :: z_dimid, lon_dimid, lat_dimid, t_dimid
      integer :: z_varid, lon_varid, lat_varid, t_varid

      call check( nf90_create(ncname, nf90_clobber, ncid) )

      call check( nf90_def_dim(ncid, latname, ny, lat_dimid) )
      call check( nf90_def_dim(ncid, lonname, nx, lon_dimid) )
      call check( nf90_def_dim(ncid, zname,   nz, z_dimid) )
      call check( nf90_def_dim(ncid, tname,   1,  t_dimid) )

      call check( nf90_def_var(ncid, latname, NF90_REAL, lat_dimid, lat_varid) )
      call check( nf90_def_var(ncid, lonname, NF90_REAL, lon_dimid, lon_varid) )
      call check( nf90_def_var(ncid, zname,   NF90_REAL, z_dimid,   z_varid) )
      call check( nf90_def_var(ncid, tname,   NF90_REAL, t_dimid,   t_varid) )

      call check( nf90_put_att(ncid, lat_varid, units, lat_units) )
      call check( nf90_put_att(ncid, lon_varid, units, lon_units) )
      call check( nf90_put_att(ncid, z_varid,   units, z_units) )
      call check( nf90_put_att(ncid, t_varid,   units, t_units) )
      call check( nf90_put_att(ncid, lat_varid, long_name, lat_longname) )
      call check( nf90_put_att(ncid, lon_varid, long_name, lon_longname) )
      call check( nf90_put_att(ncid, z_varid,   long_name, z_longname) )
      call check( nf90_put_att(ncid, t_varid,   long_name, t_longname) )

      call check( nf90_enddef(ncid) )

      call check( nf90_put_var(ncid, lat_varid, lat) )
      call check( nf90_put_var(ncid, lon_varid, lon) )
      call check( nf90_put_var(ncid, z_varid,   z) )
      call check( nf90_put_var(ncid, t_varid,   month) )

      call check( nf90_close(ncid) )
      print *,"*** SUCCESS writing coordinate information in file ",ncname,"!"

      return
      end subroutine netcdf_write_cor

      !--------------------------------------------------------
      subroutine netcdf_write_var2d(ncname,varname,var,nx,ny,var_units,var_longname,missvalue)
      !--------------------------------------------------------
      use netcdf
      implicit none
      character (len = *), parameter :: latname = "lat"
      character (len = *), parameter :: lonname = "lon"
      character (len = *), parameter :: tname = "t"
      character (len = *), parameter :: units = "units"
      character (len = *), parameter :: long_name = "long_name"
      character (len = *), parameter :: miss_name = "missing_value"
      character (len = *) :: ncname
      character (len = *) :: varname
      character (len = *) :: var_units
      character (len = *) :: var_longname
      integer :: lat_dimid, lon_dimid, t_dimid, varid
      integer :: ncid
      integer :: dimids(3)
      integer :: nx,ny
      real :: var(nx,ny)
      real :: missvalue

      call check( (nf90_open(ncname, nf90_write, ncid)) )
      call check( (nf90_redef(ncid)) )

      call check( nf90_inq_dimid(ncid, latname, lat_dimid) )
      call check( nf90_inq_dimid(ncid, lonname, lon_dimid) )
      call check( nf90_inq_dimid(ncid, tname,   t_dimid) )
      dimids = (/ lon_dimid, lat_dimid, t_dimid/)

      call check( nf90_def_var(ncid, varname, NF90_REAL, dimids, varid) )
      call check( nf90_put_att(ncid, varid, units, var_units) )
      call check( nf90_put_att(ncid, varid, long_name, var_longname) )
      call check( nf90_put_att(ncid, varid, miss_name, missvalue) )
      call check( nf90_enddef(ncid) )

      call check( nf90_put_var(ncid, varid, var) )

      call check( nf90_close(ncid) )
      print *,"*** SUCCESS writing var ",varname," in file ",ncname,"!"

      return
      end subroutine netcdf_write_var2d


      !--------------------------------------------------------
      subroutine netcdf_write_var(ncname,varname,var,nx,ny,nz,var_units,var_longname,missvalue)
      !--------------------------------------------------------
      use netcdf
      implicit none
      character (len = *), parameter :: latname = "lat"
      character (len = *), parameter :: lonname = "lon"
      character (len = *), parameter :: zname = "z"
      character (len = *), parameter :: tname = "t"
      character (len = *), parameter :: units = "units"
      character (len = *), parameter :: long_name = "long_name"
      character (len = *), parameter :: miss_name = "missing_value"
      character (len = *) :: ncname
      character (len = *) :: varname
      character (len = *) :: var_units
      character (len = *) :: var_longname
      integer :: lat_dimid, lon_dimid, z_dimid, t_dimid, varid
      integer :: ncid
      integer :: dimids(4)
      integer :: nx,ny,nz
      real :: var(nx,ny,nz)
      real :: missvalue

      call check( (nf90_open(ncname, nf90_write, ncid)) )
      call check( (nf90_redef(ncid)) )

      call check( nf90_inq_dimid(ncid, latname, lat_dimid) )
      call check( nf90_inq_dimid(ncid, lonname, lon_dimid) )
      call check( nf90_inq_dimid(ncid, zname,   z_dimid) )
      call check( nf90_inq_dimid(ncid, tname,   t_dimid) )
      dimids = (/ lon_dimid, lat_dimid, z_dimid, t_dimid/)

      call check( nf90_def_var(ncid, varname, NF90_REAL, dimids, varid) )
      call check( nf90_put_att(ncid, varid, units, var_units) )
      call check( nf90_put_att(ncid, varid, long_name, var_longname) )
      call check( nf90_put_att(ncid, varid, miss_name, missvalue) )
      call check( nf90_enddef(ncid) )

      call check( nf90_put_var(ncid, varid, var) )

      call check( nf90_close(ncid) )
      print *,"*** SUCCESS writing var ",varname," in file ",ncname,"!"

      return
      end subroutine netcdf_write_var

      subroutine check(status)
      use netcdf
      integer, intent ( in) :: status

      if(status /= nf90_noerr) then
        print *, trim(nf90_strerror(status))
        stop 2
      end if
      return
      end subroutine check
