
! Description: basic NetCDF input/output interface
!
!      Author: OU Yuyuan <ouyuyuan@lasg.iap.ac.cn>
!     Created: 2015-03-06 10:38:13 BJT
! Last Change: 2016-04-05 19:38:37 BJT

module mod_io !{{{1 
!-------------------------------------------------------{{{1
!imported variables-------------------------------------{{{1
  use netcdf
  use mod_kind, only: sglp, wp

  use mod_param, only: &
    glo_ni, glo_nj, nk, names, units, &
    missing_int, missing_float

  implicit none
  private

!public subroutines-------------------------------------{{{1

  public &
    io_create,        &
    io_create_grdvar, &
    io_get_dim_len,   &
    io_write,         &
    io_quick_output,  &
    io_read

!internal variables-------------------------------------{{{1

  integer :: ncid, varid, ndim1, ndim2, ndim3

!interfaces---------------------------------------------{{{1
!-----------------------------------------------------------
  interface io_write
    module procedure write_scalar
    module procedure write_r3d
    module procedure write_r2d
    module procedure write_r2d_rec
    module procedure write_r1d
    module procedure write_i3d
    module procedure write_i2d
  end interface

  interface io_quick_output
    module procedure quick_output_r3d
    module procedure quick_output_r2d
    module procedure quick_output_r1d
    module procedure quick_output_i3d
    module procedure quick_output_i2d
  end interface

  interface io_read
    module procedure read_r3d
    module procedure read_r1d
    module procedure read_i3d
    module procedure read_i2d
    module procedure read_scalar
  end interface

contains !{{{1
!-------------------------------------------------------{{{1

subroutine io_create (ncname) !{{{1

  character (len=*), intent(in) :: ncname

  integer :: dimid1, dimid2, dimid3, dimid4, i
  integer :: dimids(4), dimids_2d(3)

  call check ( nf90_create (ncname, NF90_CLOBBER, ncid)  )

  !def dim. {{{2
  call check ( nf90_def_dim (ncid, 'lon', glo_ni, dimid1) )
  call check ( nf90_def_dim (ncid, 'lat', glo_nj, dimid2) )
  call check ( nf90_def_dim (ncid, 'z',   nk, dimid3) )
  call check ( nf90_def_dim (ncid, 't', NF90_UNLIMITED, dimid4) )

  !def global attr. {{{2
  call check ( nf90_put_att (ncid, NF90_GLOBAL, & 
    'created', "by subroutine io_create in module mod_io") )

  ! def vars  !{{{2
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

  ! define 2d variables !{{{2
  dimids_2d = (/dimid1, dimid2, dimid4/)
  ! sea surface height
  call check ( nf90_def_var (ncid, names%ssh, nf90_float, dimids_2d, varid) )
  call check ( nf90_put_att (ncid, varid, '_FillValue', missing_float) )
  call check ( nf90_put_att (ncid, varid, 'units', units%ssh) )

  ! define 3d variables !{{{2
  dimids = (/dimid1, dimid2, dimid3, dimid4/)

  ! potential temperature
  call check ( nf90_def_var (ncid, names%pt, nf90_float, &
    dimids, varid) )
  call check ( nf90_put_att (ncid, varid, '_FillValue', &
    missing_float) )
  call check ( nf90_put_att (ncid, varid, 'units', &
    units%pt) )

  ! salinity
  call check ( nf90_def_var (ncid, names%sa, nf90_float, &
    dimids, varid) )
  call check ( nf90_put_att (ncid, varid, '_FillValue', &
    missing_float) )
  call check ( nf90_put_att (ncid, varid, 'units', &
    units%sa) )

  ! zonal velocity
  call check ( nf90_def_var (ncid, names%u, nf90_float, dimids, varid) )
  call check ( nf90_put_att (ncid, varid, '_FillValue', missing_float) )
  call check ( nf90_put_att (ncid, varid, 'units', units%u) )

  ! meridional velocity
  call check ( nf90_def_var (ncid, names%v, nf90_float, dimids, varid) )
  call check ( nf90_put_att (ncid, varid, '_FillValue', missing_float) )
  call check ( nf90_put_att (ncid, varid, 'units', units%v) )

  ! vertical velocity
  call check ( nf90_def_var (ncid, names%w, nf90_float, dimids, varid) )
  call check ( nf90_put_att (ncid, varid, '_FillValue', missing_float) )
  call check ( nf90_put_att (ncid, varid, 'units', units%w) )

  !end def {{{2
  call check (nf90_enddef(ncid) )

  call check (nf90_close(ncid) )

end subroutine io_create

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
  call check ( nf90_def_var (ncid, 'pb', nf90_float, &
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

subroutine write_r3d(ncname, varname, var, nrec) !{{{1

  character (len = *), intent(in) :: ncname, varname
  real (kind=wp), intent(in) :: var(:,:,:)

  integer :: stt(4), cnt(4), nrec

  write(*,'(a,i3,a)') '*** Output '//varname//' to file: '//ncname//' for the ', &
    nrec, 'th record ~ ~ ~' 

  ndim1 = size(var, 1)
  ndim2 = size(var, 2)
  ndim3 = size(var, 3)

  stt  = (/1, 1, 1, nrec/)
  cnt  = (/ndim1, ndim2, ndim3, 1/)

  call check (nf90_open(ncname, nf90_write, ncid)  )

  call check (nf90_inq_varid(ncid, varname, varid) )

  call check (nf90_put_var(ncid, varid, var, &
    start = stt, count = cnt) )

  call check (nf90_close(ncid) )

end subroutine write_r3d

subroutine write_r2d_rec(ncname, varname, var, nrec) !{{{1

  character (len = *), intent(in) :: ncname, varname
  real (kind=wp), intent(in) :: var(:,:)

  integer :: stt(3), cnt(3), nrec

  write(*,'(a,i3,a)') '*** Output '//varname//' to file: '//ncname//' for the ', &
    nrec, 'th record ~ ~ ~' 

  ndim1 = size(var, 1)
  ndim2 = size(var, 2)

  stt  = (/1, 1, nrec/)
  cnt  = (/ndim1, ndim2, 1/)

  call check (nf90_open(ncname, nf90_write, ncid)  )

  call check (nf90_inq_varid(ncid, varname, varid) )

  call check (nf90_put_var(ncid, varid, var, &
    start = stt, count = cnt) )

  call check (nf90_close(ncid) )

end subroutine write_r2d_rec

subroutine write_r2d(ncname, varname, var) !{{{1

  character (len = *), intent(in) :: ncname, varname
  real (kind=wp), intent(in) :: var(:,:)

  integer :: stt(2), cnt(2)

  write(*,'(a,i3,a)') '*** Output '//varname//' to file: '//ncname

  ndim1 = size(var, 1)
  ndim2 = size(var, 2)

  stt  = (/1, 1/)
  cnt  = (/ndim1, ndim2/)

  call check (nf90_open(ncname, nf90_write, ncid)  )

  call check (nf90_inq_varid(ncid, varname, varid) )

  call check (nf90_put_var(ncid, varid, var, &
    start = stt, count = cnt) )

  call check (nf90_close(ncid) )

end subroutine write_r2d

subroutine write_r1d(ncname, varname, var) !{{{1

  character (len = *), intent(in) :: ncname, varname
  real (kind=wp), intent(in) :: var(:)

  integer :: is

  call check (nf90_open(ncname, nf90_write, ncid)  )

  call check (nf90_inq_varid(ncid, varname, varid) )

  call check (nf90_put_var(ncid, varid, var) )

  call check (nf90_close(ncid) )

  write(*,*) '*** SUCCESS output '//varname//' to file: '//ncname

end subroutine write_r1d

subroutine write_i3d(ncname, varname, var) !{{{1

  character (len = *), intent(in) :: ncname, varname
  integer, intent(in) :: var(:,:,:)

  integer :: stt(3), cnt(3)

  write(*,'(a,i3,a)') '*** Output '//varname//' to file: '//ncname

  ndim1 = size(var, 1)
  ndim2 = size(var, 2)
  ndim3 = size(var, 3)

  stt  = (/1, 1, 1/)
  cnt  = (/ndim1, ndim2, ndim3/)

  call check (nf90_open(ncname, nf90_write, ncid)  )

  call check (nf90_inq_varid(ncid, varname, varid) )

  call check (nf90_put_var(ncid, varid, var, &
    start = stt, count = cnt) )

  call check (nf90_close(ncid) )

end subroutine write_i3d

subroutine write_i2d(ncname, varname, var) !{{{1

  character (len = *), intent(in) :: ncname, varname
  integer, intent(in) :: var(:,:)

  integer :: stt(2), cnt(2)

  write(*,'(a,i3,a)') '*** Output '//varname//' to file: '//ncname

  ndim1 = size(var, 1)
  ndim2 = size(var, 2)

  stt  = (/1, 1/)
  cnt  = (/ndim1, ndim2/)

  call check (nf90_open(ncname, nf90_write, ncid)  )

  call check (nf90_inq_varid(ncid, varname, varid) )

  call check (nf90_put_var(ncid, varid, var, &
    start = stt, count = cnt) )

  call check (nf90_close(ncid) )

end subroutine write_i2d

subroutine write_scalar(ncname, varname, var) !{{{1

  character (len = *), intent(in) :: ncname, varname
  real (kind=wp), intent(in) :: var

  call check (nf90_open(ncname, nf90_write, ncid)  )

  call check (nf90_inq_varid(ncid, varname, varid) )

  call check (nf90_put_var(ncid, varid, var) )

  call check (nf90_close(ncid) )

  write(*,*) '*** SUCCESS output '//varname//' to file: '//ncname

end subroutine write_scalar

subroutine quick_output_r3d ( ncname, varname, var, lon, lat, z ) !{{{1
  ! creat a new file and output a 3d real variable

  character (len = *), intent(in) :: ncname, varname
  real (kind=wp),      intent(in) :: var(:,:,:)
  real (kind=wp), dimension(:), intent(in) :: lon, lat, z

  integer :: dimids(4)
  integer :: nrec, stt(4), cnt(4), &
             dimid1, dimid2, dimid3, time_dimid, &
             lonid, latid, zid

  ndim1 = size(var, 1)
  ndim2 = size(var, 2)
  ndim3 = size(var, 3)

  nrec = 1
  stt  = (/1, 1, 1, nrec/)
  cnt  = (/ndim1, ndim2, ndim3, 1/)

  call check (nf90_create(ncname, NF90_CLOBBER, ncid)  )

  ! def dim. {{{2
  call check (nf90_def_dim(ncid, 'lon', ndim1, dimid1) )
  call check (nf90_def_dim(ncid, 'lat', ndim2, dimid2) )
  call check (nf90_def_dim(ncid, 'z',   ndim3, dimid3) )
  call check (nf90_def_dim(ncid, 'time', NF90_UNLIMITED, &
    time_dimid) )

  ! def vars !{{{2
  call check (nf90_def_var(ncid, 'lon', nf90_float, dimid1, lonid) )
  call check (nf90_def_var(ncid, 'lat', nf90_float, dimid2, latid) )
  call check (nf90_def_var(ncid, 'z',   nf90_float, dimid3, zid) )

  dimids =  (/ dimid1, dimid2, dimid3, time_dimid /)
  call check (nf90_def_var(ncid, varname, nf90_float, dimids, &
    varid) )
  call check (nf90_put_att(ncid, varid, '_FillValue', &
    missing_float) )

  call check (nf90_enddef(ncid) )

  ! write vars !{{{2
  call check (nf90_put_var(ncid, lonid, lon))
  call check (nf90_put_var(ncid, latid, lat))
  call check (nf90_put_var(ncid, zid, z))

  call check (nf90_put_var(ncid, varid, var, &
    start = stt, count = cnt) )

  ! close file !{{{2
  call check (nf90_close(ncid) )

  write(*,*) '*** SUCCESS output '//varname//' to file: '//ncname

end subroutine quick_output_r3d

subroutine quick_output_r2d ( ncname, varname, var, lon, lat ) !{{{1
  ! creat a new file and output a 2d real variable

  character (len = *), intent(in) :: ncname, varname
  real (kind=wp),      intent(in) :: var(:,:)
  real (kind=wp), dimension(:), intent(in) :: lon, lat

  integer :: dimids(3)
  integer :: nrec, stt(3), cnt(3), &
             dimid1, dimid2, time_dimid, &
             lonid, latid

  ndim1 = size(var, 1)
  ndim2 = size(var, 2)

  nrec = 1
  stt  = (/1, 1, nrec/)
  cnt  = (/ndim1, ndim2, 1/)

  call check (nf90_create(ncname, NF90_CLOBBER, ncid)  )

  ! def dim. {{{2
  call check (nf90_def_dim(ncid, 'lon', ndim1, dimid1) )
  call check (nf90_def_dim(ncid, 'lat', ndim2, dimid2) )
  call check (nf90_def_dim(ncid, 'time', NF90_UNLIMITED, &
    time_dimid) )

  ! def vars !{{{2
  dimids =  (/ dimid1, dimid2, time_dimid /)
  call check (nf90_def_var(ncid, 'lon', nf90_float, dimid1, lonid) )
  call check (nf90_def_var(ncid, 'lat', nf90_float, dimid2, latid) )

  call check (nf90_def_var(ncid, varname, nf90_float, dimids, &
    varid) )
  call check (nf90_put_att(ncid, varid, '_FillValue', &
    missing_float) )

  call check (nf90_enddef(ncid) )

  ! write var !{{{2
  call check (nf90_put_var(ncid, lonid, lon))
  call check (nf90_put_var(ncid, latid, lat))

  call check (nf90_put_var(ncid, varid, var, &
    start = stt, count = cnt) )

  ! close file !{{{2
  call check (nf90_close(ncid) )

  write(*,*) '*** SUCCESS output '//varname//' to file: '//ncname

end subroutine quick_output_r2d

subroutine quick_output_r1d ( ncname, varname, var, coor ) !{{{1
  ! creat a new file and output a 1d real variable

  character (len = *), intent(in) :: ncname, varname
  real (kind=wp), dimension(:), intent(in) :: var, coor

  integer :: dimid1, time_dimid, stt(2), cnt(2), &
    coorid, dimids(2), nrec

  ndim1 = size(var)

  nrec = 1
  stt  = (/1, nrec/)
  cnt  = (/ndim1, 1/)

  call check (nf90_create(ncname, NF90_CLOBBER, ncid)  )

  ! def dim. {{{2
  call check (nf90_def_dim(ncid, 'coor', ndim1, dimid1) )
  call check (nf90_def_var(ncid, 'coor', nf90_float, dimid1, coorid))

  call check (nf90_def_dim(ncid, 'time', NF90_UNLIMITED, &
    time_dimid) )

  ! def var !{{{2
  dimids =  (/ dimid1, time_dimid /)
  call check (nf90_def_var(ncid, varname, nf90_float, dimids, &
    varid) )
  call check (nf90_put_att(ncid, varid, '_FillValue', &
    missing_float) )

  call check (nf90_enddef(ncid) )

  ! write vars !{{{2
  call check (nf90_put_var(ncid, coorid, coor))

  call check (nf90_put_var(ncid, varid, var, &
    start = stt, count = cnt) )

  ! close file !{{{2
  call check (nf90_close(ncid) )

  write(*,*) '*** SUCCESS output '//varname//' to file: '//ncname

end subroutine quick_output_r1d

subroutine quick_output_i3d ( ncname, varname, var, lon, lat, z ) !{{{1
  ! creat a new file and output a 3d real variable

  character (len = *), intent(in) :: ncname, varname
  integer,             intent(in) :: var(:,:,:)
  real (kind=wp), dimension(:), intent(in) :: lon, lat, z

  integer :: dimids(4)
  integer :: nrec, stt(4), cnt(4), &
             dimid1, dimid2, dimid3, time_dimid, &
             lonid, latid, zid

  ndim1 = size(var, 1)
  ndim2 = size(var, 2)
  ndim3 = size(var, 3)

  nrec = 1
  stt  = (/1, 1, 1, nrec/)
  cnt  = (/ndim1, ndim2, ndim3, 1/)

  call check (nf90_create(ncname, NF90_CLOBBER, ncid)  )

  ! def dim. {{{2
  call check (nf90_def_dim(ncid, 'lon', ndim1, dimid1) )
  call check (nf90_def_dim(ncid, 'lat', ndim2, dimid2) )
  call check (nf90_def_dim(ncid, 'z',   ndim3, dimid3) )
  call check (nf90_def_dim(ncid, 'time', NF90_UNLIMITED, &
    time_dimid) )

  ! def vars !{{{2
  call check (nf90_def_var(ncid, 'lon', nf90_float, dimid1, lonid) )
  call check (nf90_def_var(ncid, 'lat', nf90_float, dimid2, latid) )
  call check (nf90_def_var(ncid, 'z',   nf90_float, dimid3, zid) )

  dimids =  (/ dimid1, dimid2, dimid3, time_dimid /)
  call check (nf90_def_var(ncid, varname, nf90_int, dimids, &
    varid) )
  call check (nf90_put_att(ncid, varid, '_FillValue', &
    missing_int) )

  call check (nf90_enddef(ncid) )

  ! write vars !{{{2
  call check (nf90_put_var(ncid, lonid, lon))
  call check (nf90_put_var(ncid, latid, lat))
  call check (nf90_put_var(ncid, zid, z))

  call check (nf90_put_var(ncid, varid, var, &
    start = stt, count = cnt) )

  ! close file !{{{2
  call check (nf90_close(ncid) )

  write(*,*) '*** SUCCESS output '//varname//' to file: '//ncname

end subroutine quick_output_i3d

subroutine quick_output_i2d ( ncname, varname, var, lon, lat ) !{{{1
  ! creat a new file and output a 2d real variable

  character (len = *), intent(in) :: ncname, varname
  integer,             intent(in) :: var(:,:)
  real (kind=wp), dimension(:), intent(in) :: lon, lat

  integer :: dimids(3)
  integer :: nrec, stt(3), cnt(3), &
             dimid1, dimid2, time_dimid, &
             lonid, latid

  ndim1 = size(var, 1)
  ndim2 = size(var, 2)

  nrec = 1
  stt  = (/1, 1, nrec/)
  cnt  = (/ndim1, ndim2, 1/)

  call check (nf90_create(ncname, NF90_CLOBBER, ncid)  )

  ! def dim. {{{2
  call check (nf90_def_dim(ncid, 'lon', ndim1, dimid1) )
  call check (nf90_def_dim(ncid, 'lat', ndim2, dimid2) )
  call check (nf90_def_dim(ncid, 'time', NF90_UNLIMITED, &
    time_dimid) )

  ! def vars !{{{2
  dimids =  (/ dimid1, dimid2, time_dimid /)
  call check (nf90_def_var(ncid, 'lon', nf90_float, dimid1, lonid) )
  call check (nf90_def_var(ncid, 'lat', nf90_float, dimid2, latid) )

  call check (nf90_def_var(ncid, varname, nf90_int, dimids, &
    varid) )
  call check (nf90_put_att(ncid, varid, '_FillValue', &
    missing_int) )

  call check (nf90_enddef(ncid) )

  ! write var !{{{2
  call check (nf90_put_var(ncid, lonid, lon))
  call check (nf90_put_var(ncid, latid, lat))

  ! quit odd, it will print the first four elements in the first row 
  ! of var to the screen when execute at BCM, maybe this is a bug of netcdf
  call check (nf90_put_var(ncid, varid, var, &
    start = stt, count = cnt) )

  ! close file !{{{2
  call check (nf90_close(ncid) )

  write(*,*) '*** SUCCESS output '//varname//' to file: '//ncname

end subroutine quick_output_i2d

subroutine read_r3d(ncname, varname, var) !{{{1
  character (len=*), intent(in) :: ncname, varname
  real (kind=wp) :: var(:,:,:)

  call check (nf90_open(ncname, NF90_NOWRITE, ncid) )

  call check (nf90_inq_varid(ncid, varname, varid) )

  call check (nf90_get_var(ncid, varid, var) )

  call check (nf90_close(ncid) )

  write(*,'(a)') 'got '//trim(varname)// ' from '//trim(ncname)

end subroutine read_r3d 

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