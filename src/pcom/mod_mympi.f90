! Description: mpi related routine
!   if nm%px = 4, nm%py = 2, then id matrix is:
!     (NOTE that east-west direction is wrap up)
!
!   0 4
!   1 5
!   2 6
!   3 7
!
!      Author: OU Yuyuan <ouyuyuan@lasg.iap.ac.cn>
!     Created: 2015-09-14 14:25:29 BJT
! Last Change: 2016-04-06 16:15:29 BJT

module mod_mympi

  ! imported variables !{{{1
  use mod_io, only: &
    io_write, io_quick_output

  use mod_kind, only: wp

  use mpi

  use mod_param, only: nm, mid, myid, npro, &
    glo_nj, glo_ni, ni, nj, nk, &
    my, our, type_my, missing_float

  use mod_type, only: &
    type_mat, type_accu_gm3d, type_accu_gr3d, &
    type_accu_gr2d, &
    type_gvar_m3d, type_gvar_m2d, type_gvar_r2d, type_gvar_r3d
    
  implicit none
  private

  ! interfaces !{{{1
  public &
    mympi_div, &
    mympi_divx, &
    mympi_divy, &
    mympi_swpbnd, &
    mympi_quick_output, &
    mympi_bcast, &
    mympi_merge, &
    mympi_output

  interface mympi_bcast
    module procedure bcast_r1d
    module procedure bcast_i
    module procedure bcast_string
  end interface

  interface mympi_quick_output
    module procedure quick_output_r3d
    module procedure quick_output_r2d
    module procedure quick_output_r1d
    module procedure quick_output_i3d
    module procedure quick_output_i2d
  end interface

  interface mympi_div
    module procedure div_gvar_m3d
    module procedure div_gvar_r3d
    module procedure div_r3d_gvar_r3d
    module procedure div_m3d
    module procedure div_r3d
    module procedure div_r2d
    module procedure div_i3d
    module procedure div_i2d
  end interface

  interface mympi_swpbnd
    module procedure swpbnd_r3d
    module procedure swpbnd_r2d
    module procedure swpbnd_r2d_2
    module procedure swpbnd_r1d
    module procedure swpbnd_i3d
    module procedure swpbnd_i2d
    module procedure swpbnd_m3d
    module procedure swpbnd_m2d
    module procedure swpbnd_m2d_4
  end interface

  interface mympi_output
    module procedure merge_out_accu_gm3d
    module procedure merge_out_accu_gm3d_mask
    module procedure merge_out_accu_gr3d
    module procedure merge_out_accu_gr3d_mask
    module procedure merge_out_accu_gr2d_mask
    module procedure merge_out_accu_gr2d_weight
    module procedure merge_out_r3d
    module procedure merge_out_r2d
    module procedure merge_out_i3d
    module procedure merge_out_i2d
  end interface

  interface mympi_merge
    module procedure merge_r2d
  end interface

  ! local variables !{{{1
  integer :: is, err, msta(mpi_status_size)

contains !{{{1

subroutine bcast_r1d (var)!{{{1
  ! broadcast var from mid to all ids
  real (kind=wp), dimension(:) :: var

  call mpi_bcast(var, size(var), mpi_real8, mid, mpi_comm_world, err)

end subroutine bcast_r1d

subroutine bcast_i (var)!{{{1
  ! broadcast var from mid to all ids
  integer :: var

  call mpi_bcast(var, 1, mpi_integer, mid, mpi_comm_world, err)

end subroutine bcast_i

subroutine bcast_string (var)!{{{1
  ! broadcast var from mid to all ids
  character (len=*) :: var

  call mpi_bcast(var, len(var), mpi_byte, mid, mpi_comm_world, err)
end subroutine bcast_string

subroutine merge_out_accu_gm3d (fname, vname1, vname2, var) !{{{1
  ! merge 3d array from other domains to mid
  character (len=*) :: fname, vname1, vname2
  type (type_accu_gm3d) :: var

  var%var%x(1)%v = var%var%x(1)%v / var%n
  var%var%x(2)%v = var%var%x(2)%v / var%n

  call merge_out_r3d (fname, vname1, var%var%x(1)%v, var%nrec + 1)
  call merge_out_r3d (fname, vname2, var%var%x(2)%v, var%nrec + 1)

  var%var%x(1)%v = 0.0 ! reset accumulated value
  var%var%x(2)%v = 0.0
  var%n    = 0 ! reset counter
  var%nrec = var%nrec + 1 ! nc record have been increased

end subroutine merge_out_accu_gm3d

subroutine merge_out_accu_gm3d_mask (fname, vname1, vname2, var, mask) !{{{1
  ! merge 3d array from other domains to mid
  ! mask out the land points
  character (len=*) :: fname, vname1, vname2
  type (type_accu_gm3d) :: var
  integer, dimension(:,:,:), intent(in) :: mask

  var%var%x(1)%v = var%var%x(1)%v / var%n
  var%var%x(2)%v = var%var%x(2)%v / var%n

  where (mask == 0)
    var%var%x(1)%v = missing_float
    var%var%x(2)%v = missing_float
  end where

  call merge_out_r3d (fname, vname1, var%var%x(1)%v, var%nrec + 1)
  call merge_out_r3d (fname, vname2, var%var%x(2)%v, var%nrec + 1)

  var%var%x(1)%v = 0.0 ! reset accumulated value
  var%var%x(2)%v = 0.0
  var%n    = 0 ! reset counter
  var%nrec = var%nrec + 1 ! nc record have been increased

end subroutine merge_out_accu_gm3d_mask

subroutine merge_out_accu_gr3d (fname, vname, var) !{{{1
  ! merge 3d array from other domains to mid
  character (len=*) :: fname, vname
  type (type_accu_gr3d) :: var

  var%var%v = var%var%v / var%n

  call merge_out_r3d (fname, vname, var%var%v, var%nrec + 1)

  var%var%v = 0.0 ! reset accumulated value
  var%n     = 0 ! reset counter
  var%nrec  = var%nrec + 1 ! nc record have been increased

end subroutine merge_out_accu_gr3d

subroutine merge_out_accu_gr3d_mask (fname, vname, var, mask) !{{{1
  ! merge 3d array from other domains to mid
  ! mask out the land points
  character (len=*) :: fname, vname
  type (type_accu_gr3d) :: var
  integer, dimension(:,:,:), intent(in) :: mask

  var%var%v = var%var%v / var%n
  where (mask == 0)
    var%var%v = missing_float
  end where

  call merge_out_r3d (fname, vname, var%var%v, var%nrec + 1)

  var%var%v = 0.0 ! reset accumulated value
  var%n     = 0 ! reset counter
  var%nrec  = var%nrec + 1 ! nc record have been increased

end subroutine merge_out_accu_gr3d_mask

subroutine merge_out_accu_gr2d_mask (fname, vname, var, mask) !{{{1
  ! merge 2d array from other domains to mid
  ! mask out the land points
  character (len=*) :: fname, vname
  type (type_accu_gr2d) :: var
  integer, dimension(:,:), intent(in) :: mask

  var%var%v = var%var%v / var%n
  where (mask == 0)
    var%var%v = missing_float
  end where

  call merge_out_r2d_rec (fname, vname, var%var%v, var%nrec + 1)

  var%var%v = 0.0 ! reset accumulated value
  var%n     = 0 ! reset counter
  var%nrec  = var%nrec + 1 ! nc record have been increased

end subroutine merge_out_accu_gr2d_mask

subroutine merge_out_accu_gr2d_weight (fname, vname, var, wg, mask) !{{{1
  ! merge 2d array from other domains to mid
  ! mask out the land points, remove global mean
  character (len=*) :: fname, vname
  type (type_accu_gr2d) :: var
  real (kind=wp), dimension(:,:), intent(in) :: wg
  integer, dimension(:,:), intent(in) :: mask

  var%var%v = var%var%v / var%n
  where (mask == 0)
    var%var%v = missing_float
  end where

  call merge_out_r2d_weight (fname, vname, var%var%v, wg, var%nrec + 1)

  var%var%v = 0.0 ! reset accumulated value
  var%n     = 0 ! reset counter
  var%nrec  = var%nrec + 1 ! nc record have been increased

end subroutine merge_out_accu_gr2d_weight

subroutine merge_out_r3d (fname, varname, var, nrec) !{{{1
  ! merge 3d array from other domains to mid
  character (len=*), intent(in) :: fname, varname
  real (kind=wp), dimension(ni,nj,nk), intent(in) :: var
  integer, intent(in) :: nrec

  real (kind=wp), allocatable, dimension(:,:,:) :: glo_var

  integer, parameter :: tag = 30
  type (type_my) :: d
  integer :: n, leng

  if (myid == mid) then

    allocate( glo_var(glo_ni, glo_nj, nk), stat=is)
    call chk(is); glo_var = 0.0

    do n = 1, npro
      d = our(n)
      leng  = (d%ge-d%gw+1) * (d%gn-d%gs+1) * nk
      if ( d%id == mid ) then
        glo_var(d%gw:d%ge, d%gs:d%gn, :) = &
          var(2:d%ni-1, 2:d%nj-1, :)
      else
        call mpi_recv (glo_var(d%gw:d%ge,d%gs:d%gn,:), &
          leng, mpi_real8, d%id, tag, mpi_comm_world, msta, err)
      end if
    end do

    call io_write (trim(fname), trim(varname), glo_var, nrec)

    deallocate(glo_var)

  else
    leng  = (my%ge-my%gw+1) * (my%gn-my%gs+1) * nk
    call mpi_ssend (var(2:my%ni-1,2:my%nj-1,:), leng, &
      mpi_real8, mid, tag, mpi_comm_world, err)
  end if

end subroutine merge_out_r3d

subroutine merge_out_r2d_weight (fname, varname, var, wg, nrec) !{{{1
  ! merge 3d array from other domains to mid
  character (len=*), intent(in) :: fname, varname
  real (kind=wp), dimension(ni,nj), intent(in) :: var
  real (kind=wp), dimension(ni,nj), intent(in) :: wg
  integer, intent(in) :: nrec

  real (kind=wp), allocatable, dimension(:,:) :: glo_var
  real (kind=wp), allocatable, dimension(:,:) :: glo_wg

  real :: mean

  if (myid == mid) then
    allocate( glo_var(glo_ni, glo_nj), stat=is)
    call chk(is); glo_var = 0.0
    allocate( glo_wg(glo_ni, glo_nj), stat=is)
    call chk(is); glo_wg = 0.0
  end if

  call merge_r2d( glo_var, var )
  call merge_r2d( glo_wg, wg )

  if (myid == mid) then
    mean = sum(glo_var*glo_wg)/sum(glo_wg)
    where (glo_var /= missing_float) glo_var = glo_var - mean

    call io_write (trim(fname), trim(varname), glo_var, nrec)

    deallocate(glo_var)
    deallocate(glo_wg)
  end if

end subroutine merge_out_r2d_weight

subroutine merge_out_r2d_rec (fname, varname, var, nrec) !{{{1
  ! merge 3d array from other domains to mid
  character (len=*), intent(in) :: fname, varname
  real (kind=wp), dimension(ni,nj), intent(in) :: var
  integer, intent(in) :: nrec

  real (kind=wp), allocatable, dimension(:,:) :: glo_var

  integer, parameter :: tag = 30
  type (type_my) :: d
  integer :: n, leng

  if (myid == mid) then

    allocate( glo_var(glo_ni, glo_nj), stat=is)
    call chk(is); glo_var = 0.0

    do n = 1, npro
      d = our(n)
      leng  = (d%ge-d%gw+1) * (d%gn-d%gs+1)
      if ( d%id == mid ) then
        glo_var(d%gw:d%ge, d%gs:d%gn) = &
          var(2:d%ni-1, 2:d%nj-1)
      else
        call mpi_recv (glo_var(d%gw:d%ge,d%gs:d%gn), &
          leng, mpi_real8, d%id, tag, mpi_comm_world, msta, err)
      end if
    end do

    call io_write (trim(fname), trim(varname), glo_var, nrec)
    deallocate(glo_var)

  else
    leng  = (my%ge-my%gw+1) * (my%gn-my%gs+1)
    call mpi_ssend (var(2:my%ni-1,2:my%nj-1), leng, &
      mpi_real8, mid, tag, mpi_comm_world, err)
  end if

end subroutine merge_out_r2d_rec

subroutine merge_out_r2d (fname, varname, var) !{{{1
  ! merge 2d array from other domains to mid
  character (len=*), intent(in) :: fname, varname
  real (kind=wp), dimension(ni,nj), intent(in) :: var

  real (kind=wp), allocatable, dimension(:,:) :: glo_var

  integer, parameter :: tag = 30
  type (type_my) :: d
  integer :: n, leng

  if (myid == mid) then

    allocate( glo_var(glo_ni, glo_nj), stat=is)
    call chk(is); glo_var = 0.0

    do n = 1, npro
      d = our(n)
      leng  = (d%ge-d%gw+1) * (d%gn-d%gs+1)
      if ( d%id == mid ) then
        glo_var(d%gw:d%ge, d%gs:d%gn) = &
          var(2:d%ni-1, 2:d%nj-1)
      else
        call mpi_recv (glo_var(d%gw:d%ge,d%gs:d%gn), &
          leng, mpi_real8, d%id, tag, mpi_comm_world, msta, err)
      end if
    end do

    call io_write (trim(fname), trim(varname), glo_var)
    deallocate(glo_var)

  else
    leng  = (my%ge-my%gw+1) * (my%gn-my%gs+1)
    call mpi_ssend (var(2:my%ni-1,2:my%nj-1), leng, &
      mpi_real8, mid, tag, mpi_comm_world, err)
  end if

end subroutine merge_out_r2d

subroutine merge_out_i3d (fname, varname, var) !{{{1
  ! merge 3d array from other domains to mid
  character (len=*), intent(in) :: fname, varname
  integer, dimension(ni,nj,nk), intent(in) :: var

  integer, allocatable, dimension(:,:,:) :: glo_var

  integer, parameter :: tag = 30
  type (type_my) :: d
  integer :: n, leng

  if (myid == mid) then

    allocate( glo_var(glo_ni, glo_nj, nk), stat=is)
    call chk(is); glo_var = 0.0

    do n = 1, npro
      d = our(n)
      leng  = (d%ge-d%gw+1) * (d%gn-d%gs+1) * nk
      if ( d%id == mid ) then
        glo_var(d%gw:d%ge, d%gs:d%gn, :) = &
          var(2:d%ni-1, 2:d%nj-1, :)
      else
        call mpi_recv (glo_var(d%gw:d%ge,d%gs:d%gn,:), &
          leng, mpi_integer, d%id, tag, mpi_comm_world, msta, err)
      end if
    end do

    call io_write (trim(fname), trim(varname), glo_var)
    deallocate(glo_var)

  else
    leng  = (my%ge-my%gw+1) * (my%gn-my%gs+1) * nk
    call mpi_ssend (var(2:my%ni-1,2:my%nj-1,:), leng, &
      mpi_integer, mid, tag, mpi_comm_world, err)
  end if

end subroutine merge_out_i3d

subroutine merge_out_i2d (fname, varname, var) !{{{1
  ! merge 3d array from other domains to mid
  character (len=*), intent(in) :: fname, varname
  integer, dimension(ni,nj), intent(in) :: var

  integer, allocatable, dimension(:,:) :: glo_var

  integer, parameter :: tag = 30
  type (type_my) :: d
  integer :: n, leng

  if (myid == mid) then

    allocate( glo_var(glo_ni, glo_nj), stat=is)
    call chk(is); glo_var = 0.0

    do n = 1, npro
      d = our(n)
      leng  = (d%ge-d%gw+1) * (d%gn-d%gs+1)
      if ( d%id == mid ) then
        glo_var(d%gw:d%ge, d%gs:d%gn) = &
          var(2:d%ni-1, 2:d%nj-1)
      else
        call mpi_recv (glo_var(d%gw:d%ge,d%gs:d%gn), &
          leng, mpi_integer, d%id, tag, mpi_comm_world, msta, err)
      end if
    end do

    call io_write (trim(fname), trim(varname), glo_var)
    deallocate(glo_var)

  else
    leng  = (my%ge-my%gw+1) * (my%gn-my%gs+1)
    call mpi_ssend (var(2:my%ni-1,2:my%nj-1), leng, &
      mpi_integer, mid, tag, mpi_comm_world, err)
  end if

end subroutine merge_out_i2d

subroutine merge_r2d (glo_var, var) !{{{1
  ! merge 2d array from other domains to mid
  real (kind=wp), dimension(glo_ni,glo_nj) :: glo_var
  real (kind=wp), dimension(ni,nj), intent(in) :: var

  integer, parameter :: tag = 30
  type (type_my) :: d
  integer :: n, leng

  if (myid == mid) then

    do n = 1, npro
      d = our(n)
      leng  = (d%ge-d%gw+1) * (d%gn-d%gs+1)
      if ( d%id == mid ) then
        glo_var(d%gw:d%ge, d%gs:d%gn) = &
          var(2:d%ni-1, 2:d%nj-1)
      else
        call mpi_recv (glo_var(d%gw:d%ge,d%gs:d%gn), &
          leng, mpi_real8, d%id, tag, mpi_comm_world, msta, err)
      end if
    end do

  else
    leng  = (my%ge-my%gw+1) * (my%gn-my%gs+1)
    call mpi_ssend (var(2:my%ni-1,2:my%nj-1), leng, &
      mpi_real8, mid, tag, mpi_comm_world, err)
  end if

end subroutine merge_r2d

subroutine quick_output_r3d (fname, varname, var, lon, lat, z) !{{{1
  ! merge 2d array from other domains to mid and output
  character (len=*) :: fname, varname
  real (kind=wp), intent(in) :: var(:,:,:)
  real (kind=wp), dimension(:), intent(in) :: lon, lat, z

  real (kind=wp), allocatable :: glo_var(:,:,:)

  integer, parameter :: tag = 50
  type (type_my) :: d
  integer :: n, leng

  if (myid == mid) then

    allocate( glo_var(glo_ni, glo_nj, nk), stat=is)
    call chk(is); glo_var = 0

    do n = 1, npro
      d = our(n)
      leng = (d%ge-d%gw+1) * (d%gn-d%gs+1) * nk
      if ( d%id == mid ) then
        glo_var(d%gw:d%ge, d%gs:d%gn, :) = var(2:my%ni-1, 2:my%nj-1, :)
      else
        call mpi_recv (glo_var(d%gw:d%ge,d%gs:d%gn,:), &
          leng, mpi_real8, d%id, tag, mpi_comm_world, msta, err)
      end if
    end do

    call io_quick_output (trim(fname), trim(varname), glo_var, lon, lat, z)
    deallocate(glo_var)

  else
    leng = (my%ge-my%gw+1) * (my%gn-my%gs+1) * nk
    call mpi_ssend (var(2:my%ni-1,2:my%nj-1,:), leng, &
      mpi_real8, mid, tag, mpi_comm_world, err)
  end if

end subroutine quick_output_r3d

subroutine quick_output_r2d (fname, varname, var, lon, lat) !{{{1
  ! merge 2d array from other domains to mid and output
  character (len=*) :: fname, varname
  real (kind=wp), intent(in) :: var(ni,nj)
  real (kind=wp), dimension(:), intent(in) :: lon, lat

  real (kind=wp), allocatable, dimension(:,:) :: glo_var

  integer, parameter :: tag = 30
  type (type_my) :: d
  integer :: n, leng

  if (myid == mid) then

    allocate( glo_var(glo_ni, glo_nj), stat=is)
    call chk(is); glo_var = 0.0

    do n = 1, npro
      d = our(n)
      leng = (d%ge-d%gw+1) * (d%gn-d%gs+1)
      if ( d%id == myid ) then
        glo_var(d%gw:d%ge, d%gs:d%gn) = var(2:d%ni-1, 2:d%nj-1)
      else
        call mpi_recv (glo_var(d%gw:d%ge,d%gs:d%gn), &
          leng, mpi_real8, d%id, tag, mpi_comm_world, msta, err)
      end if
    end do

    call io_quick_output (trim(fname), trim(varname), glo_var, lon, lat)
    deallocate(glo_var)

  else
    leng = (my%ge-my%gw+1) * (my%gn-my%gs+1)
    call mpi_ssend (var(2:my%ni-1,2:my%nj-1), leng, &
      mpi_real8, mid, tag, mpi_comm_world, err)
  end if

end subroutine quick_output_r2d

subroutine quick_output_r1d (fname, varname, var, c, coor) !{{{1
  ! merge 2d array from other domains to mid and output
  character (len=*), intent(in) :: fname, varname
  character (len=1), intent(in) :: c
  real (kind=wp), intent(in) :: var(:), coor(:)

  if ( c.eq.'y' ) then
    call output_r1d_y (fname, varname, var, coor)
  else if ( c.eq.'x' ) then
    call output_r1d_x (fname, varname, var, coor)
  else
    stop 'unhandled character in quick_output_r1d in mod_mympi'
  end if

end subroutine quick_output_r1d

subroutine output_r1d_y (fname, varname, var, coor) !{{{1
  ! merge 1d array from other domains to mid and output
  character (len=*), intent(in) :: fname, varname
  real (kind=wp), intent(in) :: var(nj), coor(glo_nj)

  real (kind=wp), allocatable, dimension(:) :: glo_var

  integer, parameter :: tag = 140
  type (type_my) :: d
  integer :: n, leng

  if (myid == mid) then
    allocate( glo_var(glo_nj), stat=is)
    call chk(is); glo_var = 0.0

    do n = 1, npro
      d = our(n)
      leng = d%nj - 2
      if ( d%id == myid ) then
        glo_var(d%gs:d%gn) = var(2:d%nj-1)
      else
        call mpi_recv (glo_var(d%gs:d%gn), &
          leng, mpi_real8, d%id, tag, mpi_comm_world, msta, err)
      end if
    end do
    call io_quick_output (trim(fname), trim(varname), glo_var, coor)
    deallocate(glo_var)

  else
    leng = my%nj - 2
    call mpi_ssend (var(2:my%nj-1), leng, &
      mpi_real8, mid, tag, mpi_comm_world, err)
  end if

end subroutine output_r1d_y

subroutine output_r1d_x (fname, varname, var, coor) !{{{1
  ! merge 1d array from other domains to mid and output
  character (len=*), intent(in) :: fname, varname
  real (kind=wp), intent(in) :: var(ni), coor(glo_ni)

  real (kind=wp), allocatable, dimension(:) :: glo_var

  integer, parameter :: tag = 150
  type (type_my) :: d
  integer :: n, leng

  if (myid == mid) then
    allocate( glo_var(glo_ni), stat=is)
    call chk(is); glo_var = 0.0

    do n = 1, npro
      d = our(n)
      leng = d%ni - 2
      if ( d%id == myid ) then
        glo_var(d%gw:d%ge) = var(2:d%ni-1)
      else
        call mpi_recv (glo_var(d%gw:d%ge), &
          leng, mpi_real8, d%id, tag, mpi_comm_world, msta, err)
      end if
    end do
    call io_quick_output (trim(fname), trim(varname), glo_var, coor)
    deallocate(glo_var)

  else
    leng = my%ni - 2
    call mpi_ssend (var(2:my%ni-1), leng, &
      mpi_real8, mid, tag, mpi_comm_world, err)
  end if

end subroutine output_r1d_x

subroutine quick_output_i3d (fname, varname, var, lon, lat, z) !{{{1
  ! merge 2d array from other domains to mid and output
  character (len=*) :: fname, varname
  integer, intent(in) :: var(ni,nj,nk)
  real (kind=wp), dimension(:), intent(in) :: lon, lat, z

  integer, allocatable :: glo_var(:,:,:)

  integer, parameter :: tag = 50
  type (type_my) :: d
  integer :: n, leng

  if (myid == mid) then

    allocate( glo_var(glo_ni, glo_nj, nk), stat=is)
    call chk(is); glo_var = 0

    do n = 1, npro
      d = our(n)
      leng = (d%ge-d%gw+1) * (d%gn-d%gs+1) * nk
      if ( d%id == myid ) then
        glo_var(d%gw:d%ge, d%gs:d%gn, :) = var(2:my%ni-1, 2:my%nj-1, :)
      else
        call mpi_recv (glo_var(d%gw:d%ge,d%gs:d%gn,:), &
          leng, mpi_integer, d%id, tag, mpi_comm_world, msta, err)
      end if
    end do

    call io_quick_output (trim(fname), trim(varname), glo_var, lon, lat, z)
    deallocate(glo_var)

  else
    leng = (my%ge-my%gw+1) * (my%gn-my%gs+1) * nk
    call mpi_ssend (var(2:my%ni-1,2:my%nj-1,:), leng, &
      mpi_integer, mid, tag, mpi_comm_world, err)
  end if

end subroutine quick_output_i3d

subroutine quick_output_i2d (fname, varname, var, lon, lat) !{{{1
  ! merge 2d array from other domains to mid and output
  character (len=*) :: fname, varname
  integer, intent(in) :: var(ni,nj)
  real (kind=wp), intent(in) :: lon(glo_ni), lat(glo_nj)

  integer, allocatable :: glo_var(:,:)

  integer, parameter :: tag = 50
  type (type_my) :: d
  integer :: n, leng

  if (myid == mid) then

    allocate( glo_var(glo_ni, glo_nj), stat=is)
    call chk(is); glo_var = 0

    do n = 1, npro
      d = our(n)
      leng = (d%ge-d%gw+1) * (d%gn-d%gs+1)
      if ( d%id == myid ) then
        glo_var(d%gw:d%ge, d%gs:d%gn) = var(2:my%ni-1, 2:my%nj-1)
      else
        call mpi_recv (glo_var(d%gw:d%ge,d%gs:d%gn), &
          leng, mpi_integer, d%id, tag, mpi_comm_world, msta, err)
      end if
    end do

    call io_quick_output (trim(fname), trim(varname), glo_var, lon, lat)
    deallocate(glo_var)

  else
    leng = (my%ge-my%gw+1) * (my%gn-my%gs+1)
    call mpi_ssend (var(2:my%ni-1,2:my%nj-1), leng, &
      mpi_integer, mid, tag, mpi_comm_world, err)
  end if

end subroutine quick_output_i2d

subroutine swpbnd_r3d (var) !{{{1
  ! horizontally swap the boundary of subdomain for 3d array

  real (kind=wp) :: var(:,:,:)

  integer, parameter :: tag = 10
  integer :: leng, nd3

  nd3 = size(var, 3)
  leng = my%nj * nd3

  ! eastward
  call mpi_sendrecv (var(my%ni-1,:,:), leng, mpi_real8, &
    my%e, tag, var(1,:,:), leng, mpi_real8, &
    my%w, tag, mpi_comm_world, msta, err)

  ! westward
  call mpi_sendrecv (var(2,:,:), leng, mpi_real8, &
    my%w, tag, var(my%ni,:,:), leng, mpi_real8, &
    my%e, tag, mpi_comm_world, msta, err)

  leng = my%ni * nd3

  ! northward
  call mpi_sendrecv (var(:,my%nj-1,:), leng, mpi_real8, &
    my%n, tag, var(:,1,:), leng, mpi_real8, &
    my%s, tag, mpi_comm_world, msta, err)

  ! southward
  call mpi_sendrecv (var(:,2,:), leng, mpi_real8, &
    my%s, tag, var(:,my%nj,:), leng, mpi_real8, &
    my%n, tag, mpi_comm_world, msta, err)

end subroutine swpbnd_r3d

subroutine swpbnd_r2d_2 (v1, v2) !{{{1
  ! horizontally swap the boundary of subdomain for 2 2d variables

  real (kind=wp), dimension(:,:) :: v1, v2

  call swpbnd_r2d( v1 )
  call swpbnd_r2d( v2 )

end subroutine swpbnd_r2d_2

subroutine swpbnd_r2d (var) !{{{1
  ! horizontally swap the boundary of subdomain for 2d array

  real (kind=wp) :: var(ni,nj)

  integer, parameter :: tag = 110
  integer :: leng

  leng = my%nj

  ! eastward
  call mpi_sendrecv (var(my%ni-1,:), leng, mpi_real8, &
    my%e, tag, var(1,:), leng, mpi_real8, &
    my%w, tag, mpi_comm_world, msta, err)

  ! westward
  call mpi_sendrecv (var(2,:), leng, mpi_real8, &
    my%w, tag, var(my%ni,:), leng, mpi_real8, &
    my%e, tag, mpi_comm_world, msta, err)

  leng = my%ni

  ! northward
  call mpi_sendrecv (var(:,my%nj-1), leng, mpi_real8, &
    my%n, tag, var(:,1), leng, mpi_real8, &
    my%s, tag, mpi_comm_world, msta, err)

  ! southward
  call mpi_sendrecv (var(:,2), leng, mpi_real8, &
    my%s, tag, var(:,my%nj), leng, mpi_real8, &
    my%n, tag, mpi_comm_world, msta, err)

end subroutine swpbnd_r2d

subroutine swpbnd_r1d (var, c) !{{{1
  ! vertically swap the boundary of subdomain for 3d array

  real (kind=wp) :: var(:)
  character (len=1) :: c

  integer, parameter :: tag = 60
  integer, parameter :: leng = 1


  if ( c.eq.'x' ) then
    ! eastward
    call mpi_sendrecv (var(my%ni-1), leng, mpi_real8, &
      my%e, tag, var(1), leng, mpi_real8, &
      my%w, tag, mpi_comm_world, msta, err)

    ! westward
    call mpi_sendrecv (var(2), leng, mpi_real8, &
      my%w, tag, var(my%ni), leng, mpi_real8, &
      my%e, tag, mpi_comm_world, msta, err)

  else if ( c.eq.'y' ) then
    ! northward
    call mpi_sendrecv (var(my%nj-1), leng, mpi_real8, &
      my%n, tag, var(1), leng, mpi_real8, &
      my%s, tag, mpi_comm_world, msta, err)

    ! southward
    call mpi_sendrecv (var(2), leng, mpi_real8, &
      my%s, tag, var(my%nj), leng, mpi_real8, &
      my%n, tag, mpi_comm_world, msta, err)

  else
    stop 'unhandled direction indicator in swpbnd_r1d'
  end if

end subroutine swpbnd_r1d

subroutine swpbnd_i2d (var) !{{{1
  ! horizontally swap the boundary of subdomain for 2d array

  integer :: var(ni,nj)

  integer, parameter :: tag = 90
  integer :: leng

  leng = my%nj

  ! eastward
  call mpi_sendrecv (var(my%ni-1,:), leng, mpi_integer, &
    my%e, tag, var(1,:), leng, mpi_integer, &
    my%w, tag, mpi_comm_world, msta, err)

  ! westward
  call mpi_sendrecv (var(2,:), leng, mpi_integer, &
    my%w, tag, var(my%ni,:), leng, mpi_integer, &
    my%e, tag, mpi_comm_world, msta, err)

  leng = my%ni

  ! northward
  call mpi_sendrecv (var(:,my%nj-1), leng, mpi_integer, &
    my%n, tag, var(:,1), leng, mpi_integer, &
    my%s, tag, mpi_comm_world, msta, err)

  ! southward
  call mpi_sendrecv (var(:,2), leng, mpi_integer, &
    my%s, tag, var(:,my%nj), leng, mpi_integer, &
    my%n, tag, mpi_comm_world, msta, err)

end subroutine swpbnd_i2d

subroutine swpbnd_i3d (var) !{{{1
  ! horizontally swap the boundary of subdomain for 3d array

  integer :: var(ni,nj,nk)

  integer, parameter :: tag = 10
  integer :: leng

  leng = my%nj * nk

  ! eastward
  call mpi_sendrecv (var(my%ni-1,:,:), leng, mpi_integer, &
    my%e, tag, var(1,:,:), leng, mpi_integer, &
    my%w, tag, mpi_comm_world, msta, err)

  ! westward
  call mpi_sendrecv (var(2,:,:), leng, mpi_integer, &
    my%w, tag, var(my%ni,:,:), leng, mpi_integer, &
    my%e, tag, mpi_comm_world, msta, err)

  leng = my%ni * nk

  ! northward
  call mpi_sendrecv (var(:,my%nj-1,:), leng, mpi_integer, &
    my%n, tag, var(:,1,:), leng, mpi_integer, &
    my%s, tag, mpi_comm_world, msta, err)

  ! southward
  call mpi_sendrecv (var(:,2,:), leng, mpi_integer, &
    my%s, tag, var(:,my%nj,:), leng, mpi_integer, &
    my%n, tag, mpi_comm_world, msta, err)

end subroutine swpbnd_i3d

subroutine swpbnd_m3d (var) !{{{1
  ! horizontally swap the boundary of subdomain for 2d matrix

  type (type_mat) :: var(:,:,:)

  call swpbnd_r3d( var%x(1) )
  call swpbnd_r3d( var%x(2) )
end subroutine swpbnd_m3d

subroutine swpbnd_m2d (var) !{{{1
  ! horizontally swap the boundary of subdomain for 2d matrix

  type (type_mat), dimension(ni,nj) :: var

  call swpbnd_r2d( var%x(1) )
  call swpbnd_r2d( var%x(2) )
end subroutine swpbnd_m2d

subroutine swpbnd_m2d_4 (v1, v2, v3, v4) !{{{1
  ! horizontally swap the boundary of subdomain for 4 2d metrics

  type (type_mat), dimension(ni,nj) :: v1, v2, v3, v4

  call swpbnd_m2d( v1 )
  call swpbnd_m2d( v2 )
  call swpbnd_m2d( v3 )
  call swpbnd_m2d( v4 )
end subroutine swpbnd_m2d_4

subroutine div_gvar_m3d (va, vb) !{{{1
  ! divide va from mid to vb of all ids
  type (type_gvar_m3d), intent(in) :: va
  type (type_gvar_m3d) :: vb

  call div_r3d( va%x(1)%v, vb%x(1)%v )
  call div_r3d( va%x(2)%v, vb%x(2)%v )

end subroutine div_gvar_m3d

subroutine div_gvar_r3d (va, vb) !{{{1
  ! divide va from mid to vb of all ids
  type (type_gvar_r3d), intent(in) :: va
  type (type_gvar_r3d) :: vb

  call div_r3d( va%v, vb%v )

end subroutine div_gvar_r3d

subroutine div_r3d_gvar_r3d (va, vb) !{{{1
  ! divide va from mid to vb of all ids
  real (kind=wp), dimension(:,:,:) :: va
  type (type_gvar_r3d) :: vb

  call div_r3d( va, vb%v )

end subroutine div_r3d_gvar_r3d

subroutine div_m3d (va, vb) !{{{1
  ! divide va from mid to vb of all ids
  type (type_mat), intent(in) :: va(:,:,:)
  type (type_mat) :: vb(:,:,:)

  call div_r3d( va%x(1), vb%x(1) )
  call div_r3d( va%x(2), vb%x(2) )

end subroutine div_m3d

subroutine div_r3d (va, vb) !{{{1
  ! divide va from mid to vb of all ids
  real (kind=wp), intent(in) :: va(:,:,:)
  real (kind=wp) :: vb(:,:,:)

  integer, parameter :: tag = 20
  type (type_my) :: d
  integer :: n, leng, nd3, is

  nd3 = size(vb, 3)
  vb = 0.0

  if (myid == mid) then

    do n = 1, npro
      d = our(n)
      leng = (d%ni - 2) * (d%nj - 2) * nd3
      if ( d%id == mid ) then
        vb(2:my%ni-1, 2:my%nj-1, :) = &
          va(d%gw:d%ge, d%gs:d%gn, :)
      else
        call mpi_ssend (va(d%gw:d%ge,d%gs:d%gn,:), &
          leng, mpi_real8, d%id, tag, mpi_comm_world, err)
      end if
    end do

  else
    leng = (my%ni - 2) * (my%nj - 2) * nd3
    call mpi_recv (vb(2:my%ni-1,2:my%nj-1,:), &
      leng, mpi_real8, mid, tag, mpi_comm_world, msta, err)
  end if

  call swpbnd_r3d (vb)

end subroutine div_r3d

subroutine div_r2d (va, vb) !{{{1
  ! divide va from mid to vb of all ids

  real (kind=wp), dimension(glo_ni,glo_nj), intent(in) :: va
  real (kind=wp), dimension(ni,nj) :: vb

  integer, parameter :: tag = 90
  type (type_my) :: d
  integer :: n, leng

  vb = 0.0

  if (myid == mid) then
    do n = 1, npro
      d = our(n)
      leng = (d%ni - 2) * (d%nj - 2)
      if ( d%id == mid ) then
        vb(2:my%ni-1, 2:my%nj-1) = va(d%gw:d%ge, d%gs:d%gn)
      else
        call mpi_ssend (va(d%gw:d%ge,d%gs:d%gn), &
          leng, mpi_real8, d%id, tag, mpi_comm_world, err)
      end if
    end do

  else
    leng = (my%ni - 2) * (my%nj - 2)
    call mpi_recv (vb(2:my%ni-1,2:my%nj-1), leng, mpi_real8, &
      mid, tag, mpi_comm_world, msta, err)
  end if

  call swpbnd_r2d (vb)

end subroutine div_r2d

subroutine mympi_divy (va, vb) !{{{1
  ! divide 1d array north-south

  real (kind=wp), dimension(glo_nj), intent(in) :: va
  real (kind=wp), dimension(nj) :: vb

  integer, parameter :: tag = 50
  type (type_my) :: d
  integer :: n, leng

  vb = 0.0

  if (myid == mid) then
    do n = 1, npro
      d  = our(n)
      leng = d%nj - 2
      if ( d%id == mid ) then
        vb(2:d%nj-1) = va(d%gs:d%gn)
      else
        call mpi_ssend (va(d%gs:d%gn), leng, mpi_real8, d%id, &
          tag, mpi_comm_world, err)
      end if
    end do

  else
    leng = my%nj - 2
    call mpi_recv (vb(2:my%nj-1), leng, mpi_real8, &
      mid, tag, mpi_comm_world, msta, err)
  end if

  call swpbnd_r1d (vb, 'y')

end subroutine mympi_divy

subroutine mympi_divx (va, vb) !{{{1
  ! divide 1d array west-east

  real (kind=wp), dimension(glo_ni), intent(in) :: va
  real (kind=wp), dimension(ni) :: vb

  integer, parameter :: tag = 50
  type (type_my) :: d
  integer :: n, leng

  vb = 0.0

  if (myid == mid) then
    do n = 1, npro
      d  = our(n)
      leng = d%ni - 2
      if ( d%id == mid ) then
        vb(2:d%ni-1) = va(d%gw:d%ge)
      else
        call mpi_ssend (va(d%gw:d%ge), leng, mpi_real8, d%id, &
          tag, mpi_comm_world, err)
      end if
    end do

  else
    leng = my%ni - 2
    call mpi_recv (vb(2:my%ni-1), leng, mpi_real8, &
      mid, tag, mpi_comm_world, msta, err)
  end if

  call swpbnd_r1d (vb, 'x')

end subroutine mympi_divx

subroutine div_i3d (va, vb) !{{{1

  integer, intent(in) :: va(:,:,:)
  integer :: vb(:,:,:)

  integer, parameter :: tag = 20
  type (type_my) :: d
  integer :: n, leng

  vb = 0.0

  if (myid == mid) then

    do n = 1, npro
      d = our(n)
      leng = (d%ni - 2) * (d%nj - 2) * nk
      if ( d%id == mid ) then
        vb(2:my%ni-1, 2:my%nj-1, :) = &
          va(d%gw:d%ge, d%gs:d%gn, :)
      else
        call mpi_ssend (va(d%gw:d%ge,d%gs:d%gn,:), &
          leng, mpi_integer, d%id, tag, mpi_comm_world, err)
      end if
    end do

  else
    leng = (my%ni - 2) * (my%nj - 2) * nk
    call mpi_recv (vb(2:my%ni-1,2:my%nj-1,:), &
      leng, mpi_integer, mid, tag, mpi_comm_world, msta, err)
  end if

  call swpbnd_i3d (vb)

end subroutine div_i3d

subroutine div_i2d (va, vb) !{{{1
  ! distribute 2d array from mid to other domains

  integer, dimension(glo_ni,glo_nj), intent(in) :: va
  integer, dimension(ni,nj) :: vb

  integer, parameter :: tag = 90
  type (type_my) :: d
  integer :: n, leng

  vb = 0

  if (myid == mid) then
    do n = 1, npro
      d = our(n)
      leng = (d%ni - 2) * (d%nj - 2)
      if ( d%id == mid ) then
        vb(2:my%ni-1, 2:my%nj-1) = va(d%gw:d%ge, d%gs:d%gn)
      else
        call mpi_ssend (va(d%gw:d%ge,d%gs:d%gn), &
          leng, mpi_integer, d%id, tag, mpi_comm_world, err)
      end if
    end do

  else
    leng = (my%ni - 2) * (my%nj - 2)
    call mpi_recv (vb(2:my%ni-1,2:my%nj-1), leng, mpi_integer, &
      mid, tag, mpi_comm_world, msta, err)
  end if

  call swpbnd_i2d (vb)

end subroutine div_i2d

subroutine chk( ista ) !{{{1
  ! check state of allocate array 

  integer, intent(in) ::  ista

  if ( ista /= 0 ) then
    write(*,*) 'Allocate array failed. Stop'
    stop 2
  end if
end subroutine chk

end module mod_mympi !{{{1
!-------------------------------------------------------{{{1
! vim:fdm=marker:fdl=0:
! vim:foldtext=getline(v\:foldstart).'...'.(v\:foldend-v\:foldstart):