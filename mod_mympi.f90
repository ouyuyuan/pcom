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
! Last Change: 2018-01-22 09:33:59 BJT

module mod_mympi

  ! imported variables !{{{1
  use mod_arrays, only: equv, equvb, eqts, eqch, bnd, prho, am
  use mod_io, only: io_write, io_read, io_create_rst

  use mod_kind, only: wp

  use mpi

  use mod_param, only: mid, myid, npro, &
    glo_nj, glo_ni, ni, nj, nk, &
    my, our, missing_float, type_my, &
    vars_info, rst_info, nm, type_nm, &
    param_set_nm

  use mod_type, only: &
    type_mat, type_accu_gr2d, &
    type_var_info, type_gvar_r3d, &
    type_frc, type_bnd, &
    type_eq_ts, type_eq_uv, type_eq_uvb, type_eq_w, type_eq_ch, &
    type_rst_info, &
    type_check_date

  implicit none
  private

  ! interfaces !{{{1
  public &
    mympi_div, &
    mympi_divx, &
    mympi_divy, &
    mympi_swpbnd, &
    mympi_bcast, &
    mympi_merge, &
    mympi_input, &
    mympi_output

  interface mympi_bcast
    module procedure bcast_r1d
    module procedure bcast_i
    module procedure bcast_string
  end interface

  interface mympi_div
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
    module procedure swpbnd_m2d
  end interface

  interface mympi_output
    module procedure merge_out_rst
    module procedure merge_out_eqts
    module procedure merge_out_equvw
    module procedure merge_out_accu_gr2d_mask
    module procedure merge_out_accu_gr2d_weight
    module procedure merge_out_r3d
    module procedure merge_out_r2d_mask
    module procedure merge_out_r2d
  end interface

  interface mympi_input
    module procedure input_rst
    module procedure input_bnd
    module procedure input_eqts
    module procedure input_equv
    module procedure input_equvb
    module procedure input_eqch
    module procedure input_frc
    module procedure input_nm
    module procedure divide_in_r3d
    module procedure divide_in_r2d
    module procedure bcast_in_string
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

subroutine merge_out_rst (rst_info) ! {{{1
  ! output time-average variables of u,v eqation
  type (type_rst_info) :: rst_info

  if ( myid == mid ) call io_create_rst ( rst_info )

  call merge_out_r2d ( rst_info%chc, eqch%chc )
  call merge_out_r2d ( rst_info%chp, eqch%chp )

  call merge_out_r2d ( rst_info%bnd_taux, bnd%taux%v )
  call merge_out_r2d ( rst_info%bnd_tauy, bnd%tauy%v )
  call merge_out_r2d ( rst_info%bnd_t, bnd%t%v )
  call merge_out_r2d ( rst_info%bnd_s, bnd%s%v )
  call merge_out_r2d ( rst_info%bnd_pa, bnd%pa%v )
  call merge_out_r2d ( rst_info%bnd_fw, bnd%fw%v )

  call merge_out_r2d ( rst_info%buc, equvb%uc )
  call merge_out_r2d ( rst_info%bvc, equvb%vc )
  call merge_out_r2d ( rst_info%bup, equvb%up )
  call merge_out_r2d ( rst_info%bvp, equvb%vp )

  call merge_out_r3d ( rst_info%tc, eqts%tc )
  call merge_out_r3d ( rst_info%sc, eqts%sc )
  call merge_out_r3d ( rst_info%tp, eqts%tp )
  call merge_out_r3d ( rst_info%sp, eqts%sp )

  call merge_out_r3d ( rst_info%prpt, prho%x(1) )
  call merge_out_r3d ( rst_info%prps, prho%x(2) )

  call merge_out_r3d ( rst_info%uc, equv%uc )
  call merge_out_r3d ( rst_info%vc, equv%vc )
  call merge_out_r3d ( rst_info%up, equv%up )
  call merge_out_r3d ( rst_info%vp, equv%vp )

  call merge_out_r3d ( rst_info%auc, equv%auc )
  call merge_out_r3d ( rst_info%avc, equv%avc )
  call merge_out_r3d ( rst_info%aup, equv%aup )
  call merge_out_r3d ( rst_info%avp, equv%avp )
  call merge_out_r3d ( rst_info%aupp,equv%aupp )
  call merge_out_r3d ( rst_info%avpp,equv%avpp )

  call merge_out_r3d ( rst_info%am, am%v )

end subroutine merge_out_rst

subroutine merge_out_eqts (eqts) ! {{{1
  ! output time-average variables of t,s eqation
  type (type_eq_ts) :: eqts

  eqts%act = eqts%act / eqts%n
  eqts%acs = eqts%acs / eqts%n
  eqts%acr = eqts%acr / eqts%n

  where (eqts%g%msk == 0)
    eqts%acr = missing_float
  end where

  call merge_out_r3d (vars_info%pt, eqts%act)
  call merge_out_r3d (vars_info%sa, eqts%acs)
  call merge_out_r3d (vars_info%rho, eqts%acr)

  eqts%act = 0.0 ! reset accumulated value
  eqts%acs = 0.0
  eqts%acr = 0.0
  eqts%n   = 0   ! reset counter

end subroutine merge_out_eqts

subroutine merge_out_equvw (equv, eqw) ! {{{1
  ! output time-average variables of u,v eqation
  type (type_eq_uv) :: equv
  type (type_eq_w) :: eqw

  equv%acu = equv%acu / equv%n
  equv%acv = equv%acv / equv%n
  eqw%acw  = eqw%acw  / eqw%n

  where (equv%g%msk == 0)
    equv%acu = missing_float
    equv%acv = missing_float
    eqw%acw  = missing_float
  end where

  call merge_out_r3d (vars_info%u, equv%acu)
  call merge_out_r3d (vars_info%v, equv%acv)
  call merge_out_r3d (vars_info%w, eqw%acw)

  equv%acu = 0.0 ! reset accumulated value
  equv%acv = 0.0
  eqw%acw  = 0.0

  equv%n   = 0   ! reset counter
  eqw%n    = 0   

end subroutine merge_out_equvw

subroutine merge_out_accu_gr2d_mask (var_info, var, mask) !{{{1
  ! merge 2d array from other domains to mid
  ! mask out the land points
  type (type_var_info), intent(in) :: var_info
  type (type_accu_gr2d) :: var
  integer, dimension(:,:), intent(in) :: mask

  var%var%v = var%var%v / var%n
  where (mask == 0)
    var%var%v = missing_float
  end where

  call merge_out_r2d (var_info, var%var%v)

  var%var%v = 0.0 ! reset accumulated value
  var%n     = 0 ! reset counter

end subroutine merge_out_accu_gr2d_mask

subroutine merge_out_accu_gr2d_weight (var_info, var, wg, mask) !{{{1
  ! merge 2d array from other domains to mid
  ! mask out the land points, remove global mean
  type (type_var_info), intent(in) :: var_info
  type (type_accu_gr2d) :: var
  real (kind=wp), dimension(:,:), intent(in) :: wg
  integer, dimension(:,:), intent(in) :: mask

  var%var%v = var%var%v / var%n
  where (mask == 0)
    var%var%v = missing_float
  end where

  call merge_out_r2d_weight (var_info, var%var%v, wg)

  var%var%v = 0.0 ! reset accumulated value
  var%n     = 0 ! reset counter

end subroutine merge_out_accu_gr2d_weight

subroutine merge_out_r3d (var_info, var) !{{{1
  ! merge 3d array from other domains to mid
  type (type_var_info), intent(in) :: var_info
  real (kind=wp), dimension(ni,nj,nk), intent(in) :: var

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

    call io_write (var_info, glo_var)

    deallocate(glo_var)

  else
    leng  = (my%ge-my%gw+1) * (my%gn-my%gs+1) * nk
    call mpi_ssend (var(2:my%ni-1,2:my%nj-1,:), leng, &
      mpi_real8, mid, tag, mpi_comm_world, err)
  end if

end subroutine merge_out_r3d

subroutine merge_out_r2d_weight (var_info, var, wg) !{{{1
  ! merge 3d array from other domains to mid
  type (type_var_info), intent(in) :: var_info
  real (kind=wp), dimension(ni,nj), intent(in) :: var
  real (kind=wp), dimension(ni,nj), intent(in) :: wg

  real (kind=wp), allocatable, dimension(:,:) :: glo_var
  real (kind=wp), allocatable, dimension(:,:) :: glo_wg

  real (kind=wp) :: mean

  if (myid == mid) then
    allocate( glo_var(glo_ni, glo_nj), stat=is)
    call chk(is); glo_var = 0.0
    allocate( glo_wg(glo_ni, glo_nj), stat=is)
    call chk(is); glo_wg = 0.0

    call merge_r2d( glo_var, var )
    call merge_r2d( glo_wg, wg )

    mean = sum(glo_var*glo_wg)/sum(glo_wg)
    where (glo_var /= missing_float) glo_var = glo_var - mean

    call io_write (var_info, glo_var)

    deallocate(glo_var)
    deallocate(glo_wg)

  else

    call merge_r2d( glo_var, var )
    call merge_r2d( glo_wg, wg )

  end if

end subroutine merge_out_r2d_weight

subroutine merge_out_r2d (var_info, var) !{{{1
  ! merge 3d array from other domains to mid
  type (type_var_info), intent(in) :: var_info
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

    call io_write (var_info, glo_var)
    deallocate(glo_var)

  else
    leng  = (my%ge-my%gw+1) * (my%gn-my%gs+1)
    call mpi_ssend (var(2:my%ni-1,2:my%nj-1), leng, &
      mpi_real8, mid, tag, mpi_comm_world, err)
  end if

end subroutine merge_out_r2d

subroutine merge_out_r2d_mask(var_info, var, mask) !{{{1
  ! merge 3d array from other domains to mid
  type (type_var_info), intent(in) :: var_info
  real (kind=wp), dimension(ni,nj) :: var
  integer, dimension(:,:), intent(in) :: mask

  real (kind=wp), allocatable, dimension(:,:) :: glo_var

  integer, parameter :: tag = 30
  type (type_my) :: d
  integer :: n, leng

  where (mask == 0)
    var = missing_float
  end where

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

    call io_write (var_info, glo_var)
    deallocate(glo_var)

  else
    leng  = (my%ge-my%gw+1) * (my%gn-my%gs+1)
    call mpi_ssend (var(2:my%ni-1,2:my%nj-1), leng, &
      mpi_real8, mid, tag, mpi_comm_world, err)
  end if

end subroutine merge_out_r2d_mask


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

subroutine input_rst ( rst_info ) ! {{{1
  ! input restart variables 
  type ( type_rst_info ) :: rst_info

  call bcast_in_string (rst_info%fname, 'cdate', rst_info%cdate)
  call bcast_in_string (rst_info%fname, 'pdate', rst_info%pdate)

  call input_eqch  ( eqch )
  call input_equv  ( equv )
  call input_equvb ( equvb )
  
  call divide_in_r3d ( rst_info%fname, 'prpt', prho%x(1) )
  call divide_in_r3d ( rst_info%fname, 'prps', prho%x(2) )

  call divide_in_r3d ( rst_info%fname, 'am', am%v )

  call input_bnd ( bnd )

end subroutine input_rst

subroutine input_bnd ( bnd ) ! {{{1
  ! input upper boundary conditions for restart
  type (type_bnd) :: bnd

  call divide_in_r2d ( rst_info%fname , 'bnd_taux' , bnd%taux%v )
  call divide_in_r2d ( rst_info%fname , 'bnd_tauy' , bnd%tauy%v )
  call divide_in_r2d ( rst_info%fname , 'bnd_t'    , bnd%t%v )
  call divide_in_r2d ( rst_info%fname , 'bnd_s'    , bnd%s%v )
  call divide_in_r2d ( rst_info%fname , 'bnd_pa'   , bnd%pa%v )
  call divide_in_r2d ( rst_info%fname , 'bnd_fw'   , bnd%fw%v )

end subroutine input_bnd

subroutine input_eqts (eqts) ! {{{1
  ! input restart variables for t,s eqation
  type (type_eq_ts) :: eqts

  if (nm%rst == 1) then 
    call divide_in_r3d (rst_info%fname, 'tc', eqts%tc)
    call divide_in_r3d (rst_info%fname, 'sc', eqts%sc)
    call divide_in_r3d (rst_info%fname, 'tp', eqts%tp)
    call divide_in_r3d (rst_info%fname, 'sp', eqts%sp)
  else
    call divide_in_r3d (nm%fi, 'pt',  eqts%tc)
    call divide_in_r3d (nm%fi, 'sa',  eqts%sc)
    eqts%tp = eqts%tc
    eqts%sp = eqts%sc
  end if

end subroutine input_eqts

subroutine input_equv (equv) ! {{{1
  ! input restart variables for u,v eqations
  type (type_eq_uv) :: equv

  call divide_in_r3d (rst_info%fname, 'uc', equv%uc)
  call divide_in_r3d (rst_info%fname, 'vc', equv%vc)
  call divide_in_r3d (rst_info%fname, 'up', equv%up)
  call divide_in_r3d (rst_info%fname, 'vp', equv%vp)

  call divide_in_r3d (rst_info%fname, 'auc', equv%auc)
  call divide_in_r3d (rst_info%fname, 'avc', equv%avc)
  call divide_in_r3d (rst_info%fname, 'aup', equv%aup)
  call divide_in_r3d (rst_info%fname, 'avp', equv%avp)
  call divide_in_r3d (rst_info%fname, 'aupp', equv%aupp)
  call divide_in_r3d (rst_info%fname, 'avpp', equv%avpp)

end subroutine input_equv

subroutine input_equvb (equvb) ! {{{1
  ! input restart variables for u,v eqations
  type (type_eq_uvb) :: equvb

  call divide_in_r2d (rst_info%fname, 'buc', equvb%uc)
  call divide_in_r2d (rst_info%fname, 'bvc', equvb%vc)
  call divide_in_r2d (rst_info%fname, 'bup', equvb%up)
  call divide_in_r2d (rst_info%fname, 'bvp', equvb%vp)

end subroutine input_equvb

subroutine input_eqch (eqch) ! {{{1
  ! input restart variables for t,s eqation
  type (type_eq_ch) :: eqch

  call divide_in_r2d (rst_info%fname, 'chc', eqch%chc)
  call divide_in_r2d (rst_info%fname, 'chp', eqch%chp)

end subroutine input_eqch

subroutine input_frc (frc) ! {{{1
  ! input forcing fields
  type (type_frc) :: frc

  call divide_in_r3d (nm%ff, 'taux', frc%taux%v) 
  call divide_in_r3d (nm%ff, 'tauy', frc%tauy%v) 
  call divide_in_r3d (nm%ff, 'bct',  frc%t%v) 
  call divide_in_r3d (nm%ff, 'bcs',  frc%s%v) 
  call divide_in_r3d (nm%ff, 'pa',   frc%pa%v) 
  call divide_in_r3d (nm%ff, 'fw',   frc%fw%v) 

end subroutine input_frc

subroutine input_nm (nm) ! {{{1
  ! broad cast namelist from mid to all ids
  type (type_nm) :: nm

  ! although we can read in namelist in all processors, 
  ! but it seems the namelist file may be destroyed by 
  ! muliple-opening (happened once)
  if ( myid == mid ) then
    call param_set_nm (nm)
    call type_check_date(nm%bd)
    call type_check_date(nm%ed)
  end if

  call bcast_i (nm%px)
  call bcast_i (nm%py)

  call bcast_string (nm%bd)
  call bcast_string (nm%ed)

  call bcast_i (nm%bt)
  call bcast_i (nm%bc)

  call bcast_string (nm%fi)
  call bcast_string (nm%ff)
  call bcast_string (nm%ff_ws)
  call bcast_string (nm%od)

  call bcast_string  (nm%out_per)
  call bcast_i (nm%rst_per)
  call bcast_i (nm%rst)

end subroutine input_nm

subroutine divide_in_r3d (fname, varname, var) !{{{1
  ! read in global 3d array and divide it to other domains from mid
  character (len=*), intent(in) :: fname, varname
  real (kind=wp), dimension(:,:,:)  :: var

  real (kind=wp), allocatable, dimension(:,:,:) :: glo_var
  integer :: nd3

  nd3 = size(var, 3)

  if (myid == mid) then
    allocate( glo_var(glo_ni, glo_nj, nd3), stat=is)
    call chk(is); glo_var = 0.0
    call io_read (fname, varname, glo_var)
  end if

  call div_r3d (glo_var, var) 

  if (myid == mid) deallocate(glo_var)

end subroutine divide_in_r3d

subroutine divide_in_r2d (fname, varname, var) !{{{1
  ! read in global 3d array and divide it to other domains from mid
  character (len=*), intent(in) :: fname, varname
  real (kind=wp), dimension(ni,nj) :: var

  real (kind=wp), allocatable, dimension(:,:) :: glo_var

  if (myid == mid) then
    allocate( glo_var(glo_ni, glo_nj), stat=is)
    call chk(is); glo_var = 0.0
    call io_read (fname, varname, glo_var)
  end if

  call div_r2d (glo_var, var) 

  if (myid == mid) deallocate(glo_var)

end subroutine divide_in_r2d

subroutine bcast_in_string (fname, varname, var) !{{{1
  ! read in string and broadcast it to other domains from mid
  character (len=*), intent(in) :: fname, varname
  character (len=*) :: var

  if (myid == mid) call io_read (fname, varname, var)
  call bcast_string (var)

end subroutine bcast_in_string

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

subroutine swpbnd_m2d (var) !{{{1
  ! horizontally swap the boundary of subdomain for 2d matrix

  type (type_mat), dimension(ni,nj) :: var

  call swpbnd_r2d( var%x(1) )
  call swpbnd_r2d( var%x(2) )
end subroutine swpbnd_m2d

subroutine div_r3d (va, vb) !{{{1
  ! divide va from mid to vb of all ids
  real (kind=wp), intent(in) :: va(:,:,:)
  real (kind=wp) :: vb(:,:,:)

  integer, parameter :: tag = 20
  type (type_my) :: d
  integer :: n, leng, nd3

  nd3 = size(vb, 3)

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

!  vb = 0

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
