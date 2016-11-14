
! Description: kind parameters to control precision
!
!      Author: OU Yuyuan <ouyuyuan@lasg.iap.ac.cn>
!     Created: 2015-02-26 08:20:12 BJT
! Last Change: 2015-09-26 15:35:56 BJT

module mod_kind

  implicit none
  public

  !precision control for real number--------------------{{{1
  ! single precision
  integer, parameter :: sglp = selected_real_kind( 6,  37) 
  ! double precision
  integer, parameter :: dblp = selected_real_kind(12, 307) 
  integer, parameter :: wp = dblp ! working precision

end module mod_kind!{{{1
!-------------------------------------------------------{{{1
! vim:fdm=marker:fdl=0:
! vim:foldtext=getline(v\:foldstart).'...'.(v\:foldend-v\:foldstart):
