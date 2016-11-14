
! Description: common routines for programming
!
!      Author: OU Yuyuan <ouyuyuan@lasg.iap.ac.cn>
!     Created: 2015-09-14 17:36:48 BJT
! Last Change: 2016-01-24 10:55:30 BJT

module mod_pro

  use mod_kind, only: wp

  implicit none
  private

  public &
    pro_allo_chk, &
    pro_days_month, &
    pro_print, &
    pro_days_year
    
  interface pro_print
    module procedure print_i3d
    module procedure print_i2d
    module procedure print_r3d
    module procedure print_r2d
    module procedure print_r1d
  end interface

contains !{{{1

subroutine print_r3d (var) !{{{1
! print 3d reall array in a nice way
  
  real (kind=wp) :: var(:,:,:)
  integer :: ni, nj, nk, i, j, k

  ni = size(var, 1)
  nj = size(var, 2)
  nk = size(var, 3)

  do k = 1, nk
    do i = 1, ni
      write(*,'(100f5.1)') var(i,:,k)
    end do
    write(*,*) ''
  end do
  
end subroutine print_r3d

subroutine print_r2d (var, opt) !{{{1
! print 2d array in a nice way
  
  real (kind=wp) :: var(:,:)
  character (len=*), optional :: opt ! print boundary 

  integer :: ni, nj, i, j

  ni = size(var, 1)
  nj = size(var, 2)

  if ( present(opt) )then

    ! print east west boundary
    if ( opt .eq. 'ew' ) then
      do i = 2, ni - 1
        write(*,'(2e7.1e1, a, 2e7.1e1)') &
          var(i,1), var(i,2), ' . . . ', var(i,nj-1), var(i,nj)
      end do
    ! print north south
    else if ( opt .eq. 'ns' ) then
      write(*,'(a,i2, a, 100e7.1e1)') &
        'row = ', 1, ', ', var(1,2:nj-1)
      write(*,'(a,i2, a, 100e7.1e1)') &
      'row = ', 2, ', ', var(2,2:nj-1)
      write(*,'(a)') ' . . . . . . '
      write(*,'(a,i2, a, 100e7.1e1)') &
      'row = ', ni - 1, ', ', var(ni-1,2:nj-1)
      write(*,'(a,i2, a, 100e7.1e1)') &
      'row = ', ni, ', ', var(ni,2:nj-1)
      write(*,*) ''
    else
      write(*,*) 'unknow option '//opt//' in routine print_r2d in module mod_io'
      stop
    end if

  else

    do i = 1, ni
      write(*,'(100e7.1e1)') var(i,:)
    end do

  end if
  
end subroutine print_r2d

subroutine print_r1d (var) !{{{1
! print 1d array in a nice way
  
  real (kind=wp) :: var(:)

  write(*,'(100i5)') int(var(:))

end subroutine print_r1d

subroutine print_i3d (var) !{{{1
! print 3d integer array in a nice way
  
  integer :: var(:,:,:)
  integer :: ni, nj, nk, i, j, k

  ni = size(var, 1)
  nj = size(var, 2)
  nk = size(var, 3)

  do k = 1, nk
    do i = 1, ni
      write(*,"(100i3)") var(i,:,k)
    end do

    write (*,*) ''
  end do
  
end subroutine print_i3d

subroutine print_i2d (var, opt) !{{{1
! print 2d array in a nice way
  
  integer :: var(:,:)
  character (len=*), optional :: opt ! print boundary 

  integer :: ni, nj, i, j

  ni = size(var, 1)
  nj = size(var, 2)

  if ( present(opt) )then

    ! print east west boundary
    if ( opt .eq. 'ew' ) then
      do i = 2, ni - 1
        write(*,'(2i3, a, 2i3)') &
          var(i,1), var(i,2), ' . . . ', var(i,nj-1), var(i,nj)
      end do
    ! print north south
    else if ( opt .eq. 'ns' ) then
      write(*,'(a,i2, a, 100i3)') &
        ', row = ', 1, ', ', var(1,2:nj-1)
      write(*,'(a,i2, a, 100i3)') &
      ', row = ', 2, ', ', var(2,2:nj-1)
      write(*,'(a)') ' . . . . . . '
      write(*,'(a,i2, a, 100i3)') &
      ', row = ', ni - 1, ', ', var(ni-1,2:nj-1)
      write(*,'(a,i2, a, 100i3)') &
      ', row = ', ni, ', ', var(ni,2:nj-1)
      write(*,*) ''
    else
      write(*,*) 'unknow option '//opt//' in routine print_i2d in module mod_io'
      stop
    end if

  else

    do i = 1, ni
      write(*,'(100i3)') int(var(i,:))
    end do

  end if
  
end subroutine print_i2d


subroutine pro_allo_chk( ista ) !{{{1
  integer, intent(in) ::  ista

  if ( ista /= 0 ) then
    stop 'Allocate array failed.'
  end if
end subroutine pro_allo_chk


function pro_days_month (y, m) !{{{1
  ! calc. days of a specific month of a specific year
  integer, intent(in) :: y, m
  integer :: pro_days_month

  integer :: days(12)

!  days(:) = (/31, 28, 31, 30, &
!              31, 30, 31, 31, &
!              30, 31, 30, 31/)
  days(:) = 30

  pro_days_month = days(m)

!  if ( m == 2 ) then
!    if ( pro_days_year (y) == 366 ) pro_days_month = 29
!  end if

end function pro_days_month


function pro_days_year (y) !{{{1
  ! calc. days of a specific year
  integer, intent(in) :: y
  integer :: pro_days_year

  if (y <= 0) stop 'year should be positive in function pro_days_year'
  
  pro_days_year = 365

  if ( mod(y, 100) == 0 ) then
    if ( mod(y, 400) == 0 ) pro_days_year = 366
  else
    if ( mod(y, 4)   == 0 ) pro_days_year = 366
  end if

end function pro_days_year

end module mod_pro !{{{1
!-------------------------------------------------------{{{1
! vim:fdm=marker:fdl=0:
! vim:foldtext=getline(v\:foldstart).'...'.(v\:foldend-v\:foldstart):
