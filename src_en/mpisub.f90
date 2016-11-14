!     =======================
      subroutine swap_array_real1d(mpiarray,jmt,north,south)
!     =======================
      implicit none
      include 'mpif.h'
      integer jmt
      integer ierr,north,south,status(mpi_status_size)
      real mpiarray(jmt)
!     send data from south to north
      call mpi_sendrecv(mpiarray(jmt-1),1,mpi_real8,north,10,  &
                        mpiarray(1),1,mpi_real8,south,10,  &
                        mpi_comm_world,status,ierr)
!     send data from north to south 
      call mpi_sendrecv(mpiarray(2),1,mpi_real8,south,11,  &
                        mpiarray(jmt),1,mpi_real8,north,11,  &
                        mpi_comm_world,status,ierr)
      return
      end

!     =======================
      subroutine swap_array_int2d(mpiarray,imt,jmt,west,east,north,south)
!     =======================
      implicit none
      include 'mpif.h'
      integer imt,jmt,i,j,swapsize
      integer ierr,west,east,north,south,status(mpi_status_size)
      integer mpiarray(imt,jmt)
      integer sendtemp(jmt),recvtemp(jmt)
      integer sendtemp2(imt),recvtemp2(imt)
      swapsize=jmt
!     send data from west to east
      do j=1,jmt
      sendtemp(j)=mpiarray(imt-1,j)
      recvtemp(j)=mpiarray(1,j)
      end do
      call mpi_sendrecv(sendtemp(1),swapsize,mpi_integer,east,10,  &
                        recvtemp(1),swapsize,mpi_integer,west,10,  &
                        mpi_comm_world,status,ierr)
      do j=1,jmt
      mpiarray(1,j)=recvtemp(j)
      end do
!     send data from east to west
      do j=1,jmt
      sendtemp(j)=mpiarray(2,j)
      recvtemp(j)=mpiarray(imt,j)
      end do
      call mpi_sendrecv(sendtemp(1),swapsize,mpi_integer,west,11,  &
                        recvtemp(1),swapsize,mpi_integer,east,11,  &
                        mpi_comm_world,status,ierr)
      do j=1,jmt
      mpiarray(imt,j)=recvtemp(j)
      end do
!     send data from south to north
      swapsize=imt
      do i=1,imt
      sendtemp2(i)=mpiarray(i,jmt-1)
      recvtemp2(i)=mpiarray(i,1)
      end do
      call mpi_sendrecv(sendtemp2(1),swapsize,mpi_integer,north,12,  &
                        recvtemp2(1),swapsize,mpi_integer,south,12,  &
                        mpi_comm_world,status,ierr)
      do i=1,imt
      mpiarray(i,1)=recvtemp2(i)
      end do
!     send data from north to south 
      do i=1,imt
      sendtemp2(i)=mpiarray(i,2)
      recvtemp2(i)=mpiarray(i,jmt)
      end do
      call mpi_sendrecv(sendtemp2(1),swapsize,mpi_integer,south,13,  &
                        recvtemp2(1),swapsize,mpi_integer,north,13,  &
                        mpi_comm_world,status,ierr)
      do i=1,imt
      mpiarray(i,jmt)=recvtemp2(i)
      end do
      return
      end

!     =======================
      subroutine swap_array_real2d(mpiarray,imt,jmt,west,east,north,south)
!     =======================
      implicit none
      include 'mpif.h'
      integer imt,jmt,i,j,swapsize
      integer ierr,west,east,north,south,status(mpi_status_size)
      real mpiarray(imt,jmt)
      real sendtemp(jmt),recvtemp(jmt)
      real sendtemp2(imt),recvtemp2(imt)
      swapsize=jmt
!     send data from west to east
      do j=1,jmt
      sendtemp(j)=mpiarray(imt-1,j)
      recvtemp(j)=mpiarray(1,j)
      end do
      call mpi_sendrecv(sendtemp(1),swapsize,mpi_real8,east,10,  &
                        recvtemp(1),swapsize,mpi_real8,west,10,  &
                        mpi_comm_world,status,ierr)
      do j=1,jmt
      mpiarray(1,j)=recvtemp(j)
      end do
!     send data from east to west
      do j=1,jmt
      sendtemp(j)=mpiarray(2,j)
      recvtemp(j)=mpiarray(imt,j)
      end do
      call mpi_sendrecv(sendtemp(1),swapsize,mpi_real8,west,11,  &
                        recvtemp(1),swapsize,mpi_real8,east,11,  &
                        mpi_comm_world,status,ierr)
      do j=1,jmt
      mpiarray(imt,j)=recvtemp(j)
      end do
!     send data from south to north
      swapsize=imt
      do i=1,imt
      sendtemp2(i)=mpiarray(i,jmt-1)
      recvtemp2(i)=mpiarray(i,1)
      end do
      call mpi_sendrecv(sendtemp2(1),swapsize,mpi_real8,north,12,  &
                        recvtemp2(1),swapsize,mpi_real8,south,12,  &
                        mpi_comm_world,status,ierr)
      do i=1,imt
      mpiarray(i,1)=recvtemp2(i)
      end do
!     send data from north to south 
      do i=1,imt
      sendtemp2(i)=mpiarray(i,2)
      recvtemp2(i)=mpiarray(i,jmt)
      end do
      call mpi_sendrecv(sendtemp2(1),swapsize,mpi_real8,south,13,  &
                        recvtemp2(1),swapsize,mpi_real8,north,13,  &
                        mpi_comm_world,status,ierr)
      do i=1,imt
      mpiarray(i,jmt)=recvtemp2(i)
      end do
      return
      end

!     =======================
      subroutine swap_array_real3d(mpiarray,imt,jmt,km,west,east,north,south)
!     =======================
      implicit none
      include 'mpif.h'
      integer imt,jmt,km,i,j,k,swapsize
      integer ierr,west,east,north,south,status(mpi_status_size)
      real mpiarray(imt,jmt,km)
      real sendtemp(jmt,km),recvtemp(jmt,km)
      real sendtemp2(imt,km),recvtemp2(imt,km)
      swapsize=jmt*km
!     send data from west to east
      do k=1,km
      do j=1,jmt
      sendtemp(j,k)=mpiarray(imt-1,j,k)
      recvtemp(j,k)=mpiarray(1,j,k)
      end do
      end do
      call mpi_sendrecv(sendtemp(1,1),swapsize,mpi_real8,east,10,  &
                        recvtemp(1,1),swapsize,mpi_real8,west,10,  &
                        mpi_comm_world,status,ierr)
      do k=1,km
      do j=1,jmt
      mpiarray(1,j,k)=recvtemp(j,k)
      end do
      end do
!     send data from east to west
      do k=1,km
      do j=1,jmt
      sendtemp(j,k)=mpiarray(2,j,k)
      recvtemp(j,k)=mpiarray(imt,j,k)
      end do
      end do
      call mpi_sendrecv(sendtemp(1,1),swapsize,mpi_real8,west,11,  &
                        recvtemp(1,1),swapsize,mpi_real8,east,11,  &
                        mpi_comm_world,status,ierr)
      do k=1,km
      do j=1,jmt
      mpiarray(imt,j,k)=recvtemp(j,k)
      end do
      end do
!     send data from south to north
      swapsize=imt*km
      do k=1,km
      do i=1,imt
      sendtemp2(i,k)=mpiarray(i,jmt-1,k)
      recvtemp2(i,k)=mpiarray(i,1,k)
      end do
      end do
      call mpi_sendrecv(sendtemp2(1,1),swapsize,mpi_real8,north,12,  &
                        recvtemp2(1,1),swapsize,mpi_real8,south,12,  &
                        mpi_comm_world,status,ierr)
      do k=1,km
      do i=1,imt
      mpiarray(i,1,k)=recvtemp2(i,k)
      end do
      end do
!     send data from north to south 
      do k=1,km
      do i=1,imt
      sendtemp2(i,k)=mpiarray(i,2,k)
      recvtemp2(i,k)=mpiarray(i,jmt,k)
      end do
      end do
      call mpi_sendrecv(sendtemp2(1,1),swapsize,mpi_real8,south,13,  &
                        recvtemp2(1,1),swapsize,mpi_real8,north,13,  &
                        mpi_comm_world,status,ierr)
      do k=1,km
      do i=1,imt
      mpiarray(i,jmt,k)=recvtemp2(i,k)
      end do
      end do
      return
      end

!     =======================
      subroutine swap_isopy_real3d(mpiarray,imt,jmt,km,west,east,north,south)
!     =======================
      implicit none
      include 'mpif.h'
      integer imt,jmt,km,i,j,k,swapsize
      integer ierr,west,east,north,south,status(mpi_status_size)
      real mpiarray(imt,km,jmt)
      real sendtemp(jmt,km),recvtemp(jmt,km)
      real sendtemp2(imt,km),recvtemp2(imt,km)
      swapsize=jmt*km
!     send data from west to east
      do k=1,km
      do j=1,jmt
      sendtemp(j,k)=mpiarray(imt-1,k,j)
      recvtemp(j,k)=mpiarray(1,k,j)
      end do
      end do
      call mpi_sendrecv(sendtemp(1,1),swapsize,mpi_real8,east,10,  &
                        recvtemp(1,1),swapsize,mpi_real8,west,10,  &
                        mpi_comm_world,status,ierr)
      do k=1,km
      do j=1,jmt
      mpiarray(1,k,j)=recvtemp(j,k)
      end do
      end do
!     send data from east to west
      do k=1,km
      do j=1,jmt
      sendtemp(j,k)=mpiarray(2,k,j)
      recvtemp(j,k)=mpiarray(imt,k,j)
      end do
      end do
      call mpi_sendrecv(sendtemp(1,1),swapsize,mpi_real8,west,11,  &
                        recvtemp(1,1),swapsize,mpi_real8,east,11,  &
                        mpi_comm_world,status,ierr)
      do k=1,km
      do j=1,jmt
      mpiarray(imt,k,j)=recvtemp(j,k)
      end do
      end do
!     send data from south to north
      swapsize=imt*km
      do k=1,km
      do i=1,imt
      sendtemp2(i,k)=mpiarray(i,k,jmt-1)
      recvtemp2(i,k)=mpiarray(i,k,1)
      end do
      end do
      call mpi_sendrecv(sendtemp2(1,1),swapsize,mpi_real8,north,12,  &
                        recvtemp2(1,1),swapsize,mpi_real8,south,12,  &
                        mpi_comm_world,status,ierr)
      do k=1,km
      do i=1,imt
      mpiarray(i,k,1)=recvtemp2(i,k)
      end do
      end do
!     send data from north to south 
      do k=1,km
      do i=1,imt
      sendtemp2(i,k)=mpiarray(i,k,2)
      recvtemp2(i,k)=mpiarray(i,k,jmt)
      end do
      end do
      call mpi_sendrecv(sendtemp2(1,1),swapsize,mpi_real8,south,13,  &
                        recvtemp2(1,1),swapsize,mpi_real8,north,13,  &
                        mpi_comm_world,status,ierr)
      do k=1,km
      do i=1,imt
      mpiarray(i,k,jmt)=recvtemp2(i,k)
      end do
      end do
      return
      end

!     =======================
      subroutine swap_isopyw_real3d(mpiarray,imt,jmt,km,west,east,north,south)
!     =======================
      implicit none
      include 'mpif.h'
      integer imt,jmt,km,i,j,k,swapsize
      integer ierr,west,east,north,south,status(mpi_status_size)
      real mpiarray(imt,0:km,jmt)
      real sendtemp(jmt,km),recvtemp(jmt,km)
      real sendtemp2(imt,km),recvtemp2(imt,km)
      swapsize=jmt*km
!     send data from west to east
      do k=1,km
      do j=1,jmt
      sendtemp(j,k)=mpiarray(imt-1,k,j)
      recvtemp(j,k)=mpiarray(1,k,j)
      end do
      end do
      call mpi_sendrecv(sendtemp(1,1),swapsize,mpi_real8,east,10,  &
                        recvtemp(1,1),swapsize,mpi_real8,west,10,  &
                        mpi_comm_world,status,ierr)
      do k=1,km
      do j=1,jmt
      mpiarray(1,k,j)=recvtemp(j,k)
      end do
      end do
!     send data from east to west
      do k=1,km
      do j=1,jmt
      sendtemp(j,k)=mpiarray(2,k,j)
      recvtemp(j,k)=mpiarray(imt,k,j)
      end do
      end do
      call mpi_sendrecv(sendtemp(1,1),swapsize,mpi_real8,west,11,  &
                        recvtemp(1,1),swapsize,mpi_real8,east,11,  &
                        mpi_comm_world,status,ierr)
      do k=1,km
      do j=1,jmt
      mpiarray(imt,k,j)=recvtemp(j,k)
      end do
      end do
!     send data from south to north
      swapsize=imt*km
      do k=1,km
      do i=1,imt
      sendtemp2(i,k)=mpiarray(i,k,jmt-1)
      recvtemp2(i,k)=mpiarray(i,k,1)
      end do
      end do
      call mpi_sendrecv(sendtemp2(1,1),swapsize,mpi_real8,north,12,  &
                        recvtemp2(1,1),swapsize,mpi_real8,south,12,  &
                        mpi_comm_world,status,ierr)
      do k=1,km
      do i=1,imt
      mpiarray(i,k,1)=recvtemp2(i,k)
      end do
      end do
!     send data from north to south 
      do k=1,km
      do i=1,imt
      sendtemp2(i,k)=mpiarray(i,k,2)
      recvtemp2(i,k)=mpiarray(i,k,jmt)
      end do
      end do
      call mpi_sendrecv(sendtemp2(1,1),swapsize,mpi_real8,south,13,  &
                        recvtemp2(1,1),swapsize,mpi_real8,north,13,  &
                        mpi_comm_world,status,ierr)
      do k=1,km
      do i=1,imt
      mpiarray(i,k,jmt)=recvtemp2(i,k)
      end do
      end do
      return
      end

!     =======================
      subroutine swap_ew_real3d(mpiarray,imt,jmt,km,west,east)
!     =======================
      implicit none
      include 'mpif.h'
      integer imt,jmt,km,j,k,swapsize
      integer ierr,west,east,status(mpi_status_size)
      real mpiarray(imt,jmt,km)
      real sendtemp(jmt,km),recvtemp(jmt,km)
      swapsize=jmt*km
!     send data from west to east
      do k=1,km
      do j=1,jmt
      sendtemp(j,k)=mpiarray(imt-1,j,k)
      recvtemp(j,k)=mpiarray(1,j,k)
      end do
      end do
      call mpi_sendrecv(sendtemp(1,1),swapsize,mpi_real8,east,10,  &
                        recvtemp(1,1),swapsize,mpi_real8,west,10,  &
                        mpi_comm_world,status,ierr)
      do k=1,km
      do j=1,jmt
      mpiarray(1,j,k)=recvtemp(j,k)
      end do
      end do
!     send data from east to west
      do k=1,km
      do j=1,jmt
      sendtemp(j,k)=mpiarray(2,j,k)
      recvtemp(j,k)=mpiarray(imt,j,k)
      end do
      end do
      call mpi_sendrecv(sendtemp(1,1),swapsize,mpi_real8,west,11,  &
                        recvtemp(1,1),swapsize,mpi_real8,east,11,  &
                        mpi_comm_world,status,ierr)
      do k=1,km
      do j=1,jmt
      mpiarray(imt,j,k)=recvtemp(j,k)
      end do
      end do
      return
      end
      
!     =======================
      subroutine swap_ns_real3d(mpiarray,imt,jmt,km,north,south)
!     =======================
      implicit none
      include 'mpif.h'
      integer imt,jmt,km,i,k,swapsize
      integer ierr,north,south,status(mpi_status_size)
      real mpiarray(imt,jmt,km)
      real sendtemp2(imt,km),recvtemp2(imt,km)
!     send data from south to north
      swapsize=imt*km
      do k=1,km
      do i=1,imt
      sendtemp2(i,k)=mpiarray(i,jmt-1,k)
      recvtemp2(i,k)=mpiarray(i,1,k)
      end do
      end do
      call mpi_sendrecv(sendtemp2(1,1),swapsize,mpi_real8,north,12,  &
                        recvtemp2(1,1),swapsize,mpi_real8,south,12,  &
                        mpi_comm_world,status,ierr)
      do k=1,km
      do i=1,imt
      mpiarray(i,1,k)=recvtemp2(i,k)
      end do
      end do
!     send data from north to south 
      do k=1,km
      do i=1,imt
      sendtemp2(i,k)=mpiarray(i,2,k)
      recvtemp2(i,k)=mpiarray(i,jmt,k)
      end do
      end do
      call mpi_sendrecv(sendtemp2(1,1),swapsize,mpi_real8,south,13,  &
                        recvtemp2(1,1),swapsize,mpi_real8,north,13,  &
                        mpi_comm_world,status,ierr)
      do k=1,km
      do i=1,imt
      mpiarray(i,jmt,k)=recvtemp2(i,k)
      end do
      end do
      return
      end

!     =======================
      subroutine swap_array_real4d(mpiarray,imt,jmt,km,km2,west,east,north,south)
!     =======================
      implicit none
      include 'mpif.h'
      integer imt,jmt,km,km2,i,j,k,k2,swapsize
      integer ierr,west,east,north,south,status(mpi_status_size)
      real mpiarray(imt,jmt,km,km2)
      real sendtemp(jmt,km,km2),recvtemp(jmt,km,km2)
      real sendtemp2(imt,km,km2),recvtemp2(imt,km,km2)
      swapsize=jmt*km*km2
!     send data from west to east
      do k2=1,km2
      do k=1,km
      do j=1,jmt
      sendtemp(j,k,k2)=mpiarray(imt-1,j,k,k2)
      recvtemp(j,k,k2)=mpiarray(1,j,k,k2)
      end do
      end do
      end do
      call mpi_sendrecv(sendtemp(1,1,1),swapsize,mpi_real8,east,10,  &
                        recvtemp(1,1,1),swapsize,mpi_real8,west,10,  &
                        mpi_comm_world,status,ierr)
      do k2=1,km2
      do k=1,km
      do j=1,jmt
      mpiarray(1,j,k,k2)=recvtemp(j,k,k2)
      end do
      end do
      end do
!     send data from east to west
      do k2=1,km2
      do k=1,km
      do j=1,jmt
      sendtemp(j,k,k2)=mpiarray(2,j,k,k2)
      recvtemp(j,k,k2)=mpiarray(imt,j,k,k2)
      end do
      end do
      end do
      call mpi_sendrecv(sendtemp(1,1,1),swapsize,mpi_real8,west,11,  &
                        recvtemp(1,1,1),swapsize,mpi_real8,east,11,  &
                        mpi_comm_world,status,ierr)
      do k2=1,km2
      do k=1,km
      do j=1,jmt
      mpiarray(imt,j,k,k2)=recvtemp(j,k,k2)
      end do
      end do
      end do
!     send data from south to north
      swapsize=imt*km*km2
      do k2=1,km2
      do k=1,km
      do i=1,imt
      sendtemp2(i,k,k2)=mpiarray(i,jmt-1,k,k2)
      recvtemp2(i,k,k2)=mpiarray(i,1,k,k2)
      end do
      end do
      end do
      call mpi_sendrecv(sendtemp2(1,1,1),swapsize,mpi_real8,north,12,  &
                        recvtemp2(1,1,1),swapsize,mpi_real8,south,12,  &
                        mpi_comm_world,status,ierr)
      do k2=1,km2
      do k=1,km
      do i=1,imt
      mpiarray(i,1,k,k2)=recvtemp2(i,k,k2)
      end do
      end do
      end do
!     send data from north to south 
      do k2=1,km2
      do k=1,km
      do i=1,imt
      sendtemp2(i,k,k2)=mpiarray(i,2,k,k2)
      recvtemp2(i,k,k2)=mpiarray(i,jmt,k,k2)
      end do
      end do
      end do
      call mpi_sendrecv(sendtemp2(1,1,1),swapsize,mpi_real8,south,13,  &
                        recvtemp2(1,1,1),swapsize,mpi_real8,north,13,  &
                        mpi_comm_world,status,ierr)
      do k2=1,km2
      do k=1,km
      do i=1,imt
      mpiarray(i,jmt,k,k2)=recvtemp2(i,k,k2)
      end do
      end do
      end do
      return
      end

!     =======================
      subroutine swap_isopy_real4d(mpiarray,imt,jmt,km,km2,west,east,north,south)
!     =======================
      implicit none
      include 'mpif.h'
      integer imt,jmt,km,km2,i,j,k,k2,swapsize
      integer ierr,west,east,north,south,status(mpi_status_size)
      real mpiarray(imt,km,jmt,km2)
      real sendtemp(jmt,km,km2),recvtemp(jmt,km,km2)
      real sendtemp2(imt,km,km2),recvtemp2(imt,km,km2)
      swapsize=jmt*km*km2
!     send data from west to east
      do k2=1,km2
      do k=1,km
      do j=1,jmt
      sendtemp(j,k,k2)=mpiarray(imt-1,k,j,k2)
      recvtemp(j,k,k2)=mpiarray(1,k,j,k2)
      end do
      end do
      end do
      call mpi_sendrecv(sendtemp(1,1,1),swapsize,mpi_real8,east,10,  &
                        recvtemp(1,1,1),swapsize,mpi_real8,west,10,  &
                        mpi_comm_world,status,ierr)
      do k2=1,km2
      do k=1,km
      do j=1,jmt
      mpiarray(1,k,j,k2)=recvtemp(j,k,k2)
      end do
      end do
      end do
!     send data from east to west
      do k2=1,km2
      do k=1,km
      do j=1,jmt
      sendtemp(j,k,k2)=mpiarray(2,k,j,k2)
      recvtemp(j,k,k2)=mpiarray(imt,k,j,k2)
      end do
      end do
      end do
      call mpi_sendrecv(sendtemp(1,1,1),swapsize,mpi_real8,west,11,  &
                        recvtemp(1,1,1),swapsize,mpi_real8,east,11,  &
                        mpi_comm_world,status,ierr)
      do k2=1,km2
      do k=1,km
      do j=1,jmt
      mpiarray(imt,k,j,k2)=recvtemp(j,k,k2)
      end do
      end do
      end do
!     send data from south to north
      swapsize=imt*km*km2
      do k2=1,km2
      do k=1,km
      do i=1,imt
      sendtemp2(i,k,k2)=mpiarray(i,k,jmt-1,k2)
      recvtemp2(i,k,k2)=mpiarray(i,k,1,k2)
      end do
      end do
      end do
      call mpi_sendrecv(sendtemp2(1,1,1),swapsize,mpi_real8,north,12,  &
                        recvtemp2(1,1,1),swapsize,mpi_real8,south,12,  &
                        mpi_comm_world,status,ierr)
      do k2=1,km2
      do k=1,km
      do i=1,imt
      mpiarray(i,k,1,k2)=recvtemp2(i,k,k2)
      end do
      end do
      end do
!     send data from north to south 
      do k2=1,km2
      do k=1,km
      do i=1,imt
      sendtemp2(i,k,k2)=mpiarray(i,k,2,k2)
      recvtemp2(i,k,k2)=mpiarray(i,k,jmt,k2)
      end do
      end do
      end do
      call mpi_sendrecv(sendtemp2(1,1,1),swapsize,mpi_real8,south,13,  &
                        recvtemp2(1,1,1),swapsize,mpi_real8,north,13,  &
                        mpi_comm_world,status,ierr)
      do k2=1,km2
      do k=1,km
      do i=1,imt
      mpiarray(i,k,jmt,k2)=recvtemp2(i,k,k2)
      end do
      end do
      end do
      return
      end

!     =======================
      subroutine swap_array_real5d(mpiarray,imt,jmt,km,km2,km3,west,east,north,south)
!     =======================
      implicit none
      include 'mpif.h'
      integer imt,jmt,km,km2,km3,i,j,k,k2,k3,swapsize
      integer ierr,west,east,north,south,status(mpi_status_size)
      real mpiarray(imt,jmt,km,km2,km3)
      real sendtemp(jmt,km,km2,km3),recvtemp(jmt,km,km2,km3)
      real sendtemp2(imt,km,km2,km3),recvtemp2(imt,km,km2,km3)
      swapsize=jmt*km*km2*km3
!     send data from west to east
      do k3=1,km3
      do k2=1,km2
      do k=1,km
      do j=1,jmt
      sendtemp(j,k,k2,k3)=mpiarray(imt-1,j,k,k2,k3)
      recvtemp(j,k,k2,k3)=mpiarray(1,j,k,k2,k3)
      end do
      end do
      end do
      end do
      call mpi_sendrecv(sendtemp(1,1,1,1),swapsize,mpi_real8,east,10,  &
                        recvtemp(1,1,1,1),swapsize,mpi_real8,west,10,  &
                        mpi_comm_world,status,ierr)
      do k3=1,km3
      do k2=1,km2
      do k=1,km
      do j=1,jmt
      mpiarray(1,j,k,k2,k3)=recvtemp(j,k,k2,k3)
      end do
      end do
      end do
      end do
!     send data from east to west
      do k3=1,km3
      do k2=1,km2
      do k=1,km
      do j=1,jmt
      sendtemp(j,k,k2,k3)=mpiarray(2,j,k,k2,k3)
      recvtemp(j,k,k2,k3)=mpiarray(imt,j,k,k2,k3)
      end do
      end do
      end do
      end do
      call mpi_sendrecv(sendtemp(1,1,1,1),swapsize,mpi_real8,west,11,  &
                        recvtemp(1,1,1,1),swapsize,mpi_real8,east,11,  &
                        mpi_comm_world,status,ierr)
      do k3=1,km3
      do k2=1,km2
      do k=1,km
      do j=1,jmt
      mpiarray(imt,j,k,k2,k3)=recvtemp(j,k,k2,k3)
      end do
      end do
      end do
      end do
!     send data from south to north
      swapsize=imt*km*km2*km3
      do k3=1,km3
      do k2=1,km2
      do k=1,km
      do i=1,imt
      sendtemp2(i,k,k2,k3)=mpiarray(i,jmt-1,k,k2,k3)
      recvtemp2(i,k,k2,k3)=mpiarray(i,1,k,k2,k3)
      end do
      end do
      end do
      end do
      call mpi_sendrecv(sendtemp2(1,1,1,1),swapsize,mpi_real8,north,12,  &
                        recvtemp2(1,1,1,1),swapsize,mpi_real8,south,12,  &
                        mpi_comm_world,status,ierr)
      do k3=1,km3
      do k2=1,km2
      do k=1,km
      do i=1,imt
      mpiarray(i,1,k,k2,k3)=recvtemp2(i,k,k2,k3)
      end do
      end do
      end do
      end do
!     send data from north to south 
      do k3=1,km3
      do k2=1,km2
      do k=1,km
      do i=1,imt
      sendtemp2(i,k,k2,k3)=mpiarray(i,2,k,k2,k3)
      recvtemp2(i,k,k2,k3)=mpiarray(i,jmt,k,k2,k3)
      end do
      end do
      end do
      end do
      call mpi_sendrecv(sendtemp2(1,1,1,1),swapsize,mpi_real8,south,13,  &
                        recvtemp2(1,1,1,1),swapsize,mpi_real8,north,13,  &
                        mpi_comm_world,status,ierr)
      do k3=1,km3
      do k2=1,km2
      do k=1,km
      do i=1,imt
      mpiarray(i,jmt,k,k2,k3)=recvtemp2(i,k,k2,k3)
      end do
      end do
      end do
      end do
      return
      end


!     =======================
      subroutine dis_var_int(svar,mat_myid,ncpux,ncpuy,myid)
!     =======================
      implicit none
      include 'mpif.h'
      integer i_id,j_id,myid,ncpux,ncpuy,ierr,status(mpi_status_size)
      integer mat_myid(ncpux+2,ncpuy)
      integer svar
      if (myid==0) then
      do j_id=1,ncpuy
      do i_id=2,ncpux+1
      if (mat_myid(i_id,j_id).gt.0) then
      call mpi_ssend(svar,1,mpi_integer,mat_myid(i_id,j_id),10,mpi_comm_world,ierr)
      end if
      end do
      end do
      end if
      
      if (myid.gt.0) then
      call mpi_recv(svar,1,mpi_integer,0,10,mpi_comm_world,status,ierr)
      end if
      
      return
      end
      
!     =======================
      subroutine dis_var_real(svar,mat_myid,ncpux,ncpuy,myid)
!     =======================
      implicit none
      include 'mpif.h'
      integer i_id,j_id,myid,ncpux,ncpuy,ierr,status(mpi_status_size)
      integer mat_myid(ncpux+2,ncpuy)
      real svar
      if (myid==0) then
      do j_id=1,ncpuy
      do i_id=2,ncpux+1
      if (mat_myid(i_id,j_id).gt.0) then
      call mpi_ssend(svar,1,mpi_real8,mat_myid(i_id,j_id),10,mpi_comm_world,ierr)
      end if
      end do
      end do
      end if
      
      if (myid.gt.0) then
      call mpi_recv(svar,1,mpi_real8,0,10,mpi_comm_world,status,ierr)
      end if
      
      return
      end
      
!     =======================
      subroutine dis_var_real1d(svar,imt,mat_myid,ncpux,ncpuy,myid)
!     =======================
      implicit none
      include 'mpif.h'
      integer imt,i_id,j_id,myid,ncpux,ncpuy,ierr,status(mpi_status_size)
      integer mat_myid(ncpux+2,ncpuy)
      real svar(imt)
      if (myid==0) then
      do j_id=1,ncpuy
      do i_id=2,ncpux+1
      if (mat_myid(i_id,j_id).gt.0) then
      call mpi_ssend(svar(1),imt,mpi_real8,mat_myid(i_id,j_id),10,mpi_comm_world,ierr)
      end if
      end do
      end do
      end if
      
      if (myid.gt.0) then
      call mpi_recv(svar(1),imt,mpi_real8,0,10,mpi_comm_world,status,ierr)
      end if
      
      return
      end

!     =======================
      subroutine div_array_real1dew(sarray,mpiarray,mat_myid,ncpux,ncpuy,simt,imt,myid)
!     =======================
      implicit none
      include 'mpif.h'
      integer simt,imt,i,j_id,i_id
      integer myid,ncpux,ncpuy,ierr,status(mpi_status_size)
      integer mat_myid(ncpux+2,ncpuy)
      real sendtemp,rectemp,sarray(simt),mpiarray(imt)
      if (myid==0) then
      do i=1,imt
      do j_id=1,ncpuy
      do i_id=2,ncpux+1
      if (mat_myid(i_id,j_id).gt.0) then
      sendtemp=sarray(i+(imt-2)*(i_id-2))
      call mpi_ssend(sendtemp,1,mpi_real8,mat_myid(i_id,j_id),10,mpi_comm_world,ierr)
      end if
      end do
      end do
      end do
      
      do i=1,imt
      mpiarray(i)=sarray(i)
      end do
      end if
      
      if (myid.gt.0) then
      do i=1,imt
      call mpi_recv(rectemp,1,mpi_real8,0,10,mpi_comm_world,status,ierr)
      mpiarray(i)=rectemp
      end do
      end if
      
      return
      end


!     =======================
      subroutine div_array_real1d(sarray,mpiarray,mat_myid,ncpux,ncpuy,sjmt,jmt,myid)
!     =======================
      implicit none
      include 'mpif.h'
      integer sjmt,jmt,j,j_id,i_id
      integer myid,ncpux,ncpuy,ierr,status(mpi_status_size)
      integer mat_myid(ncpux+2,ncpuy)
      real sendtemp,rectemp,sarray(sjmt),mpiarray(jmt)
      if (myid==0) then
      do j=1,jmt
      do j_id=1,ncpuy
      do i_id=2,ncpux+1
      if (mat_myid(i_id,j_id).gt.0) then
      sendtemp=sarray(j+(jmt-2)*(j_id-1))
      call mpi_ssend(sendtemp,1,mpi_real8,mat_myid(i_id,j_id),10,mpi_comm_world,ierr)
      end if
      end do
      end do
      end do
      
      do j=1,jmt
      mpiarray(j)=sarray(j)
      end do
      end if
      
      if (myid.gt.0) then
      do j=1,jmt
      call mpi_recv(rectemp,1,mpi_real8,0,10,mpi_comm_world,status,ierr)
      mpiarray(j)=rectemp
      end do
      end if
      
      return
      end

!     =======================
      subroutine div_array_int2d(sarray,mpiarray,mat_myid,ncpux,ncpuy,simt,sjmt,  &
                                 imt,jmt,myid)
!     =======================
      implicit none
      include 'mpif.h'
      integer simt,sjmt,imt,jmt,i,j,i_id,j_id
      integer myid,ncpux,ncpuy,ierr,status(mpi_status_size)
      integer mat_myid(ncpux+2,ncpuy)
      integer sendtemp,rectemp,sarray(simt,sjmt),mpiarray(imt,jmt)
      if (myid==0) then
      do j=1,jmt
      do i=1,imt
      do j_id=1,ncpuy
      do i_id=2,ncpux+1
      if (mat_myid(i_id,j_id).gt.0) then
      sendtemp=sarray(i+(imt-2)*(i_id-2),j+(jmt-2)*(j_id-1))
      call mpi_ssend(sendtemp,1,mpi_integer,mat_myid(i_id,j_id),10,mpi_comm_world,ierr)
      end if
      end do
      end do
      end do
      end do
      
      do j=1,jmt
      do i=1,imt
      mpiarray(i,j)=sarray(i,j)
      end do
      end do
      end if
      
      if (myid.gt.0) then
      do j=1,jmt
      do i=1,imt
      call mpi_recv(rectemp,1,mpi_integer,0,10,mpi_comm_world,status,ierr)
      mpiarray(i,j)=rectemp
      end do
      end do
      end if
      
      return
      end
      
!     =======================
      subroutine div_array_real2d(sarray,mpiarray,mat_myid,ncpux,ncpuy,simt,sjmt,  &
                                  imt,jmt,myid)
!     =======================
      implicit none
      include 'mpif.h'
      integer simt,sjmt,imt,jmt,i,j,i_id,j_id
      integer myid,ncpux,ncpuy,ierr,status(mpi_status_size)
      integer mat_myid(ncpux+2,ncpuy)
      real sendtemp,rectemp,sarray(simt,sjmt),mpiarray(imt,jmt)
      if (myid==0) then
      do j=1,jmt
      do i=1,imt
      do j_id=1,ncpuy
      do i_id=2,ncpux+1
      if (mat_myid(i_id,j_id).gt.0) then
      sendtemp=sarray(i+(imt-2)*(i_id-2),j+(jmt-2)*(j_id-1))
      call mpi_ssend(sendtemp,1,mpi_real8,mat_myid(i_id,j_id),10,mpi_comm_world,ierr)
      end if
      end do
      end do
      end do
      end do
      
      do j=1,jmt
      do i=1,imt
      mpiarray(i,j)=sarray(i,j)
      end do
      end do
      end if
      
      if (myid.gt.0) then
      do j=1,jmt
      do i=1,imt
      call mpi_recv(rectemp,1,mpi_real8,0,10,mpi_comm_world,status,ierr)
      mpiarray(i,j)=rectemp
      end do
      end do
      end if
      
      return
      end
!
!
!     =======================
      subroutine div_array_real3d(sarray,mpiarray,mat_myid,ncpux,ncpuy,simt,sjmt,  &
                                  km,imt,jmt,myid)
!     =======================
      implicit none
      include 'mpif.h'
      integer simt,sjmt,km,imt,jmt,i,j,k,i_id,j_id,sendsize
      integer myid,ncpux,ncpuy,ierr,status(mpi_status_size)
      integer mat_myid(ncpux+2,ncpuy)
      real sendtemp(imt,jmt,km),rectemp(imt,jmt,km)
      real sarray(simt,sjmt,km),mpiarray(imt,jmt,km)
      
      sendsize=imt*jmt*km
      if (myid==0) then 
      do j_id=1,ncpuy
      do i_id=2,ncpux+1
      if (mat_myid(i_id,j_id).gt.0) then
      do k=1,km
      do j=1,jmt
      do i=1,imt
      sendtemp(i,j,k)=sarray(i+(imt-2)*(i_id-2),j+(jmt-2)*(j_id-1),k)
      end do
      end do
      end do
      call mpi_ssend(sendtemp(1,1,1),sendsize,mpi_real8,mat_myid(i_id,j_id),10,  &
                     mpi_comm_world,ierr)
      end if
      end do
      end do
      
      do k=1,km
      do j=1,jmt
      do i=1,imt
      mpiarray(i,j,k)=sarray(i,j,k)
      end do
      end do
      end do
      end if
      
      if (myid.gt.0) then
      call mpi_recv(rectemp(1,1,1),sendsize,mpi_real8,0,10,mpi_comm_world,status,ierr)
      do k=1,km
      do j=1,jmt
      do i=1,imt
      mpiarray(i,j,k)=rectemp(i,j,k)
      end do
      end do
      end do
      end if
      
      return
      end

!     =======================
      subroutine div_array_real4d(sarray,mpiarray,mat_myid,ncpux,ncpuy,simt,sjmt,  &
                                  km,km2,imt,jmt,myid)
!     =======================
      implicit none
      include 'mpif.h'
      integer simt,sjmt,km,km2,imt,jmt,i,j,k,k2,i_id,j_id
      integer myid,ncpux,ncpuy,ierr,status(mpi_status_size)
      integer mat_myid(ncpux+2,ncpuy)
      real sendtemp,rectemp,sarray(simt,sjmt,km,km2),mpiarray(imt,jmt,km,km2)
      if (myid==0) then
      do k2=1,km2
      do k=1,km
      do j=1,jmt
      do i=1,imt
      do j_id=1,ncpuy
      do i_id=2,ncpux+1
      if (mat_myid(i_id,j_id).gt.0) then
      sendtemp=sarray(i+(imt-2)*(i_id-2),j+(jmt-2)*(j_id-1),k,k2)
      call mpi_ssend(sendtemp,1,mpi_real8,mat_myid(i_id,j_id),10,mpi_comm_world,ierr)
      end if
      end do
      end do
      end do
      end do
      end do
      end do
      
      do k2=1,km2
      do k=1,km
      do j=1,jmt
      do i=1,imt
      mpiarray(i,j,k,k2)=sarray(i,j,k,k2)
      end do
      end do
      end do
      end do
      end if
      
      if (myid.gt.0) then
      do k2=1,km2
      do k=1,km
      do j=1,jmt
      do i=1,imt
      call mpi_recv(rectemp,1,mpi_real8,0,10,mpi_comm_world,status,ierr)
      mpiarray(i,j,k,k2)=rectemp
      end do
      end do
      end do
      end do
      end if
      
      return
      end

!     =======================
      subroutine div_array_real5d(sarray,mpiarray,mat_myid,ncpux,ncpuy,simt,sjmt,  &
                                  km,km2,km3,imt,jmt,myid)
!     =======================
      implicit none
      include 'mpif.h'
      integer simt,sjmt,km,km2,km3,imt,jmt,i,j,k,k2,k3,i_id,j_id
      integer myid,ncpux,ncpuy,ierr,status(mpi_status_size)
      integer mat_myid(ncpux+2,ncpuy)
      real sendtemp,rectemp,sarray(simt,sjmt,km,km2,km3),mpiarray(imt,jmt,km,km2,km3)
      if (myid==0) then
      do k3=1,km3
      do k2=1,km2
      do k=1,km
      do j=1,jmt
      do i=1,imt
      do j_id=1,ncpuy
      do i_id=2,ncpux+1
      if (mat_myid(i_id,j_id).gt.0) then
      sendtemp=sarray(i+(imt-2)*(i_id-2),j+(jmt-2)*(j_id-1),k,k2,k3)
      call mpi_ssend(sendtemp,1,mpi_real8,mat_myid(i_id,j_id),10,mpi_comm_world,ierr)
      end if
      end do
      end do
      end do
      end do
      end do
      end do
      end do
      
      do k3=1,km3
      do k2=1,km2
      do k=1,km
      do j=1,jmt
      do i=1,imt
      mpiarray(i,j,k,k2,k3)=sarray(i,j,k,k2,k3)
      end do
      end do
      end do
      end do
      end do
      end if
      
      if (myid.gt.0) then
      do k3=1,km3
      do k2=1,km2
      do k=1,km
      do j=1,jmt
      do i=1,imt
      call mpi_recv(rectemp,1,mpi_real8,0,10,mpi_comm_world,status,ierr)
      mpiarray(i,j,k,k2,k3)=rectemp
      end do
      end do
      end do
      end do
      end do
      end if
      
      return
      end

!     =======================
      subroutine gath_var_real(var,mat_myid,ncpux,ncpuy,myid)
!     =======================
      implicit none
      include 'mpif.h'
      integer i_id,j_id
      integer myid,ncpux,ncpuy,ierr,status(mpi_status_size)
      integer mat_myid(ncpux+2,ncpuy)
      real var,tempvar
      
      tempvar=0
      if (myid.gt.0) then
      call mpi_ssend(var,1,mpi_real8,0,10,mpi_comm_world,ierr)
      end if
      
      if (myid==0) then
      do j_id=1,ncpuy
      do i_id=2,ncpux+1
      if (mat_myid(i_id,j_id).gt.0) then
      call mpi_recv(tempvar,1,mpi_real8,mat_myid(i_id,j_id),10,mpi_comm_world,status,ierr)
      var=var+tempvar
      end if
      end do
      end do
      end if
      
      return
      end

!     =======================
      subroutine gath_array_real1d(sarray,mpiarray,mat_myid,ncpux,ncpuy,sjmt,jmt,myid)
!     =======================
      implicit none
      include 'mpif.h'
      integer sjmt,jmt,j,i_id,j_id
      integer myid,ncpux,ncpuy,ierr,status(mpi_status_size)
      integer mat_myid(ncpux+2,ncpuy)
      real sendtemp,rectemp,sarray(sjmt),mpiarray(jmt)
      
      if (myid.gt.0) then
      do j=1,jmt
      sendtemp=mpiarray(j)
      call mpi_ssend(sendtemp,1,mpi_real8,0,10,mpi_comm_world,ierr)
      end do
      end if
      
      if (myid==0) then
      do j=1,jmt
      do j_id=1,ncpuy
      do i_id=2,ncpux+1
      if (mat_myid(i_id,j_id).gt.0) then
      call mpi_recv(rectemp,1,mpi_real8,mat_myid(i_id,j_id),10,mpi_comm_world,status,ierr)
      sarray(j+(jmt-2)*(j_id-1))=rectemp
      end if
      end do
      end do
      end do
      
      do j=1,jmt
      sarray(j)=mpiarray(j)
      end do
      end if
      
      return
      end

!     =======================
      subroutine gath_array_int2d(sarray,mpiarray,mat_myid,ncpux,ncpuy,simt,sjmt, &
                                  imt,jmt,myid)
!     =======================
      implicit none
      include 'mpif.h'
      integer simt,sjmt,imt,jmt,i,j,i_id,j_id
      integer myid,ncpux,ncpuy,ierr,status(mpi_status_size)
      integer mat_myid(ncpux+2,ncpuy)
      integer sendtemp,rectemp,sarray(simt,sjmt),mpiarray(imt,jmt)
      
      if (myid.gt.0) then
      do j=1,jmt
      do i=1,imt
      sendtemp=mpiarray(i,j)
      call mpi_ssend(sendtemp,1,mpi_integer,0,10,mpi_comm_world,ierr)
      end do
      end do
      end if
      
      if (myid==0) then
      do j=1,jmt
      do i=1,imt
      do j_id=1,ncpuy
      do i_id=2,ncpux+1
      if (mat_myid(i_id,j_id).gt.0) then
      call mpi_recv(rectemp,1,mpi_integer,mat_myid(i_id,j_id),10,mpi_comm_world,status,ierr)
      sarray(i+(imt-2)*(i_id-2),j+(jmt-2)*(j_id-1))=rectemp
      end if
      end do
      end do
      end do
      end do
      
      do j=1,jmt
      do i=1,imt
      sarray(i,j)=mpiarray(i,j)
      end do
      end do
      end if
      
      return
      end

!     =======================
      subroutine gath_array_real2d(sarray,mpiarray,mat_myid,ncpux,ncpuy,simt,sjmt, &
                                   imt,jmt,myid)
!     =======================
      implicit none
      include 'mpif.h'
      integer simt,sjmt,imt,jmt,i,j,i_id,j_id,sendsize
      integer myid,ncpux,ncpuy,ierr,status(mpi_status_size)
      integer mat_myid(ncpux+2,ncpuy)
      real sendtemp(imt,jmt),rectemp(imt,jmt),sarray(simt,sjmt),mpiarray(imt,jmt)
      
      sendsize=imt*jmt
      if (myid.gt.0) then
      do j=1,jmt
      do i=1,imt
      sendtemp(i,j)=mpiarray(i,j)
      end do
      end do
      call mpi_ssend(sendtemp(1,1),sendsize,mpi_real8,0,10,mpi_comm_world,ierr)
      end if
      
      if (myid==0) then
      do j_id=1,ncpuy
      do i_id=2,ncpux+1
      if (mat_myid(i_id,j_id).gt.0) then
      call mpi_recv(rectemp(1,1),sendsize,mpi_real8,mat_myid(i_id,j_id),10,  &
                    mpi_comm_world,status,ierr)
      do j=1,jmt
      do i=1,imt
      sarray(i+(imt-2)*(i_id-2),j+(jmt-2)*(j_id-1))=rectemp(i,j)
      end do
      end do
      end if
      end do
      end do
      
      do j=1,jmt
      do i=1,imt
      sarray(i,j)=mpiarray(i,j)
      end do
      end do
      end if
      
      return
      end

!     =======================
      subroutine gath_merid_real2d(sarray,mpiarray,mat_myid,ncpux,ncpuy,sjmt, &
                                   km,jmt,myid)
!     =======================
      implicit none
      include 'mpif.h'
      integer sjmt,jmt,km,j,k,i_id,j_id,sendsize
      integer myid,ncpux,ncpuy,ierr,status(mpi_status_size)
      integer mat_myid(ncpux+2,ncpuy)
      real sendtemp(jmt,km),rectemp(jmt,km),sarray(sjmt,km),mpiarray(jmt,km)
      
      sarray=0
      sendsize=jmt*km
      if (myid.gt.0) then
      do k=1,km
      do j=1,jmt
      sendtemp(j,k)=mpiarray(j,k)
      end do
      end do
      call mpi_ssend(sendtemp(1,1),sendsize,mpi_real8,0,10,mpi_comm_world,ierr)
      end if
      
      if (myid==0) then
      do j_id=1,ncpuy
      do i_id=2,ncpux+1
      if (mat_myid(i_id,j_id).gt.0) then
      call mpi_recv(rectemp(1,1),sendsize,mpi_real8,mat_myid(i_id,j_id),10,  &
                    mpi_comm_world,status,ierr)
      do k=1,km
      do j=1,jmt
      sarray(j+(jmt-2)*(j_id-1),k)=sarray(j+(jmt-2)*(j_id-1),k)+rectemp(j,k)
      end do
      end do
      end if
      end do
      end do
      
      do k=1,km
      do j=1,jmt
      sarray(j,k)=sarray(j,k)+mpiarray(j,k)
      end do
      end do
      end if
      
      return
      end

!     =======================
      subroutine gath_array_real3d(sarray,mpiarray,mat_myid,ncpux,ncpuy,simt,sjmt,  &
                                   km,imt,jmt,myid)
!     =======================
      implicit none
      include 'mpif.h'
      integer simt,sjmt,km,imt,jmt,i,j,k,i_id,j_id,sendsize
      integer myid,ncpux,ncpuy,ierr,status(mpi_status_size)
      integer mat_myid(ncpux+2,ncpuy)
      real sarray(simt,sjmt,km),mpiarray(imt,jmt,km)
      real sendtemp(imt,jmt,km),rectemp(imt,jmt,km)
      
      sendsize=imt*jmt*km
      if (myid.gt.0) then
      do k=1,km
      do j=1,jmt
      do i=1,imt
      sendtemp(i,j,k)=mpiarray(i,j,k)
      end do
      end do
      end do
      call mpi_ssend(sendtemp(1,1,1),sendsize,mpi_real8,0,10,mpi_comm_world,ierr)
      end if
      
      if (myid==0) then
      do j_id=1,ncpuy
      do i_id=2,ncpux+1
      if (mat_myid(i_id,j_id).gt.0) then
      call mpi_recv(rectemp(1,1,1),sendsize,mpi_real8,mat_myid(i_id,j_id),10,  &
                    mpi_comm_world,status,ierr)
      do k=1,km
      do j=1,jmt
      do i=1,imt
      sarray(i+(imt-2)*(i_id-2),j+(jmt-2)*(j_id-1),k)=rectemp(i,j,k)
      end do
      end do
      end do
      end if
      end do
      end do
      
      do k=1,km
      do j=1,jmt
      do i=1,imt
      sarray(i,j,k)=mpiarray(i,j,k)
      end do
      end do
      end do
      end if
      
      return
      end

!     =======================
      subroutine gath_array_real4d(sarray,mpiarray,mat_myid,ncpux,ncpuy,simt,sjmt,  &
                                   km,km2,imt,jmt,myid)
!     =======================
      implicit none
      include 'mpif.h'
      integer simt,sjmt,km,km2,imt,jmt,i,j,k,k2,i_id,j_id,sendsize
      integer myid,ncpux,ncpuy,ierr,status(mpi_status_size)
      integer mat_myid(ncpux+2,ncpuy)
      real sarray(simt,sjmt,km,km2),mpiarray(imt,jmt,km,km2)
      real sendtemp(imt,jmt,km,km2),rectemp(imt,jmt,km,km2)
      
      sendsize=imt*jmt*km*km2
      if (myid.gt.0) then
      do k2=1,km2
      do k=1,km
      do j=1,jmt
      do i=1,imt
      sendtemp(i,j,k,k2)=mpiarray(i,j,k,k2)
      end do
      end do
      end do
      end do
      call mpi_ssend(sendtemp(1,1,1,1),sendsize,mpi_real8,0,10,mpi_comm_world,ierr)
      end if
      
      if (myid==0) then
      do j_id=1,ncpuy
      do i_id=2,ncpux+1
      rectemp=0
      if (mat_myid(i_id,j_id).gt.0) then
      call mpi_recv(rectemp(1,1,1,1),sendsize,mpi_real8,mat_myid(i_id,j_id),10,  &
                    mpi_comm_world,status,ierr)
      do k2=1,km2
      do k=1,km
      do j=1,jmt
      do i=1,imt
      sarray(i+(imt-2)*(i_id-2),j+(jmt-2)*(j_id-1),k,k2)=rectemp(i,j,k,k2)
      end do
      end do
      end do
      end do
      end if
      end do
      end do
      
      do k2=1,km2
      do k=1,km
      do j=1,jmt
      do i=1,imt
      sarray(i,j,k,k2)=mpiarray(i,j,k,k2)
      end do
      end do
      end do
      end do
      end if
      
      return
      end

!     =======================
      subroutine gath_array_real5d(sarray,mpiarray,mat_myid,ncpux,ncpuy,simt,sjmt,  &
                                   km,km2,km3,imt,jmt,myid)
!     =======================
      implicit none
      include 'mpif.h'
      integer simt,sjmt,km,km2,km3,imt,jmt,i,j,k,k2,k3,i_id,j_id,sendsize
      integer myid,ncpux,ncpuy,ierr,status(mpi_status_size)
      integer mat_myid(ncpux+2,ncpuy)
      real sarray(simt,sjmt,km,km2,km3),mpiarray(imt,jmt,km,km2,km3)
      real sendtemp(imt,jmt,km,km2,km3),rectemp(imt,jmt,km,km2,km3)
      
      sendsize=imt*jmt*km*km2*km3
      if (myid.gt.0) then
      do k3=1,km3
      do k2=1,km2
      do k=1,km
      do j=1,jmt
      do i=1,imt
      sendtemp(i,j,k,k2,k3)=mpiarray(i,j,k,k2,k3)
      end do
      end do
      end do
      end do
      end do
      call mpi_ssend(sendtemp(1,1,1,1,1),sendsize,mpi_real8,0,10,mpi_comm_world,ierr)
      end if
      
      if (myid==0) then
      do j_id=1,ncpuy
      do i_id=2,ncpux+1
      rectemp=0
      if (mat_myid(i_id,j_id).gt.0) then
      call mpi_recv(rectemp(1,1,1,1,1),sendsize,mpi_real8,mat_myid(i_id,j_id),10,  &
                    mpi_comm_world,status,ierr)
      do k3=1,km3
      do k2=1,km2
      do k=1,km
      do j=1,jmt
      do i=1,imt
      sarray(i+(imt-2)*(i_id-2),j+(jmt-2)*(j_id-1),k,k2,k3)=rectemp(i,j,k,k2,k3)
      end do
      end do
      end do
      end do
      end do
      end if
      end do
      end do
      
      do k3=1,km3
      do k2=1,km2
      do k=1,km
      do j=1,jmt
      do i=1,imt
      sarray(i,j,k,k2,k3)=mpiarray(i,j,k,k2,k3)
      end do
      end do
      end do
      end do
      end do
      end if
      
      return
      end

