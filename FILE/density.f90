!     ==================
      subroutine density(tmask,z,t,pbt,bcp,rho,decibar,imt,jmt,km,nt,imm,jmm,fixp,  &
                         west,east,north,south,unesco,boussinesq)
!     ==================
!     calculate density/rho_0(BCOM) or reciprocal of density(PCOM)
!
      implicit none
      include 'pconst.h'
      include 'mpif.h'
!
      integer imt,jmt,km,nt,imm,jmm,i,j,k,unesco,boussinesq
      real t0,s0,p0,dens,rdens,unrdens,undens,decibar
      real tmask(imt,jmt,km),z(km),pbt(imt,jmt,2),rho(imt,jmt,km),fixp(imt,jmt,km)
      real t(imt,jmt,km,nt,2)
      real bcp(imt,jmt)
      
      integer west,east,north,south
!
!     calculate density/rho_0 or reciprocal of density
!
      if (boussinesq==1) then
      
      if (unesco==1) then
      do k=1,km
      do j=2,jmm
      do i=2,imm
      if(tmask(i,j,k).gt.c0)then
      t0         = t(i,j,k,1,tau)
      s0         = t(i,j,k,2,tau)
      p0         = fixp(i,j,k)
      rho(i,j,k) = undens(t0,s0,p0)*rrho_0
      endif
      enddo
      enddo
      enddo
      else
      do k=1,km
      do j=2,jmm
      do i=2,imm
      if(tmask(i,j,k).gt.c0)then
      t0         = t(i,j,k,1,tau)
      s0         = t(i,j,k,2,tau)
      p0         = fixp(i,j,k)
      rho(i,j,k) = dens(t0,s0,p0)*rrho_0
      endif
      enddo
      enddo
      enddo
      end if
      
      else
      
      if (unesco==1) then
      do k=1,km
      do j=2,jmm
      do i=2,imm
      if(tmask(i,j,k).gt.c0)then
      t0         = t(i,j,k,1,tau)
      s0         = t(i,j,k,2,tau)
      p0         = (pbt(i,j,tau)*z(k) + bcp(i,j))*decibar
      rho(i,j,k) = unrdens(t0,s0,p0)
      endif
      enddo
      enddo
      enddo
      else
      do k=1,km
      do j=2,jmm
      do i=2,imm
      if(tmask(i,j,k).gt.c0)then
      t0         = t(i,j,k,1,tau)
      s0         = t(i,j,k,2,tau)
      p0         = (pbt(i,j,tau)*z(k) + bcp(i,j))*decibar
      rho(i,j,k) = rdens(t0,s0,p0)
      endif
      enddo
!      rho(1  ,j,k) = rho(imm,j,k)
!      rho(imt,j,k) = rho(2  ,j,k)
      enddo
      enddo
      end if
      
      end if
      
      call swap_array_real3d(rho,imt,jmt,km,west,east,north,south)
!
      return
      end
!
!
!     ==================
      subroutine rho_ref(tmask,t,z,pbt,pt,ps,deltat,deltas,rdeltat,rdeltas,decibar,fixp, &
                         imt,jmt,km,nt,imm,jmm,west,east,north,south,unesco,boussinesq)
!     ==================
!
!     pt = {partial rho} over {partial temperature}
!     ps = {partial rho} over {partial salinity   }
!
!     potential density = pt*(t-t0) + ps*(s-s0) + den0
!
!     it will compare pt*t+ps*s rather than potential density
!     itself in subroutine convect
!
      implicit none
      include 'pconst.h'
      include 'mpif.h'
!
      integer imt,jmt,km,nt,imm,jmm,i,j,k,i_id,j_id,unesco,boussinesq
      real t1,s1,p1,dens,undens
      real rhop,rhoq,decibar,deltas,rdeltat,rdeltas,deltat
      real tmask(imt,jmt,km),z(km),pbt(imt,jmt,2),pt(imt,jmt,km),ps(imt,jmt,km)
      real t(imt,jmt,km,nt,2),fixp(imt,jmt,km)
      
      integer west,east,north,south
!
      if (boussinesq==1) then
      
      if (unesco==1) then
      do k=1,km
      do j=2,jmm
      do i=2,imm
      if(tmask(i,j,k).gt.c0)then
        t1        = t(i,j,k,1,tau)
        s1        = t(i,j,k,2,tau)
        p1        = fixp(i,j,k)
        rhop      = undens(t1+deltat,s1,p1)
        rhoq      = undens(t1-deltat,s1,p1)
        pt(i,j,k) = (rhop-rhoq)*rdeltat
        rhop      = undens(t1,s1+deltas,p1)
        rhoq      = undens(t1,s1-deltas,p1)
        ps(i,j,k) = (rhop-rhoq)*rdeltas
      endif
      enddo
      enddo
      enddo
      else
      do k=1,km
      do j=2,jmm
      do i=2,imm
      if(tmask(i,j,k).gt.c0)then
        t1        = t(i,j,k,1,tau)
        s1        = t(i,j,k,2,tau)
        p1        = fixp(i,j,k)
        rhop      = dens(t1+deltat,s1,p1)
        rhoq      = dens(t1-deltat,s1,p1)
        pt(i,j,k) = (rhop-rhoq)*rdeltat
        rhop      = dens(t1,s1+deltas,p1)
        rhoq      = dens(t1,s1-deltas,p1)
        ps(i,j,k) = (rhop-rhoq)*rdeltas
      endif
      enddo
      enddo
      enddo
      end if
      
      else
      
      if (unesco==1) then
      do k=1,km
      do j=2,jmm
      do i=2,imm
      if(tmask(i,j,k).gt.c0)then
        t1        = t(i,j,k,1,tau)
        s1        = t(i,j,k,2,tau)
        p1        = pbt(i,j,tau)*z(k)*decibar
        rhop      = undens(t1+deltat,s1,p1)
        rhoq      = undens(t1-deltat,s1,p1)
        pt(i,j,k) = (rhop-rhoq)*rdeltat
        rhop      = undens(t1,s1+deltas,p1)
        rhoq      = undens(t1,s1-deltas,p1)
        ps(i,j,k) = (rhop-rhoq)*rdeltas
      endif
      enddo
      enddo
      enddo
      else
      do k=1,km
      do j=2,jmm
      do i=2,imm
      if(tmask(i,j,k).gt.c0)then
        t1        = t(i,j,k,1,tau)
        s1        = t(i,j,k,2,tau)
        p1        = pbt(i,j,tau)*z(k)*decibar
        rhop      = dens(t1+deltat,s1,p1)
        rhoq      = dens(t1-deltat,s1,p1)
        pt(i,j,k) = (rhop-rhoq)*rdeltat
        rhop      = dens(t1,s1+deltas,p1)
        rhoq      = dens(t1,s1-deltas,p1)
        ps(i,j,k) = (rhop-rhoq)*rdeltas
      endif
      enddo
!      pt(1  ,j,k) = pt(imm,j,k)
!      pt(imt,j,k) = pt(2  ,j,k)
!      ps(1  ,j,k) = ps(imm,j,k)
!      ps(imt,j,k) = ps(2  ,j,k)
      enddo
      enddo
      end if
      
      end if
      
      call swap_array_real3d(pt,imt,jmt,km,west,east,north,south)
      call swap_array_real3d(ps,imt,jmt,km,west,east,north,south)
!
      return
      end
!
!
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

!     ====================
      function undens(t,s,p0)
!     ====================
!     this function calculates the density of seawater using the
!     standard equation of state recommended by unesco(1981).
!
!     t = potential temperature, degrees centigrade
!     s = salinity, practical salinity units
!     p = meters of depth
!
!     output  dens: gram per cubic cm
!
!     references:
!	
!	Coefficients of K is according to Jackett and Mcdougail
!	J. Atmos. & Ocean. Tech.    1995 Apr., P381-389
!
!	Coefficients of rho0 is given by Millero and A.Poisson
!	Deep-Sea Res., 28A, 625-629
!
      p=p0*1.0d-1 + 1.013	!! standard Pa
      p2=p*p
      t2=t*t
      t3=t2*t
      t4=t3*t
      t5=t4*t
      s32=s**1.5d0
      s2=s*s
      rw   =   9.99842594d2 + 6.793952d-2*t - 9.09529d-3*t2  &
             + 1.001685d-4*t3 - 1.120083d-6*t4 + 6.536336d-9*t5  &
             + (8.24493d-1 - 4.0899d-3*t + 7.6438d-5*t2  &
             - 8.2467d-7*t3 + 5.3875d-9*t4) * s  &
             + (-5.72466d-3 + 1.0227d-4*t - 1.6546d-6*t2) * s32  &
             + 4.8314d-4 * s2
      rk   =   1.965933d4   +  1.444304d2*t - 1.706103d0*t2  &
             + 9.648704d-3*t3 -4.190253d-5*t4  &
             + (5.284855d1 - 3.101089d-1*t + 6.283263d-3*t2  &
             - 5.084188d-5*t3 ) *s  &
             + (3.886640d-1 + 9.085835d-3*t - 4.619924d-4*t2)*s32  &
             + (3.186519d0 + 2.212276d-2*t -2.984642d-4*t2  &
             + 1.956415d-6*t3)*p  &
             + ((6.704388d-3 -1.847318d-4*t + 2.059331d-7*t2)*s  &
             + 1.480266d-4*s32)*p  &
             + (2.102898d-4 -1.202016d-5*t+1.394680d-7*t2  &
             - 2.040237d-6*s +6.128773d-8*s*t +6.207323d-10*s*t2)*p2
      undens = rw/(1.0d0-p/rk)*1.0d-3
      return
      end
!
!
!     ======================
      function unrdens(t,s,p0)
!     ======================
!     this function calculates the reciprocal of density
!
!     t = potential temperature, degrees centigrade
!     s = salinity, practical salinity units
!     p = meters of depth
!
!     output  dens: gram per cubic cm
!
!     references:
!	
!	Coefficients of K is according to Jackett and Mcdougail
!	J. Atmos. & Ocean. Tech.    1995 Apr., P381-389
!
!	Coefficients of rho0 is given by Millero and A.Poisson
!	Deep-Sea Res., 28A, 625-629
!
      p=p0*1.0d-1 + 1.013	!! standard Pa
      p2=p*p
      t2=t*t
      t3=t2*t
      t4=t3*t
      t5=t4*t
      s32=s**1.5d0
      s2=s*s
      rw   =   9.99842594d2 + 6.793952d-2*t - 9.09529d-3*t2  &
             + 1.001685d-4*t3 - 1.120083d-6*t4 + 6.536336d-9*t5  &
             + (8.24493d-1 - 4.0899d-3*t + 7.6438d-5*t2  &
             - 8.2467d-7*t3 + 5.3875d-9*t4) * s  &
             + (-5.72466d-3 + 1.0227d-4*t - 1.6546d-6*t2) * s32  &
             + 4.8314d-4 * s2
      rk   =   1.965933d4   +  1.444304d2*t - 1.706103d0*t2  &
             + 9.648704d-3*t3 -4.190253d-5*t4  &
             + (5.284855d1 - 3.101089d-1*t + 6.283263d-3*t2  &
             - 5.084188d-5*t3 ) *s  &
             + (3.886640d-1 + 9.085835d-3*t - 4.619924d-4*t2)*s32  &
             + (3.186519d0 + 2.212276d-2*t -2.984642d-4*t2  &
             + 1.956415d-6*t3)*p  &
             + ((6.704388d-3 -1.847318d-4*t + 2.059331d-7*t2)*s  &
             + 1.480266d-4*s32)*p  &
             + (2.102898d-4 -1.202016d-5*t+1.394680d-7*t2  &
             - 2.040237d-6*s +6.128773d-8*s*t +6.207323d-10*s*t2)*p2
!cc   dens  = rw/(1.0d0-p/rk)*1.0d-3
      unrdens = (1.0d0-p/rk)/rw*1.0d3
      return
      end
      
!     =======================
      function dens(t,s,p)
!     =======================
!
!     computer the density of seawater by using a 1st order polynomial
!     fit to the equation of state recommended by unesco(1981)
!
!     t = in-situ temperature, degrees centigrade
!     s = salinity, practical salinity units
!     p = pressure, decibars, approx. as meters of depth
!
!     d0,t0,s0,p0 are the refercence density, temperature, salinity
!     and pressure, respectively.
!
!     reference:  density/1st-order.f
!
      t0 = 4.0
      s0 = 35.0
      p0 = 2000.0
      d0 = 1.036903125227324
!
      a0 = - 0.1523022015184097*1.0e-3
      b0 =   0.7807833053448121*1.0e-3
      c0 =   4.4622974778576470*1.0e-6
!     density is only the function of T&S
      c0 =   0.0d0
!
      dens = d0 + a0*(t-t0) + b0*(s-s0) + c0*(p-p0)
!
!c       r0= 1028.106
!c       dens=(r0 + 0.7948*(s-35.0) - 0.05968*t - 0.0063*t*t
!c     *   	+ 3.7315*1.0e-5 *t*t*t)/1000.0
      return
      end
!
!
!     =======================
      function rdens(t,s,p)
!     =======================
!
!     computer the density of seawater by using a 1st order polynomial
!     fit to the equation of state recommended by unesco(1981)
!
!     t = in-situ temperature, degrees centigrade
!     s = salinity, practical salinity units
!     p = pressure, decibars, approx. as meters of depth
!
!     d0,t0,s0,p0 are the refercence density, temperature, salinity
!     and pressure, respectively.
!
!     reference:  density/1st-order.f
!
      t0 = 4.0
      s0 = 35.0
      p0 = 2000.0
      d0 = 1.036903125227324
!
      a0 = - 0.1523022015184097*1.0e-3
      b0 =   0.7807833053448121*1.0e-3
      c0 =   4.4622974778576470*1.0e-6
!     density is only the function of T&S
      c0 =   0.0d0
!
      rdens = 1.0d0/(d0 + a0*(t-t0) + b0*(s-s0) + c0*(p-p0))
!
!c       r0= 1028.106
!c       rdens=1.0d+3/(r0 + 0.7948*(s-35.0) - 0.05968*t - 0.0063*t*t
!c     *   	+ 3.7315*1.0e-5 *t*t*t)
      return
      end
