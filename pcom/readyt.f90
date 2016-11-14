!
!
!     =================
!BOP
!
! !MODULE: readyt.f90
! !DESCRIPTION: \input{sections/code-readyt}
!
! !INTERFACE:
!
      subroutine readyt(tmask,z,dz,rzu,t,pbt,spbt,pmtm,pmtp,bcp,rho,rhodp,umm,vmm,   &
                        ump,vmp,up,vp,cosu,rdxt,rdy,pax,pay,itn,ivn,pbxn,pbxs,pbye,  &
                        pbyw,pcxn,pcxs,pcye,pcyw,pdxn,pdxs,pdye,pdyw,imt,jmt,km,nt,  &
                        imm,jmm,kmp1,decibar,west,east,north,south, &
                        fixp)
!EOP
!-------------------------------------------------------------------------------
!     =================
!     calculate the variables depended on stratification
!
      implicit none
      include 'pconst.h'
      include 'mpif.h'
!
      integer imt,jmt,km,nt,imm,jmm,kmp1,i,j,k
      integer itn(imt,jmt),ivn(imt,jmt)
      real a,t1,t2,t3,rdens,decibar
      real tmask(imt,jmt,km)
      real t(imt,jmt,km,nt,2),up(imt,jmt,km,2),vp(imt,jmt,km,2)
      real ump(imt,jmt,km),umm(imt,jmt,km),vmp(imt,jmt,km),vmm(imt,jmt,km)
      real bcp(imt,jmt)
      real dz(km),rzu(imt,jmt),cosu(jmt)
      real pbt(imt,jmt,2),spbt(imt,jmt),pmtp(imt,jmt),pmtm(imt,jmt)
      real rhodp(imt,jmt,km),rho(imt,jmt,km)
      real rdxt(jmt),rdy,z(km),fixp(imt,jmt,km)
      real pax(imt,jmt),pay(imt,jmt)
      real pbxn(imt,jmt),pbxs(imt,jmt)
      real pcxn(imt,jmt),pcxs(imt,jmt)
      real pdxn(imt,jmt),pdxs(imt,jmt)
      real pbye(imt,jmt),pbyw(imt,jmt)
      real pcye(imt,jmt),pcyw(imt,jmt)
      real pdye(imt,jmt),pdyw(imt,jmt)
      
      integer west,east,north,south
!
!-----------------------------------------------------------------------
!     calculate density/rho_0 (BCOM) or 1/density (PCOM)
!-----------------------------------------------------------------------
      call density(tmask,z,t,pbt,bcp,rho,decibar,imt,jmt,km,nt,imm,jmm,fixp,  &
                   west,east,north,south)
!
!
!-----------------------------------------------------------------------
!     initialize the time-averaged pbt & up vp
!-----------------------------------------------------------------------
      do j=1,jmt
      do i=1,imt
      pmtm(i,j) = pmtp(i,j)
      pmtp(i,j) = pbt (i,j,tau)
      enddo
      enddo
!
      do j=2,jmm
      do i=2,imm
      spbt(i,j) = p5*sqrt(pbt(i,j  ,tau) + pbt(i+1,j  ,tau) +  &
                          pbt(i,j+1,tau) + pbt(i+1,j+1,tau))
      enddo
      enddo
      
      call swap_array_real2d(spbt,imt,jmt,west,east,north,south)
!
      do k=1,km
      do j=1,jmt
      do i=1,imt
      umm(i,j,k) = ump(i,j,k)
      vmm(i,j,k) = vmp(i,j,k)
      ump(i,j,k) = up(i,j,k,tau)*spbt(i,j)
      vmp(i,j,k) = vp(i,j,k,tau)*spbt(i,j)*cosu(j)
      enddo
      enddo
      enddo
!
!
!-----------------------------------------------------------------------
!     calculate pressure gradients related to stratification
!-----------------------------------------------------------------------
      call prsgrd(bcp,rdxt,rdy,pax,pay,itn,ivn,rho,rhodp,z,dz,rzu,pbxn,pbxs, &
                  pbye,pbyw,pcxn,pcxs,pcye,pcyw,pdxn,pdxs,pdye,pdyw,imt,jmt,  &
                  km,imm,jmm,kmp1,west,east,north,south)
!
!
      return
      end
