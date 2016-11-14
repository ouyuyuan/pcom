!
!     =================
      subroutine readyt_st(tmask,z,dz,rzu,t,pbt_st,bcp,rho,rhodp,   &
                    up,vp,cosu,rdxt,rdy,pax,pay,itn,ivn,pbxn,pbxs,pbye,  &
                    pbyw,pcxn,pcxs,pcye,pcyw,pdxn,pdxs,pdye,pdyw,imt,jmt,km,nt,  &
                    imm,jmm,kmp1,decibar,west,east,north,south,unesco,boussinesq, &
                    fixp,energydiag)
!     =================
!     calculate the variables depended on stratification
!
      implicit none
      include 'pconst.h'
      include 'mpif.h'
!
      integer imt,jmt,km,nt,imm,jmm,kmp1,i,j,k,k2,n,unesco,boussinesq,energydiag
      integer itn(imt,jmt),ivn(imt,jmt)
      real a,t1,t2,t3,rdens,decibar
      real tmask(imt,jmt,km)
      real t(imt,jmt,km,nt,2),up(imt,jmt,km,2),vp(imt,jmt,km,2)
      real bcp(imt,jmt)
      real dz(km),rzu(imt,jmt),cosu(jmt)
      real pbt(imt,jmt,2),pbt_st(imt,jmt,4),spbt(imt,jmt)
      real rhodp(imt,jmt,km),rho(imt,jmt,km)
      real rdxt(jmt),rdy,z(km),fixp(imt,jmt,km)
      real pax(imt,jmt),pay(imt,jmt)
      real pbxn(imt,jmt),pbxs(imt,jmt)
      real pcxn(imt,jmt),pcxs(imt,jmt)
      real pdxn(imt,jmt),pdxs(imt,jmt)
      real pbye(imt,jmt),pbyw(imt,jmt)
      real pcye(imt,jmt),pcyw(imt,jmt)
      real pdye(imt,jmt),pdyw(imt,jmt)
      
      real t_ba(imt,jmt,km,nt,2),pbt_ba(imt,jmt,2)
      
      integer west,east,north,south
!
!-----------------------------------------------------------------------
!     calculate density/rho_0 (BCOM) or 1/density (PCOM)
!-----------------------------------------------------------------------
!------BUG ????-----------------
      do k=1,4
      do j=1,jmt
         do i=1,imt
            if (pbt_st(i,j,k).eq.c0) then
               pbt_st(i,j,k)=c1
            end if
         end do
      end do
      end do
!-----BUG ????------------------

      do j=1,jmt
         do i=1,imt
            pbt(i,j,tau)=pbt_st(i,j,1)
         end do
      end do
      
      call density(tmask,z,t,pbt,bcp,rho,decibar,imt,jmt,km,nt,imm,jmm,fixp,  &
                   west,east,north,south,unesco,boussinesq)
!
!-----------------------------------------------------------------------
!     calculate pressure gradients related to stratification
!-----------------------------------------------------------------------
      call prsgrd(bcp,rdxt,rdy,pax,pay,itn,ivn,rho,rhodp,z,dz,rzu,pbxn,pbxs, &
                  pbye,pbyw,pcxn,pcxs,pcye,pcyw,pdxn,pdxs,pdye,pdyw,imt,jmt,  &
                  km,imm,jmm,kmp1,west,east,north,south,boussinesq)
!     
      return
      end subroutine readyt_st
