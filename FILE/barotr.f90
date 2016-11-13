!
!     =================
      subroutine barotr(ivn,itn,upb,vpb,r1c,r1d,sdxu,am,dub,dvb,spbt,pbt,phibx,    &
                        phiby,pbxn,pbxs,pbye,pbyw,pcxn,pcxs,pcye,pcyw,ff,ebla,     &
                        eblb,ebea,ebeb,pn,zu,cosu,rdxt,rdyt,leapfrog_b,euler_back, &
                        dtsf,c2dtsf,afb1,afb2,pmup,pmtp,nbb,imt,jmt,km,imm,jmm,    &
                        west,east,north,south,asselin_b,snbc,emp,jstn,jedn,jsts,   &
                        jeds,smth,umask,boussinesq,phib,pdxn,pdxs,pdye,pdyw)
!     =================
!     compute pbt, upb & vpb at "tau+1" time level
!
      implicit none
      include 'pconst.h'
      include 'mpif.h'
!
      integer mode_b,nbb,asselin_b,snbc,boussinesq
      integer imt,jmt,km,imm,jmm,i,j,k
      integer jstn,jedn,jsts,jeds,smth
      logical euler_back,leapfrog_b
      integer ivn(imt,jmt),itn(imt,jmt)
      real t1,t2,t3,am(imt,jmt,km),afb1,afb2,dtsf,c2dtsf
      real a(imt,jmt),b(imt,jmt),dp(imt,jmt)
      real adu(imt,jmt),adv(imt,jmt),etax(imt,jmt),etay(imt,jmt)
      real rdxt(jmt),rdyt(jmt),sdxu(jmt),umask(imt,jmt,km)
      real r1c(jmt),r1d(jmt),cosu(jmt),ff(jmt),zu(imt,jmt),pn(imt,jmt)
      real du2(imt,jmt),dv2(imt,jmt)
      real upb(imt,jmt,2),vpb(imt,jmt,2)
      real pbxn(imt,jmt),pbxs(imt,jmt)
      real pcxn(imt,jmt),pcxs(imt,jmt)
      real pbye(imt,jmt),pbyw(imt,jmt)
      real pcye(imt,jmt),pcyw(imt,jmt)
      real pdxn(imt,jmt),pdxs(imt,jmt)
      real pdye(imt,jmt),pdyw(imt,jmt),phib(imt,jmt)
      real ebea(jmt),ebeb(jmt),ebla(jmt),eblb(jmt)
      real pmup(imt,jmt),pmtp(imt,jmt),pbt(imt,jmt,2),spbt(imt,jmt)
      real dub(imt,jmt),dvb(imt,jmt)
      real phibx(imt,jmt),phiby(imt,jmt),emp(imt,jmt)
      
      integer west,east,north,south
!
!
      do 1000 mode_b = 1,nbb
!
!-----------------------------------------------------------------------
!     compute the "artificial" horizontal viscosity
!-----------------------------------------------------------------------
!
      do j=2,jmm
      do i=2,imm
      if(ivn(i,j).gt.0)then
      adu(i,j) = am(i,j,1)*( r1c(j)*(upb(i,j+1,taum)-upb(i,j,taum)) -  &
                      r1d(j)*(upb(i,j,taum)-upb(i,j-1,taum)) +   &
                     sdxu(j)*(upb(i+1,j,taum)-c2*upb(i,j,taum)+ &
                              upb(i-1,j,taum)) )*c2
      adv(i,j) = am(i,j,1)*( r1c(j)*(vpb(i,j+1,taum)-vpb(i,j,taum)) -  &
                      r1d(j)*(vpb(i,j,taum)-vpb(i,j-1,taum)) +  &
                     sdxu(j)*(vpb(i+1,j,taum)-c2*vpb(i,j,taum)+  &
                              vpb(i-1,j,taum)) )*c2
      endif
      enddo
      enddo
!
      if(mode_b.eq.1) then
        do j=2,jmm
        do i=2,imm
        dub(i,j) = dub(i,j) - adu(i,j)
        dvb(i,j) = dvb(i,j) - adv(i,j)
        enddo
        enddo
      end if
!
!
1001  continue
!


!-----------------------------------------------------------------------
!     square root of pbt at the U cell
!-----------------------------------------------------------------------
      do j=2,jmm
      do i=2,imm
      if(ivn(i,j).gt.0)then
      spbt(i,j) = p5*sqrt(pbt(i,j  ,tau) + pbt(i+1,j  ,tau) + &
                          pbt(i,j+1,tau) + pbt(i+1,j+1,tau))
      endif
!      spbt(1  ,j) = spbt(imm,j)
!      spbt(imt,j) = spbt(2  ,j)
      enddo
      enddo
      call swap_array_real2d(spbt,imt,jmt,west,east,north,south)
!
!
!-----------------------------------------------------------------------
!     calculate the pressure gradients due to change of free surface
!-----------------------------------------------------------------------
!
!
      if (boussinesq==1) then
      do j=1,jmt
      do i=1,imt
      a(i,j) = phib(i,j)*(c1-pbt(i,j,tau))
      enddo
      enddo
      end if
      
      if (boussinesq==1) then
      do j=2,jmm
      do i=2,imm
      if(ivn(i,j).gt.0)then
!
      etax(i,j) = spbt(i,j)*p5*(  &
              (a(i+1,j  )-a(i,j  ))*pdxn(i,j) &
            +(a(i+1,j+1)-a(i,j+1))*pdxs(i,j+1) + &
             (pbt(i+1,j  ,tau)-pbt(i,j  ,tau))*pbxn(i,j)   &
            +(pbt(i+1,j+1,tau)-pbt(i,j+1,tau))*pbxs(i,j+1) +  &
             (pbt(i+1,j  ,tau)+pbt(i,j  ,tau))*pcxn(i,j)    &
            +(pbt(i+1,j+1,tau)+pbt(i,j+1,tau))*pcxs(i,j+1) )
!
!
      etay(i,j) = spbt(i,j)*p5*(   &
              (a(i  ,j+1)-a(i  ,j))*pdye(i,j)   &
            +(a(i+1,j+1)-a(i+1,j))*pdyw(i+1,j) +  &
             (pbt(i  ,j+1,tau)-pbt(i,  j,tau))*pbye(i  ,j)  &
            +(pbt(i+1,j+1,tau)-pbt(i+1,j,tau))*pbyw(i+1,j) +  &
             (pbt(i  ,j+1,tau)+pbt(i,  j,tau))*pcye(i  ,j)  &
            +(pbt(i+1,j+1,tau)+pbt(i+1,j,tau))*pcyw(i+1,j) )  
!
      endif
      enddo
      enddo
      
      else
      
      do j=2,jmm
      do i=2,imm
      if(ivn(i,j).gt.0)then
!
      etax(i,j) = spbt(i,j)*p5*(  &
              phibx(i,j) + &
             (pbt(i+1,j  ,tau)-pbt(i,j  ,tau))*pbxn(i,j)   &
            +(pbt(i+1,j+1,tau)-pbt(i,j+1,tau))*pbxs(i,j+1) +  &
             (pbt(i+1,j  ,tau)+pbt(i,j  ,tau))*pcxn(i,j)    &
            +(pbt(i+1,j+1,tau)+pbt(i,j+1,tau))*pcxs(i,j+1) )
!
!
      etay(i,j) = spbt(i,j)*p5*(   &
              phiby(i,j) +  &
             (pbt(i  ,j+1,tau)-pbt(i,  j,tau))*pbye(i  ,j)  &
            +(pbt(i+1,j+1,tau)-pbt(i+1,j,tau))*pbyw(i+1,j) +  &
             (pbt(i  ,j+1,tau)+pbt(i,  j,tau))*pcye(i  ,j)  &
            +(pbt(i+1,j+1,tau)+pbt(i+1,j,tau))*pcyw(i+1,j) )  
!
      endif
      enddo
      enddo
      end if
!
!-----------------------------------------------------------------------
!     advection + viscosiy + pressure gradient + coriolis
!-----------------------------------------------------------------------
!
      do j=2,jmm
      do i=2,imm
      if(ivn(i,j).gt.0)then
        a(i,j) = dub(i,j) + adu(i,j) + ff(j)*vpb(i,j,tau) - etax(i,j)
        b(i,j) = dvb(i,j) + adv(i,j) - ff(j)*upb(i,j,tau) - etay(i,j)
      else
        a(i,j) = c0
        b(i,j) = c0
      endif
      enddo
      enddo
!
!
!-----------------------------------------------------------------------
!     coriolis adjustment
!-----------------------------------------------------------------------
!
      if(leapfrog_b)then
        do j=2,jmm
        do i=2,imm
        du2(i,j) = ebla(j)*a(i,j) + eblb(j)*b(i,j)
        dv2(i,j) = ebla(j)*b(i,j) - eblb(j)*a(i,j)
        enddo
        enddo
      else
        do j=2,jmm
        do i=2,imm
        du2(i,j) = ebea(j)*a(i,j) + ebeb(j)*b(i,j)
        dv2(i,j) = ebea(j)*b(i,j) - ebeb(j)*a(i,j)
        enddo
        enddo
      endif
!
!
!-----------------------------------------------------------------------
!     calculate the change of Pbt
!-----------------------------------------------------------------------
!
      do j=1,jmt
      do i=1,imt
      a(i,j) = zu(i,j)*spbt(i,j)*upb(i,j,tau)
      b(i,j) = zu(i,j)*spbt(i,j)*vpb(i,j,tau)*cosu(j)
      enddo
      enddo
!
!
      if (snbc==1) then
      do j=2,jmm
      do i=2,imm
      if(itn(i,j).gt.0)then
       dp(i,j) = ( - p5*(rdxt(j)*(a(i,j)+a(i,j-1)-a(i-1,j)-a(i-1,j-1)) &
                       +rdyt(j)*(b(i,j)+b(i-1,j)-b(i,j-1)-b(i-1,j-1))) &
                  - emp(i,j) &
                  )/pn(i,j)
      else
       dp(i,j) = c0
      endif
      enddo
      enddo
      else
      do j=2,jmm
      do i=2,imm
      if(itn(i,j).gt.0)then
       dp(i,j) = ( - p5*(rdxt(j)*(a(i,j)+a(i,j-1)-a(i-1,j)-a(i-1,j-1)) &
                       +rdyt(j)*(b(i,j)+b(i-1,j)-b(i,j-1)-b(i-1,j-1))) &
                  )/pn(i,j)
      else
       dp(i,j) = c0
      endif
      enddo
      enddo
      end if
!
!
      if (smth==1) then
      call swap_array_real2d(du2,imt,jmt,west,east,north,south)
      call swap_array_real2d(dv2,imt,jmt,west,east,north,south)
      call smths(du2,umask,jsts,jeds,imt,jmt,km)
      call smths(dv2,umask,jsts,jeds,imt,jmt,km)
      call smths(du2,umask,jstn,jedn,imt,jmt,km)
      call smths(dv2,umask,jstn,jedn,imt,jmt,km)
      call swap_array_real2d(du2,imt,jmt,west,east,north,south)
      call swap_array_real2d(dv2,imt,jmt,west,east,north,south)
      end if
!
!
!
!-----------------------------------------------------------------------
!     compute pbt, upb & vpb at tau+1 time level
!-----------------------------------------------------------------------
!
      if(leapfrog_b) go to 120
!
       do j=2,jmm
       do i=2,imm
       upb(i,j,tau) = upb(i,j,taum) + du2(i,j)*dtsf
       vpb(i,j,tau) = vpb(i,j,taum) + dv2(i,j)*dtsf
       pbt(i,j,tau) = pbt(i,j,taum) + dp(i,j) *dtsf
       enddo
       enddo
!
       if(euler_back) then
        euler_back = .false.
        go to 1001
       endif
!
       leapfrog_b = .true.
!
       go to 150
!
!
120   continue
!
!
      if (asselin_b==1) then
      do j=2,jmm
      do i=2,imm
      t1            = pbt(i,j,taum) + dp (i,j)*c2dtsf
      t2            = upb(i,j,taum) + du2(i,j)*c2dtsf
      t3            = vpb(i,j,taum) + dv2(i,j)*c2dtsf
      pbt(i,j,taum) = afb2*pbt(i,j,tau)+afb1*(pbt(i,j,taum)+t1)
      upb(i,j,taum) = afb2*upb(i,j,tau)+afb1*(upb(i,j,taum)+t2)
      vpb(i,j,taum) = afb2*vpb(i,j,tau)+afb1*(vpb(i,j,taum)+t3)
      pbt(i,j,tau)  = t1
      upb(i,j,tau)  = t2
      vpb(i,j,tau)  = t3
      enddo
      enddo
      else
      do j=2,jmm
      do i=2,imm
      t1            = pbt(i,j,taum) + dp (i,j)*c2dtsf
      t2            = upb(i,j,taum) + du2(i,j)*c2dtsf
      t3            = vpb(i,j,taum) + dv2(i,j)*c2dtsf
      pbt(i,j,taum) = pbt(i,j,tau)
      upb(i,j,taum) = upb(i,j,tau)
      vpb(i,j,taum) = vpb(i,j,tau)
      pbt(i,j,tau)  = t1
      upb(i,j,tau)  = t2
      vpb(i,j,tau)  = t3
      enddo
      enddo
      end if
!
150   continue
!
!
!      do j=2,jmm
!      upb(1  ,j,tau) = upb(imm,j,tau)
!      upb(imt,j,tau) = upb(2  ,j,tau)
!      vpb(1  ,j,tau) = vpb(imm,j,tau)
!      vpb(imt,j,tau) = vpb(2  ,j,tau)
!      pbt(1  ,j,tau) = pbt(imm,j,tau)
!      pbt(imt,j,tau) = pbt(2  ,j,tau)
!      upb(1  ,j,taum) = upb(imm,j,taum)
!      upb(imt,j,taum) = upb(2  ,j,taum)
!      vpb(1  ,j,taum) = vpb(imm,j,taum)
!      vpb(imt,j,taum) = vpb(2  ,j,taum)
!      pbt(1  ,j,taum) = pbt(imm,j,taum)
!      pbt(imt,j,taum) = pbt(2  ,j,taum)
!      enddo
      call swap_array_real3d(upb,imt,jmt,2,west,east,north,south)
      call swap_array_real3d(vpb,imt,jmt,2,west,east,north,south)
      call swap_array_real3d(pbt,imt,jmt,2,west,east,north,south)
      
!
!
!-----------------------------------------------------------------------
!     constract time-averaged pbt for solving baroclinc and tracer Eqs.
!-----------------------------------------------------------------------
      do j=1,jmt
      do i=1,imt
      pmup(i,j) = pmup(i,j) + pbt(i,j,tau)
      pmtp(i,j) = pmtp(i,j) + pbt(i,j,tau)
      enddo
      enddo
!
!
1000  continue
!
      return
      end
