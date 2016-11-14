!
!     =================
      subroutine barotr_st(ivn,itn,upb,vpb,r1c,r1d,sdxu,am,dub,dvb,pbt_st,phibx,    &
                    phiby,pbxn,pbxs,pbye,pbyw,pcxn,pcxs,pcye,pcyw,ff,ebla,     &
                    eblb,ebea,ebeb,pn,zu,cosu,rdxt,rdyt, &
                    dtsf,nbb,imt,jmt,km,imm,jmm,    &
                    myid,west,east,north,south,snbc,emp,jstn,jedn,jsts,   &
                    jeds,smtha,fcof,umask,boussinesq,phib,pdxn,pdxs,pdye,pdyw,  &
                    lat,lon,energydiag)
!     =================
!     compute pbt, upb & vpb at "tau+1" time level
!
      implicit none
      include 'pconst.h'
      include 'mpif.h'
!
      integer mode_b,nbb,snbc,boussinesq,nannum,energydiag
      integer imt,jmt,km,imm,jmm,i,j,k,k2,n
      integer jstn,jedn,jsts,jeds,smtha
      integer ivn(imt,jmt),itn(imt,jmt)
      real t1,t2,t3,am(imt,jmt,km),dtsf,fcof
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
      real pbt_st(imt,jmt,4),pbt(imt,jmt,2),spbt(imt,jmt)
      real pbt_st2(imt,jmt),pbt_st3(imt,jmt)
      real upb_st(imt,jmt),vpb_st(imt,jmt)
      real dub(imt,jmt),dvb(imt,jmt)
      real phibx(imt,jmt),phiby(imt,jmt),emp(imt,jmt)
      real lat(jmt),lon(imt)
      
      integer west,east,north,south,myid
!
      real gamma_b
      
      gamma_b=0.2
      
      do j=1,jmt
         do i=1,imt
            pbt(i,j,taum)=pbt_st(i,j,1)
            pbt_st2(i,j)=pbt(i,j,taum)
            pbt_st3(i,j)=pbt(i,j,taum)
            upb_st(i,j)=upb(i,j,tau)
            vpb_st(i,j)=vpb(i,j,tau)
         end do
      end do

      do 1000 mode_b = 1,nbb*c2
      
!-----------------------------------------------------------------------
!     calculate the change of Pbt at "predicts step"
!-----------------------------------------------------------------------
!
      do j=2,jmm
      do i=2,imm
      if(ivn(i,j).gt.0)then
      spbt(i,j) = p5*sqrt(pbt(i,j  ,taum) + pbt(i+1,j  ,taum) + &
                          pbt(i,j+1,taum) + pbt(i+1,j+1,taum))
      endif
      enddo
      enddo
      
      call swap_array_real2d(spbt,imt,jmt,west,east,north,south)

!-----check-------
	  nannum=0
      do j=1,jmt
        do i=1,imt
          if (spbt(i,j).ne.spbt(i,j)) then
            print "(a24,f7.2,a5,f7.2,a9,i4,a6,i4,a8,g4.1,4f7.2,2i4)", "Something Wrong at lat ", &
            lat(j),"lon",lon(i),"mode_b =",mode_b,"pro =",myid,"spbt =",spbt(i,j), &
            pbt(i,j,tau),pbt(i+1,j,tau),pbt(i,j+1,tau),pbt(i+1,j+1,tau),i,j
            nannum=nannum+1
          end if
        end do
      end do
      
      if (nannum.gt.0) then
         stop
      end if
!-----------------

      do j=1,jmt
      do i=1,imt
      a(i,j) = zu(i,j)*spbt(i,j)*upb(i,j,taum)
      b(i,j) = zu(i,j)*spbt(i,j)*vpb(i,j,taum)*cosu(j)
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
      do j=2,jmm
         do i=2,imm
            pbt(i,j,tau) = pbt(i,j,taum) + dp(i,j) *dtsf *gamma_b
         end do
      end do
      
      call swap_array_real3d(pbt,imt,jmt,2,west,east,north,south)
      
      do j=2,jmm
         do i=2,imm
            if(ivn(i,j).gt.0)then
               spbt(i,j) = p5*sqrt(pbt(i,j  ,tau) + pbt(i+1,j  ,tau) + &
                                   pbt(i,j+1,tau) + pbt(i+1,j+1,tau))
            end if
         end do
      end do
      
      call swap_array_real2d(spbt,imt,jmt,west,east,north,south)
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
        a(i,j) = dub(i,j) + adu(i,j) + p5*ff(j)*vpb(i,j,taum) - etax(i,j)
        b(i,j) = dvb(i,j) + adv(i,j) - p5*ff(j)*upb(i,j,taum) - etay(i,j)
      else
        a(i,j) = c0
        b(i,j) = c0
      endif
      enddo
      enddo
!
      do j=2,jmm
         do i=2,imm
            a(i,j) = upb(i,j,taum) + a(i,j)*dtsf
            b(i,j) = vpb(i,j,taum) + b(i,j)*dtsf
         end do
      end do
!
!-----------------------------------------------------------------------
!     coriolis adjustment
!-----------------------------------------------------------------------
!
      do j=2,jmm
         do i=2,imm
            upb(i,j,tau) = (a(i,j)+(p5*ff(j)*dtsf)*b(i,j))/(c1+(p5*ff(j)*dtsf)**c2)
            vpb(i,j,tau) = (b(i,j)-(p5*ff(j)*dtsf)*a(i,j))/(c1+(p5*ff(j)*dtsf)**c2)
         end do
      end do
      
      do j=2,jmm
         do i=2,imm
            upb(i,j,tau) = p5*(upb(i,j,taum)+upb(i,j,tau))
            vpb(i,j,tau) = p5*(vpb(i,j,taum)+vpb(i,j,tau))
         end do
      end do
      
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
!-----------------------------------------------------------------------
!     compute upb & vpb at tau+1 time level
!-----------------------------------------------------------------------      
      do j=2,jmm
         do i=2,imm
            upb(i,j,tau) = upb(i,j,taum) + a(i,j)*dtsf
            vpb(i,j,tau) = vpb(i,j,taum) + b(i,j)*dtsf
            upb(i,j,taum)= upb(i,j,tau)
            vpb(i,j,taum)= vpb(i,j,tau)
         end do
      end do

      call swap_array_real3d(upb,imt,jmt,2,west,east,north,south)
      call swap_array_real3d(vpb,imt,jmt,2,west,east,north,south)
!
!-----------------------------------------------------------------------
!     calculate the change of Pbt at "corrector step (tau+1)"
!-----------------------------------------------------------------------

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
      do j=2,jmm
         do i=2,imm
            pbt(i,j,tau) = pbt(i,j,taum) + dp(i,j) *dtsf
            pbt(i,j,taum)= pbt(i,j,tau)
         end do
      end do
      
      call swap_array_real3d(pbt,imt,jmt,2,west,east,north,south)
      
      do j=1,jmt
         do i=1,imt
            pbt_st3(i,j)=pbt_st3(i,j)+pbt(i,j,tau)
            upb_st(i,j)=upb_st(i,j)+upb(i,j,tau)
            vpb_st(i,j)=vpb_st(i,j)+vpb(i,j,tau)
         end do
      end do
      
      if (mode_b.le.nbb) then
         do j=1,jmt
            do i=1,imt
               pbt_st2(i,j)=pbt_st2(i,j)+pbt(i,j,tau)
            end do
         end do
      end if
!
1000  continue
!
      do j=1,jmt
         do i=1,imt
            pbt_st(i,j,2)=pbt_st2(i,j)/(nbb+c1)
            pbt_st(i,j,3)=pbt_st3(i,j)/(nbb*c2+c1)
            upb(i,j,tau)=upb_st(i,j)/(nbb*c2+c1)
            vpb(i,j,tau)=vpb_st(i,j)/(nbb*c2+c1)
            upb(i,j,taum)=upb(i,j,tau)
            vpb(i,j,taum)=vpb(i,j,tau)
         end do
      end do

      return
      end subroutine barotr_st
