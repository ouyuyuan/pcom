!
!     ==================
      subroutine convect(t,pt,ps,itn,dz,imt,jmt,km,nt,imm,jmm,myid,west,east,north,south,  &
                   energydiag,dpo_con,din_con,p,pbar,mass_up,bottom_h,total_in,total_po)
!     ==================
!
!     a full convective adjustment scheme, based on GFDL's version
!
!     ---------------------------------------------------------------
!     kcon = maximum number of levels at this location
!     lcon = counts levels down
!     lcona = upper layer of a convective part of water column
!     lconb = lower layer of a convective part of water column
!     rhoup = density anomaly referenced to same level
!     rholo = density anomaly referenced to level below
!     dztsum = sum of layer thicknesses
!     trasum = sum of layer tracer values
!     tramix = mixed tracer value after convection
!     lcven = number of levels ventilated (convection to surface)
!     ---------------------------------------------------------------
!
      implicit none
      include 'pconst.h'
      include 'mpif.h'
!
      integer imt,jmt,km,nt,imm,jmm,i,j,k,kk,n
      real     rhoup(km),rholo(km),trasum(2)
      real     tup,sup,tlo,slo,dztsum,tramix
      integer  kcon,lcven,l1,l,lcon,lcona,lconb,lmix
      
      real t(imt,jmt,km,nt,2),pt(imt,jmt,km),ps(imt,jmt,km),dz(km)
      integer itn(imt,jmt)
      
      integer myid,west,east,north,south
      
      integer energydiag
      real dpo_con(imt,jmt,km),din_con(imt,jmt,km)
      real p(imt,jmt,km),pbar(imt,jmt),gr_h(imt,jmt,km),gr_h2(imt,jmt,km),mass_up(imt,jmt,km)
      real total_in(imt,jmt,km),total_po(imt,jmt,km),bottom_h(imt,jmt)
      real rho_temp,t0,s0,p0,dz_temp,unrdens,gsw_internal_energy
      
      real din_test0,din_test1,din_test2,din_test3,din_test4
!
!
      if (energydiag==1) then
         do j=2,jmm
            do i=2,imm
               do k=1,itn(i,j)
                  t0=t(i,j,k,1,tau)
                  s0=t(i,j,k,2,tau)
                  p0=p(i,j,k)
                  rho_temp=unrdens(t0,s0,p0)
                  gr_h(i,j,k) = dz(k)*pbar(i,j)*rho_temp/grav
                  din_con(i,j,k)=total_in(i,j,k)
               end do
            end do
         end do
      end if

      do j=2,jmm
      do i=2,imm
!
        kcon = itn(i,j)
        if (kcon .eq. 0) goto 1310
!
          lcven = 1
          lcon  = 0
!
          do l=1,kcon-1
            l1        = l+1
            tup       = t(i,j,l1,1,tau)
            sup       = t(i,j,l1,2,tau)
            tlo       = t(i,j, l,1,tau)
            slo       = t(i,j, l,2,tau)
            rhoup(l1) = pt(i,j,l1)*tup + ps(i,j,l1)*sup
            rholo(l)  = pt(i,j,l1)*tlo + ps(i,j,l1)*slo
          enddo
!
!
!         1. initial search for uppermost unstable pair; if none is
!            found, move on to next column
!
          do k=kcon-1,1,-1
            if (rholo(k) .gt. rhoup(k+1)) lcon = k
          enddo
!
          if (lcon .eq. 0) goto 1310
!
1319      lcona = lcon
          lconb = lcon + 1
!
!         2. mix the first two unstable layers
!
          dztsum               = dz(lcona) + dz(lconb)
          do n=1,2
            trasum(n)          = t(i,j,lcona,n,tau)*dz(lcona) +  &
                                 t(i,j,lconb,n,tau)*dz(lconb)
            tramix             = trasum(n) / dztsum
            t(i,j,lcona,n,tau) = tramix
            t(i,j,lconb,n,tau) = tramix
          enddo
!
!         3. test layer below lconb
!
1306      continue
          if (lconb .eq. kcon) goto 1308
!
          l1 = lconb + 1
          rholo(lconb) = pt(i,j,l1)*t(i,j,lconb,1,tau) +  &
                         ps(i,j,l1)*t(i,j,lconb,2,tau)
!
          if (rholo(lconb) .gt. rhoup(l1)) then
            lconb  = lconb+1
            dztsum = dztsum + dz(lconb)
            do n=1,2
              trasum(n) = trasum(n) + t(i,j,lconb,n,tau)*dz(lconb)
              tramix = trasum(n) / dztsum
              do lmix=lcona,lconb
                t(i,j,lmix,n,tau) = tramix
              enddo
            enddo
            goto 1306
          end if
!
!         4. test layer above lcona
!
1308      continue
          if (lcona .gt. 1) then
            l1 = lcona-1
            rholo(l1)    = pt(i,j,lcona)*t(i,j,l1,1,tau)  +  &
                           ps(i,j,lcona)*t(i,j,l1,2,tau)
            rhoup(lcona) = pt(i,j,lcona)*t(i,j,lcona,1,tau)  +  &
                           ps(i,j,lcona)*t(i,j,lcona,2,tau)
            if (rholo(lcona-1) .gt. rhoup(lcona)) then
              lcona = lcona-1
              dztsum = dztsum + dz(lcona)
              do n=1,2
                trasum(n) = trasum(n) + t(i,j,lcona,n,tau)*dz(lcona)
                tramix    = trasum(n) / dztsum 
                do lmix=lcona,lconb
                  t(i,j,lmix,n,tau) = tramix
                enddo
              enddo
              goto 1306
            end if
          end if
!
!         5. remember the total number of levels mixed by convection
!            in this water column, as well as the ventilated column
!
          if (lcona .eq. 1) lcven = lconb - lcona + 1
!
!         6. resume search if step 3. and 4. have been passed and this
!            unstable part of the water column has thus been removed,
!            i.e. find further unstable areas further down the column
!
          if (lconb .eq. kcon) goto 1310
          lcon = lconb
!
1302      continue
          lcon = lcon + 1
!
          if (lcon .eq. kcon) goto 1310
!
          if (rholo(lcon) .le. rhoup(lcon+1)) goto 1302
!
          goto 1319

1310  continue
!
      enddo
      enddo
      
      do n=1,nt
      do k=1,km
      do j=2,jmm
      do i=2,imm
      t(i,j,k,n,taum) = t(i,j,k,n,tau)
      enddo
      enddo
      enddo
      end do
!
      call swap_array_real5d(t,imt,jmt,km,nt,2,west,east,north,south)
!
      if (energydiag==1) then
         do j=2,jmm
            do i=2,imm
               do k=1,itn(i,j)
                  t0=t(i,j,k,1,tau)
                  s0=t(i,j,k,2,tau)
                  p0=p(i,j,k)
                  rho_temp=unrdens(t0,s0,p0)
                  gr_h2(i,j,k)=dz(k)*pbar(i,j)*rho_temp/grav
                  dz_temp=gr_h2(i,j,k)-gr_h(i,j,k)
                  dpo_con(i,j,k)=dz_temp*mass_up(i,j,k)*grav
                  total_in(i,j,k)=gsw_internal_energy(s0,t0,p0)
                  din_con(i,j,k)=total_in(i,j,k)-din_con(i,j,k)
               end do
            end do
         end do
         do j=2,jmm
            do i=2,imm
               do k=1,itn(i,j)
                  dz_temp=p5*gr_h2(i,j,k)
                  do kk=k+1,itn(i,j)
                     dz_temp=dz_temp+gr_h2(i,j,kk)
                  end do
                  dz_temp=dz_temp+bottom_h(i,j)
                  total_po(i,j,k)=dz_temp*grav
               end do
            end do
         end do
      end if
      

      return
      end
