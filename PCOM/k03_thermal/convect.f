c
c     ==================
      subroutine convect
c     ==================
c
c     a full convective adjustment scheme, based on GFDL's version
c
c     ---------------------------------------------------------------
c     kcon = maximum number of levels at this location
c     lcon = counts levels down
c     lcona = upper layer of a convective part of water column
c     lconb = lower layer of a convective part of water column
c     rhoup = density anomaly referenced to same level
c     rholo = density anomaly referenced to level below
c     dztsum = sum of layer thicknesses
c     trasum = sum of layer tracer values
c     tramix = mixed tracer value after convection
c     lcven = number of levels ventilated (convection to surface)
c     ---------------------------------------------------------------
c
      implicit none
      include 'param.h'
      include 'grdvar.h'
      include 'prog.h'
c
      real     rhoup(km),rholo(km),trasum(2)
      real     tup,sup,tlo,slo,dztsum,tramix
      integer  kcon,lcven,l1,l,lcon,lcona,lconb,lmix
c
c
      do j=2,jmm
      do i=2,imm
c
        kcon = itn(i,j)
        if (kcon .eq. 0) goto 1310
c
          lcven = 1
          lcon  = 0
c
          do l=1,kcon-1
            l1        = l+1
            tup       = t(i,j,l1,1,tau)
            sup       = t(i,j,l1,2,tau)
            tlo       = t(i,j, l,1,tau)
            slo       = t(i,j, l,2,tau)
            rhoup(l1) = pt(i,j,l1)*tup + ps(i,j,l1)*sup
            rholo(l)  = pt(i,j,l1)*tlo + ps(i,j,l1)*slo
          enddo
c
c
c         1. initial search for uppermost unstable pair; if none is
c            found, move on to next column
c
          do k=kcon-1,1,-1
            if (rholo(k) .gt. rhoup(k+1)) lcon = k
          enddo
c
          if (lcon .eq. 0) goto 1310
c
1319      lcona = lcon
          lconb = lcon + 1
c
c         2. mix the first two unstable layers
c
          dztsum               = dz(lcona) + dz(lconb)
          do n=1,2
            trasum(n)          = t(i,j,lcona,n,tau)*dz(lcona) + 
     &                           t(i,j,lconb,n,tau)*dz(lconb)
            tramix             = trasum(n) / dztsum
            t(i,j,lcona,n,tau) = tramix
            t(i,j,lconb,n,tau) = tramix
          enddo
c
c         3. test layer below lconb
c
1306      continue
          if (lconb .eq. kcon) goto 1308
c
          l1 = lconb + 1
          rholo(lconb) = pt(i,j,l1)*t(i,j,lconb,1,tau) +
     &                   ps(i,j,l1)*t(i,j,lconb,2,tau)
c
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
c
c         4. test layer above lcona
c
1308      continue
          if (lcona .gt. 1) then
            l1 = lcona-1
            rholo(l1)    = pt(i,j,lcona)*t(i,j,l1,1,tau)  +
     &                     ps(i,j,lcona)*t(i,j,l1,2,tau)
            rhoup(lcona) = pt(i,j,lcona)*t(i,j,lcona,1,tau)  +
     &                     ps(i,j,lcona)*t(i,j,lcona,2,tau)
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
c
c         5. remember the total number of levels mixed by convection
c            in this water column, as well as the ventilated column
c
          if (lcona .eq. 1) lcven = lconb - lcona + 1
c
c         6. resume search if step 3. and 4. have been passed and this
c            unstable part of the water column has thus been removed,
c            i.e. find further unstable areas further down the column
c
          if (lconb .eq. kcon) goto 1310
          lcon = lconb
c
1302      continue
          lcon = lcon + 1
c
          if (lcon .eq. kcon) goto 1310
c
          if (rholo(lcon) .le. rhoup(lcon+1)) goto 1302
c
          goto 1319

1310  continue
c
      enddo
      enddo
c
      do n=1,2
      do k=1,km
      do j=2,jmm
      t(1  ,j,k,n,tau) = t(imm,j,k,n,tau)
      t(imt,j,k,n,tau) = t(2  ,j,k,n,tau)
      enddo
      enddo
      enddo
c
      return
      end
