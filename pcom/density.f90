!     ==================
!BOP
!
! !MODULE: density.f90
! !DESCRIPTION: \input{sections/code-density}
!
! !INTERFACE:
!
      subroutine density(tmask,z,t,pbt,bcp,rho,decibar,imt,jmt,km,nt,imm,jmm,fixp,  &
                         west,east,north,south)
!EOP
!-------------------------------------------------------------------------------
!     ==================
!     calculate density/rho_0(BCOM) or reciprocal of density(PCOM)
!
      implicit none
      include 'pconst.h'
      include 'mpif.h'
!
      integer imt,jmt,km,nt,imm,jmm,i,j,k
      real t0,s0,p0,dens,rdens,unrdens,undens,decibar
      real tmask(imt,jmt,km),z(km),pbt(imt,jmt,2),rho(imt,jmt,km),fixp(imt,jmt,km)
      real t(imt,jmt,km,nt,2)
      real bcp(imt,jmt)
      
      integer west,east,north,south
!
!     calculate density/rho_0 or reciprocal of density
!
      
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
      
      call swap_array_real3d(rho,imt,jmt,km,west,east,north,south)
!
      return
      end
!
!
!     ==================
!BOP
!
! !IROUTINE: rho_ref
! !DESCRIPTION: \input{sections/code-rho_ref}
!
! !INTERFACE:
!
      subroutine rho_ref(tmask,t,z,pbt,pt,ps,deltat,deltas,rdeltat,rdeltas,decibar,fixp, &
                         imt,jmt,km,nt,imm,jmm,west,east,north,south)
!EOP
!-------------------------------------------------------------------------------
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
      integer imt,jmt,km,nt,imm,jmm,i,j,k,i_id,j_id
      real t1,s1,p1,dens,undens
      real rhop,rhoq,decibar,deltas,rdeltat,rdeltas,deltat
      real tmask(imt,jmt,km),z(km),pbt(imt,jmt,2),pt(imt,jmt,km),ps(imt,jmt,km)
      real t(imt,jmt,km,nt,2),fixp(imt,jmt,km)
      
      integer west,east,north,south
!
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
      
      call swap_array_real3d(pt,imt,jmt,km,west,east,north,south)
      call swap_array_real3d(ps,imt,jmt,km,west,east,north,south)
!
      return
      end subroutine rho_ref
!
!
!     ==================
      subroutine rho_ref_st(tmask,t,z,pbt_st,pt,ps,deltat,deltas,rdeltat,rdeltas,decibar,fixp, &
                         imt,jmt,km,nt,imm,jmm,west,east,north,south)
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
      integer imt,jmt,km,nt,imm,jmm,i,j,k,i_id,j_id
      real t1,s1,p1,dens,undens
      real rhop,rhoq,decibar,deltas,rdeltat,rdeltas,deltat
      real tmask(imt,jmt,km),z(km),pbt_st(imt,jmt,4),pt(imt,jmt,km),ps(imt,jmt,km)
      real t(imt,jmt,km,nt,2),fixp(imt,jmt,km)
      
      integer west,east,north,south
!
      do k=1,km
      do j=2,jmm
      do i=2,imm
      if(tmask(i,j,k).gt.c0)then
        t1        = t(i,j,k,1,tau)
        s1        = t(i,j,k,2,tau)
        p1        = pbt_st(i,j,4)*z(k)*decibar
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
      
      call swap_array_real3d(pt,imt,jmt,km,west,east,north,south)
      call swap_array_real3d(ps,imt,jmt,km,west,east,north,south)
!
      return
      end subroutine rho_ref_st

!     ====================
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!BOP
!
! !IROUTINE: undens
! !DESCRIPTION: \input{sections/code-undens}
!
! !INTERFACE:
!
      function undens(t,s,p0)
!EOP
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
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
