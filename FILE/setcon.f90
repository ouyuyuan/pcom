!     =================
      subroutine setcon(fam,fah,fkm,fkh,am_c,ah_c,km_c,kh_c,am,ah,kappa_m,kappa_h,    &
                        gravr,decibar,deltap,deltat,deltas,rdeltap,rdeltat,rdeltas,   &
                        gamma_t,gamma_s,imt,jmt,km,boussinesq,myid,ncpux,ncpuy,simt,  &
                        sjmt,mat_myid)
                        
!     =================
!     set scalar quantities
!
      implicit none
      include 'pconst.h'
      include 'mpif.h'
!
!---------------------------------------------------------------------
!     calendar
!---------------------------------------------------------------------
      integer imt,jmt,km,i,j,k,boussinesq
      integer fam,fah,fkm,fkh
      real am_c,ah_c,km_c,kh_c
      real gravr,decibar,deltap,rdeltap,deltat,deltas
      real rdeltat,rdeltas,gamma_t,gamma_s
      real kappa_h(imt,jmt,km),kappa_m(imt,jmt,km),am(imt,jmt,km),ah(imt,jmt,km)
      
      integer myid,ncpux,ncpuy,simt,sjmt
      integer mat_myid(ncpux+2,ncpuy)
      real skappa_h(simt,sjmt,km),skappa_m(simt,sjmt,km)
      real sam(simt,sjmt,km),sah(simt,sjmt,km)

!      data monname /'jan','feb','mar','apr','may','jun', &
!                    'jul','aug','sep','oct','nov','dec'/
!     data daymd   /16,14,16,15,16,15,16,16,15,16,15,16/
!     data daypm   /31,28,31,30,31,30,31,31,30,31,30,31/
!      data daymd   /15,15,15,15,15,15,15,15,15,15,15,15/
!      data daypm   /30,30,30,30,30,30,30,30,30,30,30,30/
!      data kappa_h /10,5,1,0.5,0.1,
!     *              0.1,0.1,0.1,0.1,0.1,
!     *              0.1,0.1,0.1,0.1,0.1,
!     *              0.1,0.1,0.1,0.1,0.1,
!     *              0.1,0.1,0.1,0.1,0.1,
!     *              0.1,0.1,0.1,0.1,0.1,
!     *              0.1,0.1,0.1,0.1,0.1,
!     *              0.1,0.1,0.1,0.1,0.1,
!     *              0.1,0.1,0.1,0.1,0.1,
!     *              0.1,0.1,0.1,0.1,0.1,
!     *              0.1,0.1,0.1,0.1,0.1,
!     *              0.1,0.1,0.1,0.1,0.1/
!
!      open (15,file='vprofile.data',status='old',form='formatted')
     if (myid==0) then
     
     if (fam==1) then
     open(15,file='am.data',form='unformatted',access='direct', &
          recl=simt*sjmt*km*8,status='old')
     read(15,rec=1) sam
     close(15)
     else
     do k=1,km
     do j=1,sjmt
     do i=1,simt
     sam(i,j,k)=am_c
     end do
     end do
     end do
     end if
     
     if (fah==1) then
     open(15,file='ah.data',form='unformatted',access='direct', &
          recl=simt*sjmt*km*8,status='old')
     read(15,rec=1) sah
     close(15)
     else
     do k=1,km
     do j=1,sjmt
     do i=1,simt
     sah(i,j,k)=ah_c
     end do
     end do
     end do
     end if
     
     if (fkm==1) then
     open(15,file='km.data',form='unformatted',access='direct', &
          recl=simt*sjmt*km*8,status='old')
     read(15,rec=1) skappa_m
     close(15)
     else
     do k=1,km
     do j=1,sjmt
     do i=1,simt
     skappa_m(i,j,k)=km_c
     end do
     end do
     end do
     end if
     
     if (fkh==1) then
     open(15,file='kh.data',form='unformatted',access='direct', &
          recl=simt*sjmt*km*8,status='old')
     read(15,rec=1) skappa_h
     close(15)
     else
     do k=1,km
     do j=1,sjmt
     do i=1,simt
     skappa_h(i,j,k)=kh_c
     end do
     end do
     end do
     end if
     
     end if

     call div_array_real3d(sam,am,mat_myid,ncpux,ncpuy,simt,sjmt,  &
                           km,imt,jmt,myid)
     call div_array_real3d(sah,ah,mat_myid,ncpux,ncpuy,simt,sjmt,  &
                           km,imt,jmt,myid)
     call div_array_real3d(skappa_m,kappa_m,mat_myid,ncpux,ncpuy,simt,sjmt,  &
                           km,imt,jmt,myid)
     call div_array_real3d(skappa_h,kappa_h,mat_myid,ncpux,ncpuy,simt,sjmt,  &
                           km,imt,jmt,myid)
!     do k=1,km
!      do j=1,jmt
!       do i=1,imt
!         kappa_h(i,j,k)=kapp_h
!       enddo
!         write (6,97) j,k,kappa_h(j,k)
!     enddo
!     enddo
! 97      format (2i5,f12.6)
!
!---------------------------------------------------------------------
!     physical parameters
!---------------------------------------------------------------------
!     am      = constant lateral viscosity coeff for momentum
!     ah      = constant lateral diffusion coeff for tracers
!     kappa_m = constant vertical viscosity coefficient (cm**2/sec) 
!     kappa_h = constant vertical diffusion coefficient (cm**2/sec)
!     cdbot   = parameter for bottom drag
!
!      am      = 3.0e8
!      ah      = 1.0e7
!      kappa_m = 1.0
!cccccccccccccccccccccccccccccccccccccccc reset to 0.1      kappa_h = 0.3
!      do k=1,km
!         kappa_h(k) = 0.1
!      enddo
!      cdbot   = 2.6e-3
!
      if (boussinesq==1) then
      gravr = grav
      else
      gravr = grav*rho_0
      end if
!
!     change the unit of pressure from dynes/cm**2 to decibars
!     1 dynes/cm**2 = 0.1 N/m**2 = 1.0e-3 mbar = 1.0e-5 decibars
!
      decibar = 1.0e-5
!
!     deltax is used to compute {partial rho}/{partial x}, x=p,t,s
      deltap  = 1.0e-2
      deltat  = 1.0e-2
      deltas  = 1.0e-3
      rdeltap = p5/deltap
      rdeltat = p5/deltat
      rdeltas = p5/deltas
!
!
!---------------------------------------------------------------------
!     newtonia restoring coefficient for heat flux
!---------------------------------------------------------------------
!     gamma_t = grav * ddd/cp
!     ddd     = 40W/m**2/K
!     cp      = cpecific heat capacity of sea water, 3901 J/kg/K
!     1.0e-1  = unit change
!
      gamma_t = grav*40.0/3901.0*1.0e-1
!
!
!---------------------------------------------------------------------
!     newtonia restoring time scale for salinity
!---------------------------------------------------------------------
      gamma_s = r120*secday
!
!
!---------------------------------------------------------------------
!     asselin temporal filter parameter
!---------------------------------------------------------------------
!      afb1 = 0.25
!      afc1 = 0.25
!      aft1 = 0.25
!
!
!---------------------------------------------------------------------
!     time-step of integration
!---------------------------------------------------------------------
!     dtsf:   time step for adv_1                 (input unit: min)
!     dtuv:   time step for adv_2                 (input unit: hour)
!     dtts:   time step for tracer                (input unit: hour)
!
!      dtsf = 2.0
!      dtuv = 1.0
!      dtts = 24.0
!
!
!---------------------------------------------------------------------
!     parameters for model's output
!---------------------------------------------------------------------
!     io_tsuvp  = interval for output time average variables (unit year)
!     io_restr  = interval for output restr.data             (unit year)
!
!      io_tsuvp  = 1
!      io_restr  = 1
!
!
      return
      end
