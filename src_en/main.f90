!    +                                                               +
!    +===============================================================+
!    +             PCOM and BCOM in eta-coordinates                  +
!    +===============================================================+
!    +                                                               +
program main
      implicit none
      include 'param.h'
      include 'pconst.h'
      include 'mpif.h'
!
      integer  mode_c,mode_t,reczhy,test_mpi
      real rdy,area,startime,endtime
      character* (mpi_max_processor_name) processor_name
      integer myid,numprocs,namelen,rc,ierr,istatus,ncpux,ncpuy
      integer west,east,north,south,status(mpi_status_size)
!-----------------------------------------------------------------------------
!     define variable-sized array
!-----------------------------------------------------------------------------
      integer, allocatable :: itn (:,:)
      integer, allocatable :: ivn (:,:)
      integer, allocatable :: mat_myid(:,:)
      real, allocatable :: t      (:,:,:,:,:)
      real, allocatable :: up     (:,:,:,:)
      real, allocatable :: vp     (:,:,:,:)
      real, allocatable :: pbt    (:,:,:)
      real, allocatable :: pbt_st (:,:,:)
      real, allocatable :: upb    (:,:,:)
      real, allocatable :: vpb    (:,:,:)
      real, allocatable :: spbt   (:,:)
      real, allocatable :: rho    (:,:,:)
      real, allocatable :: w      (:,:,:)
      real, allocatable :: du     (:,:,:)
      real, allocatable :: dv     (:,:,:)
      real, allocatable :: adv_u  (:,:,:,:)
      real, allocatable :: adv_v  (:,:,:,:)
      real, allocatable :: dub    (:,:)
      real, allocatable :: dvb    (:,:)
      real, allocatable :: diffu  (:,:,:)
      real, allocatable :: diffv  (:,:,:)
      real, allocatable :: pt     (:,:,:)
      real, allocatable :: ps     (:,:,:)
      real, allocatable :: pmup   (:,:)
      real, allocatable :: pmum   (:,:)
      real, allocatable :: pmtp   (:,:)
      real, allocatable :: pmtm   (:,:)
      real, allocatable :: ump    (:,:,:)
      real, allocatable :: umm    (:,:,:)
      real, allocatable :: vmp    (:,:,:)
      real, allocatable :: vmm    (:,:,:)
      real, allocatable :: rhodp  (:,:,:)
      real, allocatable :: pax    (:,:)
      real, allocatable :: pay    (:,:)
      real, allocatable :: pbxn   (:,:)
      real, allocatable :: pbxs   (:,:)
      real, allocatable :: pcxn   (:,:)
      real, allocatable :: pcxs   (:,:)
      real, allocatable :: pdxn   (:,:)
      real, allocatable :: pdxs   (:,:)
      real, allocatable :: pbye   (:,:)
      real, allocatable :: pbyw   (:,:)
      real, allocatable :: pcye   (:,:)
      real, allocatable :: pcyw   (:,:)
      real, allocatable :: pdye   (:,:)
      real, allocatable :: pdyw   (:,:)
      real, allocatable :: phibx  (:,:)
      real, allocatable :: phiby  (:,:)
      real, allocatable :: fixp   (:,:,:)
      real, allocatable :: kappa_h(:,:,:)
      real, allocatable :: kappa_m(:,:,:)
      real, allocatable :: am     (:,:,:)
      real, allocatable :: acfl   (:)
      real, allocatable :: ah     (:,:,:)
      real, allocatable :: ahisop (:,:,:)
      real, allocatable :: athkdf (:,:,:)
      real, allocatable :: phib   (:,:)
      real, allocatable :: pn     (:,:)
      real, allocatable :: rzu    (:,:)
      real, allocatable :: zu     (:,:)
      real, allocatable :: z0     (:)
      real, allocatable :: dz0    (:)
      real, allocatable :: z      (:)
      real, allocatable :: dz     (:)
      real, allocatable :: rdz    (:)
      real, allocatable :: rdzw   (:)
      real, allocatable :: tmask  (:,:,:)
      real, allocatable :: umask  (:,:,:)
      real, allocatable :: lat    (:)
      real, allocatable :: lon    (:)
      real, allocatable :: cost   (:)
      real, allocatable :: cosu   (:)
      real, allocatable :: ff     (:)
      real, allocatable :: rdxt   (:)
      real, allocatable :: rdxu   (:)
      real, allocatable :: rdyt   (:)
      real, allocatable :: rdyu   (:)
      real, allocatable :: dxdyt  (:)
      real, allocatable :: dxdyu  (:)
      real, allocatable :: cv1    (:)
      real, allocatable :: cv2    (:)
      real, allocatable :: sdxt   (:)
      real, allocatable :: sdxu   (:)
      real, allocatable :: r1a    (:)
      real, allocatable :: r1b    (:)
      real, allocatable :: r1c    (:)
      real, allocatable :: r1d    (:)
      real, allocatable :: ebea   (:)
      real, allocatable :: ebeb   (:)
      real, allocatable :: ebla   (:)
      real, allocatable :: eblb   (:)
      real, allocatable :: epea   (:)
      real, allocatable :: epeb   (:)
      real, allocatable :: epla   (:)
      real, allocatable :: eplb   (:)
      real, allocatable :: bcf    (:,:,:,:)
      real, allocatable :: bcu    (:,:)
      real, allocatable :: bcv    (:,:)
      real, allocatable :: bct    (:,:)
      real, allocatable :: bcp    (:,:)
      real, allocatable :: bcs    (:,:)
      real, allocatable :: emp    (:,:)
      real, allocatable :: ddd    (:,:)
      real, allocatable :: pmn    (:,:)
      real, allocatable :: umn    (:,:,:)
      real, allocatable :: vmn    (:,:,:)
      real, allocatable :: wmn    (:,:,:)
      real, allocatable :: tmn_mo (:,:,:)
      real, allocatable :: smn_mo (:,:,:)
      real, allocatable :: pmn_mo (:,:)
      real, allocatable :: umn_mo (:,:,:)
      real, allocatable :: vmn_mo (:,:,:)
      real, allocatable :: wmn_mo (:,:,:)
      real, allocatable :: rhoi   (:,:,:,:)
      real, allocatable :: e      (:,:,:,:)
      real, allocatable :: K1     (:,:,:,:)
      real, allocatable :: K2     (:,:,:,:)
      real, allocatable :: K3     (:,:,:,:)
      real, allocatable :: adv_vetiso(:,:,:)
      real, allocatable :: adv_vntiso(:,:,:)
      real, allocatable :: adv_vbtiso(:,:,:)
      real, allocatable :: kref   (:)
      real, allocatable :: rdz0   (:)
      real, allocatable :: fzisop (:)
      real, allocatable :: gr_mass    (:,:,:)
      real, allocatable :: gr_mass2   (:,:,:)
      real, allocatable :: gr_mass_mo (:,:,:)
      real, allocatable :: wind_en_mo (:,:)
      real, allocatable :: bottom_h   (:,:)
      real, allocatable :: dpo_adv    (:,:,:)
      real, allocatable :: dpo_hdif   (:,:,:)
      real, allocatable :: dpo_vdif   (:,:,:)
      real, allocatable :: dpo_bc     (:,:)
      real, allocatable :: dpo_imp    (:,:,:)
      real, allocatable :: dpo_ice    (:,:,:)
      real, allocatable :: dpo_ast    (:,:,:)
      real, allocatable :: dpo_con    (:,:,:)
      real, allocatable :: dpo_bar    (:,:,:)
      real, allocatable :: dpo_bcp    (:,:,:)
      real, allocatable :: din_adv    (:,:,:)
      real, allocatable :: din_hdif   (:,:,:)
      real, allocatable :: din_vdif   (:,:,:)
      real, allocatable :: din_bc     (:,:)
      real, allocatable :: din_imp    (:,:,:)
      real, allocatable :: din_ice    (:,:,:)
      real, allocatable :: din_ast    (:,:,:)
      real, allocatable :: dsa_ast    (:,:,:)
      real, allocatable :: din_con    (:,:,:)
      real, allocatable :: din_bar    (:,:,:)
      real, allocatable :: dpo_adv_mo (:,:,:)
      real, allocatable :: dpo_hdif_mo(:,:,:)
      real, allocatable :: dpo_vdif_mo(:,:,:)
      real, allocatable :: dpo_bc_mo  (:,:)
      real, allocatable :: dpo_imp_mo (:,:,:)
      real, allocatable :: dpo_ice_mo (:,:,:)
      real, allocatable :: dpo_con_mo (:,:,:)
      real, allocatable :: dpo_bar_mo (:,:,:)
      real, allocatable :: din_adv_mo (:,:,:)
      real, allocatable :: din_hdif_mo(:,:,:)
      real, allocatable :: din_vdif_mo(:,:,:)
      real, allocatable :: din_bc_mo  (:,:)
      real, allocatable :: din_imp_mo (:,:,:)
      real, allocatable :: din_ice_mo (:,:,:)
      real, allocatable :: din_ast_mo (:,:,:)
      real, allocatable :: din_con_mo (:,:,:)
      real, allocatable :: din_bar_mo (:,:,:)
      real, allocatable :: total_in   (:,:,:)
      real, allocatable :: total_po   (:,:,:)
      real, allocatable :: p          (:,:,:)
      real, allocatable :: pbar       (:,:)
      real, allocatable :: mass_up    (:,:,:)
      real, allocatable :: dke_bcf    (:,:,:)
      real, allocatable :: dke_ape    (:,:,:,:)
      real, allocatable :: dke_fri    (:,:,:,:)
      real, allocatable :: dke_pre    (:,:,:,:)
      real, allocatable :: dke_adv    (:,:,:,:)
      real, allocatable :: dke_bar    (:,:,:,:)
      real, allocatable :: dke_cor    (:,:,:,:)
      real, allocatable :: dke_bcf_mo (:,:)
      real, allocatable :: dke_ape_mo (:,:,:)
      real, allocatable :: dke_adv_mo (:,:,:)
      real, allocatable :: dke_fri_mo (:,:,:)
      real, allocatable :: dke_pre_mo (:,:,:)
      real, allocatable :: dke_bar_mo (:,:,:)
      real, allocatable :: dke_cor_mo (:,:,:)
!-----------------------------------------------------------------------------
!     initialization mpi environment
!-----------------------------------------------------------------------------      
      call mpi_init(ierr)
      call mpi_comm_rank(mpi_comm_world,myid,ierr)
      call mpi_comm_size(mpi_comm_world,numprocs,ierr)
      call mpi_get_processor_name(processor_name,namelen,ierr)
      write(*,"('hello world! process',i4,'of',i4,' on ',20a)")  &
            myid,numprocs,processor_name
!-----------------------------------------------------------------------------
!     read and distribute base parameters
!-----------------------------------------------------------------------------
      if (myid == 0) then
         namelist /contrl/ imt,jmt,km,nt,runlen_mon,runlen_day,runlen_sec,restrt,  &
                           snbc,tnbc,gm90,implicitvmix,asselin_b,asselin_c,asselin_t,  &
                           smtha,smthb,smthc,smth_start_nlat,smth_start_slat,fcof,unesco, &
                           boussinesq,monloop,yearloop,monlong,daylong,hourlong
         namelist /mpicontrl/ ncpux,ncpuy
         namelist /tsteps/ stager_t,dtts,dtuv,dtsf
         namelist /mixing/ fam,fah,fkm,fkh,kh_max,am_c,sma_c,ah_c,km_c,kh_c,athkdf_c,cdbot
         namelist /filter/ afb1,afc1,aft1
         namelist /io/     uvtshdiag,energydiag,io_month,io_time,restr_mon,restr_time,infotime
   !
   !     read in model's control parameters
         open(15,file='namelist')
         read(15,contrl)
         read(15,mpicontrl)
         read(15,tsteps)
         read(15,mixing)
         read(15,filter)
         read(15,io)
         close(15)
   !-----------------------------------------------------------------------------
   !     set mpi parameters and distribute base parameters
   !-----------------------------------------------------------------------------
         allocate(mat_myid(ncpux+2,ncpuy),stat=istatus)
         simt=imt
         sjmt=jmt
         k=1
         do j=1,ncpuy
         do i=2,ncpux+1
           mat_myid(i,j)=k-1
           k=k+1
         end do
           mat_myid(1,j)=mat_myid(ncpux+1,j)
           mat_myid(ncpux+2,j)=mat_myid(2,j)
         end do
         
         do j=1,ncpuy
         do i=2,ncpux+1
         if (mat_myid(i,j).gt.0) then      
         call mpi_ssend(restrt,1,mpi_logical,mat_myid(i,j),10,mpi_comm_world,ierr)      
         call mpi_ssend(ncpux,1,mpi_integer,mat_myid(i,j),11,mpi_comm_world,ierr)
         call mpi_ssend(ncpuy,1,mpi_integer,mat_myid(i,j),12,mpi_comm_world,ierr)      
         end if
         end do
         end do
         
      end if
      
      if (myid.gt.0) then      
         call mpi_recv(restrt,1,mpi_logical,0,10,mpi_comm_world,status,ierr)
         call mpi_recv(ncpux,1,mpi_integer,0,11,mpi_comm_world,status,ierr)
         call mpi_recv(ncpuy,1,mpi_integer,0,12,mpi_comm_world,status,ierr)
         
         allocate(mat_myid(ncpux+2,ncpuy),stat=istatus)      
         k=1
         do j=1,ncpuy
         do i=2,ncpux+1
           mat_myid(i,j)=k-1
           k=k+1
         end do
           mat_myid(1,j)=mat_myid(ncpux+1,j)
           mat_myid(ncpux+2,j)=mat_myid(2,j)
         end do
      end if

      call dis_var_int(imt,mat_myid,ncpux,ncpuy,myid)      
      call dis_var_int(jmt,mat_myid,ncpux,ncpuy,myid)
      call dis_var_int(km,mat_myid,ncpux,ncpuy,myid)
      call dis_var_int(nt,mat_myid,ncpux,ncpuy,myid)
      call dis_var_int(restr_mon,mat_myid,ncpux,ncpuy,myid)
      call dis_var_real(restr_time,mat_myid,ncpux,ncpuy,myid)
      call dis_var_int(uvtshdiag,mat_myid,ncpux,ncpuy,myid)
      call dis_var_int(energydiag,mat_myid,ncpux,ncpuy,myid)
      call dis_var_int(io_month,mat_myid,ncpux,ncpuy,myid)
      call dis_var_real(io_time,mat_myid,ncpux,ncpuy,myid)
      call dis_var_real(infotime,mat_myid,ncpux,ncpuy,myid)
      call dis_var_int(snbc,mat_myid,ncpux,ncpuy,myid)
      call dis_var_int(tnbc,mat_myid,ncpux,ncpuy,myid)
      call dis_var_int(gm90,mat_myid,ncpux,ncpuy,myid)
      call dis_var_int(implicitvmix,mat_myid,ncpux,ncpuy,myid)
      call dis_var_int(asselin_b,mat_myid,ncpux,ncpuy,myid)
      call dis_var_int(asselin_c,mat_myid,ncpux,ncpuy,myid)
      call dis_var_int(asselin_t,mat_myid,ncpux,ncpuy,myid)
      call dis_var_int(smtha,mat_myid,ncpux,ncpuy,myid)
      call dis_var_int(smthb,mat_myid,ncpux,ncpuy,myid)
      call dis_var_int(smthc,mat_myid,ncpux,ncpuy,myid)
      call dis_var_int(unesco,mat_myid,ncpux,ncpuy,myid)
      call dis_var_int(boussinesq,mat_myid,ncpux,ncpuy,myid)
      call dis_var_int(monloop,mat_myid,ncpux,ncpuy,myid)
      call dis_var_int(yearloop,mat_myid,ncpux,ncpuy,myid)
      call dis_var_int(monlong,mat_myid,ncpux,ncpuy,myid)
      call dis_var_int(daylong,mat_myid,ncpux,ncpuy,myid)
      call dis_var_int(hourlong,mat_myid,ncpux,ncpuy,myid)
      call dis_var_int(stager_t,mat_myid,ncpux,ncpuy,myid)
      call dis_var_real(runlen_mon,mat_myid,ncpux,ncpuy,myid)
      call dis_var_real(runlen_day,mat_myid,ncpux,ncpuy,myid)
      call dis_var_real(runlen_sec,mat_myid,ncpux,ncpuy,myid)
      call dis_var_real(dtts,mat_myid,ncpux,ncpuy,myid)
      call dis_var_real(dtuv,mat_myid,ncpux,ncpuy,myid)
      call dis_var_real(dtsf,mat_myid,ncpux,ncpuy,myid)
      call dis_var_real(cdbot,mat_myid,ncpux,ncpuy,myid)
      call dis_var_real(afb1,mat_myid,ncpux,ncpuy,myid)
      call dis_var_real(afc1,mat_myid,ncpux,ncpuy,myid)
      call dis_var_real(aft1,mat_myid,ncpux,ncpuy,myid)
      call dis_var_real(smth_start_nlat,mat_myid,ncpux,ncpuy,myid)
      call dis_var_real(smth_start_slat,mat_myid,ncpux,ncpuy,myid)
      call dis_var_real(fcof,mat_myid,ncpux,ncpuy,myid)
      call dis_var_real(sma_c,mat_myid,ncpux,ncpuy,myid)
      simt=imt
      sjmt=jmt      
!-----compute base base parameters in all process-------------------------------
      do j=1,ncpuy
      do i=2,ncpux+1
        if (myid.eq.mat_myid(i,j)) then
          east=mat_myid(i+1,j)
          west=mat_myid(i-1,j)
          if (ncpuy.gt.1) then
          if (j.eq.1) then
            north=mat_myid(i,j+1)
            south=mpi_proc_null
          else if (j.eq.ncpuy) then
            north=mpi_proc_null
            south=mat_myid(i,j-1)
          else
            north=mat_myid(i,j+1)
            south=mat_myid(i,j-1)
          end if
          else
          north=mpi_proc_null
          south=mpi_proc_null
          end if
        end if
      end do
      end do
      
      imt=(simt-2)/ncpux+2
      jmt=(sjmt-2)/ncpuy+2
      imm=imt-1
      jmm=jmt-1
      kmp1=km+1
      kmm1=km-1
      test_mpi=1
      reczhy=1
      
      allocate(z0     (km)              ,stat=istatus)
      allocate(dz0    (km)              ,stat=istatus)
      allocate(z      (km)              ,stat=istatus)
      allocate(dz     (km)              ,stat=istatus)
      allocate(rdz    (km)              ,stat=istatus)
      allocate(rdzw   (km)              ,stat=istatus)
      dz0    =0
      z      =0
      dz     =0
      rdz    =0
      rdzw   =0
      runtime=0
!-----------------------------------------------------------------------------
!     open input and output file at myid=0 process
!-----------------------------------------------------------------------------     
      if (myid == 0) then
         open(66,file="diaginfo.txt",status='replace')
         
         write(66,*) 'Allocating arrays'      
      end if
!-----------------------------------------------------------------------------
!     initialization all mpi array in every process
!-----------------------------------------------------------------------------
      allocate(t      (imt,jmt,km,nt,2) ,stat=istatus)
      allocate(up     (imt,jmt,km,2)    ,stat=istatus)
      allocate(vp     (imt,jmt,km,2)    ,stat=istatus)
      allocate(pbt    (imt,jmt,2)       ,stat=istatus)
      allocate(pbt_st (imt,jmt,4)       ,stat=istatus)
      allocate(upb    (imt,jmt,2)       ,stat=istatus)
      allocate(vpb    (imt,jmt,2)       ,stat=istatus)
      allocate(spbt   (imt,jmt)         ,stat=istatus)
      allocate(rho    (imt,jmt,km)      ,stat=istatus)
      allocate(w      (imt,jmt,kmp1)    ,stat=istatus)
      allocate(du     (imt,jmt,km)      ,stat=istatus)
      allocate(dv     (imt,jmt,km)      ,stat=istatus)
      allocate(adv_u  (imt,jmt,km,3)    ,stat=istatus)
      allocate(adv_v  (imt,jmt,km,3)    ,stat=istatus)
      allocate(dub    (imt,jmt)         ,stat=istatus)
      allocate(dvb    (imt,jmt)         ,stat=istatus)
      allocate(diffu  (imt,jmt,km)      ,stat=istatus)
      allocate(diffv  (imt,jmt,km)      ,stat=istatus)
      allocate(pt     (imt,jmt,km)      ,stat=istatus)
      allocate(ps     (imt,jmt,km)      ,stat=istatus)
      allocate(pmup   (imt,jmt)         ,stat=istatus)
      allocate(pmum   (imt,jmt)         ,stat=istatus)
      allocate(pmtp   (imt,jmt)         ,stat=istatus)
      allocate(pmtm   (imt,jmt)         ,stat=istatus)
      allocate(ump    (imt,jmt,km)      ,stat=istatus)
      allocate(umm    (imt,jmt,km)      ,stat=istatus)
      allocate(vmp    (imt,jmt,km)      ,stat=istatus)
      allocate(vmm    (imt,jmt,km)      ,stat=istatus)
      allocate(rhodp  (imt,jmt,km)      ,stat=istatus)
      allocate(pax    (imt,jmt)         ,stat=istatus)
      allocate(pay    (imt,jmt)         ,stat=istatus)
      allocate(pbxn   (imt,jmt)         ,stat=istatus)
      allocate(pbxs   (imt,jmt)         ,stat=istatus)
      allocate(pcxn   (imt,jmt)         ,stat=istatus)
      allocate(pcxs   (imt,jmt)         ,stat=istatus)
      allocate(pdxn   (imt,jmt)         ,stat=istatus)
      allocate(pdxs   (imt,jmt)         ,stat=istatus)
      allocate(pbye   (imt,jmt)         ,stat=istatus)
      allocate(pbyw   (imt,jmt)         ,stat=istatus)
      allocate(pcye   (imt,jmt)         ,stat=istatus)
      allocate(pcyw   (imt,jmt)         ,stat=istatus)
      allocate(pdye   (imt,jmt)         ,stat=istatus)
      allocate(pdyw   (imt,jmt)         ,stat=istatus)
      allocate(phibx  (imt,jmt)         ,stat=istatus)
      allocate(phiby  (imt,jmt)         ,stat=istatus)
      allocate(fixp   (imt,jmt,km)      ,stat=istatus)
      allocate(kappa_h(imt,jmt,km)      ,stat=istatus)
      allocate(kappa_m(imt,jmt,km)      ,stat=istatus)
      allocate(ah     (imt,jmt,km)      ,stat=istatus)
      allocate(ahisop (imt,jmt,km)      ,stat=istatus)
      allocate(athkdf (imt,jmt,km)      ,stat=istatus)
      allocate(am     (imt,jmt,km)      ,stat=istatus)
      allocate(acfl   (jmt)             ,stat=istatus)
      allocate(phib   (imt,jmt)         ,stat=istatus)
      allocate(pn     (imt,jmt)         ,stat=istatus)
      allocate(rzu    (imt,jmt)         ,stat=istatus)
      allocate(zu     (imt,jmt)         ,stat=istatus)
      allocate(itn    (imt,jmt)         ,stat=istatus)
      allocate(ivn    (imt,jmt)         ,stat=istatus)
      allocate(tmask  (imt,jmt,km)      ,stat=istatus)
      allocate(umask  (imt,jmt,km)      ,stat=istatus)
      allocate(lat    (jmt)             ,stat=istatus)
      allocate(lon    (imt)             ,stat=istatus)
      allocate(cost   (jmt)             ,stat=istatus)
      allocate(cosu   (jmt)             ,stat=istatus)
      allocate(ff     (jmt)             ,stat=istatus)
      allocate(rdxt   (jmt)             ,stat=istatus)
      allocate(rdxu   (jmt)             ,stat=istatus)
      allocate(rdyt   (jmt)             ,stat=istatus)
      allocate(rdyu   (jmt)             ,stat=istatus)
      allocate(dxdyt  (jmt)             ,stat=istatus)
      allocate(dxdyu  (jmt)             ,stat=istatus)
      allocate(cv1    (jmt)             ,stat=istatus)
      allocate(cv2    (jmt)             ,stat=istatus)
      allocate(sdxt   (jmt)             ,stat=istatus)
      allocate(sdxu   (jmt)             ,stat=istatus)
      allocate(r1a    (jmt)             ,stat=istatus)
      allocate(r1b    (jmt)             ,stat=istatus)
      allocate(r1c    (jmt)             ,stat=istatus)
      allocate(r1d    (jmt)             ,stat=istatus)
      allocate(ebea   (jmt)             ,stat=istatus)
      allocate(ebeb   (jmt)             ,stat=istatus)
      allocate(ebla   (jmt)             ,stat=istatus)
      allocate(eblb   (jmt)             ,stat=istatus)
      allocate(epea   (jmt)             ,stat=istatus)
      allocate(epeb   (jmt)             ,stat=istatus)
      allocate(epla   (jmt)             ,stat=istatus)
      allocate(eplb   (jmt)             ,stat=istatus)
      allocate(bcf    (imt,jmt,12,7)    ,stat=istatus)
      allocate(bcu    (imt,jmt)         ,stat=istatus)
      allocate(bcv    (imt,jmt)         ,stat=istatus)
      allocate(bct    (imt,jmt)         ,stat=istatus)
      allocate(bcp    (imt,jmt)         ,stat=istatus)
      allocate(bcs    (imt,jmt)         ,stat=istatus)
      allocate(emp    (imt,jmt)         ,stat=istatus)
      allocate(ddd    (imt,jmt)         ,stat=istatus)
      allocate(pmn    (imt,jmt)         ,stat=istatus)
      allocate(umn    (imt,jmt,km)      ,stat=istatus)
      allocate(vmn    (imt,jmt,km)      ,stat=istatus)
      allocate(wmn    (imt,jmt,km)      ,stat=istatus)
      allocate(tmn_mo (imt,jmt,km)      ,stat=istatus)
      allocate(smn_mo (imt,jmt,km)      ,stat=istatus)
      allocate(pmn_mo (imt,jmt)         ,stat=istatus)
      allocate(umn_mo (imt,jmt,km)      ,stat=istatus)
      allocate(vmn_mo (imt,jmt,km)      ,stat=istatus)
      allocate(wmn_mo (imt,jmt,km)      ,stat=istatus)
      allocate(rhoi   (imt,km,jmt,1:3)  ,stat=istatus)
      allocate(e      (imt,kmp1,jmt,3)  ,stat=istatus)
      allocate(K1     (imt,km,jmt,3:3)  ,stat=istatus)
      allocate(K2     (imt,km,jmt,3:3)  ,stat=istatus)
      allocate(K3     (imt,km,jmt,1:3)  ,stat=istatus)
      allocate(adv_vetiso   (imt,km,jmt),stat=istatus)
      allocate(adv_vntiso   (imt,km,jmt),stat=istatus)
      allocate(adv_vbtiso (imt,0:km,jmt),stat=istatus)
      allocate(kref                 (km),stat=istatus)
      allocate(rdz0                 (km),stat=istatus)
      allocate(fzisop               (km),stat=istatus)
      allocate(gr_mass      (imt,jmt,km),stat=istatus)
      allocate(gr_mass2     (imt,jmt,km),stat=istatus)
      allocate(gr_mass_mo   (imt,jmt,km),stat=istatus)
      allocate(wind_en_mo   (imt,jmt)   ,stat=istatus)
      allocate(bottom_h     (imt,jmt)   ,stat=istatus)
      allocate(dpo_adv      (imt,jmt,km),stat=istatus)
      allocate(dpo_hdif     (imt,jmt,km),stat=istatus)
      allocate(dpo_vdif     (imt,jmt,km),stat=istatus)
      allocate(dpo_bc       (imt,jmt)   ,stat=istatus)
      allocate(dpo_imp      (imt,jmt,km),stat=istatus)
      allocate(dpo_ice      (imt,jmt,km),stat=istatus)
      allocate(dpo_ast      (imt,jmt,km),stat=istatus)
      allocate(dpo_con      (imt,jmt,km),stat=istatus)
      allocate(dpo_bar      (imt,jmt,km),stat=istatus)
      allocate(dpo_bcp      (imt,jmt,km),stat=istatus)
      allocate(din_adv      (imt,jmt,km),stat=istatus)
      allocate(din_hdif     (imt,jmt,km),stat=istatus)
      allocate(din_vdif     (imt,jmt,km),stat=istatus)
      allocate(din_bc       (imt,jmt)   ,stat=istatus)
      allocate(din_imp      (imt,jmt,km),stat=istatus)
      allocate(din_ice      (imt,jmt,km),stat=istatus)
      allocate(din_ast      (imt,jmt,km),stat=istatus)
      allocate(dsa_ast      (imt,jmt,km),stat=istatus)
      allocate(din_con      (imt,jmt,km),stat=istatus)
      allocate(din_bar      (imt,jmt,km),stat=istatus)
      allocate(dpo_adv_mo   (imt,jmt,km),stat=istatus)
      allocate(dpo_hdif_mo  (imt,jmt,km),stat=istatus)
      allocate(dpo_vdif_mo  (imt,jmt,km),stat=istatus)
      allocate(dpo_bc_mo    (imt,jmt)   ,stat=istatus)
      allocate(dpo_imp_mo   (imt,jmt,km),stat=istatus)
      allocate(dpo_ice_mo   (imt,jmt,km),stat=istatus)
      allocate(dpo_con_mo   (imt,jmt,km),stat=istatus)
      allocate(dpo_bar_mo   (imt,jmt,km),stat=istatus)
      allocate(din_adv_mo   (imt,jmt,km),stat=istatus)
      allocate(din_hdif_mo  (imt,jmt,km),stat=istatus)
      allocate(din_vdif_mo  (imt,jmt,km),stat=istatus)
      allocate(din_bc_mo    (imt,jmt)   ,stat=istatus)
      allocate(din_imp_mo   (imt,jmt,km),stat=istatus)
      allocate(din_ice_mo   (imt,jmt,km),stat=istatus)
      allocate(din_ast_mo   (imt,jmt,km),stat=istatus)
      allocate(din_con_mo   (imt,jmt,km),stat=istatus)
      allocate(din_bar_mo   (imt,jmt,km),stat=istatus)
      allocate(total_in     (imt,jmt,km),stat=istatus)
      allocate(total_po     (imt,jmt,km),stat=istatus)
      allocate(p            (imt,jmt,km),stat=istatus)
      allocate(mass_up      (imt,jmt,km),stat=istatus)
      allocate(pbar         (imt,jmt)   ,stat=istatus)
      allocate(dke_bcf      (imt,jmt,2)   ,stat=istatus)
      allocate(dke_ape      (imt,jmt,km,2),stat=istatus)
      allocate(dke_fri      (imt,jmt,km,2),stat=istatus)
      allocate(dke_pre      (imt,jmt,km,2),stat=istatus)
      allocate(dke_adv      (imt,jmt,km,2),stat=istatus)
      allocate(dke_bar      (imt,jmt,km,2),stat=istatus)
      allocate(dke_cor      (imt,jmt,km,2),stat=istatus)
      allocate(dke_bcf_mo   (imt,jmt)   ,stat=istatus)
      allocate(dke_adv_mo   (imt,jmt,km),stat=istatus)
      allocate(dke_ape_mo   (imt,jmt,km),stat=istatus)
      allocate(dke_fri_mo   (imt,jmt,km),stat=istatus)
      allocate(dke_pre_mo   (imt,jmt,km),stat=istatus)
      allocate(dke_bar_mo   (imt,jmt,km),stat=istatus)
      allocate(dke_cor_mo   (imt,jmt,km),stat=istatus)

      tmn_mo      = c0
      smn_mo      = c0
      pmn_mo      = c0
      umn_mo      = c0
      vmn_mo      = c0
      wmn_mo      = c0
      gr_mass     = c0
      gr_mass2    = c0
      gr_mass_mo  = c0
      wind_en_mo  = c0
      dpo_adv     = c0
      dpo_hdif    = c0
      dpo_vdif    = c0
      dpo_bc      = c0
      dpo_imp     = c0
      dpo_ice     = c0
      dpo_ast     = c0
      dpo_con     = c0
      dpo_bar     = c0
      din_adv     = c0
      din_hdif    = c0
      din_vdif    = c0
      din_bc      = c0
      din_imp     = c0
      din_ice     = c0
      din_ast     = c0
      dsa_ast     = c0
      din_con     = c0
      din_bar     = c0
      dpo_bcp     = c0
      dpo_adv_mo  = c0
      dpo_hdif_mo = c0
      dpo_vdif_mo = c0
      dpo_bc_mo   = c0
      dpo_imp_mo  = c0
      dpo_ice_mo  = c0
      dpo_con_mo  = c0
      dpo_bar_mo  = c0
      din_adv_mo  = c0
      din_hdif_mo = c0
      din_vdif_mo = c0
      din_bc_mo   = c0
      din_imp_mo  = c0
      din_ice_mo  = c0
      din_ast_mo  = c0
      din_con_mo  = c0
      din_bar_mo  = c0
      total_in    = c0
      total_po    = c0
      p           = c0
      dke_bcf     = c0
      dke_fri     = c0
      dke_ape     = c0
      dke_pre     = c0
      dke_adv     = c0
      dke_bar     = c0
      dke_cor     = c0
      dke_bcf_mo  = c0
      dke_adv_mo  = c0
      dke_fri_mo  = c0
      dke_ape_mo  = c0
      dke_pre_mo  = c0
      dke_bar_mo  = c0
      dke_cor_mo  = c0
      sa_first    = 1
      t_stepu     = 1
      
      if (myid == 0) then
         write(66,*) 'Done allocating arrays'
      end if
!-----------------------------------------------------------------------------
!     initialization basic integral arrays
!-----------------------------------------------------------------------------
!     set scalar quantities
      call setcon(fam,fah,fkm,fkh,kh_max,am_c,ah_c,km_c,kh_c,am,ah,kappa_m,kappa_h,   &
                  athkdf_c,athkdf,gravr,decibar,deltap,deltat,deltas,rdeltap,  &
                  rdeltat,rdeltas,gamma_t,gamma_s,imt,jmt,km,boussinesq,myid,  &
                  ncpux,ncpuy,simt,sjmt,mat_myid)
!
!
!     set resolution, t/u mask and the j-depended parameters
      call grdvar(z,z0,dz0,pn,itn,ivn,tmask,umask,phib,dz,rdz,rdzw,zu,rzu,   &
                  jstn,jedn,jeds,jsts,lat,lon,cost,cosu,ff,rdxt,rdxu,rdyt,rdyu,  &
                  sdxt,sdxu,r1a,r1b,r1c,r1d,cv1,cv2,dxdyt,dxdyu,area,rdy,    &
                  decibar,imt,jmt,km,dlam,dphi,imm,jmm,kmp1,kmm1,unesco, &
                  myid,ncpux,ncpuy,west,east,north,south,mat_myid,simt,sjmt, &
                  smth_start_nlat,smth_start_slat,boussinesq,bottom_h,acfl,dtts)
!
!     initialization
      call inirun(afb1,afc1,aft1,afb2,afc2,aft2,dtts,dtuv,dtsf,nss,ncc,nbb,  &
                  onbb,oncc,onbc,c2dtsf,c2dtuv,c2dtts,epea,epeb,epla,eplb,   &
                  ebea,ebeb,ebla,eblb,bcf,pt,ps,t,pbt,pbt_st,up,vp,upb,vpb,  &
                  spbt,w,dub,dvb,du,dv,diffu,diffv,pmup,pmtp,pmum,pmtm,ump, &
                  vmp,umm,vmm,pax,pay,pbxn,pcxn,pdxn,pbxs,pcxs,pdxs,pbye,  &
                  pcye,pdye,pbyw,pcyw,pdyw,rhodp,phibx,phiby,phib,rdxt,rdy,   &
                  ff,month,restrt,rho,fixp,itn,imt,jmt,km,nt,imm,jmm,kmp1, &
                  dz,decibar,myid,ncpux,ncpuy,west,east,north,south,mat_myid, &
                  simt,sjmt,unesco,boussinesq,monloop,yearloop,t_stepu,stager_t,  &
                  adv_u,adv_v,am)
!
      if (gm90==1) then
        call isopyi(imt,jmt,km,kmp1,slmxr,ahisop,ah,fzisop,kref,rdz0,   &
                    dz0,z0,K1,K2,K3,adv_vetiso,adv_vntiso,rhoi,e,adv_vbtiso)
      end if
!----------------------------------------------------------------------
!     if run time less than 1 month or 1 day, then reset integral cycle
!-----------------------------------------------------------------------
      if (runlen_mon==0.0) then
        daypm(1)=int(runlen_day)
        if (runlen_day==0.0) then
          daypm(1)=1
          nss=int(runlen_sec/dtts)
        end if
      end if
      
      if (myid==0) then
         write(66,801) nss,ncc,nbb
   801   format(/1x,'nss=',i7,2x,'ncc=',i3,2x,'nbb=',i5/)
      end if
      
!-----------------------------------------------------------------------
!     monthly cycle
!-----------------------------------------------------------------------
!

10    continue
      year = 1 + (month-1)/12
      mth  = month - (month-1)/12*12
!
!
!     set euler forward/backward scheme at the beginning of every month
!
      leapfrog_t = .false.
      leapfrog_c = .false.
      leapfrog_b = .false.
      euler_back = .true.
!
      if (monlong==1) then
         call dismonbcf(month,runlen_mon,bcf,imt,jmt,simt,sjmt,myid,ncpux,ncpuy,mat_myid)
      end if
!
!-----------------------------------------------------------------------
!     daily cycle
!-----------------------------------------------------------------------
!
      do 20 day=1,daypm(mth)
      startime=mpi_wtime()
!
!     compute coefficients for calculation of potential density anomaly
!      if (stager_t.eq.1) then
!      call rho_ref_st(tmask,t,z,pbt_st,pt,ps,deltat,deltas,rdeltat,rdeltas,decibar,fixp,  &
!                   imt,jmt,km,nt,imm,jmm,west,east,north,south,unesco,boussinesq)
!      else
      call rho_ref(tmask,t,z,pbt,pt,ps,deltat,deltas,rdeltat,rdeltas,decibar,fixp,  &
                   imt,jmt,km,nt,imm,jmm,west,east,north,south,unesco,boussinesq)
!      end if
!
!     interpolate the observed monthly mean data
      call interp(day,daymd,daypm,mth,bcf,bcu,bcv,bct,bcp,bcs,emp,ddd,imt,jmt,  &
                  simt,sjmt,myid,ncpux,ncpuy,mat_myid,monloop,yearloop,monlong)
!
!-----------------------------------------------------------------------
!     for thermal cycle
!-----------------------------------------------------------------------
      do 30 mode_t = 1,nss
!
!
!     calculate  baroclinic pressure and the relavant variables
      if (stager_t.eq.1) then
      call readyt_st(tmask,z,dz,rzu,t,pbt_st,bcp,rho,rhodp,   &
                    up,vp,cosu,rdxt,rdy,pax,pay,itn,ivn,pbxn,pbxs,pbye,  &
                    pbyw,pcxn,pcxs,pcye,pcyw,pdxn,pdxs,pdye,pdyw,imt,jmt,km,nt,  &
                    imm,jmm,kmp1,decibar,west,east,north,south,unesco,boussinesq, &
                    fixp,energydiag)
      else
      call readyt(tmask,z,dz,rzu,t,pbt,spbt,pmtm,pmtp,bcp,rho,rhodp,umm,vmm,   &
                  ump,vmp,up,vp,cosu,rdxt,rdy,pax,pay,itn,ivn,pbxn,pbxs,pbye,  &
                  pbyw,pcxn,pcxs,pcye,pcyw,pdxn,pdxs,pdye,pdyw,imt,jmt,km,nt,  &
                  imm,jmm,kmp1,decibar,west,east,north,south,unesco,boussinesq, &
                  fixp,energydiag)
      end if
!
!-----------------------------------------------------------------------
!     for baroclinic & barotropic cycle
!-----------------------------------------------------------------------
      do 40 mode_c = 1,ncc
!     calculate momentum advection, diffusion & their vertical integrals
      if (stager_t.eq.1) then
      call readyc_st(umask,tmask,ivn,pbt_st,du,dv,adv_u,adv_v,dub,  &
                    dvb,up,vp,cosu,rdxt,rdxu,rdyu,rdyt,sdxu,r1c,r1d,cv1,cv2,dz,  &
                    rdz,rdzw,rzu,pn,w,pax,pay,diffu,diffv,am,kappa_m,gravr,cdbot,  &
                    bcu,bcv,imt,jmt,km,imm,jmm,kmp1,west,east,north,  &
                    south,snbc,emp,t_stepu,energydiag,dke_bcf,dke_fri)
      else
      call readyc(umask,tmask,ivn,pmum,pmup,pbt,spbt,du,dv,dub,dvb,up,vp,cosu,  &
                  rdxt,rdxu,rdyu,rdyt,sdxu,r1c,r1d,cv1,cv2,dz,rdz,rdzw,rzu,pn,  &
                  w,pax,pay,diffu,diffv,am,kappa_m,gravr,cdbot,leapfrog_c,bcu,  &
                  bcv,imt,jmt,km,imm,jmm,kmp1,west,east,north,south,snbc,emp, &
                  energydiag)
      end if
!
      if (stager_t.eq.1) then
      call barotr_st(ivn,itn,upb,vpb,r1c,r1d,sdxu,am,dub,dvb,pbt_st,phibx,    &
                    phiby,pbxn,pbxs,pbye,pbyw,pcxn,pcxs,pcye,pcyw,ff,ebla,     &
                    eblb,ebea,ebeb,pn,zu,cosu,rdxt,rdyt, &
                    dtsf,nbb,imt,jmt,km,imm,jmm,    &
                    myid,west,east,north,south,snbc,emp,jstn,jedn,jsts,   &
                    jeds,smtha,fcof,umask,boussinesq,phib,pdxn,pdxs,pdye,pdyw,  &
                    lat,lon,energydiag)
      else
      call barotr(ivn,itn,upb,vpb,r1c,r1d,sdxu,am,dub,dvb,spbt,pbt,phibx,    &
                  phiby,pbxn,pbxs,pbye,pbyw,pcxn,pcxs,pcye,pcyw,ff,ebla,     &
                  eblb,ebea,ebeb,pn,zu,cosu,rdxt,rdyt,leapfrog_b,euler_back, &
                  dtsf,c2dtsf,afb1,afb2,pmup,pmtp,nbb,imt,jmt,km,imm,jmm,    &
                  myid,west,east,north,south,asselin_b,snbc,emp,jstn,jedn,jsts,   &
                  jeds,smtha,fcof,umask,boussinesq,phib,pdxn,pdxs,pdye,pdyw,  &
                  lat,lon,energydiag)
      end if
!
!     prediction of baroclinic mode
      if (stager_t.eq.1) then
      call bclinc_st(phib,spbt,pbt,pbt_st,rhodp,rho,rdxt,rdy,ump,vmp,upb,  &
                     vpb,up,vp,du,dv,adv_u,adv_v,diffu,diffv,epla,eplb,epea,epeb,  &
                     pax,pay,ff,umask,ivn,itn,z,dz,rzu,onbb,dtuv,  &
                     cosu,imt,jmt,km,imm,jmm,west,east,north,south,boussinesq,energydiag,  &
                     dke_pre,dke_adv,dke_ape,dke_bar,dke_bcf,dke_fri,dke_cor,myid)
      else
      call bclinc(pmup,pmum,phib,spbt,pbt,rhodp,rho,rdxt,rdy,ump,vmp,upb,    &
                  vpb,up,vp,du,dv,epla,eplb,epea,epeb,ff,umask,ivn,z,dz,rzu, &
                  onbb,leapfrog_c,dtuv,c2dtuv,afc1,afc2,cosu,imt,jmt,km,imm, &
                  jmm,west,east,north,south,asselin_c,boussinesq,itn,energydiag)
      end if
!
      if (smthb==1) then
        call swap_array_real3d(upb,imt,jmt,2,west,east,north,south)
        call swap_array_real3d(upb,imt,jmt,2,west,east,north,south)
        call smth2(upb,umask,jstn,jedn,jsts,jeds,imt,jmt,km,fcof)
        call smth2(vpb,umask,jstn,jedn,jsts,jeds,imt,jmt,km,fcof)
        call swap_array_real3d(upb,imt,jmt,2,west,east,north,south)
        call swap_array_real3d(upb,imt,jmt,2,west,east,north,south)
      end if

      t_stepu=t_stepu+1
40    continue
!
!
!     prediction of temperature and salinity
      if (stager_t.eq.1) then
      call tracer_st(imt,jmt,km,nt,imm,jmm,kmp1,dtts,t,w,ump,vmp,sdxt,rdxt,  &
                    rdxu,rdyt,rdy,rdz,rdzw,r1a,r1b,lat,cosu,tmask,itn,dz,dz0,z0,z,  &
                    pn,bct,bcs,ddd,emp,bcp,gamma_t,gamma_s,du,dv,ah,am,sma_c,acfl,kappa_h,  &
                    kappa_m,gravr,west,east,north,south,myid,gm90,implicitvmix,snbc,  &
                    tnbc,K1,K2,K3,adv_vetiso,adv_vntiso,adv_vbtiso,rhoi,e,slmxr,ahisop,  &
                    athkdf,fzisop,kref,rdz0,unesco,umn,vmn,wmn,pmn,pbt_st,rho,phib, &
                    gr_mass,energydiag,dpo_adv,dpo_hdif,dpo_vdif,dpo_bc,dpo_imp,dpo_ice,  &
                    dpo_bar,dpo_bcp,din_adv,din_hdif,din_vdif,din_bc,din_imp,din_ice,  &
                    din_bar,p,pbar,mass_up,total_in,gr_mass2)
      else
      call tracer(t,w,pmtp,pmtm,umm,vmm,ump,vmp,onbc,oncc,sdxt,rdxt,rdyt,    &
                  rdz,rdzw,r1a,r1b,tmask,itn,dz,pn,bct,bcs,ddd,gamma_t,gamma_s,  &
                  du,dv,ah,am,sma_c,acfl,kappa_h,gravr,leapfrog_t,c2dtts,dtts,  &
                  imt,jmt,km,nt,imm,jmm,kmp1,west,east,north,south,myid,asselin_t,  &
                  gm90,implicitvmix,snbc,tnbc,aft1,aft2,emp,dz0,z0,rdy,cosu,K1,   &
                  K2,K3,adv_vetiso,adv_vntiso,adv_vbtiso,rhoi,e,slmxr,ahisop,    &
                  athkdf,fzisop,kref,rdz0,unesco,umn,vmn,wmn,pmn,rho,phib,gr_mass, &
                  bottom_h,ncpux,ncpuy,mat_myid,z,nss,rdxu,lat,bcp,  &
                  energydiag,dpo_adv,dpo_hdif,dpo_vdif,dpo_bc,din_adv,din_hdif,  &
                  din_vdif,din_bc,dpo_imp,din_imp,dpo_ice,din_ice,dpo_ast,din_ast, &
                  dsa_ast,p,pbar,mass_up)
      end if
!     compute sea ice
!     call seaice
!
!
!     convective adjustment if unstable stratification ocurs
      call convect(t,pt,ps,itn,dz,imt,jmt,km,nt,imm,jmm,myid,west,east,north,south,  &
                   energydiag,dpo_con,din_con,p,pbar,mass_up,bottom_h,total_in,total_po)
!
!
      call diag(imt,jmt,simt,sjmt,km,nt,imm,jmm,kmp1,nss,dtts,mode_t,year,month,  &
                mth,day,daypm,dxdyt,dxdyu,tmask,umask,itn,area,runtime,io_time,io_month, &
                restr_time,restr_mon,infotime,myid,ncpux,ncpuy,west,east,north,south,  &
                mat_myid,stager_t,t_stepu,pbt,pbt_st,up,vp,upb,vpb,pmtp,ump,vmp,adv_u,  &
                adv_v,uvtshdiag,pmn,umn,vmn,wmn,t,tmn_mo,smn_mo,pmn_mo,umn_mo,vmn_mo,  &
                wmn_mo,am,energydiag,bcu,bcv,sa_ini,ma_ini,ke_ini,in_ini,po_ini,gr_mass, &
                dpo_adv   ,dpo_hdif   ,dpo_vdif   ,dpo_imp   ,dpo_con   ,dpo_bc    ,  &
                dpo_adv_mo,dpo_hdif_mo,dpo_vdif_mo,dpo_imp_mo,dpo_con_mo,dpo_bc_mo ,  &
                din_adv   ,din_hdif   ,din_vdif   ,din_imp   ,din_con   ,din_bc    ,  &
                din_adv_mo,din_hdif_mo,din_vdif_mo,din_imp_mo,din_con_mo,din_bc_mo ,  &
                dpo_ice   ,dpo_ast    ,dpo_bar    ,dke_bcf   ,dke_fri   ,dke_adv   ,  &
                dpo_ice_mo,dsa_ast    ,dpo_bar_mo ,dke_bcf_mo,dke_fri_mo,dke_adv_mo,  &
                din_ice   ,din_ast    ,din_bar    ,dke_ape   ,dke_pre   ,dke_bar   ,  &
                din_ice_mo,din_ast_mo ,din_bar_mo ,dke_ape_mo,dke_pre_mo,dke_bar_mo,  &
                dke_cor,dke_cor_mo,dpo_bcp,  &
                total_in,total_po,wind_en_mo,gr_mass2,gr_mass_mo,sa_first,dtuv)

30    continue
!
      if (smthc==1) then
        call swap_array_real4d(up,imt,jmt,km,2,west,east,north,south)
        call swap_array_real4d(up,imt,jmt,km,2,west,east,north,south)
        call smth3(up,umask,jstn,jedn,jsts,jeds,imt,jmt,km,fcof)
        call smth3(vp,umask,jstn,jedn,jsts,jeds,imt,jmt,km,fcof)
        call swap_array_real4d(up,imt,jmt,km,2,west,east,north,south)
        call swap_array_real4d(up,imt,jmt,km,2,west,east,north,south)
      end if
      
      endtime=mpi_wtime()
      write(*,"('month ',i5,' day ',i3,' cost ',f9.5,' seconds on proc ',i2)")   &
            month,day,endtime-startime,myid
!
20    continue
!
      month  = month  + 1
      runlen_mon = runlen_mon - 1
!
      if(runlen_mon.gt.0) goto 10
      
      if (myid==0) then
      close(66)
      
      end if
      call mpi_finalize(rc)
end program main
