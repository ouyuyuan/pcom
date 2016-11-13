!    +                                                               +
!    +===============================================================+
!    +             PCOM and BCOM in eta-coordinates                  +
!    +===============================================================+
!    +                                                               +
      implicit none
      include 'param.h'
      include 'pconst.h'
      include 'mpif.h'
!      include 'scalar.h'
!      include 'calendar.h'
!      include 'prog.h'
!      include 'grdvar.h'
!
      integer  mode_c,mode_t,reczhy,test_mpi
      real rdy,area,seadp(10000),startime,endtime
      character* (mpi_max_processor_name) processor_name
      integer myid,numprocs,namelen,rc,ierr,istatus,ncpux,ncpuy
      integer west,east,north,south,status(mpi_status_size)
      real, allocatable :: t      (:,:,:,:,:)
      real, allocatable :: up     (:,:,:,:)
      real, allocatable :: vp     (:,:,:,:)
      real, allocatable :: pbt    (:,:,:)
      real, allocatable :: upb    (:,:,:)
      real, allocatable :: vpb    (:,:,:)
      real, allocatable :: spbt   (:,:)
      real, allocatable :: rho    (:,:,:)
      real, allocatable :: w      (:,:,:)
      real, allocatable :: du     (:,:,:)
      real, allocatable :: dv     (:,:,:)
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
      real, allocatable :: am(:,:,:)
      real, allocatable :: ah(:,:,:)
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
      integer, allocatable :: itn    (:,:)
      integer, allocatable :: ivn    (:,:)
      real, allocatable :: tmask  (:,:,:)
      real, allocatable :: umask  (:,:,:)
      real, allocatable :: lat    (:)
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
      real, allocatable :: tmn    (:,:,:)
      real, allocatable :: smn    (:,:,:)
      real, allocatable :: pmn    (:,:,:)
      real, allocatable :: umn    (:,:,:)
      real, allocatable :: vmn    (:,:,:)
      real, allocatable :: wmn    (:,:,:)
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
      integer, allocatable :: mat_myid(:,:)
      
      
      call mpi_init(ierr)
      call mpi_comm_rank(mpi_comm_world,myid,ierr)
      call mpi_comm_size(mpi_comm_world,numprocs,ierr)
      call mpi_get_processor_name(processor_name,namelen,ierr)
      write(*,"('hello world! process',i2,'of',i1,' on ',20a)")  &
            myid,numprocs,processor_name

!-----read and distribute base parameters------------------------------- 
      if (myid == 0) then
      namelist /contrl/ imt,jmt,km,nt,dlam,dphi,phis,runlen,restrt,snbc,gm90,  &
                        implicitvmix,asselin_b,asselin_c,asselin_t,smth,       &
                        smth_start_nlat,smth_start_slat,unesco,boussinesq,seadp
      namelist /mpicontrl/ ncpux,ncpuy
      namelist /tsteps/ dtts,dtuv,dtsf
      namelist /mixing/ fam,fah,fkm,fkh,am_c,ah_c,km_c,kh_c,cdbot
      namelist /filter/ afb1,afc1,aft1
      namelist /io/     io_tsuvp,io_restr
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
      call dis_var_int(io_tsuvp,mat_myid,ncpux,ncpuy,myid)
      call dis_var_int(io_restr,mat_myid,ncpux,ncpuy,myid)
      call dis_var_int(snbc,mat_myid,ncpux,ncpuy,myid)
      call dis_var_int(gm90,mat_myid,ncpux,ncpuy,myid)
      call dis_var_int(implicitvmix,mat_myid,ncpux,ncpuy,myid)
      call dis_var_int(asselin_b,mat_myid,ncpux,ncpuy,myid)
      call dis_var_int(asselin_c,mat_myid,ncpux,ncpuy,myid)
      call dis_var_int(asselin_t,mat_myid,ncpux,ncpuy,myid)
      call dis_var_int(smth,mat_myid,ncpux,ncpuy,myid)
      call dis_var_int(unesco,mat_myid,ncpux,ncpuy,myid)
      call dis_var_int(boussinesq,mat_myid,ncpux,ncpuy,myid)
      call dis_var_real(dlam,mat_myid,ncpux,ncpuy,myid)
      call dis_var_real(dphi,mat_myid,ncpux,ncpuy,myid)
      call dis_var_real(phis,mat_myid,ncpux,ncpuy,myid)
      call dis_var_real(runlen,mat_myid,ncpux,ncpuy,myid)
      call dis_var_real(dtts,mat_myid,ncpux,ncpuy,myid)
      call dis_var_real(dtuv,mat_myid,ncpux,ncpuy,myid)
      call dis_var_real(dtsf,mat_myid,ncpux,ncpuy,myid)
      call dis_var_real(cdbot,mat_myid,ncpux,ncpuy,myid)
      call dis_var_real(afb1,mat_myid,ncpux,ncpuy,myid)
      call dis_var_real(afc1,mat_myid,ncpux,ncpuy,myid)
      call dis_var_real(aft1,mat_myid,ncpux,ncpuy,myid)
      call dis_var_real(smth_start_nlat,mat_myid,ncpux,ncpuy,myid)
      call dis_var_real(smth_start_slat,mat_myid,ncpux,ncpuy,myid)
      call dis_var_real1d(seadp,10000,mat_myid,ncpux,ncpuy,myid)
      simt=imt
      sjmt=jmt
!^^^^^read and distribute base parameters^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^      
      
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
      z0     =seadp(1:km)
      dz0    =0
      z      =0
      dz     =0
      rdz    =0
      rdzw   =0
      
!^^^^^compute base base parameters in all process^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!-----open info and result file in 0 process---------------------------------      
      if (myid == 0) then
      open(80,file='sbcf.data',form='unformatted',access='direct', &
           recl=simt*sjmt*7*8,status='old')
      bcfrec=1
      open(66,file="diaginfo.txt",status='replace')
      open(17,file="sfcuvts.grd",form='unformatted',access='direct', &
           recl=(simt-2)*sjmt*(km*5+1)*4,status='replace')
      write(66,*) 'Allocating arrays'      
      end if
!^^^^^open info and result file in 0 process^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!-----initialization all mpi array in every process--------------------------
      allocate(t      (imt,jmt,km,nt,2) ,stat=istatus)
      allocate(up     (imt,jmt,km,2)    ,stat=istatus)
      allocate(vp     (imt,jmt,km,2)    ,stat=istatus)
      allocate(pbt    (imt,jmt,2)       ,stat=istatus)
      allocate(upb    (imt,jmt,2)       ,stat=istatus)
      allocate(vpb    (imt,jmt,2)       ,stat=istatus)
      allocate(spbt   (imt,jmt)         ,stat=istatus)
      allocate(rho    (imt,jmt,km)      ,stat=istatus)
      allocate(w      (imt,jmt,kmp1)    ,stat=istatus)
      allocate(du     (imt,jmt,km)      ,stat=istatus)
      allocate(dv     (imt,jmt,km)      ,stat=istatus)
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
      allocate(am     (imt,jmt,km)      ,stat=istatus)
      allocate(phib   (imt,jmt)         ,stat=istatus)
      allocate(pn     (imt,jmt)         ,stat=istatus)
      allocate(rzu    (imt,jmt)         ,stat=istatus)
      allocate(zu     (imt,jmt)         ,stat=istatus)
      allocate(itn    (imt,jmt)         ,stat=istatus)
      allocate(ivn    (imt,jmt)         ,stat=istatus)
      allocate(tmask  (imt,jmt,km)      ,stat=istatus)
      allocate(umask  (imt,jmt,km)      ,stat=istatus)
      allocate(lat    (jmt)             ,stat=istatus)
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
      allocate(tmn    (imt,jmt,km)      ,stat=istatus)
      allocate(smn    (imt,jmt,km)      ,stat=istatus)
      allocate(pmn    (imt,jmt,km)      ,stat=istatus)
      allocate(umn    (imt,jmt,km)      ,stat=istatus)
      allocate(vmn    (imt,jmt,km)      ,stat=istatus)
      allocate(wmn    (imt,jmt,km)      ,stat=istatus)
      allocate(rhoi   (imt,km,jmt,1:3)  ,stat=istatus)
      allocate(e      (imt,kmp1,jmt,3)  ,stat=istatus)
      allocate(K1     (imt,km,jmt,3:3)  ,stat=istatus)
      allocate(K2     (imt,km,jmt,3:3)  ,stat=istatus)
      allocate(K3     (imt,km,jmt,1:3)  ,stat=istatus)
      allocate(adv_vetiso(imt,km,jmt)   ,stat=istatus)
      allocate(adv_vntiso(imt,km,jmt)   ,stat=istatus)
      allocate(adv_vbtiso(imt,0:km,jmt) ,stat=istatus)
      allocate(kref   (km)              ,stat=istatus)
      allocate(rdz0   (km)              ,stat=istatus)
      allocate(fzisop (km)              ,stat=istatus)
      
      t      =0
      up     =0
      vp     =0
      pbt    =0
      upb    =0
      vpb    =0
      spbt   =0
      rho    =0
      w      =0
      du     =0
      dv     =0
      dub    =0
      dvb    =0
      diffu  =0
      diffv  =0
      pt     =0
      ps     =0
      pmup   =0
      pmum   =0
      pmtp   =0
      pmtm   =0
      ump    =0
      umm    =0
      vmp    =0
      vmm    =0
      rhodp  =0
      pax    =0
      pay    =0
      pbxn   =0
      pbxs   =0
      pcxn   =0
      pcxs   =0
      pdxn   =0
      pdxs   =0
      pbye   =0
      pbyw   =0
      pcye   =0
      pcyw   =0
      pdye   =0
      pdyw   =0
      phibx  =0
      phiby  =0
      kappa_h=0
      kappa_m=0
      am     =0
      ah     =0
      phib   =0
      pn     =0
      rzu    =0
      zu     =0
      itn    =0
      ivn    =0
      tmask  =0
      umask  =0
      lat    =0
      cost   =0
      cosu   =0
      ff     =0
      rdxt   =0
      rdxu   =0
      rdyt   =0
      rdyu   =0
      dxdyt  =0
      dxdyu  =0
      cv1    =0
      cv2    =0
      sdxt   =0
      sdxu   =0
      r1a    =0
      r1b    =0
      r1c    =0
      r1d    =0
      ebea   =0
      ebeb   =0
      ebla   =0
      eblb   =0
      epea   =0
      epeb   =0
      epla   =0
      eplb   =0
      bcf    =0
      bcu    =0
      bcv    =0
      bct    =0
      bcp    =0
      bcs    =0
      emp    =0
      ddd    =0
!^^^^^initialization all mpi array in every process^^^^^^^^^^^^^^^^^^^^^^^
      if (myid == 0) then
      write(66,*) 'Done allocating arrays'
      end if
!     set scalar quantities
      call setcon(fam,fah,fkm,fkh,am_c,ah_c,km_c,kh_c,am,ah,kappa_m,kappa_h,    &
                  gravr,decibar,deltap,deltat,deltas,rdeltap,rdeltat,rdeltas,   &
                  gamma_t,gamma_s,imt,jmt,km,boussinesq,myid,ncpux,ncpuy,simt,  &
                  sjmt,mat_myid)
!
!
!     set resolution, t/u mask and the j-depended parameters
      call grdvar(z,z0,dz0,pn,itn,ivn,tmask,umask,phib,dz,rdz,rdzw,zu,rzu,   &
                  jstn,jedn,jeds,jsts,lat,cost,cosu,ff,rdxt,rdxu,rdyt,rdyu,  &
                  sdxt,sdxu,r1a,r1b,r1c,r1d,cv1,cv2,dxdyt,dxdyu,area,rdy,    &
                  decibar,imt,jmt,km,dlam,dphi,phis,imm,jmm,kmp1,kmm1,unesco, &
                  myid,ncpux,ncpuy,west,east,north,south,mat_myid,simt,sjmt, &
                  smth_start_nlat,smth_start_slat,boussinesq)
!
!     initialization
      call inirun(afb1,afc1,aft1,afb2,afc2,aft2,dtts,dtuv,dtsf,nss,ncc,nbb,  &
                  onbb,oncc,onbc,c2dtsf,c2dtuv,c2dtts,epea,epeb,epla,eplb,   &
                  ebea,ebeb,ebla,eblb,bcf,pt,ps,t,pbt,up,vp,upb,vpb,spbt,w,  &
                  dub,dvb,du,dv,diffu,diffv,pmup,pmtp,pmum,pmtm,ump,vmp,umm, &
                  vmm,pax,pay,pbxn,pcxn,pdxn,pbxs,pcxs,pdxs,pbye,pcye,pdye,  &
                  pbyw,pcyw,pdyw,rhodp,phibx,phiby,phib,rdxt,rdy,ff,month,   &
                  restrt,rho,itn,imt,jmt,km,nt,phis,imm,jmm,kmp1,dz,decibar, &
                  myid,ncpux,ncpuy,west,east,north,south,mat_myid,simt,sjmt, &
                  unesco,boussinesq,fixp)
!
      if (gm90==1) then
      call isopyi(imt,jmt,km,kmp1,slmxr,ahisop,athkdf,fzisop,kref,rdz0,   &
                  dz0,z0,K1,K2,K3,adv_vetiso,adv_vntiso,rhoi,e,adv_vbtiso)
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
      tmn    =c0
      smn    =c0
      pmn    =c0
      umn    =c0
      vmn    =c0
      wmn    =c0
!
!-----------------------------------------------------------------------
!     daily cycle
!-----------------------------------------------------------------------
!
      do 20 day=1,daypm(mth)
      startime=mpi_wtime()
!
!     compute coefficients for calculation of potential density anomaly
      call rho_ref(tmask,t,z,pbt,pt,ps,deltat,deltas,rdeltat,rdeltas,decibar,fixp,  &
                   imt,jmt,km,nt,imm,jmm,west,east,north,south,unesco,boussinesq)
!
!
!     interpolate the observed monthly mean data
      call interp(day,daymd,daypm,mth,bcf,bcu,bcv,bct,bcp,bcs,emp,ddd,imt,jmt,  &
                  bcfrec,simt,sjmt,myid,ncpux,ncpuy,mat_myid)

!
!
!-----------------------------------------------------------------------
!     for thermal cycle
!-----------------------------------------------------------------------
      do 30 mode_t = 1,nss
!
!
!     calculate  baroclinic pressure and the relavant variables
      call readyt(tmask,z,dz,rzu,t,pbt,spbt,pmtm,pmtp,bcp,rho,rhodp,umm,vmm,   &
                  ump,vmp,up,vp,cosu,rdxt,rdy,pax,pay,itn,ivn,pbxn,pbxs,pbye,  &
                  pbyw,pcxn,pcxs,pcye,pcyw,pdxn,pdxs,pdye,pdyw,imt,jmt,km,nt,  &
                  imm,jmm,kmp1,decibar,west,east,north,south,unesco,boussinesq, &
                  fixp)
!
!-----------------------------------------------------------------------
!     for baroclinic & barotropic cycle
!-----------------------------------------------------------------------
      do 40 mode_c = 1,ncc
!
!     calculate momentum advection, diffusion & their vertical integrals
      call readyc(umask,tmask,ivn,pmum,pmup,pbt,spbt,du,dv,dub,dvb,up,vp,    &
                  cosu,rdxt,rdxu,rdyu,rdyt,sdxu,r1c,r1d,cv1,cv2,dz,rdz,      &
                  rdzw,rzu,pn,w,pax,pay,diffu,diffv,am,kappa_m,gravr,        &
                  leapfrog_c,cdbot,bcu,bcv,imt,jmt,km,imm,jmm,kmp1,          &
                  west,east,north,south,snbc,emp)
!
!     prediction of barotropic mode
      call barotr(ivn,itn,upb,vpb,r1c,r1d,sdxu,am,dub,dvb,spbt,pbt,phibx,    &
                  phiby,pbxn,pbxs,pbye,pbyw,pcxn,pcxs,pcye,pcyw,ff,ebla,     &
                  eblb,ebea,ebeb,pn,zu,cosu,rdxt,rdyt,leapfrog_b,euler_back, &
                  dtsf,c2dtsf,afb1,afb2,pmup,pmtp,nbb,imt,jmt,km,imm,jmm,    &
                  west,east,north,south,asselin_b,snbc,emp,jstn,jedn,jsts,   &
                  jeds,smth,umask,boussinesq,phib,pdxn,pdxs,pdye,pdyw)
!
!     prediction of baroclinic mode
      call bclinc(pmup,pmum,phib,spbt,pbt,rhodp,rho,rdxt,rdy,ump,vmp,upb,    &
                  vpb,up,vp,du,dv,epla,eplb,epea,epeb,ff,umask,ivn,z,dz,rzu, &
                  onbb,leapfrog_c,dtuv,c2dtuv,afc1,afc2,cosu,imt,jmt,km,imm, &
                  jmm,west,east,north,south,asselin_c,boussinesq,itn)
!
      if (smth==1) then
      call swap_array_real3d(upb,imt,jmt,2,west,east,north,south)
      call swap_array_real3d(upb,imt,jmt,2,west,east,north,south)
      call smth2(upb,umask,jstn,jedn,jsts,jeds,imt,jmt,km)
      call smth2(vpb,umask,jstn,jedn,jsts,jeds,imt,jmt,km)
      call swap_array_real3d(upb,imt,jmt,2,west,east,north,south)
      call swap_array_real3d(upb,imt,jmt,2,west,east,north,south)
      end if

40    continue
!
!
!     prediction of temperature and salinity
      call tracer(t,w,pmtp,pmtm,umm,vmm,ump,vmp,onbc,oncc,sdxt,rdxt,rdyt,    &
                  rdz,rdzw,r1a,r1b,tmask,itn,dz,pn,bct,bcs,gamma_t,gamma_s,  &
                  du,dv,ah,kappa_h,gravr,leapfrog_t,c2dtts,dtts,imt,jmt,km,  &
                  nt,imm,jmm,kmp1,west,east,north,south,myid,asselin_t,gm90, &
                  implicitvmix,snbc,aft1,aft2,emp,dz0,z0,rdy,cosu,K1,K2,K3,  &
                  adv_vetiso,adv_vntiso,adv_vbtiso,rhoi,e,slmxr,ahisop,    &
                  athkdf,fzisop,kref,rdz0,unesco,umn,vmn,wmn,pmn,rho,phib)
!
!     compute sea ice
!     call seaice
!
!
!     convective adjustment if unstable stratification ocurs
      call convect(t,pt,ps,itn,dz,imt,jmt,km,nt,imm,jmm,west,east,north,south)
!
      do k=1,km
      do j=1,jmt
      do i=1,imt
      tmn(i,j,k)=tmn(i,j,k)+t(i,j,k,1,tau)
      smn(i,j,k)=smn(i,j,k)+t(i,j,k,2,tau)
      end do
      end do
      end do
!
30    continue
!
      if (smth==1) then
      call swap_array_real4d(up,imt,jmt,km,2,west,east,north,south)
      call swap_array_real4d(up,imt,jmt,km,2,west,east,north,south)
      call smth3(up,umask,jstn,jedn,jsts,jeds,imt,jmt,km)
      call smth3(vp,umask,jstn,jedn,jsts,jeds,imt,jmt,km)
      call swap_array_real4d(up,imt,jmt,km,2,west,east,north,south)
      call swap_array_real4d(up,imt,jmt,km,2,west,east,north,south)
      end if
      
      call diag(tmn,smn,pmn,umn,vmn,wmn,t,pbt,spbt,up,vp,upb,vpb,pmtp,ump,vmp,   &
                phib,rho,year,month,mth,day,daypm,io_tsuvp,io_restr,pn,dz,rdxt,  &
                dxdyt,dxdyu,tmask,umask,itn,area,imt,jmt,km,nt,imm,jmm,kmp1,nss,  &
                reczhy,myid,ncpux,ncpuy,west,east,north,south,mat_myid,simt,sjmt)
      endtime=mpi_wtime()
      write(*,"('month ',i5,' day ',i3,' cost ',f9.5,' seconds on proc ',i2)")   &
            month,day,endtime-startime,myid
!
20    continue
!
      month  = month  + 1
      runlen = runlen - 1
!
      if(runlen.gt.0) goto 10
      
      if (myid==0) then
      close(66)
      close(17)
      close(80)
      end if
      call mpi_finalize(rc)
      end
