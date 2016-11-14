!
!     ===============
      subroutine diag(imt,jmt,simt,sjmt,km,nt,imm,jmm,kmp1,nss,dtts,mode_t,year,month,  &
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
                
!     ===============
!
      implicit none
      include 'pconst.h'
      include 'mpif.h'
! basic var & array
      integer imt,jmt,km,nt,imm,jmm,kmp1,i,j,k,i2,j2,nss,mode_t,stager_t,t_stepu
      character fname*15,ftail*10,normalname*15,energyname*15,time_units*25
      integer   daypm(12)
      integer   month,year,mth,day,restr_mon,io_month
      integer uvtshdiag,energydiag,diag_out,restr_out
      real runtime,io_time,restr_time,infotime,factor,dtts,dtuv,time_id,fillvaule
      integer itn(imt,jmt)
      real tmask(imt,jmt,km),umask(imt,jmt,km),stmask(simt,sjmt,km),sumask(simt,sjmt,km)
      real dxdyt(jmt),dxdyu(jmt),sdxdyt(sjmt),spsi(simt,sjmt),area
      real lat(sjmt),lon(simt-2),z(km)
      integer myid,ncpux,ncpuy,west,east,north,south,simt,sjmt
      integer mat_myid(ncpux+2,ncpuy)
! normal output array
      real umn(imt,jmt,km),umn_mo(imt,jmt,km),sumn(simt,sjmt,km),uzhy(simt-2,sjmt,km)
      real vmn(imt,jmt,km),vmn_mo(imt,jmt,km),svmn(simt,sjmt,km),vzhy(simt-2,sjmt,km)
      real wmn(imt,jmt,km),wmn_mo(imt,jmt,km),swmn(simt,sjmt,km),wzhy(simt-2,sjmt,km)
      real                 tmn_mo(imt,jmt,km),stmn(simt,sjmt,km),tzhy(simt-2,sjmt,km)
      real                 smn_mo(imt,jmt,km),ssmn(simt,sjmt,km),szhy(simt-2,sjmt,km)
      real pmn(imt,jmt)   ,pmn_mo(imt,jmt)   ,spmn(simt,sjmt)   ,pzhy(simt-2,sjmt)
      real am (imt,jmt,km),sam(simt,sjmt,km),amzhy(simt-2,sjmt,km)
! energy output array
      real dpo_adv (imt,jmt,km),dpo_adv_mo (imt,jmt,km),sdpo_adv (simt,sjmt,km),dpo_adv_zhy (simt-2,sjmt,km)
      real din_adv (imt,jmt,km),din_adv_mo (imt,jmt,km),sdin_adv (simt,sjmt,km),din_adv_zhy (simt-2,sjmt,km)
      real dpo_hdif(imt,jmt,km),dpo_hdif_mo(imt,jmt,km),sdpo_hdif(simt,sjmt,km),dpo_hdif_zhy(simt-2,sjmt,km)
      real din_hdif(imt,jmt,km),din_hdif_mo(imt,jmt,km),sdin_hdif(simt,sjmt,km),din_hdif_zhy(simt-2,sjmt,km)
      real dpo_vdif(imt,jmt,km),dpo_vdif_mo(imt,jmt,km),sdpo_vdif(simt,sjmt,km),dpo_vdif_zhy(simt-2,sjmt,km)
      real din_vdif(imt,jmt,km),din_vdif_mo(imt,jmt,km),sdin_vdif(simt,sjmt,km),din_vdif_zhy(simt-2,sjmt,km)
      real dpo_bc  (imt,jmt)   ,dpo_bc_mo  (imt,jmt)   ,sdpo_bc  (simt,sjmt)   ,dpo_bc_zhy  (simt-2,sjmt)
      real din_bc  (imt,jmt)   ,din_bc_mo  (imt,jmt)   ,sdin_bc  (simt,sjmt)   ,din_bc_zhy  (simt-2,sjmt)
      real dpo_imp (imt,jmt,km),dpo_imp_mo (imt,jmt,km),sdpo_imp (simt,sjmt,km),dpo_imp_zhy (simt-2,sjmt,km)
      real din_imp (imt,jmt,km),din_imp_mo (imt,jmt,km),sdin_imp (simt,sjmt,km),din_imp_zhy (simt-2,sjmt,km)
      real dpo_ice (imt,jmt,km),dpo_ice_mo (imt,jmt,km),sdpo_ice (simt,sjmt,km),dpo_ice_zhy (simt-2,sjmt,km)
      real din_ice (imt,jmt,km),din_ice_mo (imt,jmt,km),sdin_ice (simt,sjmt,km),din_ice_zhy (simt-2,sjmt,km)
      real din_ast (imt,jmt,km),din_ast_mo (imt,jmt,km),sdin_ast (simt,sjmt,km),din_ast_zhy (simt-2,sjmt,km)
      real dpo_con (imt,jmt,km),dpo_con_mo (imt,jmt,km),sdpo_con (simt,sjmt,km),dpo_con_zhy (simt-2,sjmt,km)
      real din_con (imt,jmt,km),din_con_mo (imt,jmt,km),sdin_con (simt,sjmt,km),din_con_zhy (simt-2,sjmt,km)
      real wind_en (imt,jmt)   ,wind_en_mo (imt,jmt)   ,swind_en (simt,sjmt)   ,wind_en_zhy (simt-2,sjmt)
      real uwind_en(imt,jmt)   ,vwind_en   (imt,jmt)
      real dke_bcf (imt,jmt,2) ,dke_bcf_mo (imt,jmt)   ,sdke_bcf (simt,sjmt)   ,dke_bcf_zhy (simt-2,sjmt)
      real dke_adv (imt,jmt,km,2),dke_adv_mo(imt,jmt,km),sdke_adv(simt,sjmt,km),dke_adv_zhy (simt-2,sjmt,km)
      real dke_ape (imt,jmt,km,2),dke_ape_mo(imt,jmt,km),sdke_ape(simt,sjmt,km),dke_ape_zhy (simt-2,sjmt,km)
      real dke_fri (imt,jmt,km,2),dke_fri_mo(imt,jmt,km),sdke_fri(simt,sjmt,km),dke_fri_zhy (simt-2,sjmt,km)
      real dke_pre (imt,jmt,km,2),dke_pre_mo(imt,jmt,km),sdke_pre(simt,sjmt,km),dke_pre_zhy (simt-2,sjmt,km)
      real dke_bar (imt,jmt,km,2),dke_bar_mo(imt,jmt,km),sdke_bar(simt,sjmt,km),dke_bar_zhy (simt-2,sjmt,km)
      real dke_cor (imt,jmt,km,2),dke_cor_mo(imt,jmt,km),sdke_cor(simt,sjmt,km),dke_cor_zhy (simt-2,sjmt,km)
      real dpo_bar (imt,jmt,km),dpo_bar_mo (imt,jmt,km),sdpo_bar (simt,sjmt,km),dpo_bar_zhy (simt-2,sjmt,km)
      real din_bar (imt,jmt,km),din_bar_mo (imt,jmt,km),sdin_bar (simt,sjmt,km),din_bar_zhy (simt-2,sjmt,km)
      real gr_massu(imt,jmt,km),gr_masst(imt,jmt,km)
      real sgr_mass(simt,sjmt,km),gr_mass_mo(imt,jmt,km),gr_mass_zhy(simt-2,sjmt,km)
! information output var & array
      real sa_ini,ma_ini,ke_ini,in_ini,po_ini
      real sa_tol,ma_tol,ke_tol,in_tol,po_tol
      real sa_chg,ma_chg,ke_chg,in_chg,po_chg
      real dpoa_ts,dpoh_ts,dpov_ts,dpob_ts,dpoi_ts,dpoic_ts,dpoas_ts,dpoc_ts,dpob_ar,dpob_cp
      real dina_ts,dinh_ts,dinv_ts,dinb_ts,dini_ts,dinic_ts,dinas_ts,dinc_ts,dinb_ar
      real dkea_uv,dkef_uv,dkep_uv,dkeb_uv,dkeb_cf,dkea_pe ,dkec_or
      real dsaas_ts,dkebc,dkebcu,dkebcv,ev,eku,ekv,ekw
      real eku_en(imt,jmt,km),ekv_en(imt,jmt,km),ekw_en(imt,jmt,km)
      real total_in(imt,jmt,km),total_po(imt,jmt,km)
      real dpo_ast(imt,jmt,km),dsa_ast(imt,jmt,km),dpo_bcp(imt,jmt,km)
      integer sa_first
! restart output array
      real  pbt_st( imt, jmt,4), upb( imt, jmt,2), vpb( imt, jmt,2)
      real pbts_st(simt,sjmt,4),supb(simt,sjmt,2),svpb(simt,sjmt,2)
      real  up( imt, jmt,km,2), vp( imt, jmt,km,2)
      real sup(simt,sjmt,km,2),svp(simt,sjmt,km,2)
      real  adv_u( imt, jmt,km,3), adv_v( imt, jmt,km,3)
      real sadv_u(simt,sjmt,km,3),sadv_v(simt,sjmt,km,3)
      real t(imt,jmt,km,nt,2),st(simt,sjmt,km,nt,2)
      real pbt ( imt, jmt,2), ump( imt, jmt,km), vmp( imt, jmt,km), pmtp( imt, jmt)
      real pbts(simt,sjmt,2),sump(simt,sjmt,km),svmp(simt,sjmt,km),spmtp(simt,sjmt)
! other
      real bcu(imt,jmt),bcv(imt,jmt)
      real gr_mass(imt,jmt,km),gr_mass2(imt,jmt,km)
!
      call gath_array_real3d(stmask,tmask,mat_myid,ncpux,ncpuy,simt,sjmt,  &
                             km,imt,jmt,myid)
      call gath_array_real3d(sumask,umask,mat_myid,ncpux,ncpuy,simt,sjmt,  &
                             km,imt,jmt,myid)
      call gath_array_real1d(sdxdyt,dxdyt,mat_myid,ncpux,ncpuy,sjmt,jmt,myid)

!-----------------------------------------------------------------------
!     set mean parameter (factor)
!-----------------------------------------------------------------------
      runtime=runtime+dtts
      diag_out=0
      factor=0.0
      if (io_month==1) then
        factor=c1/(daypm(mth)*nss)
        if (day.eq.daypm(mth).and.mode_t.eq.nss) then
          diag_out=1
          encode(9,'(i9)',ftail) month+100000000
          normalname='N'//ftail(2:9)//'.nc'
          energyname='E'//ftail(2:9)//'.nc'
          time_units='months since 0001-01-15'
          time_id=month
        end if
      end if
      if (io_time>0) then
        factor=c1/(io_time/dtts)
        if (mod(runtime,io_time)==0) then
          diag_out=1
          encode(9,'(i9)',ftail) int(runtime+100000000)
          normalname='N'//ftail(2:9)//'.nc'
          energyname='E'//ftail(2:9)//'.nc'
          time_units='seconds since 0001-01-01-00.00.00'
          time_id=runtime
        end if
      end if
      
      restr_out=0
      if (restr_mon>0) then
        if (mod(month,restr_mon)==0) then
        if (day.eq.daypm(mth).and.mode_t.eq.nss) then
          restr_out=1
          encode(9,'(i9)',ftail) month+100000000
          fname='S'//ftail(2:9)
        end if
        end if
      end if
      if (restr_time>0) then
        if (mod(runtime,restr_time)==0) then
          restr_out=1
          encode(9,'(i9)',ftail) int(runtime+100000000)
          fname='S'//ftail(2:9)
        end if
      end if
!-----------------------------------------------------------------------
!     calculate mean array for output
!-----------------------------------------------------------------------

      do i=1,imt
        do j=1,jmt
          do k=1,km
            tmn_mo(i,j,k) = tmn_mo(i,j,k) + t(i,j,k,1,tau) * factor
            smn_mo(i,j,k) = smn_mo(i,j,k) + t(i,j,k,2,tau) * factor
            umn_mo(i,j,k) = umn_mo(i,j,k) + umn(i,j,k) * factor / 100.0d0
            vmn_mo(i,j,k) = vmn_mo(i,j,k) + vmn(i,j,k) * factor / 100.0d0
            wmn_mo(i,j,k) = wmn_mo(i,j,k) + wmn(i,j,k) * factor / 100.0d0
          end do
          pmn_mo(i,j) = pmn_mo(i,j) + pmn(i,j) * factor
        end do
      end do
      
!     calculate instantaneous information
      do i=1,imt
        do j=1,jmt
          do k=1,km
            gr_massu(i,j,k) = gr_mass(i,j,k)*dxdyu(j)/1000.0d0
            gr_masst(i,j,k) = gr_mass(i,j,k)*dxdyt(j)/1000.0d0
            gr_mass2(i,j,k) = gr_mass2(i,j,k)*dxdyt(j)/1000.0d0
            total_in(i,j,k) = gr_masst(i,j,k)*total_in(i,j,k)
            total_po(i,j,k) = gr_masst(i,j,k)*total_po(i,j,k)/10000.0d0
            eku_en(i,j,k)   = p5*gr_massu(i,j,k)*(umn(i,j,k)**2)/10000.0d0
            ekv_en(i,j,k)   = p5*gr_massu(i,j,k)*(vmn(i,j,k)**2)/10000.0d0
            ekw_en(i,j,k)   = p5*gr_massu(i,j,k)*(wmn(i,j,k)**2)/10000.0d0
          end do
          uwind_en(i,j) = bcu(i,j) * umn(i,j,1) * dxdyu(j) / 10000000.0d0
          vwind_en(i,j) = bcv(i,j) * vmn(i,j,1) * dxdyu(j) / 10000000.0d0
        end do
      end do
      
!     calculate energy
      if (energydiag.eq.1) then
      do j=2,jmm
         do i=2,imm
            do k=1,itn(i,j)
               dpo_adv(i,j,k)  = dpo_adv(i,j,k)  * dxdyt(j)/10000000.0d0/dtts
               dpo_hdif(i,j,k) = dpo_hdif(i,j,k) * dxdyt(j)/10000000.0d0/dtts
               dpo_vdif(i,j,k) = dpo_vdif(i,j,k) * dxdyt(j)/10000000.0d0/dtts
               dpo_imp(i,j,k)  = dpo_imp(i,j,k)  * dxdyt(j)/10000000.0d0/dtts
               dpo_ice(i,j,k)  = dpo_ice(i,j,k)  * dxdyt(j)/10000000.0d0/dtts
               dpo_ast(i,j,k)  = dpo_ast(i,j,k)  * dxdyt(j)/10000000.0d0/dtts
               dpo_con(i,j,k)  = dpo_con(i,j,k)  * dxdyt(j)/10000000.0d0/dtts
               dpo_bar(i,j,k)  = dpo_bar(i,j,k)  * dxdyt(j)/10000000.0d0/dtts
               dpo_bcp(i,j,k)  = dpo_bcp(i,j,k)  * dxdyt(j)/10000000.0d0/dtts
               
               din_adv(i,j,k)  = din_adv(i,j,k)  * gr_masst(i,j,k)/dtts
               din_hdif(i,j,k) = din_hdif(i,j,k) * gr_masst(i,j,k)/dtts
               din_vdif(i,j,k) = din_vdif(i,j,k) * gr_masst(i,j,k)/dtts
               din_imp(i,j,k)  = din_imp(i,j,k)  * gr_masst(i,j,k)/dtts
               din_ice(i,j,k)  = din_ice(i,j,k)  * gr_masst(i,j,k)/dtts
               din_ast(i,j,k)  = din_ast(i,j,k)  * gr_masst(i,j,k)/dtts
               din_con(i,j,k)  = din_con(i,j,k)  * gr_masst(i,j,k)/dtts
               din_bar(i,j,k)  = din_bar(i,j,k)  * gr_masst(i,j,k)/dtts
               
               dke_adv(i,j,k,1)= p5*(dke_adv(i,j,k,1)+dke_adv(i,j,k,2))  * gr_massu(i,j,k)/10000.0d0/dtuv
               dke_fri(i,j,k,1)= p5*(dke_fri(i,j,k,1)+dke_fri(i,j,k,2))  * gr_massu(i,j,k)/10000.0d0/dtuv
               dke_pre(i,j,k,1)= p5*(dke_pre(i,j,k,1)+dke_pre(i,j,k,2))  * gr_massu(i,j,k)/10000.0d0/dtuv
               dke_bar(i,j,k,1)= p5*(dke_bar(i,j,k,1)+dke_bar(i,j,k,2))  * gr_massu(i,j,k)/10000.0d0/dtuv
               dke_ape(i,j,k,1)= p5*(dke_ape(i,j,k,1)+dke_ape(i,j,k,2))  * gr_massu(i,j,k)/10000.0d0/dtuv
               dke_cor(i,j,k,1)= p5*(dke_cor(i,j,k,1)+dke_cor(i,j,k,2))  * gr_massu(i,j,k)/10000.0d0/dtuv
               
               dpo_adv_mo(i,j,k)  = dpo_adv_mo(i,j,k)  + dpo_adv(i,j,k)  * factor*10000.0d0/dxdyt(j)
               dpo_hdif_mo(i,j,k) = dpo_hdif_mo(i,j,k) + dpo_hdif(i,j,k) * factor*10000.0d0/dxdyt(j)
               dpo_vdif_mo(i,j,k) = dpo_vdif_mo(i,j,k) + dpo_vdif(i,j,k) * factor*10000.0d0/dxdyt(j)
               dpo_imp_mo(i,j,k)  = dpo_imp_mo(i,j,k)  + dpo_imp(i,j,k)  * factor*10000.0d0/dxdyt(j)
               dpo_ice_mo(i,j,k)  = dpo_ice_mo(i,j,k)  + dpo_ice(i,j,k)  * factor*10000.0d0/dxdyt(j)
               dpo_con_mo(i,j,k)  = dpo_con_mo(i,j,k)  + dpo_con(i,j,k)  * factor*10000.0d0/dxdyt(j)
               dpo_bar_mo(i,j,k)  = dpo_bar_mo(i,j,k)  + dpo_bar(i,j,k)  * factor*10000.0d0/dxdyt(j)
               
               din_adv_mo(i,j,k)  = din_adv_mo(i,j,k)  + din_adv(i,j,k)  * factor*10000.0d0/dxdyt(j)
               din_hdif_mo(i,j,k) = din_hdif_mo(i,j,k) + din_hdif(i,j,k) * factor*10000.0d0/dxdyt(j)
               din_vdif_mo(i,j,k) = din_vdif_mo(i,j,k) + din_vdif(i,j,k) * factor*10000.0d0/dxdyt(j)
               din_imp_mo(i,j,k)  = din_imp_mo(i,j,k)  + din_imp(i,j,k)  * factor*10000.0d0/dxdyt(j)
               din_ice_mo(i,j,k)  = din_ice_mo(i,j,k)  + din_ice(i,j,k)  * factor*10000.0d0/dxdyt(j)
               din_ast_mo(i,j,k)  = din_ast_mo(i,j,k)  + din_ast(i,j,k)  * factor*10000.0d0/dxdyt(j)
               din_con_mo(i,j,k)  = din_con_mo(i,j,k)  + din_con(i,j,k)  * factor*10000.0d0/dxdyt(j)
               din_bar_mo(i,j,k)  = din_bar_mo(i,j,k)  + din_bar(i,j,k)  * factor*10000.0d0/dxdyt(j)
               
               dke_adv_mo(i,j,k)  = dke_adv_mo(i,j,k)  + dke_adv(i,j,k,1)  * factor*10000.0d0/dxdyu(j)
               dke_fri_mo(i,j,k)  = dke_fri_mo(i,j,k)  + dke_fri(i,j,k,1)  * factor*10000.0d0/dxdyu(j)
               dke_pre_mo(i,j,k)  = dke_pre_mo(i,j,k)  + dke_pre(i,j,k,1)  * factor*10000.0d0/dxdyu(j)
               dke_bar_mo(i,j,k)  = dke_bar_mo(i,j,k)  + dke_bar(i,j,k,1)  * factor*10000.0d0/dxdyu(j)
               dke_ape_mo(i,j,k)  = dke_ape_mo(i,j,k)  + dke_ape(i,j,k,1)  * factor*10000.0d0/dxdyu(j)
               dke_cor_mo(i,j,k)  = dke_cor_mo(i,j,k)  + dke_cor(i,j,k,1)  * factor*10000.0d0/dxdyu(j)
               gr_mass_mo(i,j,k)  = gr_mass_mo(i,j,k)  + gr_masst(i,j,k) * factor
               
               
            end do
            dpo_bc(i,j) = dpo_bc(i,j) * dxdyt(j)/10000000.0d0/dtts
            din_bc(i,j) = din_bc(i,j) * gr_masst(i,j,1)/dtts

            dke_bcf(i,j,1)  = p5*(dke_bcf(i,j,1)+dke_bcf(i,j,2))  * gr_massu(i,j,1)/10000.0d0/dtuv
            
            dpo_bc_mo(i,j) = dpo_bc_mo(i,j) + dpo_bc(i,j) * factor*10000.0d0/dxdyt(j)
            din_bc_mo(i,j) = din_bc_mo(i,j) + din_bc(i,j) * factor*10000.0d0/dxdyt(j)
            
            dke_bcf_mo(i,j)=dke_bcf_mo(i,j)+dke_bcf(i,j,1) * factor*10000.0d0/dxdyu(j)
            wind_en_mo(i,j)=wind_en_mo(i,j)+(uwind_en(i,j)+vwind_en(i,j))*factor*10000d0/dxdyu(j)
         end do
      end do
      end if
!-----------------------------------------------------------------------
!     output
!-----------------------------------------------------------------------
      if (diag_out==1) then
        fillvaule=9.99e30
        if (uvtshdiag==1) then
          call gath_array_real3d(sumn,umn_mo,mat_myid,ncpux,ncpuy,simt,sjmt,km,imt,jmt,myid)
          call gath_array_real3d(svmn,vmn_mo,mat_myid,ncpux,ncpuy,simt,sjmt,km,imt,jmt,myid)
          call gath_array_real3d(swmn,wmn_mo,mat_myid,ncpux,ncpuy,simt,sjmt,km,imt,jmt,myid)
          call gath_array_real3d(stmn,tmn_mo,mat_myid,ncpux,ncpuy,simt,sjmt,km,imt,jmt,myid)
          call gath_array_real3d(ssmn,smn_mo,mat_myid,ncpux,ncpuy,simt,sjmt,km,imt,jmt,myid)
          call gath_array_real2d(spsi,pmn_mo,mat_myid,ncpux,ncpuy,simt,sjmt,imt,jmt,myid)
          !call gath_array_real3d(sam,am,mat_myid,ncpux,ncpuy,simt,sjmt,km,imt,jmt,myid)
          if (myid==0) then
            ev = c0
            do j=2,sjmt-1
              do i=2,simt-1
                ev = ev + spsi(i,j)*sdxdyt(j)*stmask(i,j,1)/area
              end do
            end do
            do j=2,sjmt-1
              do i=2,simt-1
                spsi(i,j) = spsi(i,j) - ev
              end do
            end do
            uzhy=fillvaule
            vzhy=fillvaule
            wzhy=fillvaule
            tzhy=fillvaule
            szhy=fillvaule
            pzhy=fillvaule
            amzhy=fillvaule
            do i=2,simt-1
              do j=1,sjmt
                do k=1,km
                  if (sumask(i,j,k).gt.0) then
                    uzhy(i-1,j,k)=sumn(i,j,k)
                    vzhy(i-1,j,k)=svmn(i,j,k)
                    wzhy(i-1,j,k)=swmn(i,j,k)
                  end if
                  if (stmask(i,j,k).gt.0) then
                    tzhy(i-1,j,k)=stmn(i,j,k)
                    szhy(i-1,j,k)=ssmn(i,j,k)
                    amzhy(i-1,j,k)=sam(i,j,k)
                  end if
                end do
                if (stmask(i,j,1).gt.0) then
                  pzhy(i-1,j)=spsi(i,j)/100d0
                end if
              end do
            end do
            call netcdf_read_cor(ncname,lat,lon,z,simt-2,sjmt,km)
            call netcdf_write_cor(normalname,lat,lon,z,simt-2,sjmt,km,time_id,time_units)
            call netcdf_write_var(normalname,uname,uzhy,simt-2,sjmt,km,u_units,u_longname,fillvaule)
            call netcdf_write_var(normalname,vname,vzhy,simt-2,sjmt,km,v_units,v_longname,fillvaule)
            call netcdf_write_var(normalname,wname,wzhy,simt-2,sjmt,km,w_units,w_longname,fillvaule)
            call netcdf_write_var(normalname,tname,tzhy,simt-2,sjmt,km,t_units,t_longname,fillvaule)
            call netcdf_write_var(normalname,sname,szhy,simt-2,sjmt,km,s_units,s_longname,fillvaule)
!            call netcdf_write_var(normalname,amname,amzhy,simt-2,sjmt,km,am_units,am_longname,fillvaule)
            call netcdf_write_var2d(normalname,pname,pzhy,simt-2,sjmt,p_units,p_longname,fillvaule)
          end if !(myid==0)
          tmn_mo = c0
          smn_mo = c0
          pmn_mo = c0
          umn_mo = c0
          vmn_mo = c0
          wmn_mo = c0
        end if !(uvtshdiag==0)
        
        if (energydiag.eq.1) then
          call swap_array_real3d(dpo_adv_mo ,imt,jmt,km,west,east,north,south)
          call swap_array_real3d(dpo_hdif_mo,imt,jmt,km,west,east,north,south)
          call swap_array_real3d(dpo_vdif_mo,imt,jmt,km,west,east,north,south)
          call swap_array_real3d(dpo_imp_mo ,imt,jmt,km,west,east,north,south)
          call swap_array_real3d(dpo_ice_mo ,imt,jmt,km,west,east,north,south)
          call swap_array_real3d(dpo_con_mo ,imt,jmt,km,west,east,north,south)
          call swap_array_real3d(dpo_bar_mo ,imt,jmt,km,west,east,north,south)
          
          call swap_array_real3d(din_adv_mo ,imt,jmt,km,west,east,north,south)
          call swap_array_real3d(din_hdif_mo,imt,jmt,km,west,east,north,south)
          call swap_array_real3d(din_vdif_mo,imt,jmt,km,west,east,north,south)
          call swap_array_real3d(din_imp_mo ,imt,jmt,km,west,east,north,south)
          call swap_array_real3d(din_ice_mo ,imt,jmt,km,west,east,north,south)
          call swap_array_real3d(din_ast_mo ,imt,jmt,km,west,east,north,south)
          call swap_array_real3d(din_con_mo ,imt,jmt,km,west,east,north,south)
          call swap_array_real3d(din_bar_mo ,imt,jmt,km,west,east,north,south)
          
          call swap_array_real3d(dke_adv_mo ,imt,jmt,km,west,east,north,south)
          call swap_array_real3d(dke_pre_mo ,imt,jmt,km,west,east,north,south)
          call swap_array_real3d(dke_fri_mo ,imt,jmt,km,west,east,north,south)
          call swap_array_real3d(dke_bar_mo ,imt,jmt,km,west,east,north,south)
          call swap_array_real3d(dke_ape_mo ,imt,jmt,km,west,east,north,south)
          call swap_array_real3d(dke_cor_mo ,imt,jmt,km,west,east,north,south)
          
          call swap_array_real2d(dpo_bc_mo ,imt,jmt,west,east,north,south)
          call swap_array_real2d(din_bc_mo ,imt,jmt,west,east,north,south)
          call swap_array_real2d(wind_en_mo,imt,jmt,west,east,north,south)
          call swap_array_real2d(dke_bcf_mo,imt,jmt,west,east,north,south)
          
          call gath_array_real3d(sdpo_adv ,dpo_adv_mo ,mat_myid,ncpux,ncpuy,simt,sjmt,km,imt,jmt,myid)
          call gath_array_real3d(sdpo_hdif,dpo_hdif_mo,mat_myid,ncpux,ncpuy,simt,sjmt,km,imt,jmt,myid)
          call gath_array_real3d(sdpo_vdif,dpo_vdif_mo,mat_myid,ncpux,ncpuy,simt,sjmt,km,imt,jmt,myid)
          call gath_array_real3d(sdpo_imp ,dpo_imp_mo ,mat_myid,ncpux,ncpuy,simt,sjmt,km,imt,jmt,myid)
          call gath_array_real3d(sdpo_ice ,dpo_ice_mo ,mat_myid,ncpux,ncpuy,simt,sjmt,km,imt,jmt,myid)
          call gath_array_real3d(sdpo_con ,dpo_con_mo ,mat_myid,ncpux,ncpuy,simt,sjmt,km,imt,jmt,myid)
          call gath_array_real3d(sdpo_bar ,dpo_bar_mo ,mat_myid,ncpux,ncpuy,simt,sjmt,km,imt,jmt,myid)
          
          call gath_array_real3d(sdin_adv ,din_adv_mo ,mat_myid,ncpux,ncpuy,simt,sjmt,km,imt,jmt,myid)
          call gath_array_real3d(sdin_hdif,din_hdif_mo,mat_myid,ncpux,ncpuy,simt,sjmt,km,imt,jmt,myid)
          call gath_array_real3d(sdin_vdif,din_vdif_mo,mat_myid,ncpux,ncpuy,simt,sjmt,km,imt,jmt,myid)
          call gath_array_real3d(sdin_imp ,din_imp_mo ,mat_myid,ncpux,ncpuy,simt,sjmt,km,imt,jmt,myid)
          call gath_array_real3d(sdin_ice ,din_ice_mo ,mat_myid,ncpux,ncpuy,simt,sjmt,km,imt,jmt,myid)
          call gath_array_real3d(sdin_ast ,din_ast_mo ,mat_myid,ncpux,ncpuy,simt,sjmt,km,imt,jmt,myid)
          call gath_array_real3d(sdin_con ,din_con_mo ,mat_myid,ncpux,ncpuy,simt,sjmt,km,imt,jmt,myid)
          call gath_array_real3d(sdin_bar ,din_bar_mo ,mat_myid,ncpux,ncpuy,simt,sjmt,km,imt,jmt,myid)
          
          call gath_array_real3d(sdke_adv ,dke_adv_mo ,mat_myid,ncpux,ncpuy,simt,sjmt,km,imt,jmt,myid)
          call gath_array_real3d(sdke_fri ,dke_fri_mo ,mat_myid,ncpux,ncpuy,simt,sjmt,km,imt,jmt,myid)
          call gath_array_real3d(sdke_pre ,dke_pre_mo ,mat_myid,ncpux,ncpuy,simt,sjmt,km,imt,jmt,myid)
          call gath_array_real3d(sdke_bar ,dke_bar_mo ,mat_myid,ncpux,ncpuy,simt,sjmt,km,imt,jmt,myid)
          call gath_array_real3d(sdke_ape ,dke_ape_mo ,mat_myid,ncpux,ncpuy,simt,sjmt,km,imt,jmt,myid)
          call gath_array_real3d(sdke_cor ,dke_cor_mo ,mat_myid,ncpux,ncpuy,simt,sjmt,km,imt,jmt,myid)
          
          call gath_array_real2d(sdpo_bc ,dpo_bc_mo ,mat_myid,ncpux,ncpuy,simt,sjmt,imt,jmt,myid)
          call gath_array_real2d(sdin_bc ,din_bc_mo ,mat_myid,ncpux,ncpuy,simt,sjmt,imt,jmt,myid)
          call gath_array_real2d(swind_en,wind_en_mo,mat_myid,ncpux,ncpuy,simt,sjmt,imt,jmt,myid)
          call gath_array_real2d(sdke_bcf,dke_bcf_mo,mat_myid,ncpux,ncpuy,simt,sjmt,imt,jmt,myid)
          
          dpo_adv_zhy =fillvaule
          dpo_hdif_zhy=fillvaule
          dpo_vdif_zhy=fillvaule
          dpo_imp_zhy =fillvaule
          dpo_ice_zhy =fillvaule
          dpo_con_zhy =fillvaule
          dpo_bar_zhy =fillvaule
          
          din_adv_zhy =fillvaule
          din_hdif_zhy=fillvaule
          din_vdif_zhy=fillvaule
          din_imp_zhy =fillvaule
          din_ice_zhy =fillvaule
          din_ast_zhy =fillvaule
          din_con_zhy =fillvaule
          din_bar_zhy =fillvaule
          
          dke_adv_zhy =fillvaule
          dke_fri_zhy =fillvaule
          dke_pre_zhy =fillvaule
          dke_bar_zhy =fillvaule
          dke_ape_zhy =fillvaule
          dke_cor_zhy =fillvaule
          
          dpo_bc_zhy =fillvaule
          din_bc_zhy =fillvaule
          wind_en_zhy=fillvaule
          dke_bcf_zhy=fillvaule
          
          if (myid==0) then
          do j=1,sjmt
             do i=2,simt-1
                do k=1,km
                   if (stmask(i,j,k).gt.0) then
                      dpo_adv_zhy(i-1,j,k)  = sdpo_adv(i,j,k)
                      dpo_hdif_zhy(i-1,j,k) = sdpo_hdif(i,j,k)
                      dpo_vdif_zhy(i-1,j,k) = sdpo_vdif(i,j,k)
                      dpo_imp_zhy(i-1,j,k)  = sdpo_imp(i,j,k)
                      dpo_ice_zhy(i-1,j,k)  = sdpo_ice(i,j,k)
                      dpo_con_zhy(i-1,j,k)  = sdpo_con(i,j,k)
                      dpo_bar_zhy(i-1,j,k)  = sdpo_bar(i,j,k)
                      
                      din_adv_zhy(i-1,j,k)  = sdin_adv(i,j,k)
                      din_hdif_zhy(i-1,j,k) = sdin_hdif(i,j,k)
                      din_vdif_zhy(i-1,j,k) = sdin_vdif(i,j,k)
                      din_imp_zhy(i-1,j,k)  = sdin_imp(i,j,k)
                      din_ice_zhy(i-1,j,k)  = sdin_ice(i,j,k)
                      din_ast_zhy(i-1,j,k)  = sdin_ast(i,j,k)
                      din_con_zhy(i-1,j,k)  = sdin_con(i,j,k)
                      din_bar_zhy(i-1,j,k)  = sdin_bar(i,j,k)
                   end if
                   if (sumask(i,j,k).gt.0) then
                      dke_adv_zhy(i-1,j,k)  = sdke_adv(i,j,k)
                      dke_fri_zhy(i-1,j,k)  = sdke_fri(i,j,k)
                      dke_pre_zhy(i-1,j,k)  = sdke_pre(i,j,k)
                      dke_bar_zhy(i-1,j,k)  = sdke_bar(i,j,k)
                      dke_ape_zhy(i-1,j,k)  = sdke_ape(i,j,k)
                      dke_cor_zhy(i-1,j,k)  = sdke_cor(i,j,k)
                   end if
                end do
                if (stmask(i,j,1).gt.0) then
                   dpo_bc_zhy(i-1,j)=sdpo_bc(i,j)
                   din_bc_zhy(i-1,j)=sdin_bc(i,j)
                end if
                if (sumask(i,j,1).gt.0) then
                   wind_en_zhy(i-1,j)=swind_en(i,j)
                   dke_bcf_zhy(i-1,j)=sdke_bcf(i,j)
                end if
             end do
          end do
          call netcdf_read_cor(ncname,lat,lon,z,simt-2,sjmt,km)
          call netcdf_write_cor(energyname,lat,lon,z,simt-2,sjmt,km,time_id,time_units)
          
          call netcdf_write_var(energyname,dpo_adv_name ,dpo_adv_zhy ,simt-2,sjmt,km,dpo_adv_units ,dpo_adv_longname ,fillvaule)
          call netcdf_write_var(energyname,dpo_hdif_name,dpo_hdif_zhy,simt-2,sjmt,km,dpo_hdif_units,dpo_hdif_longname,fillvaule)
          call netcdf_write_var(energyname,dpo_vdif_name,dpo_vdif_zhy,simt-2,sjmt,km,dpo_vdif_units,dpo_vdif_longname,fillvaule)
          call netcdf_write_var(energyname,dpo_imp_name ,dpo_imp_zhy ,simt-2,sjmt,km,dpo_imp_units ,dpo_imp_longname ,fillvaule)
          call netcdf_write_var(energyname,dpo_ice_name ,dpo_ice_zhy ,simt-2,sjmt,km,dpo_ice_units ,dpo_ice_longname ,fillvaule)
          call netcdf_write_var(energyname,dpo_con_name ,dpo_con_zhy ,simt-2,sjmt,km,dpo_con_units ,dpo_con_longname ,fillvaule)
          call netcdf_write_var(energyname,dpo_bar_name ,dpo_bar_zhy ,simt-2,sjmt,km,dpo_bar_units ,dpo_bar_longname ,fillvaule)
          
          call netcdf_write_var(energyname,din_adv_name ,din_adv_zhy ,simt-2,sjmt,km,din_adv_units ,din_adv_longname ,fillvaule)
          call netcdf_write_var(energyname,din_hdif_name,din_hdif_zhy,simt-2,sjmt,km,din_hdif_units,din_hdif_longname,fillvaule)
          call netcdf_write_var(energyname,din_vdif_name,din_vdif_zhy,simt-2,sjmt,km,din_vdif_units,din_vdif_longname,fillvaule)
          call netcdf_write_var(energyname,din_imp_name ,din_imp_zhy ,simt-2,sjmt,km,din_imp_units ,din_imp_longname ,fillvaule)
          call netcdf_write_var(energyname,din_ice_name ,din_ice_zhy ,simt-2,sjmt,km,din_ice_units ,din_ice_longname ,fillvaule)
          !call netcdf_write_var(energyname,din_ast_name ,din_ast_zhy ,simt-2,sjmt,km,din_ast_units ,din_ast_longname ,fillvaule)
          call netcdf_write_var(energyname,din_con_name ,din_con_zhy ,simt-2,sjmt,km,din_con_units ,din_con_longname ,fillvaule)
          call netcdf_write_var(energyname,din_bar_name ,din_bar_zhy ,simt-2,sjmt,km,din_bar_units ,din_bar_longname ,fillvaule)
          
          call netcdf_write_var(energyname,dke_adv_name ,dke_adv_zhy ,simt-2,sjmt,km,dke_adv_units ,dke_adv_longname ,fillvaule)
          call netcdf_write_var(energyname,dke_fri_name ,dke_fri_zhy ,simt-2,sjmt,km,dke_fri_units ,dke_fri_longname ,fillvaule)
          call netcdf_write_var(energyname,dke_pre_name ,dke_pre_zhy ,simt-2,sjmt,km,dke_pre_units ,dke_pre_longname ,fillvaule)
          call netcdf_write_var(energyname,dke_bar_name ,dke_bar_zhy ,simt-2,sjmt,km,dke_bar_units ,dke_bar_longname ,fillvaule)
          call netcdf_write_var(energyname,dke_ape_name ,dke_ape_zhy ,simt-2,sjmt,km,dke_ape_units ,dke_ape_longname ,fillvaule)
          call netcdf_write_var(energyname,dke_cor_name ,dke_cor_zhy ,simt-2,sjmt,km,dke_cor_units ,dke_cor_longname ,fillvaule)
          
          call netcdf_write_var2d(energyname,dpo_bc_name ,dpo_bc_zhy ,simt-2,sjmt,dpo_bc_units ,dpo_bc_longname ,fillvaule)
          call netcdf_write_var2d(energyname,din_bc_name ,din_bc_zhy ,simt-2,sjmt,din_bc_units ,din_bc_longname ,fillvaule)          
          call netcdf_write_var2d(energyname,wind_en_name,wind_en_zhy,simt-2,sjmt,wind_en_units,wind_en_longname,fillvaule)
          call netcdf_write_var2d(energyname,dke_bcf_name,dke_bcf_zhy,simt-2,sjmt,dke_bcf_units,dke_bcf_longname,fillvaule)
          
          end if!(myid==0)
          dpo_adv_mo    = c0
          dpo_hdif_mo   = c0
          dpo_vdif_mo   = c0
          dpo_imp_mo    = c0
          dpo_ice_mo    = c0
          dpo_con_mo    = c0
          dpo_bar_mo    = c0
          
          din_adv_mo    = c0
          din_hdif_mo   = c0
          din_vdif_mo   = c0
          din_imp_mo    = c0
          din_ice_mo    = c0
          din_ast_mo    = c0
          din_con_mo    = c0
          din_bar_mo    = c0
          
          dke_adv_mo    = c0
          dke_fri_mo    = c0
          dke_pre_mo    = c0
          dke_bar_mo    = c0
          dke_ape_mo    = c0
          dke_cor_mo    = c0
          
          dpo_bc_mo     = c0
          din_bc_mo     = c0
          wind_en_mo    = c0
          dke_bcf_mo    = c0
        end if!(energydiag==0)
      end if !(diag_out==0)
!-----------------------------------------------------------------------
!     calculate instantaneous information
!-----------------------------------------------------------------------
      if (infotime>0) then
      if (mod(runtime,infotime)==0) then
      eku   = c0
      ekv   = c0
      ekw   = c0
      ma_tol= c0
      in_tol= c0
      sa_tol= c0
      po_tol= c0
      
      dpoa_ts = c0
      dpoh_ts = c0
      dpov_ts = c0
      dpoi_ts = c0
      dpoic_ts= c0
      dpoas_ts= c0
      dpoc_ts = c0
      dpob_ar = c0
      dpob_cp = c0
      
      dina_ts = c0
      dinh_ts = c0
      dinv_ts = c0
      dini_ts = c0
      dinic_ts= c0
      dinas_ts= c0
      dsaas_ts= c0
      dinc_ts = c0
      dinb_ar = c0
      
      dkebcu = c0
      dkebcv = c0
      dpob_ts = c0
      dinb_ts = c0
      
      dkea_uv = c0
      dkef_uv = c0
      dkep_uv = c0
      dkeb_uv = c0
      dkea_pe = c0
      dkec_or = c0
      dkeb_cf = c0
      
      if (energydiag.eq.1) then
!     advance information output
      do j=2,jmm
        do i=2,imm
          do k=1,km
            if (umask(i,j,k).gt.0) then
              eku = eku + eku_en(i,j,k)
              ekv = ekv + ekv_en(i,j,k)
              ekw = ekw + ekw_en(i,j,k)
              dkea_uv = dkea_uv + dke_adv(i,j,k,1)
              dkef_uv = dkef_uv + dke_fri(i,j,k,1)
              dkep_uv = dkep_uv + dke_pre(i,j,k,1)
              dkeb_uv = dkeb_uv + dke_bar(i,j,k,1)
              dkea_pe = dkea_pe + dke_ape(i,j,k,1)
              dkec_or = dkec_or + dke_cor(i,j,k,1)
            end if
            if (tmask(i,j,k).gt.0) then
              in_tol = in_tol + total_in(i,j,k)
              po_tol = po_tol + total_po(i,j,k)
              ma_tol = ma_tol + gr_mass2(i,j,k)
              sa_tol = sa_tol + t(i,j,k,2,tau)*gr_masst(i,j,k)/1000.0d0
              dsaas_ts=dsaas_ts + dsa_ast(i,j,k)*gr_masst(i,j,k)/1000.0d0
              
              dpoa_ts = dpoa_ts + dpo_adv(i,j,k)
              dpoh_ts = dpoh_ts + dpo_hdif(i,j,k)
              dpov_ts = dpov_ts + dpo_vdif(i,j,k)
              dpoi_ts = dpoi_ts + dpo_imp(i,j,k)
              dpoic_ts= dpoic_ts+ dpo_ice(i,j,k)
              dpoas_ts= dpoas_ts+ dpo_ast(i,j,k)
              dpoc_ts = dpoc_ts + dpo_con(i,j,k)
              dpob_ar = dpob_ar + dpo_bar(i,j,k)
              dpob_cp = dpob_cp + dpo_bcp(i,j,k)
              
              dina_ts = dina_ts + din_adv(i,j,k)
              dinh_ts = dinh_ts + din_hdif(i,j,k)
              dinv_ts = dinv_ts + din_vdif(i,j,k)
              dini_ts = dini_ts + din_imp(i,j,k)
              dinic_ts= dinic_ts+ din_ice(i,j,k)
              dinas_ts= dinas_ts+ din_ast(i,j,k)
              dinc_ts = dinc_ts + din_con(i,j,k)
              dinb_ar = dinb_ar + din_bar(i,j,k)
            end if
          end do
          if (umask(i,j,1).gt.0) then
             dkebcu = dkebcu + uwind_en(i,j)
             dkebcv = dkebcv + vwind_en(i,j)
             dkeb_cf = dkeb_cf + dke_bcf(i,j,1)
          end if
          if (tmask(i,j,1).gt.0) then
             dpob_ts = dpob_ts + dpo_bc(i,j)
             dinb_ts = dinb_ts + din_bc(i,j)
          end if
        end do
      end do
      
      call gath_var_real(eku,mat_myid,ncpux,ncpuy,myid)
      call gath_var_real(ekv,mat_myid,ncpux,ncpuy,myid)
      call gath_var_real(ekw,mat_myid,ncpux,ncpuy,myid)
      call gath_var_real(in_tol,mat_myid,ncpux,ncpuy,myid)
      call gath_var_real(po_tol,mat_myid,ncpux,ncpuy,myid)
      call gath_var_real(ma_tol,mat_myid,ncpux,ncpuy,myid)
      call gath_var_real(sa_tol,mat_myid,ncpux,ncpuy,myid)
      
      call gath_var_real(dpoa_ts,mat_myid,ncpux,ncpuy,myid)
      call gath_var_real(dpoh_ts,mat_myid,ncpux,ncpuy,myid)
      call gath_var_real(dpov_ts,mat_myid,ncpux,ncpuy,myid)
      call gath_var_real(dpoi_ts,mat_myid,ncpux,ncpuy,myid)
      call gath_var_real(dpoic_ts,mat_myid,ncpux,ncpuy,myid)
      call gath_var_real(dpoas_ts,mat_myid,ncpux,ncpuy,myid)
      call gath_var_real(dpoc_ts,mat_myid,ncpux,ncpuy,myid)
      call gath_var_real(dpob_ts,mat_myid,ncpux,ncpuy,myid)
      call gath_var_real(dpob_ar,mat_myid,ncpux,ncpuy,myid)
      call gath_var_real(dpob_cp,mat_myid,ncpux,ncpuy,myid)
      
      call gath_var_real(dina_ts,mat_myid,ncpux,ncpuy,myid)
      call gath_var_real(dinh_ts,mat_myid,ncpux,ncpuy,myid)
      call gath_var_real(dinv_ts,mat_myid,ncpux,ncpuy,myid)
      call gath_var_real(dini_ts,mat_myid,ncpux,ncpuy,myid)
      call gath_var_real(dinic_ts,mat_myid,ncpux,ncpuy,myid)
      call gath_var_real(dinas_ts,mat_myid,ncpux,ncpuy,myid)
      call gath_var_real(dsaas_ts,mat_myid,ncpux,ncpuy,myid)
      call gath_var_real(dinc_ts,mat_myid,ncpux,ncpuy,myid)
      call gath_var_real(dinb_ts,mat_myid,ncpux,ncpuy,myid)
      call gath_var_real(dinb_ar,mat_myid,ncpux,ncpuy,myid)
      
      call gath_var_real(dkea_uv,mat_myid,ncpux,ncpuy,myid)
      call gath_var_real(dkef_uv,mat_myid,ncpux,ncpuy,myid)
      call gath_var_real(dkep_uv,mat_myid,ncpux,ncpuy,myid)
      call gath_var_real(dkeb_uv,mat_myid,ncpux,ncpuy,myid)
      call gath_var_real(dkea_pe,mat_myid,ncpux,ncpuy,myid)
      call gath_var_real(dkec_or,mat_myid,ncpux,ncpuy,myid)
      
      call gath_var_real(dkebcu,mat_myid,ncpux,ncpuy,myid)
      call gath_var_real(dkebcv,mat_myid,ncpux,ncpuy,myid)
      call gath_var_real(dkeb_cf,mat_myid,ncpux,ncpuy,myid)
      
      ke_tol=eku+ekv+ekw
      if (myid==0) then
      if (sa_first==1) then
         ke_ini=ke_tol
         in_ini=in_tol
         po_ini=po_tol
         ma_ini=ma_tol
         sa_ini=sa_tol
         write(66,"(5e15.6)") ke_ini,in_ini,po_ini,ma_ini,sa_ini
         sa_first=0
      end if
      end if
      
      dkebc=dkebcu+dkebcv
      sa_chg=sa_tol-sa_ini
      ma_chg=ma_tol-ma_ini
      ke_chg=ke_tol-ke_ini
      in_chg=in_tol-in_ini
      po_chg=po_tol-po_ini
      
      sa_ini=sa_tol
      ma_ini=ma_tol
      ke_ini=ke_tol
      in_ini=in_tol
      po_ini=po_tol
      
      call gath_array_real2d(spsi,pmn,mat_myid,ncpux,ncpuy,simt,sjmt, &
                             imt,jmt,myid)      
      
      if (myid==0) then
      ev = c0
      do j=2,sjmt-1
      do i=2,simt-1
      ev = ev + spsi(i,j)*sdxdyt(j)*stmask(i,j,1)/area
      enddo
      enddo
      !write(*,"(a5,2e15.6)") "GPE",po_chg,(dpoa_ts+dpoh_ts+dpov_ts+dpoi_ts+dpob_ts+dpoic_ts+dpoc_ts+dpob_ar)*dtts
      !write(*,"(a5,2e15.6)") "INE",in_chg,(dina_ts+dinh_ts+dinv_ts+dini_ts+dinb_ts+dinic_ts+dinc_ts+dinb_ar)*dtts
      !write(*,"(a5,2e15.6)") "TKE",ke_chg,(dkea_uv+dkef_uv+dkep_uv+dkeb_uv+dkeb_cf+dkea_pe+dkec_or)*dtuv
      write(66,"(i5,i3,f15.2,38e15.6,f15.4)") month,day,runtime,  &
      ke_tol , in_tol , po_tol , ma_tol , sa_tol ,  &
      dpoa_ts, dpoh_ts, dpov_ts, dpoi_ts, dpob_ts , dpoic_ts, dpoas_ts, dpoc_ts, dpob_ar, dpob_cp,  &
      dina_ts, dinh_ts, dinv_ts, dini_ts, dinb_ts , dinic_ts, dinas_ts, dinc_ts, dinb_ar,  &
      ke_chg , in_chg , po_chg , ma_chg , sa_chg  , dsaas_ts, &
      dkea_uv, dkef_uv, dkep_uv, dkeb_uv, dkeb_cf , dkea_pe , dkec_or , dkebc  , ev
      end if
      
      else
!     basic information output
      do k=1,km
        do j=2,jmm
          do i=2,imm
            if (umask(i,j,k).gt.0) then
              eku = eku + eku_en(i,j,k)
              ekv = ekv + ekv_en(i,j,k)
              ekw = ekw + ekw_en(i,j,k)
            end if
          end do
        end do
      end do
      call gath_var_real(eku,mat_myid,ncpux,ncpuy,myid)
      call gath_var_real(ekv,mat_myid,ncpux,ncpuy,myid)
      call gath_var_real(ekw,mat_myid,ncpux,ncpuy,myid)
      ke_tol= eku+ekv+ekw
      call gath_array_real2d(spsi,pmn,mat_myid,ncpux,ncpuy,simt,sjmt, &
                             imt,jmt,myid)
      if (myid==0) then
      ev = c0
      do j=2,sjmt-1
      do i=2,simt-1
      ev = ev + spsi(i,j)*sdxdyt(j)*stmask(i,j,1)/area
      enddo
      enddo
      write(66,"(i5,i3,f15.2,e15.6,f15.4)") month,day,runtime,ke_tol,ev
      end if
      
      end if !energydiag.eq.1
      end if !(mod(runtime,infotime)==0)
      end if !infotime>0
!-----------------------------------------------------------------------
!     save datafile for restarting integration
!-----------------------------------------------------------------------
      if (restr_out==1) then
      if (stager_t.eq.1) then
        call swap_array_real4d(adv_u,imt,jmt,km,3,west,east,north,south)
        call swap_array_real4d(adv_v,imt,jmt,km,3,west,east,north,south)
        
        call gath_array_real3d(pbts_st,pbt_st,mat_myid,ncpux,ncpuy,simt,sjmt,4      ,imt,jmt,myid)
        call gath_array_real3d(supb   ,upb   ,mat_myid,ncpux,ncpuy,simt,sjmt,2      ,imt,jmt,myid)
        call gath_array_real3d(svpb   ,vpb   ,mat_myid,ncpux,ncpuy,simt,sjmt,2      ,imt,jmt,myid)
        call gath_array_real3d(sam    ,am    ,mat_myid,ncpux,ncpuy,simt,sjmt,km     ,imt,jmt,myid)
        call gath_array_real4d(sup    ,up    ,mat_myid,ncpux,ncpuy,simt,sjmt,km,2   ,imt,jmt,myid)
        call gath_array_real4d(svp    ,vp    ,mat_myid,ncpux,ncpuy,simt,sjmt,km,2   ,imt,jmt,myid)
        call gath_array_real4d(sadv_u ,adv_u ,mat_myid,ncpux,ncpuy,simt,sjmt,km,3   ,imt,jmt,myid)
        call gath_array_real4d(sadv_v ,adv_v ,mat_myid,ncpux,ncpuy,simt,sjmt,km,3   ,imt,jmt,myid)
        call gath_array_real5d(st     ,t     ,mat_myid,ncpux,ncpuy,simt,sjmt,km,nt,2,imt,jmt,myid)
        if (myid==0) then
          open(87,file=fname,form='unformatted',status='unknown')
          write(87) month,t_stepu,st,pbts_st,sup,svp,supb,svpb,sadv_u,sadv_v,sam
          print *,"write restart field at month ",month," for stager scheme"
          close(87)
        end if
      else
        call gath_array_real2d(spmtp,pmtp,mat_myid,ncpux,ncpuy,simt,sjmt,imt,jmt,myid)
        call gath_array_real3d(pbts ,pbt ,mat_myid,ncpux,ncpuy,simt,sjmt,2,imt,jmt,myid)
        call gath_array_real3d(supb ,upb ,mat_myid,ncpux,ncpuy,simt,sjmt,2,imt,jmt,myid)
        call gath_array_real3d(svpb ,vpb ,mat_myid,ncpux,ncpuy,simt,sjmt,2,imt,jmt,myid)
        call gath_array_real3d(sump ,ump ,mat_myid,ncpux,ncpuy,simt,sjmt,km,imt,jmt,myid)
        call gath_array_real3d(svmp ,vmp ,mat_myid,ncpux,ncpuy,simt,sjmt,km,imt,jmt,myid)
        call gath_array_real4d(sup  ,up  ,mat_myid,ncpux,ncpuy,simt,sjmt,km,2,imt,jmt,myid)
        call gath_array_real4d(svp  ,vp  ,mat_myid,ncpux,ncpuy,simt,sjmt,km,2,imt,jmt,myid)
        call gath_array_real5d(st   ,t   ,mat_myid,ncpux,ncpuy,simt,sjmt,km,nt,2,imt,jmt,myid)
        call gath_array_real3d(sam  ,am  ,mat_myid,ncpux,ncpuy,simt,sjmt,km,imt,jmt,myid)
        if (myid==0) then
          open(87,file=fname,form='unformatted',status='unknown')
          write(87)month,st,pbts,sup,svp,supb,svpb,spmtp,sump,svmp,sam
          print *,"write restart field at month ",month
          close(87)
        end if
      end if
      end if
!
      return
      end
