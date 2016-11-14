!====================== include file "pconst.h" ========================
!
!
!     rules for parameter constants
!
!     use prefix of "c" for whole real numbers (eg: c57 for 57.0)
!     use "m" after prefix to designate negative values (minus sign)
!       (eg: cm7 for -7.0)
!     use prefix of "p" for non repeating fractions (eg: p5 for 0.5)
!     use prefix of "r" for reciprocals (eg: r3 for 1/3.0)
!     combine use of prefix above and "e" for scientific notation, with
!       (eg: c5e4 for 5.0e4, c1em10 for 1.0e-10)
!
      real c0,c1,c2,c3,c4,c5,c6,c25,c35
      parameter (c0=0.0d0,c1=1.0d0,c2=2.0d0,c3=3.0d0,c4=4.0d0)
      parameter (c5=5.0d0,c6=6.0d0,c25=25.0d0,c35=35.0d0)
!
      real p125,p25,p5,p75
      parameter (p125=0.125d0,p25=0.25d0, p5=0.5d0, p75=0.75d0)
!
      real c60,c1440,c3600
      parameter (c60=60.0d0, c1440=1440.0d0, c3600=3600.0d0)
!
      real c1e3,c1em4,c1em12
      parameter (c1e3=1.0d3,c1em4=1.0d-4,c1em12=1.0d-12)
!
      real r120,r150,r180,r365
      parameter (r120=c1/120.0d0,r150=c1/150.0d0)
      parameter (r180=c1/180.0d0,r365=c1/365.0d0)
!
      real secday
      parameter (secday=c1/(c60*c1440))
!
      real rho_0,rrho_0
      parameter (rho_0=1.029d0,rrho_0=c1/rho_0)
!
      real tbice
      parameter (tbice=-1.5d0)
!
!      real fcof
!      parameter (fcof=0.05d0)
!
      real pi,torad,radius,omega,grav,cwater,kelvin
      parameter (pi     = 3.141592653589793)
      parameter (torad  = pi*r180    )
      parameter (radius = 6370.0d5   )
      parameter (omega  = pi/43082.0 )
      parameter (grav   = 980.0)
      parameter (cwater = 3901.0/1000.0)
      parameter (kelvin = 273.16)
!
!     grav = earth's gravitational acceleration (cm/sec**2)
!
      integer tau,taum
      parameter (tau=1,taum=2)
      
      
      character (len = *), parameter :: ncname = "pcom_ini.nc"
      character (len = *), parameter :: ctname = "ct"
      character (len = *), parameter :: saname = "sa"
      character (len = *), parameter :: t30name = "t30"
      character (len = *), parameter :: s30name = "s30"
            
      character (len = *), parameter :: uname = "u"
      character (len = *), parameter :: vname = "v"
      character (len = *), parameter :: wname = "w"
      character (len = *), parameter :: tname = "pt"
      character (len = *), parameter :: sname = "sa"
      character (len = *), parameter :: pname = "ssh"
      character (len = *), parameter :: amname = "am"
      character (len = *), parameter :: u_units = "m/s"
      character (len = *), parameter :: v_units = "m/s"
      character (len = *), parameter :: w_units = "m/s"
      character (len = *), parameter :: t_units = "C"
      character (len = *), parameter :: s_units = "g/kg"
      character (len = *), parameter :: p_units = "cm"
      character (len = *), parameter :: am_units = "cm**2/s"
      character (len = *), parameter :: u_longname = "zonal velocity"
      character (len = *), parameter :: v_longname = "meridional velocity"
      character (len = *), parameter :: w_longname = "vertical velocity"
      character (len = *), parameter :: t_longname = "potential temperature"
      character (len = *), parameter :: s_longname = "absolute salinity"
      character (len = *), parameter :: p_longname = "sea surface level"
      character (len = *), parameter :: am_longname = "horizontal turbulent dissipation coefficient"

      character (len = *), parameter :: dpo_adv_name  = "dpo_adv"
      character (len = *), parameter :: dpo_hdif_name = "dpo_hdif"
      character (len = *), parameter :: dpo_vdif_name = "dpo_vdif"
      character (len = *), parameter :: dpo_imp_name  = "dpo_imp"
      character (len = *), parameter :: dpo_ice_name  = "dpo_ice"
      character (len = *), parameter :: dpo_con_name  = "dpo_con"
      character (len = *), parameter :: dpo_bc_name   = "dpo_bc"
      character (len = *), parameter :: dpo_bar_name  = "dpo_bar"
      
      character (len = *), parameter :: din_adv_name  = "din_adv"
      character (len = *), parameter :: din_hdif_name = "din_hdif"
      character (len = *), parameter :: din_vdif_name = "din_vdif"
      character (len = *), parameter :: din_imp_name  = "din_imp"
      character (len = *), parameter :: din_ice_name  = "din_ice"
      character (len = *), parameter :: din_ast_name  = "din_ast"
      character (len = *), parameter :: din_con_name  = "din_con"
      character (len = *), parameter :: din_bc_name   = "din_bc"
      character (len = *), parameter :: din_bar_name  = "din_bar"
      
      character (len = *), parameter :: dke_adv_name  = "dke_adv"
      character (len = *), parameter :: dke_fri_name  = "dke_fri"
      character (len = *), parameter :: dke_pre_name  = "dke_pre"
      character (len = *), parameter :: dke_bar_name  = "dke_bar"
      character (len = *), parameter :: dke_bcf_name  = "dke_bcf"
      character (len = *), parameter :: dke_ape_name  = "dke_ape"
      character (len = *), parameter :: dke_cor_name  = "dke_cor"
      
      character (len = *), parameter :: wind_en_name  = "wind_en"
      
      character (len = *), parameter :: dpo_adv_units  = "W/m**2"
      character (len = *), parameter :: dpo_hdif_units = "W/m**2"
      character (len = *), parameter :: dpo_vdif_units = "W/m**2"
      character (len = *), parameter :: dpo_imp_units  = "W/m**2"
      character (len = *), parameter :: dpo_ice_units  = "W/m**2"
      character (len = *), parameter :: dpo_con_units  = "W/m**2"
      character (len = *), parameter :: dpo_bc_units   = "W/m**2"
      character (len = *), parameter :: dpo_bar_units  = "W/m**2"
      
      character (len = *), parameter :: din_adv_units  = "W/m**2"
      character (len = *), parameter :: din_hdif_units = "W/m**2"
      character (len = *), parameter :: din_vdif_units = "W/m**2"
      character (len = *), parameter :: din_imp_units  = "W/m**2"
      character (len = *), parameter :: din_ice_units  = "W/m**2"
      character (len = *), parameter :: din_ast_units  = "W/m**2"
      character (len = *), parameter :: din_con_units  = "W/m**2"
      character (len = *), parameter :: din_bc_units   = "W/m**2"
      character (len = *), parameter :: din_bar_units  = "W/m**2"
      
      character (len = *), parameter :: dke_adv_units  = "W/m**2"
      character (len = *), parameter :: dke_fri_units  = "W/m**2"
      character (len = *), parameter :: dke_pre_units  = "W/m**2"
      character (len = *), parameter :: dke_bar_units  = "W/m**2"
      character (len = *), parameter :: dke_bcf_units  = "W/m**2"
      character (len = *), parameter :: dke_ape_units  = "W/m**2"
      character (len = *), parameter :: dke_cor_units  = "W/m**2"
      
      character (len = *), parameter :: wind_en_units  = "W/m**2"
      
      character (len = *), parameter :: dpo_adv_longname  = "GPE change due to advection"
      character (len = *), parameter :: dpo_hdif_longname = "GPE change due to horizontal mixing"
      character (len = *), parameter :: dpo_vdif_longname = "GPE change due to vertical mixing"
      character (len = *), parameter :: dpo_imp_longname  = "GPE change due to implicit diffusion"
      character (len = *), parameter :: dpo_ice_longname  = "GPE change due to sea ice model"
      character (len = *), parameter :: dpo_con_longname  = "GPE change due to convective adjustment"
      character (len = *), parameter :: dpo_bc_longname   = "GPE change due to T-S force"
      character (len = *), parameter :: dpo_bar_longname  = "GPE change due to barotropy mass advection"
      
      character (len = *), parameter :: din_adv_longname  = "internal energy change due to advection"
      character (len = *), parameter :: din_hdif_longname = "internal energy change due to horizontal mixing"
      character (len = *), parameter :: din_vdif_longname = "internal energy change due to vertical mixing"
      character (len = *), parameter :: din_imp_longname  = "internal energy change due to implicit diffusion"
      character (len = *), parameter :: din_ice_longname  = "internal energy change due to sea ice model"
      character (len = *), parameter :: din_ast_longname  = "internal energy change due to Robert-Asselin time filter"
      character (len = *), parameter :: din_con_longname  = "internal energy change due to convective adjustment"
      character (len = *), parameter :: din_bc_longname   = "internal energy change due to T-S force"
      character (len = *), parameter :: din_bar_longname  = "internal energy change due to barotropy mass advection"
      
      character (len = *), parameter :: dke_adv_longname  = "kinetic energy change due to advection"
      character (len = *), parameter :: dke_fri_longname  = "kinetic energy change due to friction"
      character (len = *), parameter :: dke_pre_longname  = "kinetic energy change due to pressure gradient"
      character (len = *), parameter :: dke_bar_longname  = "kinetic energy change due to barotropy mass advection"
      character (len = *), parameter :: dke_bcf_longname  = "kinetic energy change due to wind stress"
      character (len = *), parameter :: dke_ape_longname  = "kinetic energy change due to atmospheric pressure"
      character (len = *), parameter :: dke_cor_longname  = "kinetic energy change due to Coriolis adjust"
      
      character (len = *), parameter :: wind_en_longname  = "wind force input energy"
      
      character (len = *), parameter :: bcfname = "pcom_bcf.nc"
      character (len = *), parameter :: bcuname = "bcu"
      character (len = *), parameter :: bcvname = "bcv"
      character (len = *), parameter :: bctname = "bct"
      character (len = *), parameter :: bcpname = "bcp"
      character (len = *), parameter :: bcsname = "bcs"
      character (len = *), parameter :: empname = "emp"
      character (len = *), parameter :: dddname = "ddd"
      
      character (len = *), parameter :: vmncname = "vmix.nc"
      character (len = *), parameter :: vmarname = "vmix"
      