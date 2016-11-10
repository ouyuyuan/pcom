c====================== include file "param.h" =========================
c
c     main parameter file which sets ocean characteristics:
c
c                                                                       
c     imt = number of grid points in longitudinal direction          
c     jmt = number of grid points in latitudinal  direction           
c     km  = number of the sigma levels
c     nt  = number of tracers (temperature, salinity, ...)
c                                                                       
      integer imt,jmt,km,nt,imm,jmm,kmp1,kmm1,i,j,k,n,m
c
      parameter (imt=32,jmt=32,km=30)
      parameter (nt=2)
      parameter (imm=imt-1,jmm=jmt-1,kmp1=km+1,kmm1=km-1)
