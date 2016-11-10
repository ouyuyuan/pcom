c====================== include file "cvbc.h" ==========================
c
c     vertical boundary condition variables:
c
c     bcf = monthly mean surface forcing fileds. use linear
c           interpolation to produce:
c
c     bcu = bcf(1) :   sea surface zonal windstres       (dynes/cm**2)
c     bcv = bcf(2) :   sea surface meridional windstres  (dynes/cm**2)
c     bct = bcf(3) :   sea surface air temperature       (celsius)
c     bcp = bcf(4) :   sea surface air presure           (dynes/cm**2)
c     bcs = bcf(5) :   sea surface salinity              (model unit)
c     emp = bcf(6) :   rate of evaporation minus precipitation (cm/s)
c     ddd = bcf(7) :   coefficient for calculation of HF (w/m2/c)
c
      real*4 bcf,bcu,bcv,bct,bcp,bcs,emp,ddd
c
      common /cvbc/ bcf(imt,jmt,12,7)
      common /cvbc/ bcu(imt,jmt),bcv(imt,jmt)
      common /cvbc/ bct(imt,jmt),bcp(imt,jmt),bcs(imt,jmt)
      common /cvbc/ emp(imt,jmt),ddd(imt,jmt)
