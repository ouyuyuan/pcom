!====================== include file "cvbc.h" ==========================
!
!     vertical boundary condition variables:
!
!     bcf = monthly mean surface forcing fileds. use linear
!           interpolation to produce:
!
!     bcu = bcf(1) :   sea surface zonal windstres       (dynes/cm**2)
!     bcv = bcf(2) :   sea surface meridional windstres  (dynes/cm**2)
!     bct = bcf(3) :   sea surface air temperature       (celsius)
!     bcp = bcf(4) :   sea surface air presure           (dynes/cm**2)
!     bcs = bcf(5) :   sea surface salinity              (model unit)
!     emp = bcf(6) :   rate of evaporation minus precipitation (cm/s)
!     ddd = bcf(7) :   coefficient for calculation of HF (w/m2/c)
!
      real*4 bcf,bcu,bcv,bct,bcp,bcs,emp,ddd
!
      common /cvbc/ bcf(imt,jmt,12,7)
      common /cvbc/ bcu(imt,jmt),bcv(imt,jmt)
      common /cvbc/ bct(imt,jmt),bcp(imt,jmt),bcs(imt,jmt)
      common /cvbc/ emp(imt,jmt),ddd(imt,jmt)
