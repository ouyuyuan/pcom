c====================== include file "diag.h" ==========================
c
c     tmn    = monthly/annual averaged potential temperature (c)
c     smn    = monthly/annual averaged salinity            (model unit)
c     pmn    = monthly/annual averaged pbt                 (dynes/cm**2)
c     umn    = monthly/annual averaged up
c     vmn    = monthly/annual averaged vp
c
      real tmn,smn,pmn,umn,vmn
c
      common /cdiag/ tmn(imt,jmt,km),smn(imt,jmt,km)
      common /cdiag/ umn(imt,jmt,km),vmn(imt,jmt,km)
      common /cdiag/ pmn(imt,jmt)
