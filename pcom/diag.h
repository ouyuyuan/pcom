!====================== include file "diag.h" ==========================
!
!     tmn    = monthly/annual averaged potential temperature (c)
!     smn    = monthly/annual averaged salinity            (model unit)
!     pmn    = monthly/annual averaged pbt                 (dynes/cm**2)
!     umn    = monthly/annual averaged up
!     vmn    = monthly/annual averaged vp
!
      real tmn,smn,pmn,umn,vmn
!
      common /cdiag/ tmn(imt,jmt,km),smn(imt,jmt,km)
      common /cdiag/ umn(imt,jmt,km),vmn(imt,jmt,km)
      common /cdiag/ pmn(imt,jmt)
