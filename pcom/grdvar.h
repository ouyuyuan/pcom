!====================== include file "grdvar.h" ========================
!
!     variables which are functions of the grid
!
!     z0    = depth at the center of "t" and "u" grid cells
!     dz0   = layer-thickness in z-coordinate
!     z     = i.e. "eat",  at the center of "t" and "u" grid cells
!     dz    = thickness of "t" and "u" grid cells
!     rdz   = reciprocal of dz
!     rdzw  = reciprocal of thickness of "w" cell (in cm)
!     phib  = potential high at bottom
!     pn    = constant bottom pressure (PCOM); (-1)*phib (BCOM)
!     zu    = int(dz) of "u" grid cell
!     rzu   = reciprocal of zu
!
      real z0,dz0,z,dz,rdz,rdzw,phib,pn,rzu,zu
      common /grdv/ z0(km),dz0(km),z(km),dz(km),rdz(km),rdzw(km)
      common /grdv/ phib(imt,jmt),pn(imt,jmt),rzu(imt,jmt),zu(imt,jmt)
!
!
!     tmask = tracer cell land/sea mask   = (0.0, 1.0) on (land, sea)
!     umask = velocity cell land/sea mask = (0.0, 1.0) on (land, sea)
      real tmask,umask
      integer itn,ivn
      common /grdv/ tmask(imt,jmt,km),umask(imt,jmt,km)
      common /grdv/ itn(imt,jmt),ivn(imt,jmt)
!
!
!     cost  = cosine of "t" grid point latitude
!     cosu  = cosine of "u,v" grid point latitude
!     ff    = 2*omega*sine(j)
!     rdxt  = reciprocal of longitudinal width of "t" grid box (in cm)
!     rdxu  = reciprocal of longitudinal width of "u" grid box (in cm)
!     rdy   = reciprocal of latitudinal height
!     rdyt  = reciprocal of (latitudinal height * cost)
!     rdyu  = reciprocal of (latitudinal height * cosu)
!     dxdyt = dxt * dy
!     dxdyu = dxu * dy
!     area  = global ocean area
!     cv1   = square of sinu/cosu/a
!     cv2   = sinu/cosu * rdxu/radius
!     sdxt  = 1/2 square of rdxt
!     sdxu  = 1/2 square of rdxu
!     r1#   = coefficients for calculation of diffusion
!     ep#   = coefficients for semi-implicitly handle of Coriolis term
!
      real lat,cost,cosu,ff
      real rdxt,rdxu,rdyt,rdyu,rdy
      real dxdyt,dxdyu,area
      real cv1,cv2,sdxt,sdxu,r1a,r1b,r1c,r1d
      real ebea,ebeb,ebla,eblb,epea,epeb,epla,eplb
!
      common /grdv/ lat(jmt),cost(jmt),cosu(jmt),ff(jmt)
      common /grdv/ rdxt(jmt),rdxu(jmt),rdyt(jmt),rdyu(jmt),rdy
      common /grdv/ dxdyt(jmt),dxdyu(jmt),area
      common /grdv/ cv1(jmt),cv2(jmt),sdxt(jmt),sdxu(jmt)
      common /grdv/ r1a(jmt),r1b(jmt),r1c(jmt),r1d(jmt)
      common /grdv/ ebea(jmt),ebeb(jmt),ebla(jmt),eblb(jmt)
      common /grdv/ epea(jmt),epeb(jmt),epla(jmt),eplb(jmt)
!
