c====================== include file "grdvar.h" ========================
c
c     variables which are functions of the grid
c
c     z0    = depth at the center of "t" and "u" grid cells
c     dz0   = layer-thickness in z-coordinate
c     z     = i.e. "eat",  at the center of "t" and "u" grid cells
c     dz    = thickness of "t" and "u" grid cells
c     rdz   = reciprocal of dz
c     rdzw  = reciprocal of thickness of "w" cell (in cm)
c     phib  = potential high at bottom
c     pn    = constant bottom pressure (PCOM); (-1)*phib (BCOM)
c     zu    = int(dz) of "u" grid cell
c     rzu   = reciprocal of zu
c
      real z0,dz0,z,dz,rdz,rdzw,phib,pn,rzu,zu
      common /grdv/ z0(km),dz0(km),z(km),dz(km),rdz(km),rdzw(km)
      common /grdv/ phib(imt,jmt),pn(imt,jmt),rzu(imt,jmt),zu(imt,jmt)
c
c
c     tmask = tracer cell land/sea mask   = (0.0, 1.0) on (land, sea)
c     umask = velocity cell land/sea mask = (0.0, 1.0) on (land, sea)
      real tmask,umask
      integer itn,ivn
      common /grdv/ tmask(imt,jmt,km),umask(imt,jmt,km)
      common /grdv/ itn(imt,jmt),ivn(imt,jmt)
c
c
c     cost  = cosine of "t" grid point latitude
c     cosu  = cosine of "u,v" grid point latitude
c     ff    = 2*omega*sine(j)
c     rdxt  = reciprocal of longitudinal width of "t" grid box (in cm)
c     rdxu  = reciprocal of longitudinal width of "u" grid box (in cm)
c     rdy   = reciprocal of latitudinal height
c     rdyt  = reciprocal of (latitudinal height * cost)
c     rdyu  = reciprocal of (latitudinal height * cosu)
c     dxdyt = dxt * dy
c     dxdyu = dxu * dy
c     area  = global ocean area
c     cv1   = square of sinu/cosu/a
c     cv2   = sinu/cosu * rdxu/radius
c     sdxt  = 1/2 square of rdxt
c     sdxu  = 1/2 square of rdxu
c     r1#   = coefficients for calculation of diffusion
c     ep#   = coefficients for semi-implicitly handle of Coriolis term
c
      real lat,cost,cosu,ff
      real rdxt,rdxu,rdyt,rdyu,rdy
      real dxdyt,dxdyu,area
      real cv1,cv2,sdxt,sdxu,r1a,r1b,r1c,r1d
      real ebea,ebeb,ebla,eblb,epea,epeb,epla,eplb
c
      common /grdv/ lat(jmt),cost(jmt),cosu(jmt),ff(jmt)
      common /grdv/ rdxt(jmt),rdxu(jmt),rdyt(jmt),rdyu(jmt),rdy
      common /grdv/ dxdyt(jmt),dxdyu(jmt),area
      common /grdv/ cv1(jmt),cv2(jmt),sdxt(jmt),sdxu(jmt)
      common /grdv/ r1a(jmt),r1b(jmt),r1c(jmt),r1d(jmt)
      common /grdv/ ebea(jmt),ebeb(jmt),ebla(jmt),eblb(jmt)
      common /grdv/ epea(jmt),epeb(jmt),epla(jmt),eplb(jmt)
c
