c====================== include file "pconst.h" ========================
c
c
c     rules for parameter constants
c
c     use prefix of "c" for whole real numbers (eg: c57 for 57.0)
c     use "m" after prefix to designate negative values (minus sign)
c       (eg: cm7 for -7.0)
c     use prefix of "p" for non repeating fractions (eg: p5 for 0.5)
c     use prefix of "r" for reciprocals (eg: r3 for 1/3.0)
c     combine use of prefix above and "e" for scientific notation, with
c       (eg: c5e4 for 5.0e4, c1em10 for 1.0e-10)
c
      real c0,c1,c2,c3,c4,c5,c6,c25,c35
      parameter (c0=0.0,c1=1.0,c2=2.0,c3=3.0,c4=4.0)
      parameter (c5=5.0,c6=6.0,c25=25.0,c35=35.0)
c
      real p125,p25,p5
      parameter (p125=0.125,p25=0.25, p5=0.5)
c
      real c60,c1440,c3600
      parameter (c60=60.0, c1440=1440.0, c3600=3600.0)
c
      real c1e3,c1em4,c1em12
      parameter (c1e3=1.0e3,c1em4=1.0e-4,c1em12=1.0e-12)
c
      real r120,r150,r180,r365
      parameter (r120=c1/120.0,r150=c1/150.0)
      parameter (r180=c1/180.0,r365=c1/365.0)
c
      real secday
      parameter (secday=c1/(c60*c1440))
c
      real rho_0,rrho_0
      parameter (rho_0=1.029,rrho_0=c1/rho_0)
c
      real tbice
      parameter (tbice=-1.5)
c
c
      real pi,torad,radius,omega,grav
      parameter (pi     = 3.141592653589793)
      parameter (torad  = pi*r180    )
      parameter (radius = 6370.0d5   )
      parameter (omega  = pi/43082.0 )
      parameter (grav   = 980.0)
c
c     grav = earth's gravitational acceleration (cm/sec**2)
c
