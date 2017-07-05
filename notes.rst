
**********************
 potential bugs in v2
**********************

op_vint, op_vint_ns, op_vint_ew, divided prh, but prh can be zero.

also, they don't deal with missing values.

*****************
 Changes to v1.0
*****************

density function
================

p0 not add 1-atm pressure, according to
\citet{Millero1981}, \citet[P.20]{Fofonoff1983}.

The check value in \citet{Jackett1995} is 
rho (3.0 deg., 35.5 psu, 300e5 Pa) = 1041.83267 kg / m^3

If add 1-atm pressure, the calculate result will be 
0.04 kg/m^3 greater than the check value. And if not add
1-atm pressue, the difference between the calculated value
and the check value is within 3.6e-7.

*********
 Caveats
*********

The intent(out) attribute problem
=================================

the grid variables in pcom contain an allocatable array, 
if a function or subroutine has a dumb parameter of type
grid variable, and the parameter attrubute is set to
intent(out), which is a typical case in the reload of
assignment oporator, then it will cause problem. Because
maybe the compile create a new variable of type grid, but
have not allocate memory for it, thus the result is that the
parameter passed in seems to 'lose' it memory.

I'm not sure whether this is a bug of ifort compiler. An
turn around solution is always use intent(inout) instead of
intent(out).


******************
unsolved questions
******************

pbt in int_btrop
================

upb is set as the average value in the whole cycle of
barotropic integration at the end of that cycle, but why not
pbt the same? ( Note that in batrop.f90, upb and pbt is not
set to the average values at the end of barotropic cycle)

pressure gradient of barotropic mode
====================================

in int_pgra of mod_int.f90, why use the odd pbxn, pbxs,
pcxn, pcxs?

initial sea bottom geopotential height gradient
===============================================

In calculating graphib, phib is on grid 1, graphib need be
on grid 3, why difference quotient first then do grid
interpolation, why not do grid interpolation first then do
difference quotient?


***************************
 remainning potential bugs 
***************************

  The follow issues are bugs by my understanding, but I'm no
so sure, since the v1.0 code can produce reasonably good
simulations. So I list them here for future reference.

output fluctuation of ssh
==========================
why only output fluctuation? the steric effects is hidden
then.

odd calculation of density calculating sea surface height
=========================================================
why use pbt%bc instead of pbt%tc? why t and s of previous
time step?
acw lies on g3 is ad hoc
========================
as the title, see mod_arrays

mask out in output
==================

when output u,v variables, pcom mask out the land points in
mod_mympi, but T,s variables does not do so. Why u,v is
special?

output variables
================

all the output variables are output to the same grid
withouth any interpolation.

frictin operator
================
In calculating the cv2* term, again, why use the center
difference scheme instead of the normal scheme? just for no
need to interpolation?

advection operator
==================
in op_adv of mod_op, when calc. divergence , it doesn't
adopt the standard op_div, but quite oddly (first multiply
cosine factor to v, then interpolate to grid, then multiply
u)

Smagorinsky scheme
==================

in calculating partial difference, why use center scheme, 
not the usual sheme, then interporate to the original grid.

calculation of u (ump)
======================

in int_trop, when calc. u, it is divided by 
pbt%bc2, which is on the T-grid, why not interpolate on
U-grid first?

days of per month
=================

For comparison with v1.0, v2.0 is month day is set to 30 for all 12
month. In future, this should be set back to correct days
(modify pro_days_month in mod_pro.f90)

Corrios adjustment
==================

In the deviation of appendix D, the terms need to adjust
is the tendency term, but in int_btrop of mod_int, the terms
to adjust is the velocity terms. Also, the sign of the
numerator is opposite as in the code.
(in fact, the deviation of appendix D match barotr.f90
exactly well, but not barotr_st.f90 of v1.0, is this a but
of Zhang Yu's code?)

corilios force in the momentum buget
====================================

in int_btrop, the Coriolis force is calulated as
0.5*2*omega*sin(phi), but I think the 0.5 factor should be
deleted. In the Coriolis adjustment step, however, the 0.5
factor disappears. Also, according to the continum equation, f should
be f* (including a metric term), why omit it in the code?

vertical viscocity
==================

in fri_gr3d of mod_gvar, when calculating friction
contributed by vertical shear, why use p U/ p z instead of
p U / p p, this is not consistent with the pressure
coordinate model.

calc cv2
=========
use the direct tan function to calc cv2 and cv1, instead of
the indrect method as v1.0 does

the order of interpolate and multiplication
===========================================

in the op_adv in mod_op, when calculating the nolinear
terms, why first interpolate, then muliply, why not multiply
first, then interpolate? (uU, vU, etc)

vertically interpolate horizontall velocity
===========================================

in op_adv of mod_op (vter_r3d), when calculate verticall advection,
need to interpolate u, v from vg2 to vg1, but the code
simpley use average, in the case of unuiform vertical
resolution, this may be not correct. Do we have to change to
layer weighted interpolation ?

upwelling.f90
=============

v1.0 of upwelling.f90, only calculated w(2:km), not
calculate w(kmp1), is that mean to set w(kmp1( to zero?

this problem also remain in upwelling of mod_int in v2.0

average density
===============

in mod_con.f90, set rho to 1029, more approprate value is 1035

averge and differentiate
========================

It is common to average first before differentiate, but in
main.f90 when calculate graphib, it fifferentiate first,
then average to the desire grid, this is not consistent with
many other variables derivative.

missing values of input file
============================

The missing values of the input file (intial.nc and
forcing.nc) is not at and only at lands. This can be
demonstrote by ::

  !  where (spread(itn, 3, 12) == 0)
  where (frc%tau%x(1) == missing_float)

of inistat in main.f90. This problem may be cause by
hand-modify the topography by Zhangyu.

days of month
=============

in pro_days_month of mod_pro, the calculating is not
correct, just for consistent with pcom1.0 . Also, it didnot
check for leap year (has been comment out)

vertical coordinates
====================

the pressure coordinates in vg1 and vg2 is proportional to
z, and independent of time. Thus, the model is accutually
barotropic.

parameters
==========

a in mod_con, should be 6371e3

equation of state
=================

in den_rho of mod_den.f90, the pressure add 1.013 before
calculation, this is wrong.

pressure gradient force of barotropic mode
==========================================

in int_pgra, when calc. pcne, it integrate grapc at u grid,
but acturally, grapc is at grid 2 and grid 4, not grid 3 (u
grid). And similar situation happens with pbne, pbsw.

density in calculation atmospheric pressure gradient
====================================================

  In int_pgra of mod_int.f90, PCOM calculate grapa with
constant density rho_0, not the surface density rho(z=0),
this is not correct, as pointed out by Griffies(2004) P.60.

pressure calculate from s
=========================

in den_rrho of mod_den, the pressure for the current layer
is calculate from the normalized pressure s, but s is only
valid for the bottom, not the current layer. Why all layers
use the same s?

barotropic vsu, vsv
===================

the code remove vsu from pubt at the first step of
barotropic integration, but add it again to pubt at every
barotropic step, this seems unneccessary, see barotropic
subroutine in main.f90

initial bottom pressure
=======================

the initial bottom pressure ph is calculated using
the mean t and s for each layer. This is equivalence to the
situation that no horizontal gridents of salinity and
temperture exists in the initial field. However, this is not
the case. And remember that ph is used later for dianostic
the ture pressure of each layer, and the normalized sea
bottom pressure is set to 1, which means that initial bottom
pressure can effect the dianostic pressure, and not just as
a reference role. I think using 3d T and S is needed to
calculate ph.
