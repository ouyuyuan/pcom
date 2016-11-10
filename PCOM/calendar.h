c======================= include file "calendar.h"======================
c
c                       calendar specification arrays
c
c-----------------------------------------------------------------------
c
c     monname = character names of months
c     daypm   = array of month lengths in days
c     daymd   = array of the middle date in monthes
c     month   = cumulative monthes of the integration
c     year    = the year of model
c     mth     = calendar month (from 1 to 12)
c     day     = calendar day (from 1 to daypm(mth))
c
c-----------------------------------------------------------------------
      integer   daypm,daymd
      integer   month,year,mth,day
      character monname(12)*3
c
      common /calend/ monname
      common /calend/ daypm(12), daymd(12)
      common /calend/ month,year,mth,day
