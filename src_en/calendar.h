!======================= include file "calendar.h"======================
!
!                       calendar specification arrays
!
!-----------------------------------------------------------------------
!
!     monname = character names of months
!     daypm   = array of month lengths in days
!     daymd   = array of the middle date in monthes
!     month   = cumulative monthes of the integration
!     year    = the year of model
!     mth     = calendar month (from 1 to 12)
!     day     = calendar day (from 1 to daypm(mth))
!
!-----------------------------------------------------------------------
      integer   daypm(12),daymd(12)
      integer   month,year,mth,day
      character monname(12)*3
!

      data monname /'jan','feb','mar','apr','may','jun', &
                    'jul','aug','sep','oct','nov','dec'/
      data daymd   /15,15,15,15,15,15,15,15,15,15,15,15/
      data daypm   /30,30,30,30,30,30,30,30,30,30,30,30/
!      common /calend/ monname
!      common /calend/ daypm(12), daymd(12)
!      common /calend/ month,year,mth,day
