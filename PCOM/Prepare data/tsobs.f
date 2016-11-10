      parameter (imt=32,jmt=32,km=30)
      real*8 t(imt,jmt,km),s(imt,jmt,km)

	do k=1,km
       do j=1,jmt
      do i=1,imt
            t(i,j,k)=10.0
            s(i,j,k)=35.0
       enddo
       enddo
       enddo

 
      open(72,file='tsobs.data',form='unformatted',status='new')
      write(72) t,s
      close(72)
      stop
      end
