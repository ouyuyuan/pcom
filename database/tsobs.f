      parameter (imt=182,jmt=78,km=30)
      real*8 t(imt,jmt,km),s(imt,jmt,km)

       do j=1,jmt
      do i=1,imt
	do k=1,km
            t(i,j,k)=-1.8
            s(i,j,k)=35.0
       enddo
       enddo
       enddo

      open(72,file='tsobs.data',form='unformatted')
      write(72) t,s
      close(72)
      stop
      end
