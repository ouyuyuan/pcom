      parameter (imt=182,jmt=78)
      real*4 bcf(imt,jmt,12,7)

       do 888 j=1,jmt

	if(j.le.jmt/2)then
		t=-1.8
	else
		t=-20
	endif
	
       do i=1,imt
       do k=1,12
            bcf(i,j,k,1) = 0.0
            bcf(i,j,k,2) = 0.0
            bcf(i,j,k,3) = t
            bcf(i,j,k,4) = 0.0
            bcf(i,j,k,5) = 35.0
            bcf(i,j,k,6) = 0.0
            bcf(i,j,k,7) = 40.0
       enddo
       enddo
888    continue

      open(72,file='sbcf.data',form='unformatted',status='new')
      write(72) bcf
      close(72)
      stop
      end
