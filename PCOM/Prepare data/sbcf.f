
      parameter (imt=32,jmt=32)

      real*4 bcf(imt,jmt,12,7)


       do 888 j=1,jmt


	tlat = -1 + (j-1)*2
	
	t = 25 * (1 - tlat/60.0)

	print*, tlat, t

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
