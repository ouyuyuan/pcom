      parameter (imt=182,jmt=78)
      dimension itn(imt,jmt)

       do j=1,jmt
      do i=1,imt
            itn(i,j)=0
       enddo
       enddo

      do j=2,jmt-1
       do i=1,imt
            itn(i,j)=30		!km=30
       enddo
       enddo

      open(72,file='topog.data',form='unformatted')
      write(72) itn
      close(72)
      stop
      end
