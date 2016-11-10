      parameter (imt=32,jmt=32)
      dimension itn(imt,jmt)

       do j=1,jmt
      do i=1,imt
            itn(i,j)=0
       enddo
       enddo

      do j=2,jmt-1
       do i=2,imt-1
            itn(i,j)=30
       enddo
       enddo


      open(72,file='topog.data',form='unformatted',status='new')
      write(72) itn
      close(72)
      stop
      end
