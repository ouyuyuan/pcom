	real*8 t(30),s(30)
	do k=1,30
	t(k)=-1.8
	s(k)=35.0
	enddo
c
	open(82,file='ts30.data',form='formatted')
        write(82,*) t,s
	close(82)
c
	end
