subroutine tracerinit(stepl) 
  USE header
  integer :: i,j,k,stepl
!     initializes tracer fields                                         
!     TRACER 1                                                          
!     =========                                                         
    it=1 
	do j=0,NJ+1
		do i=0,NI+1
			do k=NK,1,-1
      			Tr(it,i,j,k,0) = -zc(i,j,k)
			end do
		end do
	end do
end subroutine tracerinit
