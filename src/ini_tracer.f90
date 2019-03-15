      subroutine tracerinit(stepl)
!     --------------------
!     initializes tracer fields 
!     
      USE header
      integer  i,j,k,l,n,stepl
      double precision z
!     
!     Tracer 1 is a tracer in y

!     TRACERS
!     =========

!     Initialze nutrient
      do j=0,NJ+1
         do i=0,NI+1
            do k=NK,1,-1
				Tr(1,i,j,k,0) = yc(j)
            end do
         end do
      end do


!     Periodic Boundaries (periodic e-w)
      do n=0,1
         do k=0,NK+1
            do j=0,NJ+1
               do it=1,ntr
                  Tr(it,0,j,k,n)= Tr(it,NI,j,k,n)
                  Tr(it,NI+1,j,k,n)= Tr(it,1,j,k,n)
               end do
            end do
         end do
      end do
      return
      end

