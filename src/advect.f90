subroutine advect(var,uvarx,is_sink) 
!     ---------------------------------------------                     
! MODIFIED TO ADD A SINKING or SETTLING VELOCITY TO TRACER 2     (AM. 5/7/2016)                   
! is_sink dtermines if the settling velocity is to be added.  This is added in advect.f90
USE header
!     convecn(m,n,dtime) means use level m and write s,T to level n     
!     computes Cx,Cy,Cz at the cell centers.                            
!     uflx,vflx,wflx,sflx,Tflx are the fluxes defined at                
!     cell faces (using QUICK).                                         
!                                                                       

INTEGER :: i,i2,j,k,n,m,step,is_sink
! is_sink =0 if tracer is not sinking, =1 if tracer to be acted on by sinking in advect.f90


REAL(kind=rc_kind), PARAMETER  :: Eighth=0.125d0,Half=0.5
REAL(kind=rc_kind), PARAMETER :: courinv=2.d0, w0_sink=1/86400. , D_sink= 1000.     !wi_sink is in m/s

REAL(kind=rc_kind) :: absvel,left,right,ctr,fc,upos,uneg 

REAL(kind=rc_kind) :: dtimel

REAL(kind=rc_kind), dimension(    0:NI+1,0:NJ+1, 0:NK+1) :: var      ! var=(s,T,u,v,w,Tr(it,:,:,:))                    
REAL(kind=rc_kind), dimension(    0:NI+1,0:NJ+1, 0:NK+1) :: uvarx    ! uvarx  is the source term, divergence of the advective fluxes                   
REAL(kind=rc_kind), dimension(    0:NI  ,0:NJ  , 0:NK  ) :: varflx   ! varflx is a temporary variable to compute fluxes at the faces
REAL(kind=rc_kind) :: dvarx(0:NI+1),dvary(0:NJ),dvarz(0:NK+1)        ! dvarX is d(var)/dX at the faces (linear difference).
REAL(kind=rc_kind) :: wf_sink(0:NK),wfface,w_sink,zdep,wi_sink       ! wf_sink is the face sinking velocity of particles

INTEGER OMP_GET_THREAD_NUM,OMP_GET_NUM_THREADS,CHUNK,NPROC
                       
  !------------------------------------------
  ! computation of the source term in x.

  wi_sink = -w0_sink/WL

  !$OMP PARALLEL DO PRIVATE(k,j,dvarx,left,ctr,fc,i2,right) SHARED(var,uf) COLLAPSE(2) 
  do k=1,NK 
    do j=1,NJ 
      !write(6,"(a,i3,a,i3,a,i2,a,i2)") " OpenMP:",k,"   ",j,"    ,N_threads = ",OMP_GET_NUM_THREADS()," thread = ", OMP_GET_THREAD_NUM() 
      dvarx(0)=var(1,j,k)-var(NI,j,k)
      do i=1,NI-1 
        dvarx(i)= var(i+1,j,k) - var(i,j,k) 
      end do 
      dvarx(NI)=dvarx(0)
      dvarx(NI+1)=dvarx(1)

      var(NI+1,j,k)= var(1,j,k) 
      var(0,j,k)= var(NI,j,k) 
       
      do i=1,NI
        if (uf(i,j,k).gt.0.d0) then 
          left= Eighth*(dvarx(i) -dvarx(i-1)) 
          ctr= Half* (var(i+1,j,k) +var(i,j,k)) 
          fc= ctr -left 
          call ultim(var(i-1,j,k),var(i,j,k),var(i+1,j,k),dvarx(i),dvarx(i-1),courinv,fc)         
          varflx(i,j,k)= uf(i,j,k)*fc 
         else 
          if(i==NI) then                  
            i2=0;
           else
            i2=i;
          endif
          right= Eighth*(dvarx(i2+1) -dvarx(i2)) 
          ctr= Half* (var(i2+1,j,k) +var(i,j,k)) 
          fc= ctr - right 
          call ultim(var(i2+2,j,k),var(i2+1,j,k),var(i,j,k),dvarx(i),dvarx(i2+1),courinv,fc)           
          varflx(i,j,k)= uf(i,j,k)*fc 
        end if 
      end do
      varflx(0,j,k)= varflx(NI,j,k) 
    end do 
  end do 
  !$OMP END PARALLEL DO

  do i=1,NI 
    uvarx(i,1:NJ,1:NK)=varflx(i,1:NJ,1:NK) -varflx(i-1,1:NJ,1:NK) 
  end do 

  ! computation of the source term in x done.
  !------------------------------------------


                                                       
  !------------------------------------------
  ! computation of the source term in y.
                                                               
  !$OMP PARALLEL DO PRIVATE(k,i,dvary,left,ctr,fc,right) SHARED(var,vf) COLLAPSE(2) 
  do k=1,NK 
    do i=1,NI
 
      dvary(0)= 0.d0; 
      dvary(NJ)= 0.d0; 
      do j=1,NJ-1 
        dvary(j)= var(i,j+1,k) - var(i,j,k) 
      end do 

      var(i,0,k)= var(i,1,k) 
      var(i,NJ+1,k)= var(i,NJ,k) 

      do j=1,NJ-1 
        if (vf(i,j,k).ge.0.d0) then 
          left= Eighth*(dvary(j) -dvary(j-1)) 
          ctr= Half* (var(i,j+1,k) +var(i,j,k)) 
          fc= ctr -left
          call ultim(var(i,j-1,k),var(i,j,k),var(i,j+1,k),dvary(j),dvary(j-1),courinv,fc)         
          varflx(i,j,k)= vf(i,j,k)*fc 
        else 
         right= Eighth*(dvary(j+1) -dvary(j)) 
         ctr= Half* (var(i,j+1,k) +var(i,j,k)) 
         fc= ctr - right
         call ultim(var(i,j+2,k),var(i,j+1,k),var(i,j,k),dvary(j),dvary(j+1),courinv,fc)           
         varflx(i,j,k)= vf(i,j,k)*fc 
       end if 
     end do 
     varflx(i,0,k)=0.d0;
     varflx(i,NJ,k)=0.d0;
    end do 
  end do 
  !$OMP END PARALLEL DO
      

  do j=1,NJ 
    uvarx(1:NI,j,1:NK)= (varflx(1:NI,j,1:NK) -varflx(1:NI,j-1,1:NK) ) + uvarx(1:NI,j,1:NK) 
  end do 
            

  ! computation of the source term in y done.
  !------------------------------------------


  !------------------------------------------
  ! computation of the source term in z.
                                                           
  !$OMP PARALLEL DO PRIVATE(j,i,dvarz,left,ctr,fc,right) SHARED(var,wf) COLLAPSE(2)
  do j=1,NJ
    do i=1,NI 

      var(i,j,0)= var(i,j,1) 

      dvarz(0)= 0.d0 
      do k=1,NK-1 
         dvarz(k)= var(i,j,k+1) - var(i,j,k) 
      end do 
      dvarz(NK:NK+1)= 0.d0  !     Top boundary - Zero gradient                                      

      var(i,j,NK+1)= var(i,j,NK) 

      if (is_sink.eq.1) then
         wf_sink(0)= 0.
         wf_sink(NK)= 0.
         do k=1,NK-1
            ! w sinking velocity w_sink= wi_sink (1-z/D_sink)**0.67
            ! This relationship is from our sinking flux paper Omand et al. wi_sink= max sinking vel, Dsink is where w_snink= 0.
            ! wi_sink and D_sink must be defined paramters.
            zdep= -zf(i,j,k)*DL/D_sink   ! zdep is positive, and it's max value should be 1
            w_sink= wi_sink*(1. -  zdep)              ! w_sink varies between wi_sink at the surface, and 0 below D_sink.
            wf_sink(k)= 0.5*( Jac(i,j,k+1)*wz(i,j,k+1)+ Jac(i,j,k)*wz(i,j,k) )*w_sink            ! convert to face fluxes. 
            !PRINT *, w_sink, wf_sink(k), zdep
         end do
      else
         wf_sink(0:NK)= 0.d0
      end if

      do k=1,NK 
         wfface= wf(i,j,k)+wf_sink(k)
         !if(i.eq.NI/2.AND.j.eq.NJ/2.AND.is_sink.eq.1) PRINT *, wfface, wf(i,j,k), wf_sink(k)
        if (wfface.ge.0.d0) then 
          left= Eighth*(dvarz(k) -dvarz(k-1)) 
          ctr= Half* (var(i,j,k+1) +var(i,j,k)) 
          fc= ctr -left 
          call ultim(var(i,j,k-1),var(i,j,k),var(i,j,k+1),dvarz(k),dvarz(k-1),courinv,fc)         
          varflx(i,j,k)= (wf(i,j,k)+wf_sink(k))*fc
         else 
          if (k.eq.NK) then 
            varflx(i,j,k)= (wf(i,j,k)+wf_sink(k))*Half*(var(i,j,k+1) +var(i,j,k))
           else
            right= Eighth*(dvarz(k+1) -dvarz(k)) 
            ctr= Half* (var(i,j,k+1) +var(i,j,k)) 
            fc= ctr - right 
            call ultim(var(i,j,k+2),var(i,j,k+1),var(i,j,k),dvarz(k),dvarz(k+1),courinv,fc)           
            varflx(i,j,k)= (wf(i,j,k)+wf_sink(k))*fc
          end if 
        end if 
      end do 
      varflx(i,j,0)= 0.d0   !     Solid Boundary                                                    
    end do 
  end do 
  !$OMP END PARALLEL DO

  do k=1,NK 
    uvarx(1:NI,1:NJ,k)=  ( varflx(1:NI,1:NJ,k) -varflx(1:NI,1:NJ,k-1) ) + uvarx(1:NI,1:NJ,k) 
  enddo 

  ! computation of the source term in z done.
  !------------------------------------------

return 
                                                                        
END                                           
      
!----------------------------------------------------------------------
! SUBROUTINE ultim                                                                  
!     Implementation of the ULTIMATE limiter based on upstrm,dnstrm,ctr 
!     points and d(phi)/dx at the up- and dnstrm locations              
!     Returns the "ulitimate limited" face value of the variable        
!     that should be multiplied by uf,vf or wf to give the flux.        
                                                                         
SUBROUTINE ultim(up,ct,dn,gdn,gup,courinv,fc) 
USE header, only : rc_kind
                                                                  
REAL(kind=rc_kind) :: up,ct,dn,gup,gdn,courinv,fc
REAL(kind=rc_kind) :: ref,del,adel,acurv 
                                                                  
del= dn -up 
adel= abs(del) 
acurv= abs(gdn -gup) 
if (acurv.ge.adel)  then 
  fc= ct 
  return 
 else 
  ref= up + (ct -up)*courinv 
  if (del.gt.0) then 
    fc= max(ct,fc) 
    ref= min(ref,dn) 
    fc= min(fc,ref) 
  else 
    fc= min(ct,fc) 
    ref= max(ref,dn) 
    fc= max(ref,fc) 
  endif 
end if 
return 
END                                           
