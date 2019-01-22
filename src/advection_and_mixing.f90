subroutine advection_and_mixing(m,n,dtimel,step) 
!     ---------------------------------------------                     
USE header

INTEGER :: i,j,k
INTEGER :: n,m,step,is_sink

REAL(kind=rc_kind) :: dtimel

REAL(kind=rc_kind), dimension(    0:NI+1,0:NJ+1, 0:NK+1) :: var             ! var=(s,T,u,v,w,Tr(it,:,:,:)
REAL(kind=rc_kind), dimension(    0:NI+1,0:NJ+1, 0:NK+1) :: uvarx    ! uvarx  is the source term, divergence of the advective fluxes
REAL(kind=rc_kind), dimension(    1:NI  ,1:NJ  , 1:NK  ) :: vardif                ! vardif is the source term from diabatic processes
INTEGER :: iv_compute_kz

iv_compute_kz=1

do it=1,ntr
    ! Flux terms initialized to 0
    uvarx = 0.d0; vardif = 0.d0
    var=Tr(it,:,:,:,m)
    ! ---------------------------------------------------------------
    ! computation of the advective fluxes, using QUICK scheme
    is_sink = 0 ! no sinking
    CALL advect(var,uvarx,is_sink)
    ! ---------------------------------------------------------------
    ! computation of the horizontal diabatic fluxes
    vardif=0.;
    !call mixing_horizontal(var,vardif);
    call mixing_isopycnal(var,vardif,1.);!PRINT*,"VISCOUS REDI";
    uvarx(1:NI,1:NJ,1:NK)=uvarx(1:NI,1:NJ,1:NK)-vardif(1:NI,1:NJ,1:NK)
    ! ---------------------------------------------------------------
    ! computation of the vertical diabatic fluxes
    vardif=0.;
    !call mixing_vertical(var,vardif,m,step,iv_compute_kz);
    uvarx(1:NI,1:NJ,1:NK)=uvarx(1:NI,1:NJ,1:NK)-vardif(1:NI,1:NJ,1:NK)
    iv_compute_kz=0;
    ! ---------------------------------------------------------------
    ! final summation
    Tr(it,1:NI,1:NJ,1:NK,n) = Tr(it,1:NI,1:NJ,1:NK,0)-dtimel*Jacinv(1:NI,1:NJ,1:NK)*uvarx(1:NI,1:NJ,1:NK)
enddo ! it loop

return 
                                                                        
END
