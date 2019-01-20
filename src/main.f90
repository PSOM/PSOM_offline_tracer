PROGRAM main
!----------------------------------------------------
#include "cppdefs.h"

USE header
implicit none
INTEGER :: nbegin
REAL(kind=rc_kind) ::  dtim
REAL(kind=rc_kind) :: hpd,hpdinv
REAL(kind=rc_kind) :: z,dnkm1,ep,epm1,dnkm1p,zpd
REAL(kind=rc_kind) ::   xdu(0:NI+1,0:NJ+1),ydu(0:NI+1,0:NJ+1),xdv(0:NI+1,0:NJ+1),ydv(0:NI+1,0:NJ+1)

#include "main_declarations.h"
#include "ini_param.h"

nbegin = pickup_step

! 1. Initialize the tracers
CALL ini_setup(pcorr)
pfac = 3.5d0
CALL read_cdf_velocities(nbegin)
CALL sigma           ! calculates metric terms in vertical for the moving part of the grid
CALL staticsigma     ! calculates metric terms in vertical for the fixed part of the grid
                     ! It is important to call sigma before staticsigma and ini_st          
                     ! sigma needs to be called in momentum after each time advancement

do i=0,NI
	do j=0,NJ
		xdu(i,j) = (xc(i+1)-xc(i))/LEN*1.d3
		ydv(i,j) = (yc(j+1)-yc(j))/LEN*1.d3
		J2d(i,j) = xdu(i,j)*ydv(i,j) ! fix this for completeness
		xdu(NI+1,j) = (xc(NI+1)-xc(NI))/LEN*1.d3
		ydv(NI+1,j) = (yc(j+1)-yc(j))/LEN*1.d3
		J2d(NI+1,j) = xdu(NI+1,j)*ydv(NI+1,j)
	end do
	xdu(i,NJ+1) = xdu(i,NJ)
	ydv(i,NJ+1) = ydv(i,NJ)
	J2d(i,NJ+1) = xdu(i,NJ+1)*ydv(i,NJ+1)
end do

do i=0,NI+1
	do j=0,NJ+1
		do k=0,NK
			Jac(i,j,k) = J2d(i,j)/(zc(i,j,k+1)-zc(i,j,k))
			Jacinv(i,j,k) = 1/Jac(i,j,k)
		end do
		Jac(i,j,0) = Jac(i,j,1)
		Jacinv(i,j,0) = 1/Jac(i,j,0)
	end do
end do
 
CALL tracerinit(0)    !initializes tracer
! 2. advection routine
do step = pickup_step,(pickup_step+nsteps)
    ! 2a. load the velocity fields and interpolate in time
    if (mod(step,out3d_int).eq.0) then
        nbegin = step
		CALL read_cdf_velocities(nbegin)
    endif
	
    ! 2b. advect the tracer using the velocity fields and calculate reaction term
    do ivb=1,3
        if(ivb==1) then; dtim=dtf/3.d0; ivs=0;ivf=1;endif;
        if(ivb==2) then; dtim=0.5d0*dtf;ivs=1;ivf=1;endif;
        if(ivb==3) then; dtim=dtf      ;ivs=1;ivf=0;endif;
        tsp = dtim*1.d05
        CALL advection_and_mixing(ivs,ivf,dtim,step)
		CALL tracersource(ivs,ivf,dtim)
    enddo ! ivb
    ! 3. Save the tracer in a netcdf file
    if (mod(step,out3d_int).eq.0) then
            CALL write_cdf_3D(step,0)
    endif
    print*,step
enddo ! steps

END PROGRAM main
