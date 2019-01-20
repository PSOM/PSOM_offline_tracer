PROGRAM main
!----------------------------------------------------
#include "cppdefs.h"

USE header
implicit none
INTEGER :: nbegin
REAL(kind=rc_kind) ::  dtim, dstep

#include "main_declarations.h"
#include "ini_param.h"

nbegin = pickup_step

! 1. Initialize the variables and tracers
CALL ini_setup(pcorr)
CALL sigma           ! calculates metric terms in vertical for the moving part of the grid
CALL staticsigma(nbegin)     ! calculates metric terms in vertical for the fixed part of the grid
                     ! It is important to call sigma before staticsigma and ini_st          
                     ! sigma needs to be called in momentum after each time advancement
CALL tracerinit(0)    !initializes tracer
! 2. advection routine
do step = pickup_step,(pickup_step+nsteps)
    ! 2a. load the velocity fields and interpolate in time
	dstep = mod(step,out3d_int)/out3d_int
    if (dstep.eq.0) then
        nbegin = step
		CALL read_cdf_velocities(nbegin)
    endif
	uf = uf0*(1-dstep)+uf1*dstep
	vf = vf0*(1-dstep)+vf1*dstep
	wf = wf0*(1-dstep)+wf1*dstep
	h = h0*(1-dstep)+h1*dstep
	CALL sigma           ! calculates metric terms in vertical for the moving part of the grid
	
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
    if (mod(step,100).eq.0) then
            CALL write_cdf_3D(step,0)
    endif
    print*,step
enddo ! steps

END PROGRAM main
