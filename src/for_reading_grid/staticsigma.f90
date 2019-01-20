subroutine staticsigma(nstp)

#include "cppdefs.h"
  USE header
  implicit none
  #include "netcdf.inc"
  !     ------------------------                                          
  !     FOR THE REGION BELOW THE TOPMOST LAYER                            
  !     modified a for periodicew bc                                      
  !     This subroutine updates all those quantities that are a function  
  !     of time and  (time and depth). It is called  every time step.     
  !     The surface bc  wfbct is updated here.                            
!  implicit REAL(kind=rc_kind) :: (a-h,o-z) 
  integer i,j,k 
  REAL(kind=rc_kind) :: hpd,hpdinv,hu,hv,hx,hy,hxpdx,hypdy,z,temp,       &
       be2,wxk,wyk,d2,dnk,dnkm1,sig,ep,epm1,dnkm1p,zpd,           &
       g13(0:NI+1,0:NJ+1,0:NK+1),g23(0:NI+1,0:NJ+1,0:NK+1)          
  !     Ddx and Ddy are known from init                                   
  !                                                                       
  !    use grid instead of grid parameters

  integer :: nstp
  integer :: idInFile,idzSliceFile
  integer :: idxc,idyc,idzc,idh
  REAL(kind=rc_kind) ::  rcode
  character (len = 550) :: zslice_data
  integer start(3), count(3), countuf(3), countvf(3), countwf(3)
  integer start2d(2), count2d(2)

  DATA start /1, 1, 1/
  DATA start2d /1, 1/
  DATA count2d /NI, NJ/

count(1)= NI+2
count(2)= NJ+2
count(3)= NK+2

WRITE(zslice_data,'("face_",I5.5,".cdf")') nstp
print *, TRIM(dirout)//zslice_data
idzSliceFile = ncopn(TRIM(dirout)//zslice_data, NCNOWRIT,rcode)

idzc = ncvid(idzSliceFile,'zf',rcode)
call ncvgt( idzSliceFile, idzc, start, count, zf, rcode )
call ncclos(idzSliceFile, rcode)


WRITE(zslice_data,'("full_",I5.5,".cdf")') nstp
print *, TRIM(dirout)//zslice_data
idzSliceFile = ncopn(TRIM(dirout)//zslice_data, NCNOWRIT,rcode)

idxc = ncvid(idzSliceFile,'xc',rcode)
idyc = ncvid(idzSliceFile,'yc',rcode)
idh = ncvid(idzSliceFile,'h',rcode)
call ncvgt( idzSliceFile, idxc, start(1), count(1), xc, rcode )
call ncvgt( idzSliceFile, idyc, start(2), count(2), yc, rcode )
call ncvgt( idzSliceFile, idh, start, count, h, rcode)
call ncclos(idzSliceFile, rcode)



  dnk= dble(NK) 

#ifdef fixed_bottom_thickness
  dnkm1= dble(NK-1-1) 
#else
  dnkm1= dble(NK-1) 
#endif

  !dnkm1p= dnkm1/pfac 
  be2= beta*EPS*EPS 
  !epm1= exp(pfac) -1.d0 
  !ep = pfac*exp(1.d0) 
  !                                                                       
  do j=0,NJ+1 
     do i=0,NI+1 
        !     All these variables are functions of time                         
            hpd= dztop +D(i,j) 
            hpdinv= 1.d0/hpd                                     
!                                                                       
            do 20 k=0,NK-1 
               z= zc(i,j,k) 
               zpd= z +dztop 
               wz(i,j,k)= 1/(zf(i,j,k+1)-zf(i,j,k))*DL
               Jac(i,j,k)= J2d(i,j)/wz(i,j,k)                  
               sig= dble(k)-0.5d0 

			   wx(i,j,k)= wz(i,j,k)*zpd*hpdinv*Ddx(i,j) 
               wy(i,j,k)= wz(i,j,k)*zpd*hpdinv*Ddy(i,j) 
               g13(i,j,k)= ux(i,j)*wx(i,j,k) +uy(i,j)*wy(i,j,k) 
               g23(i,j,k)= vx(i,j)*wx(i,j,k) +vy(i,j)*wy(i,j,k) 
!               g33(i,j,k)= wx(i,j,k)*wx(i,j,k) +wy(i,j,k)*wy(i,j,k) +  
!     &              be2*wz(i,j)*wz(i,j)                                
   20       continue 
            do k=0,NK-1 
               sig= dble(k) 
               wt(i,j,k)= 0.d0 
               z= zc(i,j,k) 
               zpd= z +dztop 
               wzk(i,j,k)= zc(i,j,k+1)-zc(i,j,k) 
            end do 
            wzk(i,j,NK-1)= 0.5*(wz(i,j,NK) + wz(i,j,NK-1)) 
      
            do 15 k=0,NK                
               z= zf(i,j,k) 
               zpd= z +dztop 
               sig= dble(k)                         
                                                                        
               temp= (epm1*dnkm1p/(epm1*zpd +hpd))*zpd*hpdinv 
               wxk= Ddx(i,j)*temp 
               wyk= Ddy(i,j)*temp 
               gqk(i,j,k,1)= qpr*Jac(i,j,k)*(ux(i,j)*wxk +uy(i,j)*wyk) 
               gqk(i,j,k,2)= qpr*Jac(i,j,k)*(vx(i,j)*wxk +vy(i,j)*wyk) 
               gqk(i,j,k,3)= Jac(i,j,k)*(qpr*(wxk*wxk +wyk*wyk) +       &
     &              be2*wz(i,j,k)*wz(i,j,k))                            
   15       continue 
      enddo
    enddo
                                                                  
      do 19 k=1,NK 
      do 21 i=0,NI 
         do 31 j=1,NJ 
            Jifc(i,j,k)= 0.5d0*(Jac(i,j,k)+ Jac(i+1,j,k)) 
            gi(i,j,k,1)= 0.5d0*(g11(i,j) +g11(i+1,j))*Jifc(i,j,k) 
            gi(i,j,k,2)= 0.5d0*(g12(i,j) +g12(i+1,j))*Jifc(i,j,k) 
            gqi(i,j,k,1)= qpr*gi(i,j,k,1) 
            gqi(i,j,k,2)= qpr*gi(i,j,k,2) 
   31    continue 
   21 continue 
   19 continue 
      do 28 k=1,NK 
      do 22 i=1,NI 
         do 32 j=0,NJ 
            Jjfc(i,j,k)= 0.5d0*(Jac(i,j,k)+ Jac(i,j+1,k)) 
            gj(i,j,k,1)= 0.5d0*(g12(i,j) +g12(i,j+1))*Jjfc(i,j,k) 
            gj(i,j,k,2)= 0.5d0*(g22(i,j) +g22(i,j+1))*Jjfc(i,j,k) 
            gqj(i,j,k,1)= qpr*gj(i,j,k,1) 
            gqj(i,j,k,2)= qpr*gj(i,j,k,2) 
                                                                        
   32    continue 
   22 continue 
   28 continue 
!                                                                       
      do j=1,NJ 
         do i=0,NI 
            do k=1,NK 
               gi3(i,j,k)= 0.5d0*(g13(i,j,k) +g13(i+1,j,k))*Jifc(i,j,k) 
               gqi3(i,j,k)= qpr*gi3(i,j,k) 
            enddo
         enddo
      enddo

      do j=0,NJ 
         do i=1,NI 
            do k=1,NK 
               gj3(i,j,k)= 0.5d0*(g23(i,j,k) +g23(i,j+1,k))*Jjfc(i,j,k) 
               gqj3(i,j,k)= qpr*gj3(i,j,k) 

            enddo
         enddo
      enddo

 do i=0,NI
   do j=0,NJ
     do k=0,NK
       Jacinv(i,j,k)=1./Jac(i,j,k)
     enddo
   enddo
 enddo

                                                    
       return 
      END                                           
