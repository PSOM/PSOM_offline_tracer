subroutine write_cdf_3D(stepl,n)

!----------------------------------------------------------
! 3D output routine.
! It writes both center and face values to be written in .cdf files.
!----------------------------------------------------------

  USE header, only : NI,NJ,NK,xc,yc,zc,consump,Tr,nconsume,dirout,rc_kind,ntr,UL,WL,EPS
#include "netcdf.inc"

  integer i,j,k,n,stepl,nstp,it

  character(LEN=150) outname_data, outname_face

  REAL(kind=rc_kind) ::  Trwrite(0:NI+1,0:NJ+1,0:NK+1,ntr)
  REAL(kind=rc_kind) :: z(0:NI+1,0:NJ+1,0:NK+1)

  integer start(3), count(3), start2d(2), count2d(2),count2dtr(3),count4(4),start4(4),start2dtr(3),count4consump(4)
  integer countuf(3), countvf(3), countwf(3)

  integer dims(3), dims2d(2),dims2dtr(3),dims4(4),dimsconsump(4)
  integer dimuf(3), dimvf(3), dimwf(3)

  integer :: iddatfile, idigit, idudx, idudy, idudz, idvbysection, idvby, idvbz, idvcon, idvcy, idvbx
  integer :: idvcz, idvc, idvdivfeddy, idvdivfreyn, idvdx, idvdz, idvd, idvfb, idvh, idvn2bar
  integer :: idvn2, idvnsq100m, idvnsq30m, idvpe,idvpsiv,idvpsiw,idvpv, idvp,idvrhbar,idvrho,idvrnk
  integer :: idvstrain,idvstress,idvstr,idvs,idvtbar,idvtemp,idvtim,idvtr,idvt,idvu,idvvb,idvvc,idvvor,idvv
  integer :: idvwb,idvwc,idvwpv,idvw,idvy,idvzsave,idvz,idvz3,idwdx,idwdy,idwdz,iimday,ipos
  integer :: idvx,idvcon100,idFaceFile,idvuf, idvvf, idvwf, idvKf
  REAL(kind=rc_kind) ::  rcode

  DATA start /1, 1, 1/
  DATA start4 /1, 1, 1, 1/
  DATA start2d /1, 1/
  DATA start2dtr /1, 1, 1/

  count(1)= NI+2;
  count(2)= NJ+2;
  count(3)= NK+2;

  count4(1)= NI+2;
  count4(2)= NJ+2;
  count4(3)= NK+2;
  count4(4)= ntr;

  count4consump(1)= NI+2;
  count4consump(2)= NJ+2;
  count4consump(3)= NK+2;
  count4consump(4)= nconsume;

  do k=0,NK+1
    do j=0,NJ+1
      do i=0,NI+1
        z(i,j,k)= zc(i,j,k)
        do it=1,ntr
          Trwrite(i,j,k,it)= Tr(it,i,j,k,n)
        end do
      end do
    end do
  end do


  ! For the sake of better plots
  do j=0,NJ+1
    do i=0,NI+1
      z(i,j,NK+1)= 0.d0 ;   z(i,j,0)= 0.5*(z(i,j,0)+z(i,j,1));
      do it=1,ntr
        Trwrite(i,j,NK+1,it)= Trwrite(i,j,NK,it); Trwrite(i,j,0,it)= Trwrite(i,j,1,it);
      end do
    end do
  end do

  ! ---- END OF THE INITIALIZATION PART -----------------------------------------------

  ! --- START WRITING

  ! Output file names
  WRITE(outname_data,'("adv_25steps_full_",I5.5,".cdf")') stepl     ! Cell centers
  WRITE(outname_face,'("adv_face_",I5.5,".cdf")') stepl     ! Cell faces

  !---------------------------------------
  ! Write values at the cell centers:
  idDatFile =  nccre(TRIM(dirout)//outname_data,NCCLOB,rcode)

  ! 3D dimensions for most variables
  dims(1) = ncddef(idDatFile,'x',NI+2,rcode);
  dims(2) = ncddef(idDatFile,'y',NJ+2,rcode);
  dims(3) = ncddef(idDatFile,'sigma',NK+2,rcode);

  ! 4D dimensions (for tracer)
  dims4(1) = dims(1);
  dims4(2) = dims(2);
  dims4(3) = dims(3);
  dims4(4) = ncddef(idDatFile,'ntr',ntr,rcode);

  ! 3D dimensions (for consump)
  dimsconsump(1) = dims(1);
  dimsconsump(2) = dims(2);
  dimsconsump(3) = dims(3);
  dimsconsump(4) = ncddef(idDatFile,'ntrcon',nconsume,rcode);

  idvx = ncvdef(idDatFile,'xc',NCDOUBLE,1,dims(1),rcode)
  idvy = ncvdef(idDatFile,'yc',NCDOUBLE,1,dims(2),rcode)
  idvz = ncvdef(idDatFile,'zc',NCDOUBLE,3,dims,rcode)
  !idvc = ncvdef(idDatFile,'consump',NCDOUBLE,4,dimsconsump,rcode)
  idvtr = ncvdef(idDatFile,'tr',NCDOUBLE,4,dims4,rcode)

  CALL ncendf(idDatFile,rcode)
  CALL ncvpt(idDatFile,idvx, start(1), count(1), xc, rcode)
  CALL ncvpt(idDatFile,idvy, start(2), count(2), yc, rcode)
  CALL ncvpt(idDatFile,idvz, start, count, z, rcode)
  !CALL ncvpt(idDatFile,idvc,start4,count4consump,consump,rcode)
  CALL ncvpt(idDatFile,idvtr, start4, count4, Trwrite, rcode)

  CALL ncclos(idDatFile,rcode)

  ! End of writing of cell centers values.
  !---------------------------------------

  return
end
