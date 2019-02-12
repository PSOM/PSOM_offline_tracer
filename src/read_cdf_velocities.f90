subroutine read_cdf_velocities(nstp)
  !subroutine readcdfcn(Tr,uf,vf,wf,xc,yc,zc,xf,yf,zf,nbegin)
  !     ------------------------------------------------------
  !     reads in the netCDF files written by the routine outcdf.f

  !use header, only : NI,NJ,NK,ntr,rc_kind
  use header

  implicit none

#include "netcdf.inc"
  ! include 'dims.f90'

  integer :: nstp
  integer :: idInFile,idzSliceFile
  integer :: idxc,idyc,idzc,iduf,idvf,idwf,idh

  REAL(kind=rc_kind) ::  rcode

  character (len = 550) :: inname_data, zslice_data

  integer start(3), count(3), countuf(3), countvf(3), countwf(3)
  integer start2d(2), count2d(2)

  DATA start /1, 1, 1/
  DATA start2d /1, 1/
  DATA count2d /NI, NJ/
  
  count(1)= NI+2
  count(2)= NJ+2
  count(3)= NK+2

  countuf(1)= NI+1
  countuf(2)= NJ
  countuf(3)= NK

  countvf(1)= NI
  countvf(2)= NJ+1
  countvf(3)= NK

  countwf(1)= NI
  countwf(2)= NJ
  countwf(3)= NK+1

  WRITE(zslice_data,'("full_",I5.5,".cdf")') nstp
  print *, TRIM(dirout)//zslice_data
  idzSliceFile = ncopn(TRIM(dirout)//zslice_data, NCNOWRIT,rcode)

  idxc = ncvid(idzSliceFile,'xc',rcode)
  idyc = ncvid(idzSliceFile,'yc',rcode)
  idzc = ncvid(idzSliceFile,'zc',rcode)
  idh = ncvid(idzSliceFile,'h',rcode)
  call ncvgt( idzSliceFile, idxc, start(1), count(1), xc, rcode )
  call ncvgt( idzSliceFile, idyc, start(2), count(2), yc, rcode )
  call ncvgt( idzSliceFile, idzc, start, count, zc, rcode )
  call ncvgt( idzSliceFile, idh, start, count, h, rcode)
  call ncclos(idzSliceFile, rcode)
  
  WRITE(zslice_data,'("full_",I5.5,".cdf")') nstp+out3d_int
  print *, zslice_data
  idzSliceFile = ncopn(TRIM(dirout)//zslice_data, NCNOWRIT,rcode)

  idh = ncvid(idzSliceFile,'h',rcode)
  call ncvgt( idzSliceFile, idh, start, count, h1, rcode)
  call ncclos(idzSliceFile, rcode)

  WRITE(inname_data,'("face_",I5.5,".cdf")') nstp
  print*, inname_data
  idInFile = ncopn(TRIM(dirout)//inname_data, NCNOWRIT,rcode)

  iduf = ncvid(idInFile,'uf',rcode)
  idvf = ncvid(idInFile,'vf',rcode)
  idwf = ncvid(idInFile,'wf',rcode)

  call ncvgt( idInFile, iduf, start, countuf, uf, rcode )
  call ncvgt( idInFile, idvf, start, countvf, vf, rcode )
  call ncvgt( idInFile, idwf, start, countwf, wf, rcode )
  
  WRITE(inname_data,'("face_",I5.5,".cdf")') nstp+out3d_int
  print*, inname_data
  idInFile = ncopn(TRIM(dirout)//inname_data, NCNOWRIT,rcode)

  iduf = ncvid(idInFile,'uf',rcode)
  idvf = ncvid(idInFile,'vf',rcode)
  idwf = ncvid(idInFile,'wf',rcode)

  call ncvgt( idInFile, iduf, start, countuf, uf1, rcode )
  call ncvgt( idInFile, idvf, start, countvf, vf1, rcode )
  call ncvgt( idInFile, idwf, start, countwf, wf1, rcode )

  call ncclos(idInFile, rcode)
  
  uf = uf/UL
  vf = vf/UL
  wf = wf*EPS/WL

  return
end
