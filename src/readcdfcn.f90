      subroutine readcdfcn(nbegin)
		  !subroutine readcdfcn(Tr,uf,vf,wf,xc,yc,zc,xf,yf,zf,nbegin)
!     ------------------------------------------------------
!     reads in the netCDF files written by the routine outcdf.f

	  !use header, only : NI,NJ,NK,ntr,rc_kind
	  use header
	  
	  implicit none

      INCLUDE '/usr/local/include/netcdf.inc'
	  

     ! include 'dims.f90'

      integer NI2,NJ2,NK2,i,j,k,l,nbegin
      parameter ( NI2=NI+2,NJ2=NJ+2,NK2=NK+2)
	  
	  !integer :: iddatfile, idigit, idudx, idudy, idudz, idvbysection, idvby, idvbz, idvcon, idvcy, idvbx
	  !integer :: idvcz, idvc, idvdivfeddy, idvdivfreyn, idvdx, idvdz, idvd, idvfb, idvh, idvn2bar
	  !integer :: idvn2, idvnsq100m, idvnsq30m, idvpe,idvpsiv,idvpsiw,idvpv, idvp,idvrhbar,idvrho,idvrnk
	  !integer :: idvstrain,idvstress,idvstr,idvs,idvtbar,idvtemp,idvtim,idvtr,idvt,idvu,idvvb,idvvc,idvvor,idvv
	  !integer :: idvwb,idvwc,idvwpv,idvw,idvy,idvzsave,idvz,idvz3,idwdx,idwdy,idwdz,iimday,ipos
	  integer :: idvx,idvcon100,idFaceFile,idvuf,idvvf,idvwf,idvKf,idF
	  integer :: nstp
	  integer :: idT
	  REAL(kind=rc_kind) ::  rcode
!
      !double precision Tr(0:NI+1,0:NJ+1,0:NK+1,ntr), &
      !     uf(0:NI,NJ,NK),vf(NI,0:NJ,NK),wf(NI,NJ,0:NK)

      character*80 outname
      character*80 facename

      integer start(3), count(3), start2d(2), count2d(2), &
           countuf(3), countvf(3), countwf(3),count4(4),start4(4)

      DATA start /1, 1, 1/
      DATA start4 /1, 1, 1, 1/
      DATA start2d /1, 1/
      DATA count /NI2, NJ2, NK2/
      DATA count4 /NI2, NJ2, NK2, ntr/
      DATA count2d /NI2, NJ2/


      countuf(1)= NI+1
      countuf(2)= NJ
      countuf(3)= NK

      countvf(1)= NI
      countvf(2)= NJ+1
      countvf(3)= NK

      countwf(1)= NI
      countwf(2)= NJ
      countwf(3)= NK+1


!     read netCDF file
      write(6,*) 'readcdfcn, nbegin=',nbegin
!      outname = 'outxxxxx.cdf'
!      facename = 'vfcxxxxx.cdf'
      outname = 'full_xxxxx.cdf'
      facename = 'face_xxxxx.cdf'
      nstp = nbegin
      do 10 ipos = 10, 6, -1
       idigit = mod(nstp,10)
       outname(ipos:ipos) = char(ichar('0') + idigit)
       facename(ipos:ipos) = char(ichar('0') + idigit)
       nstp = nstp / 10
 10   continue

      idT = ncopn(outname, NCNOWRIT,rcode)
      idvx = ncvid(idT,'xc',rcode)
      idvy = ncvid(idT,'yc',rcode)
      idvz = ncvid(idT,'yc',rcode)

      call ncvgt( idT, idvx, start, count, xc, rcode )
      call ncvgt( idT, idvy, start, count, yc, rcode )
      call ncvgt( idT, idvz, start, count, zf, rcode )


      call ncclos(idT, rcode)
      write(6,*) 'grid read'

      return
      write(6,*) 'reading in file:', facename
      idF = ncopn(facename, NCNOWRIT,rcode)
!
      idvuf = ncvid(idF,'uf',rcode)
      idvvf = ncvid(idF,'vf',rcode)
      idvwf = ncvid(idF,'wf',rcode)
!
      call ncvgt( idF, idvuf, start, countuf, uf, rcode )
      call ncvgt( idF, idvvf, start, countvf, vf, rcode )
      call ncvgt( idF, idvwf, start, countwf, wf, rcode )
!      call ncvgt( idF, idvkf, start, countwf, Kz, rcode )

!      call ncclos(idG, rcode)
      call ncclos(idF, rcode)
      write(6,*) 'return from readcdfcn'

      return
      end
