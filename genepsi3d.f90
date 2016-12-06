program genepsi3d
implicit none

!*****************************************************************!
! 0- This program will generate all the files necessary for our
!    customize IMB based on Lagrange reconstructions 
! 1- It has to be use with ilag=1 in incompact3d.prm
! 2- To be launch before the simulation
! 3- The object is defined in the cylinder subroutine
! 4- You can add your own subroutine for your own object
! 5- The option -assume byterecl with the intel fortan compiler
! 6- To swich to double precision: real(4) --> real(8) and
!    RECL=4 --> RECL=8
! 7- Please cite the following paper if you are using this file:
! Gautier R., Laizet S. & Lamballais E., 2014, A DNS study of 
! jet control with microjets using an alterna ng direc on forcing
! strategy, Int. J. of Computa onal Fluid Dynamics, 28, 393--410
!*****************************************************************!
!
integer,parameter                  :: nx=361,ny=217,nz=48       !***! domain size (if too big then size of the sub-domain containing the object with inbig=1)
integer,parameter                  :: nclx= 2,ncly= 1,nclz= 0   !***! boundary conditions
integer,parameter                  :: nraf=10                   !***! refinement factor for the location of the solid object
integer,parameter                  :: istret=0                  !***! stretching in the y-direction (same as in incompact3d.prm)
integer,parameter                  :: nxraf=(nx-1)*nraf+1,nyraf=(ny-1)*nraf+1,nzraf=(nz)*nraf+1  !***! nx-1/ny-1/nz-1 if ncl=1/2, otherwise nx/ny/nz
integer,parameter                  :: nxbig=361,nybig=217,nzbig=48  !***! size of the full simulation (useless if inbig=0)
integer,parameter                  :: npif=2        !***! do not modify
integer,parameter                  :: izap=1        !***! do not modify
integer,parameter                  :: inbig=0       !***! equal to 1 if a sub-domain is used for the solid object otherwise 0
integer,parameter                  :: iepm=0        !***! do not modify
integer,parameter                  :: nobjmax=4     !***! number of objects that can be found in a direction
integer,dimension(ny,nz)           :: nobjx         
integer,dimension(nx,nz)           :: nobjy         
integer,dimension(nx,ny)           :: nobjz         
integer,dimension(0:nobjmax,ny,nz) :: nxipif,nxfpif 
integer,dimension(0:nobjmax,nx,nz) :: nyipif,nyfpif 
integer,dimension(0:nobjmax,nx,ny) :: nzipif,nzfpif 
real(4),dimension(  nobjmax,ny,nz) :: xi,xf         
real(4),dimension(  nobjmax,nx,nz) :: yi,yf         
real(4),dimension(  nobjmax,nx,ny) :: zi,zf         
real(4),dimension(       nx,ny,nz) :: epsi          
real(4),dimension(       nx,ny,nz) :: epsim         
real(4),dimension(             ny) :: yp            
real(4)                            :: dx,dy,dz
real(4)                            :: xlx,yly,zlz
real(4)                            :: xstart,xend
real(4)                            :: ystart,yend
real(4)                            :: zstart,zend
real(4)                            :: epm
integer                            :: i,j,k,numvis,count
integer                            :: ibig,jbig,kbig
integer                            :: istart,iend
integer                            :: jstart,jend
integer                            :: kstart,kend
integer                            :: itest,n
character(len=4) suffixe
!
!  Parameters:
   xlx=20.     !***!size of the domain or the sub-domain
   yly=12.     
   zlz=3.  
   if(nclx.eq.0)dx=xlx/nx
   if(nclx.eq.1.or.nclx.eq.2)dx=xlx/(nx-1)
   if (istret==0) then
   if(ncly.eq.0)dy=yly/ny
   if(ncly.eq.1.or.ncly.eq.2)dy=yly/(ny-1)
   do j=1,ny
      yp(j)=(j-1)*dy
   enddo
   else
      open(1,file='yp.dat',form='unformatted')
      read(1) yp
      close(1)
   endif
   if(nclz.eq.0)dz=zlz/nz
   if(nclz.eq.1.or.nclz.eq.2)dz=zlz/(nz-1)
   
   if(inbig.eq.1)then
      istart=1       
      iend=nx       
      jstart=1     
      jend=ny       
      kstart=1     
      kend=nz       
   endif
   if(inbig.eq.1)then  
      xstart=0.        
      ystart=yp(jstart)
      zstart=3.        
   endif

   if(iepm.eq.1)then
      epm=1.        
   endif

   call gene_epsi_3D(epsi,nx,ny,nz,dx,dz,xlx,yly,zlz ,&
                        nclx,ncly,nclz,nxraf,nyraf,nzraf   ,&
                        xi,xf,yi,yf,zi,zf,nobjx,nobjy,nobjz,&
                        nobjmax,yp,nraf)
   call verif_epsi(epsi,npif,izap,nx,ny,nz,nobjmax,&
                   nxipif,nxfpif,nyipif,nyfpif,nzipif,nzfpif)

   if(inbig.eq.0)then
      call write_geomcomplex(nx,ny,nz,nobjx,nobjy,nobjz,xi,xf,yi,yf,zi,zf,&
                             nxipif,nxfpif,nyipif,nyfpif,nzipif,nzfpif,nobjmax,npif,&
                             inbig,istart,jstart,kstart,iend,jend,kend,nxbig,nybig,nzbig,&
                             xstart,ystart,zstart)
      OPEN(10,FILE='epsilon.dat',FORM='UNFORMATTED',&
           ACCESS='DIRECT', RECL=4)
      COUNT = 1
      DO K=1,nz
      DO J=1,ny
      DO I=1,nx
         WRITE(10,REC=COUNT) epsi(I,J,K)
         COUNT = COUNT + 1
      ENDDO
      ENDDO
      ENDDO
      CLOSE(10)
      if(iepm.eq.1)then
         call gene_epsim(epsi,epsim,nx,ny,nz,nobjx,nobjy,nobjz,&
                         xi,xf,yi,yf,zi,zf,nobjmax,dx,dz,epm,&
                         xlx,yly,zlz,yp)
         OPEN(10,FILE='epsilonm.dat',FORM='UNFORMATTED',&
              ACCESS='DIRECT', RECL=4)
         COUNT = 1
         DO K=1,nz
         DO J=1,ny
         DO I=1,nx
            WRITE(10,REC=COUNT) epsim(I,J,K)
            COUNT = COUNT + 1
         ENDDO
         ENDDO
         ENDDO
         CLOSE(10)
      endif
   else
      call write_geomcomplex(nx,ny,nz,nobjx,nobjy,nobjz,xi,xf,yi,yf,zi,zf,&
                             nxipif,nxfpif,nyipif,nyfpif,nzipif,nzfpif,nobjmax,npif,&
                             inbig,istart,jstart,kstart,iend,jend,kend,nxbig,nybig,nzbig,&
                             xstart,ystart,zstart)
      OPEN(10,FILE='epsilon.dat',FORM='UNFORMATTED',&
           ACCESS='DIRECT', RECL=4)
      COUNT = 1
      do kbig=1,nzbig
      do jbig=1,nybig
      do ibig=1,nxbig
         if(kbig.lt.kstart.or.kbig.gt.kend.or.&
            jbig.lt.jstart.or.jbig.gt.jend.or.&
            ibig.lt.istart.or.ibig.gt.iend)then
            WRITE(10,REC=COUNT) 0.
            COUNT = COUNT + 1
         else
            k=kbig+1-kstart
            j=jbig+1-jstart
            i=ibig+1-istart
            WRITE(10,REC=COUNT) epsi(i,j,k)
            COUNT = COUNT + 1
         endif
      enddo
      enddo
      enddo
      close(10)
      if(iepm.eq.1)then
         call gene_epsim(epsi,epsim,nx,ny,nz,nobjx,nobjy,nobjz,&
                         xi,xf,yi,yf,zi,zf,nobjmax,dx,dz,epm,&
                         xlx,yly,zlz,yp)
         OPEN(10,FILE='epsilonm.dat',FORM='UNFORMATTED',&
              ACCESS='DIRECT', RECL=4)
         COUNT = 1
         do kbig=1,nzbig
         do jbig=1,nybig
         do ibig=1,nxbig
            if(kbig.lt.kstart.or.kbig.gt.kend.or.&
               jbig.lt.jstart.or.jbig.gt.jend.or.&
               ibig.lt.istart.or.ibig.gt.iend)then
               WRITE(10,REC=COUNT) 0.
               COUNT = COUNT + 1
            else
               k=kbig+1-kstart
               j=jbig+1-jstart
               i=ibig+1-istart
               WRITE(10,REC=COUNT) epsim(i,j,k)
               COUNT = COUNT + 1
            endif
         enddo
         enddo
         enddo
         close(10)
      endif
   endif
!
end
!
!***************************************************************************
!***************************************************************************
!***************************************************************************
!
subroutine gene_epsi_3D(epsi,nx,ny,nz,dx,dz,xlx,yly,zlz ,&
                        nclx,ncly,nclz,nxraf,nyraf,nzraf   ,&
                        xi,xf,yi,yf,zi,zf,nobjx,nobjy,nobjz,&
                        nobjmax,yp,nraf)
   implicit none
!
   real(4),dimension(nx,ny,nz)      :: epsi
   integer                          :: nx,ny,nz,nobjmax
   real(4)                          :: dx,dz
   real(4)                          :: xlx,yly,zlz
   integer                          :: nclx,ncly,nclz
   integer                          :: nxraf,nyraf,nzraf
   integer                          :: nraf
   integer,dimension(ny,nz)         :: nobjx,nobjxraf
   integer,dimension(nx,nz)         :: nobjy,nobjyraf
   integer,dimension(nx,ny)         :: nobjz,nobjzraf
   real(4),dimension(nobjmax,ny,nz) :: xi,xf
   real(4),dimension(nobjmax,nx,nz) :: yi,yf
   real(4),dimension(nobjmax,nx,ny) :: zi,zf
   real(4),dimension(  nxraf,ny,nz) :: xepsi
   real(4),dimension(  nx,nyraf,nz) :: yepsi
   real(4),dimension(  nx,ny,nzraf) :: zepsi
   real(4),dimension(           ny) :: yp
   real(4),dimension(        nyraf) :: ypraf
   real(4)                          :: dxraf,dzraf
   integer                          :: i,j,k
   integer                          :: ii,jj,kk
   real(4)                          :: x,y,z
   integer                          :: inum,jnum,knum
   integer                          :: ibug,jbug,kbug
   integer                          :: iobj,jobj,kobj
   integer                          :: iflu,jflu,kflu
   integer                          :: isol,jsol,ksol
   integer                          :: iraf,jraf,kraf
   integer                          :: nobjxmax ,nobjymax ,nobjzmax
   integer                          :: nobjxmaxraf,nobjymaxraf,nobjzmaxraf
   integer                          :: idebraf,jdebraf,kdebraf
   integer                          :: ifinraf,jfinraf,kfinraf
   character(len=4) suffixe
   integer                          :: numvis

   epsi(:,:,:)=0.
   call cylinder(epsi,nx,ny,nz,dx,yp,dz,1.) 
   print*,'step 1'
   if(nclx.eq.0)then
      dxraf =xlx/nxraf
   elseif(nclx.eq.1.or.nclx.eq.2)then
      dxraf =xlx/(nxraf-1)
   endif
   xepsi(:,:,:)=0.
   call cylinder(xepsi,nxraf,ny,nz,dxraf,yp,dz,1.) 
   print*,'step 2'
   do j=1,ny-1
      do jraf=1,nraf
         ypraf(jraf+nraf*(j-1))=yp(j)+(jraf-1)*(yp(j+1)-yp(j))/nraf
      enddo
   enddo
   if(ncly.ne.0)ypraf(nyraf)=yp(ny)

   yepsi(:,:,:)=0.
   call cylinder(yepsi,nx,nyraf,nz,dx,ypraf,dz,1.)
print*,'step 3'
   if(nclz.eq.0)then
      dzraf=zlz/nzraf
   elseif(nclz.eq.1.or.nclz.eq.2)then
      dzraf=zlz/(nzraf-1)
   endif
   zepsi(:,:,:)=0.
   call cylinder(zepsi,nx,ny,nzraf,dx,yp,dzraf,1.) 
   print*,'step 4'

   nobjx(:,:)=0
   nobjxmax=0
   do k=1,nz
   do j=1,ny
      inum=0
      if(epsi(1,j,k).eq.1.)then
         inum=1
         nobjx(j,k)=1
      endif
      do i=1,nx-1
         if(epsi(i,j,k).eq.0..and.epsi(i+1,j,k).eq.1.)then
            inum=inum+1
            nobjx(j,k)=nobjx(j,k)+1
         endif
      enddo
      if(inum.gt.nobjxmax)then
         nobjxmax=inum
      endif
   enddo
   enddo
   print*,'nobjxmax=',nobjxmax

   nobjxraf(:,:)=0
   ibug=0
   nobjxmaxraf=0
   inum=0
   do k=1,nz
   do j=1,ny
      inum=0
      if(xepsi(1,j,k).eq.1.)then
         inum=1
         nobjxraf(j,k)=1
      endif
      do i=1,nxraf-1
         if(xepsi(i,j,k).eq.0..and.xepsi(i+1,j,k).eq.1.)then
            inum=inum+1
            nobjxraf(j,k)=nobjxraf(j,k)+1
         endif
      enddo
      if(inum.gt.nobjxmaxraf)then
         nobjxmaxraf=inum
      endif
      if(nobjx(j,k).ne.nobjxraf(j,k))then
         ibug=ibug+1
      endif
   enddo
   enddo
   print*,'nobjxmaxraf=',nobjxmaxraf
   print*,'ibug=',ibug
   print*,'step 5'

   nobjy(:,:)=0
   nobjymax=0
   do k=1,nz
   do i=1,nx
      jnum=0
      if(epsi(i,1,k).eq.1.)then
         jnum=1
         nobjy(i,k)=1
      endif
      do j=1,ny-1
         if(epsi(i,j,k).eq.0..and.epsi(i,j+1,k).eq.1.)then
            jnum=jnum+1
            nobjy(i,k)=nobjy(i,k)+1
         endif
      enddo
      if(jnum.gt.nobjymax)then
         nobjymax=jnum
      endif
   enddo
   enddo
   print*,'nobjymax=',nobjymax

   nobjyraf(:,:)=0
   jbug=0
   nobjymaxraf=0
   jnum=0
   do k=1,nz
   do i=1,nx
      jnum=0
      if(yepsi(i,1,k).eq.1.)then
         jnum=1
         nobjyraf(i,k)=1
      endif
      do j=1,nyraf-1
         if(yepsi(i,j,k).eq.0..and.yepsi(i,j+1,k).eq.1.)then
            jnum=jnum+1
            nobjyraf(i,k)=nobjyraf(i,k)+1
         endif
      enddo
      if(jnum.gt.nobjymaxraf)then
         nobjymaxraf=jnum
      endif
      if(nobjy(i,k).ne.nobjyraf(i,k))then
         jbug=jbug+1
      endif
   enddo
   enddo
   print*,'nobjymaxraf=',nobjymaxraf
   print*,'jbug=',jbug
   print*,'step 6'

   nobjz(:,:)=0
   nobjzmax=0
   do j=1,ny
   do i=1,nx
      knum=0
      if(epsi(i,j,1).eq.1.)then
         knum=1
         nobjz(i,j)=1
      endif
      do k=1,nz-1
         if(epsi(i,j,k).eq.0..and.epsi(i,j,k+1).eq.1.)then
            knum=knum+1
            nobjz(i,j)=nobjz(i,j)+1
         endif
      enddo
      if(knum.gt.nobjzmax)then
         nobjzmax=knum
      endif
   enddo
   enddo
   print*,'nobjzmax=',nobjzmax

   nobjzraf(:,:)=0
   kbug=0
   nobjzmaxraf=0
   knum=0
   do j=1,ny
   do i=1,nx
      knum=0
      if(zepsi(i,j,1).eq.1.)then
         knum=1
         nobjzraf(i,j)=1
      endif
      do k=1,nzraf-1
         if(zepsi(i,j,k).eq.0..and.zepsi(i,j,k+1).eq.1.)then
            knum=knum+1
            nobjzraf(i,j)=nobjzraf(i,j)+1
         endif
      enddo
      if(knum.gt.nobjzmaxraf)then
         nobjzmaxraf=knum
      endif
      if(nobjz(i,j).ne.nobjzraf(i,j))then
         kbug=kbug+1
      endif
   enddo
   enddo
   print*,'nobjzmaxraf=',nobjzmaxraf
   print*,'kbug=',kbug
   print*,'step 7'

   do k=1,nz
   do j=1,ny
      inum=0
      if(xepsi(1,j,k).eq.1.)then
         inum=inum+1
         xi(inum,j,k)=-dx!-xlx
      endif
      do i=1,nxraf-1
         if(xepsi(i,j,k).eq.0..and.xepsi(i+1,j,k).eq.1.)then
            inum=inum+1
            xi(inum,j,k)=dxraf*(i-1)+dxraf/2.
         elseif(xepsi(i,j,k).eq.1..and.xepsi(i+1,j,k).eq.0.)then
            xf(inum,j,k)=dxraf*(i-1)+dxraf/2.
         endif
      enddo
      if(xepsi(nxraf,j,k).eq.1.)then
         xf(inum,j,k)=xlx+dx!2.*xlx
      endif
   enddo
   enddo

   if(ibug.ne.0)then
      do k=1,nz
      do j=1,ny
         if(nobjx(j,k).ne.nobjxraf(j,k))then
            iobj=0
            if(epsi(1,j,k).eq.1.)iobj=iobj+1
            do i=1,nx-1
               if(epsi(i,j,k).eq.0..and.epsi(i+1,j,k).eq.1.)iobj=iobj+1
               if(epsi(i,j,k).eq.0..and.epsi(i+1,j,k).eq.0.)iflu=1
               if(epsi(i,j,k).eq.1..and.epsi(i+1,j,k).eq.1.)isol=1
               do iraf=1,nraf
                  if(xepsi(iraf+nraf*(i-1)  ,j,k).eq.0..and.&
                     xepsi(iraf+nraf*(i-1)+1,j,k).eq.1.)idebraf=iraf+nraf*(i-1)+1
                  if(xepsi(iraf+nraf*(i-1)  ,j,k).eq.1..and.&
                     xepsi(iraf+nraf*(i-1)+1,j,k).eq.0.)ifinraf=iraf+nraf*(i-1)+1
               enddo
               if(idebraf.ne.0.and.ifinraf.ne.0.and.&
                  idebraf.lt.ifinraf.and.iflu.eq.1)then
                  iobj=iobj+1
                  do ii=iobj,nobjmax-1
                     xi(ii,j,k)=xi(ii+1,j,k)
                     xf(ii,j,k)=xf(ii+1,j,k)
                  enddo
                  iobj=iobj-1
               endif
               if(idebraf.ne.0.and.ifinraf.ne.0.and.&
                  idebraf.gt.ifinraf.and.isol.eq.1)then
                  iobj=iobj+1
                  do ii=iobj,nobjmax-1
                     xi(ii,j,k)=xi(ii+1,j,k)
                  enddo
                  iobj=iobj-1
                  do ii=iobj,nobjmax-1
                     xf(ii,j,k)=xf(ii+1,j,k)
                  enddo
               endif
               idebraf=0
               ifinraf=0
               iflu=0
            enddo
         endif
      enddo
      enddo
   endif
   print*,'step 8'

   do k=1,nz
   do i=1,nx
      jnum=0
      if(yepsi(i,1,k).eq.1.)then
         jnum=jnum+1
         yi(jnum,i,k)=-(yp(2)-yp(1))!-yly
      endif
      do j=1,nyraf-1
         if(yepsi(i,j,k).eq.0..and.yepsi(i,j+1,k).eq.1.)then
            jnum=jnum+1
            yi(jnum,i,k)=ypraf(j)+(ypraf(j+1)-ypraf(j))/2.!dyraf*(j-1)+dyraf/2.
         elseif(yepsi(i,j,k).eq.1..and.yepsi(i,j+1,k).eq.0.)then
            yf(jnum,i,k)=ypraf(j)+(ypraf(j+1)-ypraf(j))/2.!dyraf*(j-1)+dyraf/2.!
         endif
      enddo
      if(yepsi(i,nyraf,k).eq.1.)then
         yf(jnum,i,k)=yly+(yp(ny)-yp(ny-1))/2.!2.*yly
      endif
   enddo
   enddo

   if(jbug.ne.0)then
      do k=1,nz
      do i=1,nx
         if(nobjy(i,k).ne.nobjyraf(i,k))then
            jobj=0
            if(epsi(i,1,k).eq.1.)jobj=jobj+1
            do j=1,ny-1
               if(epsi(i,j,k).eq.0..and.epsi(i,j+1,k).eq.1.)jobj=jobj+1
               if(epsi(i,j,k).eq.0..and.epsi(i,j+1,k).eq.0.)jflu=1
               if(epsi(i,j,k).eq.1..and.epsi(i,j+1,k).eq.1.)jsol=1
               do jraf=1,nraf
                  if(yepsi(i,jraf+nraf*(j-1)  ,k).eq.0..and.&
                     yepsi(i,jraf+nraf*(j-1)+1,k).eq.1.)jdebraf=jraf+nraf*(j-1)+1
                  if(yepsi(i,jraf+nraf*(j-1)  ,k).eq.1..and.&
                     yepsi(i,jraf+nraf*(j-1)+1,k).eq.0.)jfinraf=jraf+nraf*(j-1)+1
               enddo
               if(jdebraf.ne.0.and.jfinraf.ne.0.and.&
                  jdebraf.lt.jfinraf.and.jflu.eq.1)then
                  jobj=jobj+1
                  do jj=jobj,nobjmax-1
                     yi(jj,i,k)=yi(jj+1,i,k)
                     yf(jj,i,k)=yf(jj+1,i,k)
                  enddo
                  jobj=jobj-1
               endif
               if(jdebraf.ne.0.and.jfinraf.ne.0.and.&
                  jdebraf.gt.jfinraf.and.jsol.eq.1)then
                  jobj=jobj+1
                  do jj=jobj,nobjmax-1
                     yi(jj,i,k)=yi(jj+1,i,k)
                  enddo
                  jobj=jobj-1
                  do jj=jobj,nobjmax-1
                     yf(jj,i,k)=yf(jj+1,i,k)
                  enddo
               endif
               jdebraf=0
               jfinraf=0
               jflu=0
            enddo
         endif
      enddo
      enddo
   endif
   print*,'step 9'

   do j=1,ny
   do i=1,nx
      knum=0
      if(zepsi(i,j,1).eq.1.)then
         knum=knum+1
         zi(knum,i,j)=-dz!zlz
      endif
      do k=1,nzraf-1
         if(zepsi(i,j,k).eq.0..and.zepsi(i,j,k+1).eq.1.)then
            knum=knum+1
            zi(knum,i,j)=dzraf*(k-1)+dzraf/2.
         elseif(zepsi(i,j,k).eq.1..and.zepsi(i,j,k+1).eq.0.)then
            zf(knum,i,j)=dzraf*(k-1)+dzraf/2.
         endif
      enddo
      if(zepsi(i,j,nzraf).eq.1.)then
         zf(knum,i,j)=zlz+dz!2.*zlz
      endif
   enddo
   enddo

   if(kbug.ne.0)then
      do j=1,ny
      do i=1,nx
         if(nobjz(i,j).ne.nobjzraf(i,j))then
            kobj=0
            if(epsi(i,j,1).eq.1.)kobj=kobj+1
            do k=1,nz-1
               if(epsi(i,j,k).eq.0..and.epsi(i,j,k+1).eq.1.)kobj=kobj+1
               if(epsi(i,j,k).eq.0..and.epsi(i,j,k+1).eq.0.)kflu=1
               if(epsi(i,j,k).eq.1..and.epsi(i,j,k+1).eq.1.)ksol=1
               do kraf=1,nraf
                  if(zepsi(i,j,kraf+nraf*(k-1)  ).eq.0..and.&
                     zepsi(i,j,kraf+nraf*(k-1)+1).eq.1.)kdebraf=kraf+nraf*(k-1)+1
                  if(zepsi(i,j,kraf+nraf*(k-1)  ).eq.1..and.&
                     zepsi(i,j,kraf+nraf*(k-1)+1).eq.0.)kfinraf=kraf+nraf*(k-1)+1
               enddo
               if(kdebraf.ne.0.and.kfinraf.ne.0.and.&
                  kdebraf.lt.kfinraf.and.kflu.eq.1)then
                  kobj=kobj+1
                  do kk=kobj,nobjmax-1
                     zi(kk,i,j)=zi(kk+1,i,j)
                     zf(kk,i,j)=zf(kk+1,i,j)
                  enddo
                  kobj=kobj-1
               endif
               if(kdebraf.ne.0.and.kfinraf.ne.0.and.&
                  kdebraf.gt.kfinraf.and.ksol.eq.1)then
                  kobj=kobj+1
                  do kk=kobj,nobjmax-1
                     zi(kk,i,j)=zi(kk+1,i,j)
                  enddo
                  kobj=kobj-1
                  do kk=kobj,nobjmax-1
                     zf(kk,i,j)=zf(kk+1,i,j)
                  enddo
               endif
               kdebraf=0
               kfinraf=0
               kflu=0
            enddo
         endif
      enddo
      enddo
   endif
   print*,'step 10'
!
   return
end subroutine gene_epsi_3D
!
!***************************************************************************
!***************************************************************************
!***************************************************************************
!
subroutine gene_epsim(epsi,epsim,nx,ny,nz,nobjx,nobjy,nobjz,&
                      xi,xf,yi,yf,zi,zf,nobjmax,dx,dz,epm,&
                      xlx,yly,zlz,yp)
   implicit none
!
   real(4),dimension(       nx,ny,nz) :: epsi
   real(4),dimension(       nx,ny,nz) :: epsim
   real(4),dimension(  nobjmax,ny,nz) :: xi,xf
   real(4),dimension(  nobjmax,nx,nz) :: yi,yf
   real(4),dimension(  nobjmax,nx,ny) :: zi,zf
   integer,dimension(          ny,nz) :: nobjx
   integer,dimension(          nx,nz) :: nobjy
   integer,dimension(          nx,ny) :: nobjz
   real(4),dimension(             ny) :: yp
   real(4)                            :: x,y,z
   real(4)                            :: dx,dz
   real(4)                            :: xlx,yly,zlz
   real(4)                            :: xe,ye,ze,epm
   integer                            :: i,j,k
   integer                            :: nx,ny,nz
   integer                            :: ix,jy,kz
   integer                            :: nobjmax
!
   epsim(:,:,:)=epsi(:,:,:)
   xe=epm*dx
   do k=1,nz
   do j=1,ny
   do i=1,nx
      if(epsi(i,j,k).eq.1..and.&
         nobjx (j,k).ne.0)then
         x=dx*(i-1)
         do ix=1,nobjx(j,k)
            if(x  .ge.xi(ix,j,k)   .and.& 
               x  .le.xi(ix,j,k)+xe.and.&
               0. .lt.xi(ix,j,k)   .or. & 
               x  .le.xf(ix,j,k)   .and.&
               x  .ge.xf(ix,j,k)-xe.and.&
               xlx.gt.xf(ix,j,k)   )then
               epsim(i,j,k)=0.
            endif
         enddo
      endif
   enddo
   enddo
   enddo
   do k=1,nz
   do i=1,nx
   do j=2,ny-1
      if(epsi(i,j,k).eq.1..and.&
         nobjy (i,k).gt.0)then
         y=yp(j)
         ye=epm*(yp(j+1)-yp(j-1))/2.
         do jy=1,nobjy(i,k)
            if(y  .ge.yi(jy,i,k)   .and.& 
               y  .le.yi(jy,i,k)+ye.and.&
               0. .lt.yi(jy,i,k)   .or. & 
               y  .le.yf(jy,i,k)   .and.&
               y  .ge.yf(jy,i,k)-ye.and.&
               yly.gt.yf(jy,i,k)   )then
               epsim(i,j,k)=0.
            endif
         enddo
      endif
   enddo
   enddo
   enddo
   ze=epm*dz
   do j=1,ny
   do i=1,nx
   do k=1,nz
      if(epsi(i,j,k).eq.1..and.&
         nobjz (i,j).gt.0)then
         z=dz*(k-1)
         do kz=1,nobjz(i,j)
            if(z  .ge.zi(kz,i,j)   .and.& 
               z  .le.zi(kz,i,j)+ze.and.&
               0. .lt.zi(kz,i,j)   .or. & 
               z  .le.zf(kz,i,j)   .and.&
               z  .ge.zf(kz,i,j)-ze.and.&
               zlz.gt.zf(kz,i,j)   )then
               epsim(i,j,k)=0.
            endif
         enddo
      endif
   enddo
   enddo
   enddo
!
   return
end subroutine gene_epsim
!
!***************************************************************************
!***************************************************************************
!***************************************************************************
!
subroutine verif_epsi(epsi,npif,izap,nx,ny,nz,nobjmax,&
                      nxipif,nxfpif,nyipif,nyfpif,nzipif,nzfpif)
   implicit none
!
   integer                            :: nx,ny,nz,nobjmax
   real(4),dimension(       nx,ny,nz) :: epsi
   integer,dimension(0:nobjmax,ny,nz) :: nxipif,nxfpif
   integer,dimension(0:nobjmax,nx,nz) :: nyipif,nyfpif
   integer,dimension(0:nobjmax,nx,ny) :: nzipif,nzfpif
   integer                            :: npif,izap
   integer                            :: i,j,k
   integer                            :: inum ,jnum ,knum
   integer                            :: iflu ,jflu ,kflu
   integer                            :: ising,jsing,ksing,itest

   nxipif(:,:,:)=npif
   nxfpif(:,:,:)=npif
   ising=0
   do k=1,nz
   do j=1,ny
      inum=0
      iflu=0
      if(epsi(1,j,k).eq.1.)inum=inum+1
      if(epsi(1,j,k).eq.0.)iflu=iflu+1
      do i=2,nx
         if(epsi(i  ,j,k).eq.0.)iflu=iflu+1
         if(epsi(i-1,j,k).eq.0..and.&
            epsi(i  ,j,k).eq.1.)then
            inum=inum+1
            if(inum.eq.1)then
               nxipif(inum  ,j,k)=iflu-izap
               if(iflu-izap.lt.npif)ising=ising+1
               if(iflu-izap.ge.npif)nxipif(inum  ,j,k)=npif
               iflu=0
            else
               nxipif(inum  ,j,k)=iflu-izap
               nxfpif(inum-1,j,k)=iflu-izap
               if(iflu-izap.lt.npif)ising=ising+1
               if(iflu-izap.ge.npif)nxipif(inum  ,j,k)=npif
               if(iflu-izap.ge.npif)nxfpif(inum-1,j,k)=npif
               iflu=0
            endif
         endif
         if(epsi(i,j,k).eq.1.)iflu=0
      enddo
      if(epsi(nx,j,k).eq.0.)then
         nxfpif(inum,j,k)=iflu-izap
         if(iflu-izap.lt.npif)ising=ising+1
         if(iflu-izap.lt.npif)nxfpif(inum,j,k)=npif
      endif
   enddo
   enddo
   print*,'number of points with potential problem in X :',ising
   print*,' step 11'

   nyipif(:,:,:)=npif
   nyfpif(:,:,:)=npif
   jsing=0
   do k=1,nz
   do i=1,nx
      jnum=0
      jflu=0
      if(epsi(i,1,k).eq.1.)jnum=jnum+1
      if(epsi(i,1,k).eq.0.)jflu=jflu+1
      do j=2,ny
         if(epsi(i,j  ,k).eq.0.)jflu=jflu+1
         if(epsi(i,j-1,k).eq.0..and.&
            epsi(i,j  ,k).eq.1.)then
            jnum=jnum+1
            if(jnum.eq.1)then
               nyipif(jnum  ,i,k)=jflu-izap
               if(jflu-izap.lt.npif)jsing=jsing+1
               if(jflu-izap.ge.npif)nyipif(jnum  ,i,k)=npif
               jflu=0
            else
               nyipif(jnum  ,i,k)=jflu-izap
               nyfpif(jnum-1,i,k)=jflu-izap
               if(jflu-izap.lt.npif)jsing=jsing+1
               if(jflu-izap.ge.npif)nyipif(jnum  ,i,k)=npif
               if(jflu-izap.ge.npif)nyfpif(jnum-1,i,k)=npif
               jflu=0
            endif
         endif
         if(epsi(i,j,k).eq.1.)jflu=0
      enddo
      if(epsi(i,ny,k).eq.0.)then
         nyfpif(jnum,i,k)=jflu-izap
         if(jflu-izap.lt.npif)jsing=jsing+1
         if(jflu-izap.lt.npif)nyfpif(jnum,i,k)=npif
      endif
   enddo
   enddo
   print*,'number of points with potential problem in Y :',jsing
   print*,'step 12'

   if(nz.gt.1)then
      nzipif(:,:,:)=npif
      nzfpif(:,:,:)=npif
      ksing=0
      do j=1,ny
      do i=1,nx
         knum=0
         kflu=0
         if(epsi(i,j,1).eq.1.)knum=knum+1
         if(epsi(i,j,1).eq.0.)kflu=kflu+1
         do k=2,nz
            if(epsi(i,j,k  ).eq.0.)kflu=kflu+1
            if(epsi(i,j,k-1).eq.0..and.&
               epsi(i,j,k  ).eq.1.)then
               knum=knum+1
               if(knum.eq.1)then
                  nzipif(knum  ,i,j)=kflu-izap
                  if(kflu-izap.lt.npif)ksing=ksing+1
                  if(kflu-izap.ge.npif)nzipif(knum  ,i,j)=npif
                  kflu=0
               else
                  nzipif(knum  ,i,j)=kflu-izap
                  nzfpif(knum-1,i,j)=kflu-izap
                  if(kflu-izap.lt.npif)ksing=ksing+1
                  if(kflu-izap.ge.npif)nzipif(knum  ,i,j)=npif
                  if(kflu-izap.ge.npif)nzfpif(knum-1,i,j)=npif
                  kflu=0
               endif
            endif
            if(epsi(i,j,k).eq.1.)kflu=0
         enddo
         if(epsi(i,j,nz).eq.0.)then
            nzfpif(knum,i,j)=kflu-izap
            if(kflu-izap.lt.npif)ksing=ksing+1
            if(kflu-izap.lt.npif)nzfpif(knum,i,j)=npif
         endif
      enddo
      enddo
      print*,'number of points with potential problem in Z :',ksing
   endif  
   print*,'step 13'
!
   return
end subroutine verif_epsi
!
!***************************************************************************
!***************************************************************************
!***************************************************************************
!
subroutine write_geomcomplex(nx,ny,nz,nobjx,nobjy,nobjz,xi,xf,yi,yf,zi,zf,&
                             nxipif,nxfpif,nyipif,nyfpif,nzipif,nzfpif,nobjmax,npif,&
                             inbig,istart,jstart,kstart,iend,jend,kend,nxbig,nybig,nzbig,&
                             xstart,ystart,zstart)
   implicit none
!
   integer                            :: nx,ny,nz,nobjmax
   integer,dimension(ny,nz)           :: nobjx
   integer,dimension(nx,nz)           :: nobjy
   integer,dimension(nx,ny)           :: nobjz
   real(4),dimension(  nobjmax,ny,nz) :: xi,xf
   real(4),dimension(  nobjmax,nx,nz) :: yi,yf
   real(4),dimension(  nobjmax,nx,ny) :: zi,zf
   integer,dimension(0:nobjmax,ny,nz) :: nxipif,nxfpif
   integer,dimension(0:nobjmax,nx,nz) :: nyipif,nyfpif
   integer,dimension(0:nobjmax,nx,ny) :: nzipif,nzfpif
   integer                            :: i,j,k
   integer                            :: ibig,jbig,kbig
   integer                            :: nxbig,nybig,nzbig
   integer                            :: inbig,npif
   integer                            :: istart,iend
   integer                            :: jstart,jend
   integer                            :: kstart,kend
   real(4)                            :: xstart
   real(4)                            :: ystart
   real(4)                            :: zstart
!
   if(inbig.eq.0)then
      open(11,file='nobjx.dat'  ,form='formatted')
      do k=1,nz
      do j=1,ny
         write(11,*)nobjx(j,k)
      enddo
      enddo
      close(11)
      open(12,file='nobjy.dat'  ,form='formatted')
      do k=1,nz
      do i=1,nx
         write(12,*)nobjy(i,k)
      enddo
      enddo
      close(12)
      open(13,file='nobjz.dat'  ,form='formatted')
      do j=1,ny
      do i=1,nx
         write(13,*)nobjz(i,j)
      enddo
      enddo
      close(13)
      open(21,file='nxifpif.dat',form='formatted')
      do k=1,nz
      do j=1,ny
      do i=0,nobjmax
         write(21,*)nxipif(i,j,k),nxfpif(i,j,k)
      enddo
      enddo
      enddo
      close(21)
      open(22,file='nyifpif.dat',form='formatted')
      do k=1,nz
      do i=1,nx
      do j=0,nobjmax
         write(22,*)nyipif(j,i,k),nyfpif(j,i,k)
      enddo
      enddo
      enddo
      close(22)
      open(23,file='nzifpif.dat',form='formatted')
      do j=1,ny
      do i=1,nx
      do k=0,nobjmax
         write(23,*)nzipif(k,i,j),nzfpif(k,i,j)
      enddo
      enddo
      enddo
      close(23)
      open(31,file='xixf.dat'   ,form='formatted')
      do k=1,nz
      do j=1,ny
      do i=1,nobjmax
         write(31,*)xi(i,j,k),xf(i,j,k)
      enddo
      enddo
      enddo
      close(31)
      open(32,file='yiyf.dat'   ,form='formatted')
      do k=1,nz
      do i=1,nx
      do j=1,nobjmax
         write(32,*)yi(j,i,k),yf(j,i,k)
      enddo
      enddo
      enddo
      close(32)
      open(33,file='zizf.dat'   ,form='formatted')
      do j=1,ny
      do i=1,nx
      do k=1,nobjmax
         write(33,*)zi(k,i,j),zf(k,i,j)
      enddo
      enddo
      enddo
      close(33)
   else
      open(11,file='nobjx.dat'  ,form='formatted')
      do kbig=1,nzbig
      do jbig=1,nybig
         if(kbig.lt.kstart.or.kbig.gt.kend.or.&
            jbig.lt.jstart.or.jbig.gt.jend)then
            write(11,*)0
         else
            k=kbig+1-kstart
            j=jbig+1-jstart
            write(11,*)nobjx(j,k)
         endif
      enddo
      enddo
      close(11)
      open(12,file='nobjy.dat'  ,form='formatted')
      do kbig=1,nzbig
      do ibig=1,nxbig
         if(kbig.lt.kstart.or.kbig.gt.kend.or.&
            ibig.lt.istart.or.ibig.gt.iend)then
            write(12,*)0
         else
            k=kbig+1-kstart
            i=ibig+1-istart
            write(12,*)nobjy(i,k)
         endif
      enddo
      enddo
      close(12)
      open(13,file='nobjz.dat'  ,form='formatted')
      do jbig=1,nybig
      do ibig=1,nxbig
         if(jbig.lt.jstart.or.jbig.gt.jend.or.&
            ibig.lt.istart.or.ibig.gt.iend)then
            write(13,*)0
         else
            j=jbig+1-jstart
            i=ibig+1-istart
            write(13,*)nobjz(i,j)
         endif
      enddo
      enddo
      close(13)
      open(21,file='nxifpif.dat',form='formatted')
      do kbig=1,nzbig
      do jbig=1,nybig
         if(kbig.lt.kstart.or.kbig.gt.kend.or.&
            jbig.lt.jstart.or.jbig.gt.jend)then
            do i=0,nobjmax
               write(21,*)npif,npif
            enddo
         else
            k=kbig+1-kstart
            j=jbig+1-jstart
            do i=0,nobjmax
               write(21,*)nxipif(i,j,k),nxfpif(i,j,k)
            enddo
         endif
      enddo
      enddo
      close(21)
      open(22,file='nyifpif.dat',form='formatted')
      do kbig=1,nzbig
      do ibig=1,nxbig
         if(kbig.lt.kstart.or.kbig.gt.kend.or.&
            ibig.lt.istart.or.ibig.gt.iend)then
            do j=0,nobjmax
               write(22,*)npif,npif
            enddo
         else
            k=kbig+1-kstart
            i=ibig+1-istart
            do j=0,nobjmax
               write(22,*)nyipif(j,i,k),nyfpif(j,i,k)
            enddo
         endif
      enddo
      enddo
      close(22)
      open(23,file='nzifpif.dat',form='formatted')
      do jbig=1,nybig
      do ibig=1,nxbig
         if(jbig.lt.jstart.or.jbig.gt.jend.or.&
            ibig.lt.istart.or.ibig.gt.iend)then
            do k=0,nobjmax
               write(23,*)npif,npif
            enddo
         else
            j=jbig+1-jstart
            i=ibig+1-istart
            do k=0,nobjmax
               write(23,*)nzipif(k,i,j),nzfpif(k,i,j)
            enddo
         endif
      enddo
      enddo
      close(23)
      open(31,file='xixf.dat'   ,form='formatted')
      do kbig=1,nzbig
      do jbig=1,nybig
         if(kbig.lt.kstart.or.kbig.gt.kend.or.&
            jbig.lt.jstart.or.jbig.gt.jend)then
            do i=1,nobjmax
               write(31,*)0.,0.
            enddo
         else
            k=kbig+1-kstart
            j=jbig+1-jstart
            do i=1,nobjmax
               write(31,*)xi(i,j,k)+xstart,xf(i,j,k)+xstart
            enddo
         endif
      enddo
      enddo
      close(31)
      open(32,file='yiyf.dat'   ,form='formatted')
      do kbig=1,nzbig
      do ibig=1,nxbig
         if(kbig.lt.kstart.or.kbig.gt.kend.or.&
            ibig.lt.istart.or.ibig.gt.iend)then
            do j=1,nobjmax
               write(32,*)0.,0.
            enddo
         else
            k=kbig+1-kstart
            i=ibig+1-istart
            do j=1,nobjmax
               write(32,*)yi(j,i,k)+ystart,yf(j,i,k)+ystart
            enddo
         endif
      enddo
      enddo
      close(32)
      open(33,file='zizf.dat'   ,form='formatted')
      do jbig=1,nybig
      do ibig=1,nxbig
         if(jbig.lt.jstart.or.jbig.gt.jend.or.&
            ibig.lt.istart.or.ibig.gt.iend)then
            do k=1,nobjmax
               write(33,*)0.,0.
            enddo
         else
            j=jbig+1-jstart
            i=ibig+1-istart
            do k=1,nobjmax
               write(33,*)zi(k,i,j)+zstart,zf(k,i,j)+zstart
            enddo
         endif
      enddo
      enddo
      close(33)
   endif
!
   return
end subroutine write_geomcomplex


subroutine cylinder(epsi,nx,ny,nz,dx,yp,dz,remp)
   implicit none
!
   real(4),dimension(nx,ny,nz) :: epsi
   real(4),dimension(      ny) :: yp
   integer                     :: nx,ny,nz
   real(4)                     :: dx,dz
   real(4)                     :: xc,yc,r
   real(4)                     :: remp
   integer                     :: i,j,k,icheck
   real(4)                     :: x,y,x0,y0,x1,y1,xm,ym,pi,xa,xlx,xy,x00
   real(4)                     :: alpha

!
pi=acos(-1.)

do k=1,nz
do j=1,ny
do i=1,nx
   xm=(i-1)*dx 
   ym=yp(j)
   r=sqrt((xm-5.)*(xm-5.)+(ym-6.)*(ym-6.))
   if (r-0.5 >= 0.) cycle
   epsi(i,j,k)=remp
enddo
enddo
enddo

!
   return
end subroutine cylinder
!
