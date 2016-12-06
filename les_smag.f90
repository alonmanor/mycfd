!################################################################################
!This file is part of Incompact3d.
!
!Incompact3d
!Copyright (c) 2012 Eric Lamballais and Sylvain Laizet
!eric.lamballais@univ-poitiers.fr / sylvain.laizet@gmail.com
!
!    Incompact3d is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation.
!
!    Incompact3d is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with the code.  If not, see <http://www.gnu.org/licenses/>.
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!    We kindly request that you cite Incompact3d in your publications and 
!    presentations. The following citations are suggested:
!
!    1-Laizet S. & Lamballais E., 2009, High-order compact schemes for 
!    incompressible flows: a simple and efficient method with the quasi-spectral 
!    accuracy, J. Comp. Phys.,  vol 228 (15), pp 5989-6015
!
!    2-Laizet S. & Li N., 2011, Incompact3d: a powerful tool to tackle turbulence 
!    problems with up to 0(10^5) computational cores, Int. J. of Numerical 
!    Methods in Fluids, vol 67 (11), pp 1735-1757
!################################################################################

!********************************************************************
!
subroutine les_smag(ux1,uy1,uz1,ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1,&
     ux2,uy2,uz2,ta2,tb2,tc2,td2,te2,tf2,tg2,th2,ti2,tj2,di2,&
     ux3,uy3,uz3,ta3,tb3,tc3,td3,te3,tf3,tg3,th3,ti3,di3)
! 
!********************************************************************
USE param
USE variables
USE decomp_2d


implicit none

real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1
real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ux2,uy2,uz2 
real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ta2,tb2,tc2,td2,te2,tf2,tg2,th2,ti2,tj2,di2
real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: ux3,uy3,uz3
real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: ta3,tb3,tc3,td3,te3,tf3,tg3,th3,ti3,di3

integer :: ijk,nvect1,nvect2,nvect3,i,j,k
real(mytype) :: x,y,z



nvect1=xsize(1)*xsize(2)*xsize(3)
nvect2=ysize(1)*ysize(2)*ysize(3)
nvect3=zsize(1)*zsize(2)*zsize(3)
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! LES MODULE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!*****************************************
! calculate the rate of strees tensor Sij. 
! preformed independant of specific SFS model
! X   >  Y  >  Z
!*****************************************
call derx (td1,ux1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0) !dudx -> td = Sxx
call derx (te1,uy1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1) !dvdx -> te
call derx (tf1,uz1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1) !dwdx -> tf
!Move to Y-PENCILS
call transpose_x_to_y(ux1,ux2)
call transpose_x_to_y(uy1,uy2)
call transpose_x_to_y(uz1,uz2)
call transpose_x_to_y(td1,td2)
call transpose_x_to_y(te1,te2)
call transpose_x_to_y(tf1,tf2)
call dery (tg2,ux2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1) !dudy -> tg
te2(:,:,:) = 0.5 * (te2(:,:,:) + tg2(:,:,:) ) ! te = 0.5 * (dvdx + dudy) = Sxy = Syx
!now tg is free
call dery (tg2,uy2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0) !dvdy -> tg = Syy
call dery (th2,uz2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1) !dwdy -> th
!Move to Z-PENCILS
call transpose_y_to_z(ux2,ux3)
call transpose_y_to_z(uy2,uy3)
call transpose_y_to_z(uz2,uz3)
call transpose_y_to_z(td2,td3)
call transpose_y_to_z(te2,te3)
call transpose_y_to_z(tf2,tf3)
call transpose_y_to_z(tg2,tg3)
call transpose_y_to_z(th2,th3)
call derz (ti3,ux3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1) !dudz -> ti
tf3(:,:,:) = 0.5 * (tf3(:,:,:) + ti3(:,:,:) ) ! tf = 0.5 * (dwdx + dudz) = Sxz = Szx
!now ti is free
call derz (ti3,uy3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1) !dvdz -> ti
th3(:,:,:) = 0.5 * (th3(:,:,:) + ti3(:,:,:) ) ! th = 0.5 * (dwdy + dvdz) = Syz = Szy
!ti is free again
call derz (ti3,uz3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0) !dwdz -> ti = Szz
!now all the Sij elements are calculated
! Sij = Sji, hence we need only 6 arrays: td3, te3, tf3, tg, th, ti
! Sxx   =   dudx   =      td      Sxy = 0.5*(dudy+dvdx) = te   Sxz = 0.5*(dudz+dwdx) = tf
! Syx = 0.5*(dudy+dvdx) = te      Syy =      dvdy       = tg   Syz = 0.5*(dvdz+dwdy) = th
! Szx = 0.5*(dudz+dwdx) = tf      Szy = 0.5*(dvdz+dwdy) = th   Szz =        dwdz     = ti

!*****************************************
! calculate eddy viscousity
!*****************************************
! calculate sqrt(2 |Sij Sij|) - rate of stress tensor magnitude. store in tj
! note that for most SFS models, tj has to be kept for the scalar treatment
do ijk=1,nvect3
  tj3(ijk,1,1) = 2.0*(td3(ijk,1,1)**2.0+tg3(ijk,1,1)**2.0+ti3(ijk,1,1)**2.0)
  tj3(ijk,1,1) = tj3(ijk,1,1)+te3(ijk,1,1)**2.0+tf3(ijk,1,1)**2.0+th3(ijk,1,1)**2.0
  tj3(ijk,1,1) = C_smag * sqrt( tj3(ijk,1,1) )
enddo 
do k=1,zsize(3)
do i=1,zsize(1)
do j=1,zsize(2)
  ! to get the final eddy viscousity, multiply by deltax^2.
  ! deltax is evaluated as (dx*dy*dz)^0.3333
  tj3(i,j,k)=tj3(i,j,k)*( dx * dz * (yp(j+1)-yp(j)) )**(0.666666)
  if (j.eq.zsize(2)) then
    tj3(i,j,k)=tj3(i,j,k)*( dx * dz * (yp(j)-yp(j-1)) )**(0.666666) 
  endif
enddo
enddo
enddo
!*****************************************
! calculate stress tensor tau_ij and drivate
! Z  >  Y  >  X
!*****************************************
do ijk=1,nvect3
    tl3(ijk,1,1) = -2.0 * tj3(ijk,1,1) * ti3(ijk,1,1) !tau_zz = -2*ni*S_zz
enddo
call derz (tk3,tl3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1) ! tk = dtau_zz/dz
do ijk=1,nvect3
    ti3(ijk,1,1) = -2.0 * tj3(ijk,1,1) * th3(ijk,1,1) !tau_yz = -2*ni*S_yz
enddo
   call derz (tl3,ti3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0) ! tl = dtau_yz/dz    
do ijk=1,nvect3
    tn3(ijk,1,1) = -2.0 * tj3(ijk,1,1) * tf3(ijk,1,1) !tau_xz = -2*ni*S_xz
enddo
call derz (ti3,tn3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0) ! ti = dtau_zz/dz


end subroutine les_smag

