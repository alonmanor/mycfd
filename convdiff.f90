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
subroutine convdiff(ux1,uy1,uz1,ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1,&
     ux2,uy2,uz2,ta2,tb2,tc2,td2,te2,tf2,tg2,th2,ti2,tj2,di2,&
     ux3,uy3,uz3,ta3,tb3,tc3,td3,te3,tf3,tg3,th3,ti3,di3,&
     duxdx1,duxdx2,duxdx3,duydx1,duydx2,duydx3,duzdx1,duzdx2,duzdx3,&
     duxdy1,duxdy2,duxdy3,            duydy2,duydy3,            duzdy2,duzdy3,&
     duxdz1,duxdz2,duxdz3,            duydz2,duydz3,                        duzdz3,&
     xnu_sgs1,xnu_sgs2,xnu_sgs3,les_a1,les_a2,les_a3,les_b1,les_b2,&
     div_tau_x1,div_tau_y1,div_tau_z1,&
     div_tau_x2,div_tau_y2,div_tau_z2,&
     div_tau_x3,div_tau_y3,div_tau_z3,&
     k_sgs1    ,k_sgs2   ,k_sgs3,&
     e_svm_x1  ,e_svm_y1  ,e_svm_z1,&
     e_svm_x2  ,e_svm_y2  ,e_svm_z2,&
     e_svm_x3  ,e_svm_y3  ,e_svm_z3, phi1 )
! 
!********************************************************************
USE param
USE variables
USE decomp_2d
use mymath


implicit none


real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,phi1
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1
real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ux2,uy2,uz2 
real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ta2,tb2,tc2,td2,te2,tf2,tg2,th2,ti2,tj2,di2
real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: ux3,uy3,uz3
real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: ta3,tb3,tc3,td3,te3,tf3,tg3,th3,ti3,di3
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: duxdx1,duydx1,duzdx1,duxdy1,duxdz1,les_a1,les_b1
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: div_tau_x1,div_tau_y1,div_tau_z1,xnu_sgs1
real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: duxdx2,duydx2,duzdx2,les_a2,les_b2
real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: duxdy2,duydy2,duzdy2,duxdz2,duydz2
real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: div_tau_x2,div_tau_y2,div_tau_z2,xnu_sgs2
real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: duxdx3,duydx3,duzdx3,les_a3
real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: duxdy3,duydy3,duzdy3
real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: duxdz3,duydz3,duzdz3
real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: div_tau_x3,div_tau_y3,div_tau_z3,xnu_sgs3
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: k_sgs1,e_svm_x1,e_svm_y1,e_svm_z1
real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: k_sgs2,e_svm_x2,e_svm_y2,e_svm_z2
real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: k_sgs3,e_svm_x3,e_svm_y3,e_svm_z3

integer :: ijk,nvect1,nvect2,nvect3,i,j,k
real(mytype) :: x,y,z
integer,dimension(3) :: vdebug


nvect1=xsize(1)*xsize(2)*xsize(3)
nvect2=ysize(1)*ysize(2)*ysize(3)
nvect3=zsize(1)*zsize(2)*zsize(3)

if (iskew==0) then !UROTU!
!WORK X-PENCILS
   call derx (ta1,uy1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
   call derx (tb1,uz1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
   call transpose_x_to_y(ux1,ux2)
   call transpose_x_to_y(uy1,uy2)
   call transpose_x_to_y(uz1,uz2)
   call transpose_x_to_y(ta1,ta2)
   call transpose_x_to_y(tb1,tb2)
!WORK Y-PENCILS
   call dery (tc2,ux2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1) 
   call dery (td2,uz2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1) 
   call transpose_y_to_z(ux2,ux3)
   call transpose_y_to_z(uy2,uy3)
   call transpose_y_to_z(uz2,uz3)
   call transpose_y_to_z(ta2,ta3)
   call transpose_y_to_z(tb2,tb3)
   call transpose_y_to_z(tc2,tc3)
   call transpose_y_to_z(td2,td3)
!WORK Z-PENCILS
   call derz (te3,ux3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
   call derz (tf3,uy3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
   do ijk=1,nvect3
      ta3(ijk,1,1)=uz3(ijk,1,1)*(te3(ijk,1,1)-tb3(ijk,1,1))-&
           uy3(ijk,1,1)*(ta3(ijk,1,1)-tc3(ijk,1,1))
      tb3(ijk,1,1)=ux3(ijk,1,1)*(ta3(ijk,1,1)-tc3(ijk,1,1))-&
           uz3(ijk,1,1)*(td3(ijk,1,1)-tf3(ijk,1,1))
      tc3(ijk,1,1)=uy3(ijk,1,1)*(td3(ijk,1,1)-tf3(ijk,1,1))-&
           ux3(ijk,1,1)*(te3(ijk,1,1)-tb3(ijk,1,1))
   enddo
else !SKEW!
!WORK X-PENCILS
   do ijk=1,nvect1
      ta1(ijk,1,1)=ux1(ijk,1,1)*ux1(ijk,1,1)
      tb1(ijk,1,1)=ux1(ijk,1,1)*uy1(ijk,1,1)
      tc1(ijk,1,1)=ux1(ijk,1,1)*uz1(ijk,1,1)
   enddo
   call derx (td1,ta1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
   call derx (te1,tb1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0)
   call derx (tf1,tc1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0)
   call derx (ta1,ux1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0)
   call derx (tb1,uy1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
   call derx (tc1,uz1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
   if (iles > 0) then
        do ijk=1,nvect1
            duxdx1(ijk,1,1) = ta1(ijk,1,1)
            duydx1(ijk,1,1) = tb1(ijk,1,1)
            duzdx1(ijk,1,1) = tc1(ijk,1,1)
        enddo
        if ((iles.eq.4).or.(iles.eq.5)) then
			call k_sgs_circular(ux1,uy1,uz1,k_sgs1,xsize(1),&
				xsize(2),xsize(3),e_svm_y1,1)
        endif
   endif

   do ijk=1,nvect1
      ta1(ijk,1,1)=0.5*td1(ijk,1,1)+0.5*ux1(ijk,1,1)*ta1(ijk,1,1)
      tb1(ijk,1,1)=0.5*te1(ijk,1,1)+0.5*ux1(ijk,1,1)*tb1(ijk,1,1)
      tc1(ijk,1,1)=0.5*tf1(ijk,1,1)+0.5*ux1(ijk,1,1)*tc1(ijk,1,1)      
   enddo

   call transpose_x_to_y(ux1,ux2)
   call transpose_x_to_y(uy1,uy2)
   call transpose_x_to_y(uz1,uz2)
   call transpose_x_to_y(ta1,ta2)
   call transpose_x_to_y(tb1,tb2)
   call transpose_x_to_y(tc1,tc2)
   if (iles > 0) then
       call transpose_x_to_y(duxdx1,duxdx2)
       call transpose_x_to_y(duydx1,duydx2)
       call transpose_x_to_y(duzdx1,duzdx2)
       if ((iles.eq.4).or.(iles.eq.5)) then
		call transpose_y_to_z(k_sgs1,k_sgs2)
       endif
   endif
!WORK Y-PENCILS
   do ijk=1,nvect2
      td2(ijk,1,1)=ux2(ijk,1,1)*uy2(ijk,1,1)
      te2(ijk,1,1)=uy2(ijk,1,1)*uy2(ijk,1,1)
      tf2(ijk,1,1)=uz2(ijk,1,1)*uy2(ijk,1,1)
   enddo
   call dery (tg2,td2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0) 
   call dery (th2,te2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
   call dery (ti2,tf2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0) 
   call dery (td2,ux2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1) 
   call dery (te2,uy2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0)
   call dery (tf2,uz2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)    
   if (iles > 0) then
        do ijk=1,nvect2
            duxdy2(ijk,1,1) = td2(ijk,1,1)
            duydy2(ijk,1,1) = te2(ijk,1,1)
            duzdy2(ijk,1,1) = tf2(ijk,1,1)
        enddo
   endif
   do ijk=1,nvect2
      ta2(ijk,1,1)=ta2(ijk,1,1)+0.5*tg2(ijk,1,1)+0.5*uy2(ijk,1,1)*td2(ijk,1,1)
      tb2(ijk,1,1)=tb2(ijk,1,1)+0.5*th2(ijk,1,1)+0.5*uy2(ijk,1,1)*te2(ijk,1,1)
      tc2(ijk,1,1)=tc2(ijk,1,1)+0.5*ti2(ijk,1,1)+0.5*uy2(ijk,1,1)*tf2(ijk,1,1)      
   enddo
   call transpose_y_to_z(ux2,ux3)
   call transpose_y_to_z(uy2,uy3)
   call transpose_y_to_z(uz2,uz3)
   call transpose_y_to_z(ta2,ta3)
   call transpose_y_to_z(tb2,tb3)
   call transpose_y_to_z(tc2,tc3)
   if (iles > 0) then
       call transpose_y_to_z(duxdy2,duxdy3)
       call transpose_y_to_z(duydy2,duydy3)
       call transpose_y_to_z(duzdy2,duzdy3)
       call transpose_y_to_z(duxdx2,duxdx3)
       call transpose_y_to_z(duydx2,duydx3)
       call transpose_y_to_z(duzdx2,duzdx3)
       if ((iles.eq.4).or.(iles.eq.5)) then
		call transpose_y_to_z(k_sgs2,k_sgs3)
       endif
   endif
!WORK Z-PENCILS
   do ijk=1,nvect3
      td3(ijk,1,1)=ux3(ijk,1,1)*uz3(ijk,1,1)
      te3(ijk,1,1)=uy3(ijk,1,1)*uz3(ijk,1,1)
      tf3(ijk,1,1)=uz3(ijk,1,1)*uz3(ijk,1,1)
   enddo
   call derz (tg3,td3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0)
   call derz (th3,te3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0)
   call derz (ti3,tf3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
   call derz (td3,ux3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
   call derz (te3,uy3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
   call derz (tf3,uz3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0)
   if (iles > 0) then
        do ijk=1,nvect3
            duxdz3(ijk,1,1) = td3(ijk,1,1)
            duydz3(ijk,1,1) = te3(ijk,1,1)
            duzdz3(ijk,1,1) = tf3(ijk,1,1)
        enddo
        if ((iles.eq.4).or.(iles.eq.5)) then
			call k_sgs_circular(ux3,uy3,uz3,k_sgs3,zsize(1),&
				zsize(2),zsize(3),e_svm_y3,3)
        endif
   endif
   do ijk=1,nvect3
      ta3(ijk,1,1)=ta3(ijk,1,1)+0.5*tg3(ijk,1,1)+0.5*uz3(ijk,1,1)*td3(ijk,1,1)
      tb3(ijk,1,1)=tb3(ijk,1,1)+0.5*th3(ijk,1,1)+0.5*uz3(ijk,1,1)*te3(ijk,1,1)
      tc3(ijk,1,1)=tc3(ijk,1,1)+0.5*ti3(ijk,1,1)+0.5*uz3(ijk,1,1)*tf3(ijk,1,1)   
   enddo
endif
!ALL THE CONVECTIVE TERMS ARE IN TA3, TB3 and TC3

td3(:,:,:)=ta3(:,:,:)
te3(:,:,:)=tb3(:,:,:)
tf3(:,:,:)=tc3(:,:,:)
if (iles.eq.1) then !calculate nu_SGS for Smagorinsky SGS model
	call xnu_smag(duxdx3,duydx3,duzdx3,duxdy3,duydy3,duzdy3,&
     duxdz3,duydz3,duzdz3,xnu_sgs3)
     !call showval3(xnu_sgs3,5,5,5)
endif
if (iles.eq.2) then !calculate nu_SGS for sigma SGS model
	if (itime>5) then !march a few timesteps with smagorinsky, to allow legal G_ij to form
		call xnu_sigma(duxdx3,duydx3,duzdx3,duxdy3,duydy3,duzdy3,&
		 duxdz3,duydz3,duzdz3,xnu_sgs3)
    else
		call xnu_smag(duxdx3,duydx3,duzdx3,duxdy3,duydy3,duzdy3,&
		 duxdz3,duydz3,duzdz3,xnu_sgs3)
    endif
endif
if (iles.eq.4) then !SVM using largest EV of S_ij
	call e_svm_s(duxdx3,duydx3,duzdx3,duxdy3,duydy3,duzdy3,&
		 duxdz3,duydz3,duzdz3,e_svm_x3,e_svm_y3,e_svm_z3)
endif
if (iles.eq.5) then !SVM using vorticity
	call e_svm_v(duxdx3,duydx3,duzdx3,duxdy3,duydy3,duzdy3,&
		 duxdz3,duydz3,duzdz3,e_svm_x3,e_svm_y3,e_svm_z3)
endif
if ((iles.eq.4).or.(iles.eq.5)) then
	call transpose_x_to_y(e_svm_x3,e_svm_x2)
	call transpose_x_to_y(e_svm_y3,e_svm_y2)
	call transpose_x_to_y(e_svm_z3,e_svm_z2)
	call transpose_y_to_z(e_svm_x2,e_svm_x1)
	call transpose_y_to_z(e_svm_y2,e_svm_y1)
	call transpose_y_to_z(e_svm_z2,e_svm_z1)
	!partial SGS TKE in k_sgs1
	call k_sgs_circular(ux1,uy1,uz1,k_sgs1,xsize(1),&
				xsize(2),xsize(3),e_svm_y1,1)
	call transpose_x_to_y(k_sgs1,k_sgs2)
	call transpose_y_to_z(k_sgs2,k_sgs3)
	call k_sgs_circular(ux3,uy3,uz3,k_sgs3,zsize(1),&
				zsize(2),zsize(3),e_svm_y3,3)
	!full SGS TKE in k_sgs3
endif

	
!DIFFUSIVE TERMS IN Z
call derzz (ta3,ux3,di3,sz,sfzp,sszp,swzp,zsize(1),zsize(2),zsize(3),1)
call derzz (tb3,uy3,di3,sz,sfzp,sszp,swzp,zsize(1),zsize(2),zsize(3),1)
call derzz (tc3,uz3,di3,sz,sfz ,ssz ,swz ,zsize(1),zsize(2),zsize(3),0)



if ((iles.eq.1).or.(iles.eq.2)) then
    les_a3(:,:,:) = -2.0*xnu_sgs3(:,:,:) * 0.5*(duxdz3(:,:,:) + duzdx3(:,:,:) ) ! tau_xz = -2.0*nu*Sxz 
    call derz(div_tau_x3, les_a3 ,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0) ! div_tau_x=d(tau_xz)/dz
    les_a3(:,:,:) = -2.0*xnu_sgs3(:,:,:) * 0.5*(duydz3(:,:,:) + duzdy3(:,:,:) ) ! tau_yz = -2.0*nu*Syz 
    call derz(div_tau_y3, les_a3 ,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0) ! div_tau_y=d(tau_yz)/dz
    les_a3(:,:,:) = -2.0*xnu_sgs3(:,:,:) * duzdz3(:,:,:)                                   ! tau_zz = -2.0*nu*Szz 
    call derz(div_tau_z3, les_a3 ,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1) ! div_tau_z=d(tau_zz)/dz
endif
if ((iles.eq.4).or.(iles.eq.5)) then
	les_a3(:,:,:) =  - k_sgs3(:,:,:) * (      e_svm_x3(:,:,:) * e_svm_z3(:,:,:)) !tau_xz = -K(e_x*e_z)
	call derz(div_tau_x3, les_a3 ,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0) ! div_tau_x=d(tau_xz)/dz
	les_a3(:,:,:) =  - k_sgs3(:,:,:) * (      e_svm_y3(:,:,:) * e_svm_z3(:,:,:)) !tau_yz = -K(e_y*e_z)
	call derz(div_tau_y3, les_a3 ,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0) ! div_tau_y=d(tau_yz)/dz
	les_a3(:,:,:) =    k_sgs3(:,:,:) * (1.0 - e_svm_z3(:,:,:) * e_svm_z3(:,:,:)) !tau_zz = -K(1 - e_z*e_z)
	call derz(div_tau_z3, les_a3 ,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1) ! div_tau_z=d(tau_zz)/dz
endif


!WORK Y-PENCILS
call transpose_z_to_y(ta3,ta2)
call transpose_z_to_y(tb3,tb2)
call transpose_z_to_y(tc3,tc2)
call transpose_z_to_y(td3,td2)
call transpose_z_to_y(te3,te2)
call transpose_z_to_y(tf3,tf2)

if (iles > 0) then
	call transpose_z_to_y(div_tau_x3,div_tau_x2)
    call transpose_z_to_y(div_tau_y3,div_tau_y2)
    call transpose_z_to_y(div_tau_z3,div_tau_z2)
    if ((iles.eq.4).or.(iles.eq.5)) then
		call transpose_z_to_y(k_sgs3,k_sgs2)
    endif
    if ((iles.eq.1).or.(iles.eq.2)) then
		call transpose_z_to_y(xnu_sgs3,xnu_sgs2)
		call transpose_z_to_y(duydz3,duydz2)
		call transpose_z_to_y(duxdz3,duxdz2)
    endif
endif

tg2(:,:,:)=td2(:,:,:)
th2(:,:,:)=te2(:,:,:)
ti2(:,:,:)=tf2(:,:,:)

!DIFFUSIVE TERMS IN Y
!-->for ux
if (istret.ne.0) then 
   call deryy (td2,ux2,di2,sy,sfyp,ssyp,swyp,ysize(1),ysize(2),ysize(3),1)
   call dery (te2,ux2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
   do k=1,ysize(3)
   do j=1,ysize(2)
   do i=1,ysize(1)
      td2(i,j,k)=td2(i,j,k)*pp2y(j)-pp4y(j)*te2(i,j,k)
   enddo
   enddo
   enddo
else
   call deryy (td2,ux2,di2,sy,sfyp,ssyp,swyp,ysize(1),ysize(2),ysize(3),1) 
endif
!-->for uy
if (istret.ne.0) then 
   call deryy (te2,uy2,di2,sy,sfy,ssy,swy,ysize(1),ysize(2),ysize(3),0)
   call dery (tf2,uy2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0)
   do k=1,ysize(3)
   do j=1,ysize(2)
   do i=1,ysize(1)
      te2(i,j,k)=te2(i,j,k)*pp2y(j)-pp4y(j)*tf2(i,j,k)
   enddo
   enddo
   enddo
else
   call deryy (te2,uy2,di2,sy,sfy,ssy,swy,ysize(1),ysize(2),ysize(3),0) 
endif
!-->for uz
if (istret.ne.0) then 
   call deryy (tf2,uz2,di2,sy,sfyp,ssyp,swyp,ysize(1),ysize(2),ysize(3),1)
   call dery (tj2,uz2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
   do k=1,ysize(3)
   do j=1,ysize(2)
   do i=1,ysize(1)
      tf2(i,j,k)=tf2(i,j,k)*pp2y(j)-pp4y(j)*tj2(i,j,k)
   enddo
   enddo
   enddo
else
   call deryy (tf2,uz2,di2,sy,sfyp,ssyp,swyp,ysize(1),ysize(2),ysize(3),1) 
endif

if ((iles.eq.1).or.(iles.eq.2)) then
    les_a2(:,:,:) = -2.0*xnu_sgs2(:,:,:) *0.5* (duxdy2(:,:,:) + duydx2(:,:,:) ) ! tau_xy = -2.0*nu*Sxy 
    call dery(les_b2, les_a2 ,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0) ! d(tau_xy)/dy
    div_tau_x2(:,:,:) = div_tau_x2(:,:,:) + les_b2(:,:,:) !div_tau_x = div_tau_x + d(tau_xy)/dy
    les_a2(:,:,:) = -2.0*xnu_sgs2(:,:,:) * duydy2(:,:,:)                            ! tau_yy = -2.0*nu*Syy 
    call dery(les_b2, les_a2 ,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1) ! d(tau_yy)/dy
    div_tau_y2(:,:,:) = div_tau_y2(:,:,:) + les_b2(:,:,:) !div_tau_y = div_tau_y + d(tau_yy)/dy
    les_a2(:,:,:) = -2.0*xnu_sgs2(:,:,:)*0.5*(duzdy2(:,:,:)+duydz2(:,:,:))! tau_zy = -2.0*nu*Szy 
    call dery(les_b2, les_a2 ,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0) ! d(tau_zy)/dy
    div_tau_z2(:,:,:) = div_tau_z2(:,:,:) + les_b2(:,:,:) !div_tau_z = div_tau_z + d(tau_zy)/dy
endif
if ((iles.eq.4).or.(iles.eq.5)) then
	les_a2(:,:,:) = - k_sgs2(:,:,:) * (      e_svm_x2(:,:,:) * e_svm_y2(:,:,:)) ! tau_xy = -K(e_x*e_y)
    call dery(les_b2, les_a2 ,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0) ! d(tau_xy)/dy
    div_tau_x2(:,:,:) = div_tau_x2(:,:,:) + les_b2(:,:,:) !div_tau_x = div_tau_x + d(tau_xy)/dy
    les_a2(:,:,:) =   k_sgs2(:,:,:) * (1.0 - e_svm_y2(:,:,:) * e_svm_y2(:,:,:)) !  tau_yy = K(1-e_y*e_y)
    call dery(les_b2, les_a2 ,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1) ! d(tau_yy)/dy
    div_tau_y2(:,:,:) = div_tau_y2(:,:,:) + les_b2(:,:,:) !div_tau_y = div_tau_y + d(tau_yy)/dy
    les_a2(:,:,:) =  - k_sgs2(:,:,:) * (      e_svm_z2(:,:,:) * e_svm_y2(:,:,:)) ! tau_zy = -K(e_z*e_y)
    call dery(les_b2, les_a2 ,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0) ! d(tau_zy)/dy
    div_tau_z2(:,:,:) = div_tau_z2(:,:,:) + les_b2(:,:,:) !div_tau_z = div_tau_z + d(tau_zy)/dy
endif

ta2(:,:,:)=ta2(:,:,:)+td2(:,:,:)
tb2(:,:,:)=tb2(:,:,:)+te2(:,:,:)
tc2(:,:,:)=tc2(:,:,:)+tf2(:,:,:)

!WORK X-PENCILS
call transpose_y_to_x(ta2,ta1)
call transpose_y_to_x(tb2,tb1)
call transpose_y_to_x(tc2,tc1) !diff
call transpose_y_to_x(tg2,td1)
call transpose_y_to_x(th2,te1)
call transpose_y_to_x(ti2,tf1) !conv

if (iles > 0) then
    call transpose_y_to_x(div_tau_x2,div_tau_x1)
    call transpose_y_to_x(div_tau_y2,div_tau_y1)
    call transpose_y_to_x(div_tau_z2,div_tau_z1)
    if ((iles.eq.4).or.(iles.eq.5)) then
		call transpose_y_to_x(k_sgs2,k_sgs1)
    endif
    if ((iles.eq.1).or.(iles.eq.2)) then
		call transpose_y_to_x(xnu_sgs2,xnu_sgs1)
		call transpose_y_to_x(duxdy2,duxdy1)
		call transpose_y_to_x(duxdz2,duxdz1)
    endif
endif

tg1(:,:,:)=td1(:,:,:)
th1(:,:,:)=te1(:,:,:)
ti1(:,:,:)=tf1(:,:,:)

!DIFFUSIVE TERMS IN X
call derxx (td1,ux1,di1,sx,sfx ,ssx ,swx ,xsize(1),xsize(2),xsize(3),0)
call derxx (te1,uy1,di1,sx,sfxp,ssxp,swxp,xsize(1),xsize(2),xsize(3),1)
call derxx (tf1,uz1,di1,sx,sfxp,ssxp,swxp,xsize(1),xsize(2),xsize(3),1)

ta1(:,:,:)=ta1(:,:,:)+td1(:,:,:)
tb1(:,:,:)=tb1(:,:,:)+te1(:,:,:)
tc1(:,:,:)=tc1(:,:,:)+tf1(:,:,:)



!if (nrank==1) print *,'WARNING rotating channel',itime
!tg1(:,:,:)=tg1(:,:,:)-2./18.*uy1(:,:,:)
!th1(:,:,:)=th1(:,:,:)-2./18.*ux1(:,:,:)

if ((iles.eq.1).or.(iles.eq.2)) then
    les_a1(:,:,:) = -2.0*xnu_sgs1(:,:,:) * duxdx1(:,:,:)  ! tau_xx = -2.0*nu*Sxx 
    call derx(les_b1, les_a1 ,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1) ! d(tau_xx)/dx
    div_tau_x1(:,:,:) = div_tau_x1(:,:,:) + les_b1(:,:,:) !div_tau_x = div_tau_x + d(tau_xx)/dx
    les_a1(:,:,:) = -2.0*xnu_sgs1(:,:,:) * 0.5 *(duxdy1(:,:,:)+duydx1(:,:,:))   ! tau_xy = -2.0*nu*Sxy 
    call derx(les_b1, les_a1 ,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0) ! d(tau_xy)/dx
    div_tau_y1(:,:,:) = div_tau_y1(:,:,:) + les_b1(:,:,:) !div_tau_y = div_tau_y + d(tau_xy)/dx
    les_a1(:,:,:) = -2.0*xnu_sgs1(:,:,:) * 0.5* (duzdx1(:,:,:) + duxdz1(:,:,:) )! tau_xz = -2.0*nu*Sxz 
    call derx(les_b1, les_a1 ,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0) ! d(tau_xz)/dx
    div_tau_z1(:,:,:) = div_tau_z1(:,:,:) + les_b1(:,:,:) !div_tau_z = div_tau_z + d(tau_xz)/dx
endif
if ((iles.eq.4).or.(iles.eq.5)) then
    les_a1(:,:,:) =   k_sgs1(:,:,:) * (1.0 - e_svm_y1(:,:,:) * e_svm_y1(:,:,:))   ! tau_xx = K(1-e_x*e_x)
    call derx(les_b1, les_a1 ,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1) ! d(tau_xx)/dx
    div_tau_x1(:,:,:) = div_tau_x1(:,:,:) + les_b1(:,:,:) !div_tau_x = div_tau_x + d(tau_xx)/dx
    les_a1(:,:,:) = - k_sgs1(:,:,:) * (      e_svm_x1(:,:,:) * e_svm_y1(:,:,:))   ! tau_xy = -K(e_x*e_y)
    call derx(les_b1, les_a1 ,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0) ! d(tau_xy)/dx
    div_tau_y1(:,:,:) = div_tau_y1(:,:,:) + les_b1(:,:,:) !div_tau_y = div_tau_y + d(tau_xy)/dx
    les_a1(:,:,:) = - k_sgs1(:,:,:) * (      e_svm_x1(:,:,:) * e_svm_z1(:,:,:))   ! tau_xz = -K(e_x*e_z)
    call derx(les_b1, les_a1 ,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0) ! d(tau_xz)/dx
    div_tau_z1(:,:,:) = div_tau_z1(:,:,:) + les_b1(:,:,:) !div_tau_z = div_tau_z + d(tau_xz)/dx
endif

!FINAL SUM: DIFF TERMS + CONV TERMS
ta1(:,:,:)=xnu*ta1(:,:,:)-tg1(:,:,:)
tb1(:,:,:)=xnu*tb1(:,:,:)-th1(:,:,:)
tc1(:,:,:)=xnu*tc1(:,:,:)-ti1(:,:,:)
do i=1,xsize(1)
	do j=1,xsize(2)
		do k=1,xsize(3)
			tb1(i,j,k) = tb1(i,j,k) + buoy_param*phi1(i,j,k)
		enddo
	enddo
enddo
if (iles > 0) then
	ta1(:,:,:)=-div_tau_x1(:,:,:) + ta1(:,:,:)
    tb1(:,:,:)=-div_tau_y1(:,:,:) + tb1(:,:,:)
    tc1(:,:,:)=-div_tau_z1(:,:,:) + tc1(:,:,:)
endif
if (itype.eq.10) then
    ta1(:,:,:)=ta1(:,:,:)+dpdx_drive
endif
if ((itype.eq.11).and.(damp.gt.0.0)) then !damping layer
	call damping(ta1,tb1,ux1,uy1)
endif


end subroutine convdiff


!************************************************************
!
subroutine damping(ta1,tb1,ux1,uy1)
!
!************************************************************

USE param
USE variables
USE decomp_2d

implicit none

real(mytype),dimension(xsize(1),xsize(2),xsize(3)),intent(inout) :: ta1,tb1
real(mytype),dimension(xsize(1),xsize(2),xsize(3)),intent(in) :: ux1,uy1
integer i,j,k
real(mytype) :: ux_geos, uy_geos

ux_geos = u_geos
uy_geos = 0.0
do i=1,xsize(1)
	do k=1,xsize(3)
		do j=xstart(2),xend(2)
			ta1(i,j,k) = ta1(i,j,k) + f_damp(j)*(ux_geos - ux1(i,j,k))
			tb1(i,j,k) = tb1(i,j,k) + f_damp(j)*(uy_geos - uy1(i,j,k))
		enddo
	enddo
enddo

return

end subroutine damping


!************************************************************
!
subroutine scalar(ux1,uy1,uz1,phi1,phis1,phiss1,di1,ta1,tb1,tc1,td1,&
     uy2,uz2,phi2,di2,ta2,tb2,tc2,td2,uz3,phi3,di3,ta3,tb3,epsi)
!
!************************************************************

USE param
USE variables
USE decomp_2d

implicit none

real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,phi1,phis1,&
                                              phiss1,di1,ta1,tb1,tc1,td1,epsi
real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: uy2,uz2,phi2,di2,ta2,tb2,tc2,td2
real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: uz3,phi3,di3,ta3,tb3

integer :: ijk,nvect1,nvect2,nvect3,i,j,k,nxyz
real(mytype) :: x,y,z

nvect1=xsize(1)*xsize(2)*xsize(3)
nvect2=ysize(1)*ysize(2)*ysize(3)
nvect3=zsize(1)*zsize(2)*zsize(3)

!X PENCILS
do ijk=1,nvect1
   ta1(ijk,1,1)=ux1(ijk,1,1)*phi1(ijk,1,1)
enddo
call derx (tb1,ta1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0)
call derxx (ta1,phi1,di1,sx,sfxp,ssxp,swxp,xsize(1),xsize(2),xsize(3),1)

call transpose_x_to_y(phi1,phi2)
call transpose_x_to_y(uy1,uy2)
call transpose_x_to_y(uz1,uz2)

!Y PENCILS
do ijk=1,nvect2
   ta2(ijk,1,1)=uy2(ijk,1,1)*phi2(ijk,1,1)
enddo
call dery (tb2,ta2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0)
if (istret.ne.0) then 
   call deryy (ta2,phi2,di2,sy,sfyp,ssyp,swyp,ysize(1),ysize(2),ysize(3),1)
   call dery (tc2,phi2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0)
   do k=1,ysize(3)
   do j=1,ysize(2)
   do i=1,ysize(1)
      ta2(i,j,k)=ta2(i,j,k)*pp2y(j)-pp4y(j)*tc2(i,j,k)
   enddo
   enddo
   enddo
else
   call deryy (ta2,phi2,di2,sy,sfyp,ssyp,swyp,ysize(1),ysize(2),ysize(3),1) 
endif

call transpose_y_to_z(phi2,phi3)
call transpose_y_to_z(uz2,uz3)

!Z PENCILS
do ijk=1,nvect3
   ta3(ijk,1,1)=uz3(ijk,1,1)*phi3(ijk,1,1)
enddo
call derz (tb3,ta3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0)
call derzz (ta3,phi3,di3,sz,sfzp,sszp,swzp,zsize(1),zsize(2),zsize(3),1)

call transpose_z_to_y(ta3,tc2)
call transpose_z_to_y(tb3,td2)

!Y PENCILS ADD TERMS
do ijk=1,nvect2
   tc2(ijk,1,1)=tc2(ijk,1,1)+ta2(ijk,1,1)
   td2(ijk,1,1)=td2(ijk,1,1)+tb2(ijk,1,1)
enddo

call transpose_y_to_x(tc2,tc1)
call transpose_y_to_x(td2,td1)

!X PENCILS ADD TERMS
do ijk=1,nvect1
   ta1(ijk,1,1)=ta1(ijk,1,1)+tc1(ijk,1,1) !SECOND DERIVATIVE
   tb1(ijk,1,1)=tb1(ijk,1,1)+td1(ijk,1,1) !FIRST DERIVATIVE
enddo
 
do ijk=1,nvect1
   ta1(ijk,1,1)=xnu/sc*ta1(ijk,1,1)-tb1(ijk,1,1) 
enddo

!TIME ADVANCEMENT
nxyz=xsize(1)*xsize(2)*xsize(3)  

if ((nscheme.eq.1).or.(nscheme.eq.2)) then
   if ((nscheme.eq.1.and.itime.eq.1.and.ilit.eq.0).or.&
        (nscheme.eq.2.and.itr.eq.1)) then
      do ijk=1,nxyz
         phi1(ijk,1,1)=gdt(itr)*ta1(ijk,1,1)+phi1(ijk,1,1)
         phis1(ijk,1,1)=ta1(ijk,1,1)          
      enddo
   else
      do ijk=1,nxyz
         phi1(ijk,1,1)=adt(itr)*ta1(ijk,1,1)+bdt(itr)*phis1(ijk,1,1)+phi1(ijk,1,1)
         phis1(ijk,1,1)=ta1(ijk,1,1)          
      enddo
endif
endif

if (nscheme.eq.3) then 
if (nrank==0) print *,'Not ready'
stop 
endif

if (nscheme==4) then
   if ((itime.eq.1).and.(ilit.eq.0)) then
      if (nrank==0) print *,'start with Euler',itime
      do ijk=1,nxyz !start with Euler
         phi1(ijk,1,1)=dt*ta1(ijk,1,1)+phi1(ijk,1,1)
         phis1(ijk,1,1)=ta1(ijk,1,1)          
      enddo
   else
      if  ((itime.eq.2).and.(ilit.eq.0)) then
         if (nrank==0) print *,'then with AB2',itime
         do ijk=1,nxyz
            phi1(ijk,1,1)=1.5*dt*ta1(ijk,1,1)-0.5*dt*phis1(ijk,1,1)+phi1(ijk,1,1)
            phiss1(ijk,1,1)=phis1(ijk,1,1)
            phis1(ijk,1,1)=ta1(ijk,1,1)
         enddo 
      else
         do ijk=1,nxyz
            phi1(ijk,1,1)=adt(itr)*ta1(ijk,1,1)+bdt(itr)*phis1(ijk,1,1)+&
                 cdt(itr)*phiss1(ijk,1,1)+phi1(ijk,1,1)
            phiss1(ijk,1,1)=phis1(ijk,1,1)
            phis1(ijk,1,1)=ta1(ijk,1,1)
         enddo
      endif
   endif
endif


 end subroutine scalar
 
!~  
!************************************************************
!
subroutine scalar_les_eddy_visc(ux1,uy1,uz1,phi1,phis1,phiss1,di1,ta1,tb1,tc1,td1,&
     uy2,uz2,phi2,di2,ta2,tb2,tc2,td2,uz3,phi3,di3,ta3,tb3,epsi,&
     xnu_sgs1,xnu_sgs2,xnu_sgs3,tau_phi_x1,tau_phi_y2,tau_phi_z3,&
     phimax1,phimin1,phimax2,phimin2,phimax3,phimin3)
!
!************************************************************

USE param
USE variables
USE decomp_2d
use mymath
use parfiX
use parfiY
use parfiZ

implicit none

real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,phi1,phis1,&
                                              phiss1,di1,ta1,tb1,tc1,td1,epsi,xnu_sgs1,&
                                              tau_phi_x1,phimax1,phimin1
real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: uy2,uz2,phi2,di2,ta2,tb2,&
											  tc2,td2,xnu_sgs2,tau_phi_y2,phimax2,phimin2
real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: uz3,phi3,di3,ta3,tb3,&
											  xnu_sgs3,tau_phi_z3,phimax3,phimin3

integer :: ijk,nvect1,nvect2,nvect3,i,j,k,nxyz
real(mytype) :: x,y,z

nvect1=xsize(1)*xsize(2)*xsize(3)
nvect2=ysize(1)*ysize(2)*ysize(3)
nvect3=zsize(1)*zsize(2)*zsize(3)

if (ilimitadvec.eq.1) then
	call set_scalar_minmax(phi1,phi2,phi3,phimax1,phimax2,&
	phimax3,phimin1,phimin2,phimin3)
endif

!X PENCILS
do ijk=1,nvect1
   ta1(ijk,1,1)=ux1(ijk,1,1)*phi1(ijk,1,1)
enddo
call derx (tb1,ta1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0) !d(u*phi)/dx -> tb1
call derx (ta1,phi1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1) !dphi/dx -> ta1
do ijk=1,nvect1
	!(nu_mol+nu_sgs)*(dphi/dx) -> ta1
	ta1(ijk,1,1)=(xnu_sgs1(ijk,1,1)/prtdl)*ta1(ijk,1,1)  
enddo
!d/dx((nu_mol+nu_sgs)*dphi/dx )-> tau_phi_x1 dissipation+SGS term
!phi is symmetric across the boundary, so ta1 is a-symmetric
call derx (tau_phi_x1,ta1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0)
call derxx (ta1,phi1,di1,sx,sfxp,ssxp,swxp,xsize(1),xsize(2),xsize(3),1)
tau_phi_x1 = tau_phi_x1 + xnu/sc*ta1


call transpose_x_to_y(uy1,uy2)
call transpose_x_to_y(uz1,uz2)
if (ilimitadvec.ne.1) then !if ilimitadvec.eq.1, this transposition is done in set_scalar_minmax()
	call transpose_x_to_y(phi1,phi2)
endif


!Y PENCILS


do ijk=1,nvect2
   ta2(ijk,1,1)=uy2(ijk,1,1)*phi2(ijk,1,1)
enddo
call dery (tb2,ta2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0)!d(v*phi)/dy -> tb2
call dery (ta2,phi2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)!dphi/dy -> ta2
do ijk=1,nvect2
	!(nu_mol+nu_sgs)*(dphi/dy) -> ta2
	ta2(ijk,1,1)=(xnu_sgs2(ijk,1,1)/prtdl)*ta2(ijk,1,1) 
enddo
!d/dy((nu_mol+nu_sgs)*dphi/dy ) -> tau_phi_y2 dissipation+SGS term
call dery (tau_phi_y2,ta2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0)  
if (istret.ne.0) then 
   call deryy (ta2,phi2,di2,sy,sfyp,ssyp,swyp,ysize(1),ysize(2),ysize(3),1)
   call dery (tc2,phi2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0)
   do k=1,ysize(3)
   do j=1,ysize(2)
   do i=1,ysize(1)
      ta2(i,j,k)=ta2(i,j,k)*pp2y(j)-pp4y(j)*tc2(i,j,k)
   enddo
   enddo
   enddo
else
   call deryy (ta2,phi2,di2,sy,sfyp,ssyp,swyp,ysize(1),ysize(2),ysize(3),1) 
endif
tau_phi_y2 = tau_phi_y2 + xnu/sc*ta2

!call showval2(tau_phi_y2, 1,2,1)
if (ilimitadvec.ne.1) then !if ilimitadvec.eq.1, this transposition is done in set_scalar_minmax()
	call transpose_y_to_z(phi2,phi3)
endif
call transpose_y_to_z(uz2,uz3)


!Z PENCILS

do ijk=1,nvect3
   ta3(ijk,1,1)=uz3(ijk,1,1)*phi3(ijk,1,1)
enddo
call derz (tb3,ta3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0)!d(w*phi)/dz -> tb3
call derz (ta3,phi3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),0)!dphi/dz -> ta2


do ijk=1,nvect3
	!(nu_mol+nu_sgs)*(dphi/dz) -> ta3
	ta3(ijk,1,1)=(xnu_sgs3(ijk,1,1)/prtdl)*ta3(ijk,1,1) 
enddo
!d/dz((nu_mol+nu_sgs)*dphi/dz ) -> tau_phi_z3 dissipation+SGS term
call derz (tau_phi_z3,ta3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),0)
call derzz (ta3,phi3,di3,sz,sfzp,sszp,swzp,zsize(1),zsize(2),zsize(3),1)
tau_phi_z3 = tau_phi_z3 + xnu/sc*ta3


call transpose_z_to_y(tau_phi_z3,tc2)
call transpose_z_to_y(tb3,td2)


!Y PENCILS ADD TERMS
do ijk=1,nvect2
   tau_phi_y2(ijk,1,1)=tc2(ijk,1,1)+tau_phi_y2(ijk,1,1)
   td2(ijk,1,1)=td2(ijk,1,1)+tb2(ijk,1,1)
enddo

call transpose_y_to_x(tau_phi_y2,tc1)
call transpose_y_to_x(td2,td1)


!X PENCILS ADD TERMS
do ijk=1,nvect1
   tau_phi_x1(ijk,1,1)=tau_phi_x1(ijk,1,1)+tc1(ijk,1,1) !diffusion + SGS term
   tb1(ijk,1,1)=tb1(ijk,1,1)+td1(ijk,1,1) !advection term
enddo

do ijk=1,nvect1
   ta1(ijk,1,1)=tau_phi_x1(ijk,1,1)-tb1(ijk,1,1)
enddo


!TIME ADVANCEMENT
nxyz=xsize(1)*xsize(2)*xsize(3)  

!~ 

if ((nscheme.eq.1).or.(nscheme.eq.2)) then
   if ((nscheme.eq.1.and.itime.eq.1.and.ilit.eq.0).or.&
        (nscheme.eq.2.and.itr.eq.1)) then
      do ijk=1,nxyz
         phi1(ijk,1,1)=gdt(itr)*ta1(ijk,1,1)+phi1(ijk,1,1)
         phis1(ijk,1,1)=ta1(ijk,1,1)          
      enddo
   else
      do ijk=1,nxyz
         phi1(ijk,1,1)=adt(itr)*ta1(ijk,1,1)+bdt(itr)*phis1(ijk,1,1)+phi1(ijk,1,1)
         phis1(ijk,1,1)=ta1(ijk,1,1)          
      enddo
endif
endif

if (nscheme.eq.3) then 
if (nrank==0) print *,'Not ready'
stop 
endif

if (nscheme==4) then
   if ((itime.eq.1).and.(ilit.eq.0)) then
      if (nrank==0) print *,'start with Euler',itime
      do ijk=1,nxyz !start with Euler
         phi1(ijk,1,1)=dt*ta1(ijk,1,1)+phi1(ijk,1,1)
         phis1(ijk,1,1)=ta1(ijk,1,1)          
      enddo
   else
      if  ((itime.eq.2).and.(ilit.eq.0)) then
         if (nrank==0) print *,'then with AB2',itime
         do ijk=1,nxyz
            phi1(ijk,1,1)=1.5*dt*ta1(ijk,1,1)-0.5*dt*phis1(ijk,1,1)+phi1(ijk,1,1)
            phiss1(ijk,1,1)=phis1(ijk,1,1)
            phis1(ijk,1,1)=ta1(ijk,1,1)
         enddo 
      else
         do ijk=1,nxyz
            phi1(ijk,1,1)=adt(itr)*ta1(ijk,1,1)+bdt(itr)*phis1(ijk,1,1)+&
                 cdt(itr)*phiss1(ijk,1,1)+phi1(ijk,1,1)
            phiss1(ijk,1,1)=phis1(ijk,1,1)
            phis1(ijk,1,1)=ta1(ijk,1,1)
         enddo
      endif
   endif
endif



if (ilimitadvec.eq.1) then
	call clip_to_scalar_minmax(phi1,phimax1,phimin1)
endif
end subroutine scalar_les_eddy_visc



!************************************************************




!************************************************************

!
subroutine scalar_les_svm(ux1,uy1,uz1,phi1,phis1,phiss1,di1,ta1,tb1,tc1,td1,&
     uy2,uz2,phi2,di2,ta2,tb2,tc2,td2,uz3,phi3,di3,ta3,tb3,epsi,&
     dphidx1	   ,dphidx2   ,dphidx3,&
     dphidy1	   ,dphidy2   ,dphidy3,&
     dphidz1	   ,dphidz2   ,dphidz3,&
     k_sgs1    ,k_sgs2   ,k_sgs3,&
	 e_svm_x1  ,e_svm_y1  ,e_svm_z1,&
	 e_svm_x2  ,e_svm_y2  ,e_svm_z2,&
	 e_svm_x3  ,e_svm_y3  ,e_svm_z3) 
!
!************************************************************

USE param
USE variables
USE decomp_2d

implicit none

real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,phi1,phis1,&
                                              phiss1,di1,ta1,tb1,tc1,td1,epsi,xnu_sgs1,&
                                              k_sgs1,e_svm_x1,e_svm_y1,e_svm_z1,&
                                              dphidx1,dphidy1,dphidz1
real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: uy2,uz2,phi2,di2,ta2,tb2,&
											  tc2,td2,xnu_sgs2,k_sgs2,&
											  e_svm_x2,e_svm_y2,e_svm_z2,&
											  dphidx2,dphidy2,dphidz2
real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: uz3,phi3,di3,ta3,tb3,&
											  k_sgs3,e_svm_x3,e_svm_y3,e_svm_z3,&
											  dphidx3,dphidy3,dphidz3

integer :: ijk,nvect1,nvect2,nvect3,i,j,k,nxyz
real(mytype) :: x,y,z

nvect1=xsize(1)*xsize(2)*xsize(3)
nvect2=ysize(1)*ysize(2)*ysize(3)
nvect3=zsize(1)*zsize(2)*zsize(3)

!X PENCILS
do ijk=1,nvect1
   ta1(ijk,1,1)=ux1(ijk,1,1)*phi1(ijk,1,1)
enddo
call derx (tb1,ta1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0)
call derxx (ta1,phi1,di1,sx,sfxp,ssxp,swxp,xsize(1),xsize(2),xsize(3),1)


call transpose_x_to_y(phi1,phi2)
call transpose_x_to_y(uy1,uy2)
call transpose_x_to_y(uz1,uz2)

!Y PENCILS
do ijk=1,nvect2
   ta2(ijk,1,1)=uy2(ijk,1,1)*phi2(ijk,1,1)
enddo
call dery (tb2,ta2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0)
if (istret.ne.0) then 
   call deryy (ta2,phi2,di2,sy,sfyp,ssyp,swyp,ysize(1),ysize(2),ysize(3),1)
   call dery (tc2,phi2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0)
   do k=1,ysize(3)
   do j=1,ysize(2)
   do i=1,ysize(1)
      ta2(i,j,k)=ta2(i,j,k)*pp2y(j)-pp4y(j)*tc2(i,j,k)
   enddo
   enddo
   enddo
else
   call deryy (ta2,phi2,di2,sy,sfyp,ssyp,swyp,ysize(1),ysize(2),ysize(3),1) 
endif

call transpose_y_to_z(phi2,phi3)
call transpose_y_to_z(uz2,uz3)


!Z PENCILS
do ijk=1,nvect3
   ta3(ijk,1,1)=uz3(ijk,1,1)*phi3(ijk,1,1)
enddo
call derz (tb3,ta3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0)
call derzz (ta3,phi3,di3,sz,sfzp,sszp,swzp,zsize(1),zsize(2),zsize(3),1)

call transpose_z_to_y(ta3,tc2)
call transpose_z_to_y(tb3,td2)

!Y PENCILS ADD TERMS
do ijk=1,nvect2
   tc2(ijk,1,1)=tc2(ijk,1,1)+ta2(ijk,1,1)
   td2(ijk,1,1)=td2(ijk,1,1)+tb2(ijk,1,1)
enddo

call transpose_y_to_x(tc2,tc1)
call transpose_y_to_x(td2,td1)

!X PENCILS ADD TERMS
do ijk=1,nvect1
   ta1(ijk,1,1)=ta1(ijk,1,1)+tc1(ijk,1,1) !SECOND DERIVATIVE
   tb1(ijk,1,1)=tb1(ijk,1,1)+td1(ijk,1,1) !FIRST DERIVATIVE
enddo
 
do ijk=1,nvect1
   ta1(ijk,1,1)=xnu/sc*ta1(ijk,1,1)-tb1(ijk,1,1) 
enddo

!SVM scalar treatment (pullin, 2000)

call derx (dphidx1,phi1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
call transpose_x_to_y(dphidx1,dphidx2)
call dery (dphidy2,phi2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
call transpose_y_to_z(dphidy2,dphidy3)
call transpose_y_to_z(dphidx2,dphidx3)
call derz (dphidz3,phi3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
! q_z -> ta3
ta3     = sqrt(k_sgs3(:,:,:))*(&
			  (   -e_svm_x3(:,:,:)*e_svm_z3(:,:,:))*dphidx3(:,:,:) &
			+ (   -e_svm_y3(:,:,:)*e_svm_z3(:,:,:))*dphidy3(:,:,:) &
			+ (1.0-e_svm_z3(:,:,:)*e_svm_z3(:,:,:))*dphidz3(:,:,:)  )
! d(q_z)/dz -> tb3
call derz (tb3,ta3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0)

call transpose_z_to_y(tb3,tb2)
call transpose_z_to_y(dphidz3,dphidz2)
! q_y -> ta2
ta2     = sqrt(k_sgs2(:,:,:))*(&
			  (   -e_svm_x2(:,:,:)*e_svm_y2(:,:,:))*dphidx2(:,:,:) &
			+ (1.0-e_svm_y2(:,:,:)*e_svm_y2(:,:,:))*dphidy2(:,:,:) &
			+ (   -e_svm_z2(:,:,:)*e_svm_y2(:,:,:))*dphidz2(:,:,:)  )
! d(q_y)/dy -> tc2
call dery (tc2,ta2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0)
! d(q_z)/dz + d(q_y)/dy
tb2(:,:,:) = tb2(:,:,:) + tc2(:,:,:)

call transpose_y_to_x(tb2,tb1)
call transpose_y_to_x(dphidz2,dphidz1)
call transpose_y_to_x(dphidy2,dphidy1)
! q_x -> td1
td1 = sqrt(k_sgs1(:,:,:))*(&
			  (1.0-e_svm_x1(:,:,:)*e_svm_x1(:,:,:))*dphidx1(:,:,:) &
			+ (   -e_svm_y1(:,:,:)*e_svm_x1(:,:,:))*dphidy1(:,:,:) &
			+ (   -e_svm_z1(:,:,:)*e_svm_x1(:,:,:))*dphidz1(:,:,:)  )
! d(q_x)/dx -> tc1
call derx (tc1,td1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0)
! d(q_z)/dz + d(q_y)/dy + d(q_x)/dx
tb1(:,:,:) = tb1(:,:,:) + tc1(:,:,:)

! add to the total trend
if (istret > 0) then
	do j=1,xsize(2)
		do i=1,xsize(1)
			do k=1,xsize(3)
				ta1(i,j,k) = ta1(i,j,k) + (0.25*delta_bar(zstart(2)-1+j) * tb1(i,j,k))
			enddo
		enddo
	enddo
else
	ta1(:,:,:) = ta1(:,:,:) + (0.25*delta_bar(1) * tb1(:,:,:))
endif


!TIME ADVANCEMENT
nxyz=xsize(1)*xsize(2)*xsize(3)  

if ((nscheme.eq.1).or.(nscheme.eq.2)) then
   if ((nscheme.eq.1.and.itime.eq.1.and.ilit.eq.0).or.&
        (nscheme.eq.2.and.itr.eq.1)) then
      do ijk=1,nxyz
         phi1(ijk,1,1)=gdt(itr)*ta1(ijk,1,1)+phi1(ijk,1,1)
         phis1(ijk,1,1)=ta1(ijk,1,1)          
      enddo
   else
      do ijk=1,nxyz
         phi1(ijk,1,1)=adt(itr)*ta1(ijk,1,1)+bdt(itr)*phis1(ijk,1,1)+phi1(ijk,1,1)
         phis1(ijk,1,1)=ta1(ijk,1,1)          
      enddo
endif
endif

if (nscheme.eq.3) then 
if (nrank==0) print *,'Not ready'
stop 
endif

if (nscheme==4) then
   if ((itime.eq.1).and.(ilit.eq.0)) then
      if (nrank==0) print *,'start with Euler',itime
      do ijk=1,nxyz !start with Euler
         phi1(ijk,1,1)=dt*ta1(ijk,1,1)+phi1(ijk,1,1)
         phis1(ijk,1,1)=ta1(ijk,1,1)          
      enddo
   else
      if  ((itime.eq.2).and.(ilit.eq.0)) then
         if (nrank==0) print *,'then with AB2',itime
         do ijk=1,nxyz
            phi1(ijk,1,1)=1.5*dt*ta1(ijk,1,1)-0.5*dt*phis1(ijk,1,1)+phi1(ijk,1,1)
            phiss1(ijk,1,1)=phis1(ijk,1,1)
            phis1(ijk,1,1)=ta1(ijk,1,1)
         enddo 
      else
         do ijk=1,nxyz
            phi1(ijk,1,1)=adt(itr)*ta1(ijk,1,1)+bdt(itr)*phis1(ijk,1,1)+&
                 cdt(itr)*phiss1(ijk,1,1)+phi1(ijk,1,1)
            phiss1(ijk,1,1)=phis1(ijk,1,1)
            phis1(ijk,1,1)=ta1(ijk,1,1)
         enddo
      endif
   endif
endif
if (itype.gt.9) then
!clip scalar to mix-max values
!~ do ijk=1,nxyz
!~ 	phi1(ijk,1,1) = max(phi1(ijk,1,1), min(phi_bottom,phi_top))
!~ 	phi1(ijk,1,1) = min(phi1(ijk,1,1), max(phi_bottom,phi_top))
!~ enddo
endif
 end subroutine scalar_les_svm

!************************************************************
! 
subroutine set_scalar_minmax(phi1,phi2,phi3,phimax1,phimax2,&
	phimax3,phimin1,phimin2,phimin3)
!
!************************************************************
USE param
USE variables
USE decomp_2d
USE mymath

implicit none

real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: phi1,phimax1,phimin1
real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: phi2,phimax2,phimin2
real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: phi3,phimax3,phimin3

integer :: ijk,nvect1,nvect2,nvect3,i,j,k,nxyz
real(mytype) :: x,y,z

!~ do i=1,xsize(1)
!~ 	do j=1,xsize(2)
!~ 		do k=1,xsize(3)
!~ 			phimax1(i,j,k) = phi1(i,j,k)
!~ 			phimin1(i,j,k) = phi1(i,j,k)
!~ 		enddo
!~ 	enddo
!~ enddo
 
!~ call showval1(phimin1,124,92,49)

do i=2,xsize(1)-1
	do j=1,xsize(2)
		do k=1,xsize(3)	
			phimax1(i,j,k) = max(phi1(i,j,k),max(phi1(i+1,j,k),phi1(i-1,j,k)))
			phimin1(i,j,k) = min(phi1(i,j,k),min(phi1(i+1,j,k),phi1(i-1,j,k)))
		enddo
	enddo
enddo

!boundary points
do j=1,xsize(2)
	do k=1,xsize(3)
		if (nclx.eq.0) then
			phimax1(1,j,k) = max(phimax1(1,j,k),max(phi1(2,j,k),phi1(xsize(1),j,k)))
			phimin1(1,j,k) = min(phimin1(1,j,k),min(phi1(2,j,k),phi1(xsize(1),j,k)))
			phimax1(xsize(1),j,k) = max(phimax1(xsize(1),j,k),max(phi1(1,j,k),phi1(xsize(1)-1,j,k)))
			phimin1(xsize(1),j,k) = min(phimax1(xsize(1),j,k),min(phi1(1,j,k),phi1(xsize(1)-1,j,k)))
		endif
	enddo
enddo
call transpose_x_to_y(phi1,phi2)
call transpose_x_to_y(phimax1,phimax2)
call transpose_x_to_y(phimin1,phimin2)

do i=1,ysize(1)
	do j=2,ysize(2)-1
		do k=1,ysize(3)
			phimax2(i,j,k) = max(phimax2(i,j,k),max(phi2(i,j-1,k),phi2(i,j+1,k)))
			phimin2(i,j,k) = min(phimin2(i,j,k),min(phi2(i,j-1,k),phi2(i,j+1,k)))
		enddo
	enddo
enddo

call transpose_y_to_z(phi2,phi3)
call transpose_y_to_z(phimax2,phimax3)
call transpose_y_to_z(phimin2,phimin3)

do i=1,zsize(1)
	do j=1,zsize(2)
		do k=2,zsize(3)-1
			phimax3(i,j,k) = max(phimax3(i,j,k),max(phi3(i,j,k-1),phi3(i,j,k+1)))
			phimin3(i,j,k) = min(phimin3(i,j,k),min(phi3(i,j,k-1),phi3(i,j,k+1)))
		enddo
	enddo
enddo

!boundary points
do j=1,zsize(2)
	do i=1,zsize(1)
		if (nclz.eq.0) then
			phimax3(i,j,1) = max(phimax3(i,j,1),max(phi3(i,j,2),phi3(i,j,zsize(3))))
			phimin3(i,j,1) = min(phimin3(i,j,1),min(phi3(i,j,2),phi3(i,j,zsize(3))))
			phimax3(i,j,zsize(3)) = max(phimax3(i,j,zsize(3)),max(phi3(i,j,1),phi3(i,j,zsize(3)-1)))
			phimin3(i,j,zsize(3)) = min(phimin3(i,j,zsize(3)),min(phi3(i,j,1),phi3(i,j,zsize(3)-1)))
		endif
	enddo
enddo

call transpose_z_to_y(phimax3,phimax2)
call transpose_z_to_y(phimin3,phimin2)
call transpose_y_to_x(phimax2,phimax1)
call transpose_y_to_x(phimin2,phimin1)

end subroutine set_scalar_minmax


!************************************************************
! 
subroutine clip_to_scalar_minmax(phi1,phimax1,phimin1)
!***************************
USE param
USE variables
USE decomp_2d

implicit none

real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: phi1,phimax1,phimin1
integer :: ijk,nvect1,nvect2,nvect3,i,j,k,nxyz,minclipcount,maxclipcount
real(mytype) :: x,y,z,tempphi,vpermin,vpermax

minclipcount = 0
maxclipcount = 0
vpermin = 0.0
vpermax = 0.0
do i=1,xsize(1)
do j=1,xsize(2)
do k=1,xsize(3)
tempphi = phi1(i,j,k)
phi1(i,j,k) = max(phi1(i,j,k), phimin1(i,j,k))
phi1(i,j,k) = min(phi1(i,j,k), phimax1(i,j,k))
if (tempphi>phi1(i,j,k)) then
maxclipcount = maxclipcount + 1
vpermax = vpermax + (tempphi-phi1(i,j,k))/phi1(i,j,k)
endif
if (tempphi<phi1(i,j,k)) then
minclipcount = minclipcount + 1
vpermin = vpermin + (phi1(i,j,k)-tempphi)/phi1(i,j,k)
endif
enddo
enddo
enddo
!~ print *,'max clip- thread',nrank,100*float(maxclipcount)/float(xsize(1)*xsize(2)*xsize(3)),'clipped',vpermax/float(maxclipcount)
!~ print *,'min clip- thread',nrank,100*float(minclipcount)/float(xsize(1)*xsize(2)*xsize(3)),'clipped',vpermin/float(minclipcount)

end subroutine clip_to_scalar_minmax



