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
     div_tau_x3,div_tau_y3,div_tau_z3)
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
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: duxdx1,duydx1,duzdx1,duxdy1,duxdz1,les_a1,les_b1
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: div_tau_x1,div_tau_y1,div_tau_z1,xnu_sgs1
real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: duxdx2,duydx2,duzdx2,les_a2,les_b2
real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: duxdy2,duydy2,duzdy2,duxdz2,duydz2
real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: div_tau_x2,div_tau_y2,div_tau_z2,xnu_sgs2
real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: duxdx3,duydx3,duzdx3,les_a3
real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: duxdy3,duydy3,duzdy3
real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: duxdz3,duydz3,duzdz3
real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: div_tau_x3,div_tau_y3,div_tau_z3,xnu_sgs3


integer :: ijk,nvect1,nvect2,nvect3,i,j,k
real(mytype) :: x,y,z


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

if (iles > 0) then !calculate nu_SGS=(delta*C_S)^2 * sqrt(2*Sij*Sij)
    do ijk=1,nvect3
       xnu_sgs3(ijk,1,1) = 2.0 * duxdx3(ijk,1,1) * duxdx3(ijk,1,1)
       xnu_sgs3(ijk,1,1) = 2.0 * duydy3(ijk,1,1) * duydy3(ijk,1,1) + xnu_sgs3(ijk,1,1)
       xnu_sgs3(ijk,1,1) = 2.0 * duzdz3(ijk,1,1) * duzdz3(ijk,1,1) + xnu_sgs3(ijk,1,1)
       xnu_sgs3(ijk,1,1) = ( duxdy3(ijk,1,1) + duydx3(ijk,1,1) )  ** 2.0 + xnu_sgs3(ijk,1,1)
       xnu_sgs3(ijk,1,1) = ( duxdz3(ijk,1,1) + duzdx3(ijk,1,1) )  ** 2.0 + xnu_sgs3(ijk,1,1)
       xnu_sgs3(ijk,1,1) = ( duydz3(ijk,1,1) + duzdy3(ijk,1,1) )  ** 2.0 + xnu_sgs3(ijk,1,1)
    enddo
    if (istret > 0) then
		do j=1,zsize(2)
			do i=1,zsize(1)
				do k=1,zsize(3)
					xnu_sgs3(i,j,k) = (delta_bar(zstart(2)-1+j) * Csmag)**2.0 * sqrt(xnu_sgs3(i,j,k))
				enddo
			enddo
		enddo
    else
		xnu_sgs3(:,:,:) = (delta_bar(1) * Csmag)**2.0 * sqrt(xnu_sgs3(:,:,:))
    endif
endif

!DIFFUSIVE TERMS IN Z
call derzz (ta3,ux3,di3,sz,sfzp,sszp,swzp,zsize(1),zsize(2),zsize(3),1)
call derzz (tb3,uy3,di3,sz,sfzp,sszp,swzp,zsize(1),zsize(2),zsize(3),1)
call derzz (tc3,uz3,di3,sz,sfz ,ssz ,swz ,zsize(1),zsize(2),zsize(3),0)

if (iles == 1) then
    les_a3(:,:,:) = -2.0*xnu_sgs3(:,:,:) * 0.5*(duxdz3(:,:,:) + duzdx3(:,:,:) ) ! tau_xz = -2.0*nu*Sxz 
    call derz(div_tau_x3, les_a3 ,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0) ! div_tau_x=d(tau_xz)/dz
    les_a3(:,:,:) = -2.0*xnu_sgs3(:,:,:) * 0.5*(duydz3(:,:,:) + duzdy3(:,:,:) ) ! tau_yz = -2.0*nu*Syz 
    call derz(div_tau_y3, les_a3 ,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0) ! div_tau_y=d(tau_yz)/dz
    les_a3(:,:,:) = -2.0*xnu_sgs3(:,:,:) * duzdz3(:,:,:)                                   ! tau_zz = -2.0*nu*Szz 
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
    call transpose_z_to_y(xnu_sgs3,xnu_sgs2)
    call transpose_z_to_y(duydz3,duydz2)
    call transpose_z_to_y(duxdz3,duxdz2)
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

if (iles == 1) then
    les_a2(:,:,:) = -2.0*xnu_sgs2(:,:,:) *0.5* (duxdy2(:,:,:) + duydx2(:,:,:) ) ! tau_xy = -2.0*nu*Sxy 
    call dery(les_b2, les_a3 ,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0) ! d(tau_xy)/dy
    div_tau_x2(:,:,:) = div_tau_x2(:,:,:) + les_b2(:,:,:) !div_tau_x = div_tau_x + d(tau_xy)/dy
    les_a2(:,:,:) = -2.0*xnu_sgs2(:,:,:) * duydy2(:,:,:)                            ! tau_yy = -2.0*nu*Syy 
    call dery(les_b2, les_a2 ,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1) ! d(tau_yy)/dy
    div_tau_y2(:,:,:) = div_tau_y2(:,:,:) + les_b2(:,:,:) !div_tau_y = div_tau_y + d(tau_yy)/dy
    les_a2(:,:,:) = -2.0*xnu_sgs2(:,:,:)*0.5*(duzdy2(:,:,:)+duydz2(:,:,:))! tau_zy = -2.0*nu*Szy 
    call dery(les_b2, les_a3 ,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0) ! d(tau_zy)/dy
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
    call transpose_y_to_x(xnu_sgs2,xnu_sgs1)
    call transpose_y_to_x(duxdy2,duxdy1)
    call transpose_y_to_x(duxdz2,duxdz1)
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

if (iles == 1) then
    les_a1(:,:,:) = -2.0*xnu_sgs1(:,:,:) * duxdx1(:,:,:)  ! tau_xx = -2.0*nu*Sxx 
    call derx(les_b1, les_a1 ,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1) ! d(tau_xx)/dx
    div_tau_x1(:,:,:) = div_tau_x1(:,:,:) + les_b1(:,:,:) !div_tau_x = div_tau_x + d(tau_xx)/dx
    les_a1(:,:,:) = -2.0*xnu_sgs1(:,:,:) * 0.5 *(duxdy1(:,:,:)+duydx1(:,:,:))                            ! tau_xy = -2.0*nu*Sxy 
    call derx(les_b1, les_a1 ,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0) ! d(tau_xy)/dx
    div_tau_y1(:,:,:) = div_tau_y1(:,:,:) + les_b1(:,:,:) !div_tau_y = div_tau_y + d(tau_xy)/dx
    les_a1(:,:,:) = -2.0*xnu_sgs1(:,:,:) * 0.5* (duzdx1(:,:,:) + duxdz1(:,:,:) )! tau_xz = -2.0*nu*Sxz 
    call derx(les_b1, les_a1 ,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0) ! d(tau_xz)/dx
    div_tau_z1(:,:,:) = div_tau_z1(:,:,:) + les_b1(:,:,:) !div_tau_z = div_tau_z + d(tau_xz)/dx
endif

!FINAL SUM: DIFF TERMS + CONV TERMS
ta1(:,:,:)=xnu*ta1(:,:,:)-tg1(:,:,:)
tb1(:,:,:)=xnu*tb1(:,:,:)-th1(:,:,:)
tc1(:,:,:)=xnu*tc1(:,:,:)-ti1(:,:,:)
if (iles == 1) then
	ta1(:,:,:)=div_tau_x1(:,:,:) + ta1(:,:,:)
    tb1(:,:,:)=div_tau_y1(:,:,:) + tb1(:,:,:)
    tc1(:,:,:)=div_tau_z1(:,:,:) + tc1(:,:,:)
endif

end subroutine convdiff



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
 
 
 !************************************************************
!
subroutine scalar_les(ux1,uy1,uz1,phi1,phis1,phiss1,di1,ta1,tb1,tc1,td1,&
     uy2,uz2,phi2,di2,ta2,tb2,tc2,td2,uz3,phi3,di3,ta3,tb3,epsi,&
     xnu_sgs1,xnu_sgs2,xnu_sgs3,tau_phi_x1,tau_phi_y2,tau_phi_z3)
!
!************************************************************

USE param
USE variables
USE decomp_2d

implicit none

real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,phi1,phis1,&
                                              phiss1,di1,ta1,tb1,tc1,td1,epsi,xnu_sgs1,&
                                              tau_phi_x1
real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: uy2,uz2,phi2,di2,ta2,tb2,&
											  tc2,td2,xnu_sgs2,tau_phi_y2
real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: uz3,phi3,di3,ta3,tb3,&
											  xnu_sgs3,tau_phi_z3

integer :: ijk,nvect1,nvect2,nvect3,i,j,k,nxyz
real(mytype) :: x,y,z

nvect1=xsize(1)*xsize(2)*xsize(3)
nvect2=ysize(1)*ysize(2)*ysize(3)
nvect3=zsize(1)*zsize(2)*zsize(3)

!X PENCILS
do ijk=1,nvect1
   ta1(ijk,1,1)=ux1(ijk,1,1)*phi1(ijk,1,1)
enddo
call derx (tb1,ta1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0) !d(u*phi)/dx -> tb1
call derx (ta1,phi1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1) !dphi/dx -> ta1
do ijk=1,nvect1
	!(nu_mol+nu_sgs)*(dphi/dx) -> ta1
	ta1(ijk,1,1)=(xnu_sgs1(ijk,1,1)+xnu/sc)*ta1(ijk,1,1)  
enddo
!d/dx((nu_mol+nu_sgs)*dphi/dx )-> tau_phi_x1 dissipation+SGS term
!phi is symmetric across the boundary, so ta1 is a-symmetric
call derx (tau_phi_x1,ta1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0)  

call transpose_x_to_y(phi1,phi2)
call transpose_x_to_y(uy1,uy2)
call transpose_x_to_y(uz1,uz2)

!Y PENCILS
do ijk=1,nvect2
   ta2(ijk,1,1)=uy2(ijk,1,1)*phi2(ijk,1,1)
enddo
call dery (tb2,ta2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0)!d(v*phi)/dy -> tb2
call dery (ta2,phi2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)!dphi/dy -> ta2
do ijk=1,nvect2
	!(nu_mol+nu_sgs)*(dphi/dy) -> ta2
	ta2(ijk,1,1)=(xnu_sgs2(ijk,1,1)+xnu/sc)*ta2(ijk,1,1) 
enddo
!d/dy((nu_mol+nu_sgs)*dphi/dy ) -> tau_phi_y2 dissipation+SGS term
call dery (tau_phi_y2,ta2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0)  

call transpose_y_to_z(phi2,phi3)
call transpose_y_to_z(uz2,uz3)

!Z PENCILS
do ijk=1,nvect3
   ta3(ijk,1,1)=uz3(ijk,1,1)*phi3(ijk,1,1)
enddo
call derz (tb3,ta3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0)!d(w*phi)/dz -> tb3
call derz (ta3,phi3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),0)!dphi/dz -> ta2
do ijk=1,nvect3
	!(nu_mol+nu_sgs)*(dphi/dz) -> ta3
	ta3(ijk,1,1)=(xnu_sgs3(ijk,1,1)+xnu/sc)*ta3(ijk,1,1) 
enddo
!d/dz((nu_mol+nu_sgs)*dphi/dz ) -> tau_phi_z3 dissipation+SGS term
call derz (tau_phi_z3,ta3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),0)


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
   ta1(ijk,1,1)=ta1(ijk,1,1)-tb1(ijk,1,1) 
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


 end subroutine scalar_les



