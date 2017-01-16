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

module var

use decomp_2d
USE variables
USE param

! define all major arrays here

real(mytype), save, allocatable, dimension(:,:,:) :: ux1, ux2, ux3,po3,dv3,pp3
real(mytype), save, allocatable, dimension(:,:,:) :: uy1, uy2, uy3
real(mytype), save, allocatable, dimension(:,:,:) :: uz1, uz2, uz3
real(mytype), save, allocatable, dimension(:,:,:) :: phi1, phi2, phi3
real(mytype), save, allocatable, dimension(:,:,:) :: gx1, gy1, gz1, hx1, hy1, hz1, phis1,phiss1
real(mytype), save, allocatable, dimension(:,:,:) :: px1, py1, pz1
real(mytype), save, allocatable, dimension(:,:,:) :: ep1

!arrays for statistic collection
real(mytype), save, allocatable, dimension(:,:,:) :: umean,vmean,wmean,uumean,vvmean,wwmean,uvmean,uwmean,vwmean,tmean
real(mytype), save, allocatable, dimension(:,:,:) :: phimean, phiphimean

!arrays for visualization
real(mytype), save, allocatable, dimension(:,:,:) :: uvisu

! define all work arrays here
real(mytype), save, allocatable, dimension(:,:,:) :: ta1,tb1,tc1,td1,&
     te1,tf1,tg1,th1,ti1,di1
real(mytype), save, allocatable, dimension(:,:,:) :: ta2,tb2,tc2,td2,&
     te2,tf2,tg2,th2,ti2,tj2,di2
real(mytype), save, allocatable, dimension(:,:,:) :: ta3,tb3,tc3,td3,&
     te3,tf3,tg3,th3,ti3,di3

! arrays for explicit LES 
real(mytype), save, allocatable, dimension(:,:,:) :: duxdx1, duxdx2, duxdx3 !dudx
real(mytype), save, allocatable, dimension(:,:,:) :: duydx1, duydx2, duydx3 !dvdx
real(mytype), save, allocatable, dimension(:,:,:) :: duzdx1, duzdx2, duzdx3 !dwdx
real(mytype), save, allocatable, dimension(:,:,:) :: duxdy1, duxdy2, duxdy3 !dudy
real(mytype), save, allocatable, dimension(:,:,:) :: duydy2, duydy3 !dvdy
real(mytype), save, allocatable, dimension(:,:,:) :: duzdy2, duzdy3 !dwdy
real(mytype), save, allocatable, dimension(:,:,:) :: duxdz1,duxdz2,duxdz3 !dudz
real(mytype), save, allocatable, dimension(:,:,:) :: duydz2,duydz3 !dvdz
real(mytype), save, allocatable, dimension(:,:,:) :: duzdz3 !dwdz
real(mytype), save, allocatable, dimension(:,:,:) :: xnu_sgs1, xnu_sgs2, xnu_sgs3 !nu_sgs
real(mytype), save, allocatable, dimension(:,:,:) :: div_tau_x1, div_tau_x2, div_tau_x3 !d(tau_xj)/dxj
real(mytype), save, allocatable, dimension(:,:,:) :: div_tau_y1, div_tau_y2, div_tau_y3 !d(tau_yj)/dxj
real(mytype), save, allocatable, dimension(:,:,:) :: div_tau_z1, div_tau_z2, div_tau_z3 !d(tau_zj)/dxj
real(mytype), save, allocatable, dimension(:,:,:) :: les_a1,les_a2,les_a3,les_b1,les_b2 !les working variables
real(mytype), save, allocatable, dimension(:,:,:) :: tau_phi_x1,tau_phi_y2,tau_phi_z3   ! nu*d(phi)/dx_i
!variables for stretched vortex SGS model
real(mytype), save, allocatable, dimension(:,:,:) :: k_sgs1,e_svm_x1,e_svm_y1,e_svm_z1  
real(mytype), save, allocatable, dimension(:,:,:) :: k_sgs2,e_svm_x2,e_svm_y2,e_svm_z2
real(mytype), save, allocatable, dimension(:,:,:) :: k_sgs3,e_svm_x3,e_svm_y3,e_svm_z3
real(mytype), save, allocatable, dimension(:,:,:) :: dphidx1,dphidy1,dphidz1
real(mytype), save, allocatable, dimension(:,:,:) :: dphidx2,dphidy2,dphidz2
real(mytype), save, allocatable, dimension(:,:,:) :: dphidx3,dphidy3,dphidz3

! 
integer, save :: nxmsize, nymsize, nzmsize 


contains

  
  subroutine init_variables

    TYPE(DECOMP_INFO), save :: ph  ! decomposition object

    if (nclx==0) then
       nxmsize = xsize(1)
    else
       nxmsize = xsize(1) -1
    endif
    if (ncly==0) then
       nymsize = ysize(2)
    else
       nymsize = ysize(2) -1
    endif
    if (nclz==0) then
       nzmsize = zsize(3)
    else
       nzmsize = zsize(3) -1
    endif
    call decomp_info_init(nxmsize, nymsize, nzmsize, ph)
    

!X PENCILS
    call alloc_x(ux1, opt_global=.true.)
    call alloc_x(uy1, opt_global=.true.)
#ifndef TWOD
    call alloc_x(uz1, opt_global=.true.)
    call alloc_x(pz1, opt_global=.true.)
#else
    allocate (uz1(1,1,1))
    allocate (pz1(1,1,1))
#endif
    call alloc_x(px1, opt_global=.true.)
    call alloc_x(py1, opt_global=.true.)
    call alloc_x(phi1, opt_global=.true.)
    call alloc_x(gx1);call alloc_x(gy1);call alloc_x(gz1);call alloc_x(phis1) 
    call alloc_x(hx1);call alloc_x(hy1);call alloc_x(hz1);call alloc_x(phiss1)
    call alloc_x(ta1);call alloc_x(tb1);call alloc_x(tc1)
    call alloc_x(td1);call alloc_x(te1);call alloc_x(tf1)
    call alloc_x(tg1);call alloc_x(th1);call alloc_x(ti1)
    call alloc_x(di1);call alloc_x(ep1)
    if (iles > 0) then
		! allocation for explicit LES 
		call alloc_x(duxdx1);  call alloc_x(duydx1); call alloc_x(duzdx1); 
		call alloc_x(duxdy1);  call alloc_x(duxdz1); 
		call alloc_x(div_tau_x1);  call alloc_x(div_tau_y1); call alloc_x(div_tau_z1);
		call alloc_x(xnu_sgs1);  call alloc_x(les_a1); call alloc_x(les_b1); 
		if (iscalar.eq.1) then 
			call alloc_x(tau_phi_x1); 
		endif
		if ((iles.eq.4).or.(iles.eq.5)) then
			call alloc_x(k_sgs1);  call alloc_x(e_svm_x1); 
			call alloc_x(e_svm_y1);  call alloc_x(e_svm_z1);
			call alloc_x(dphidx1);  call alloc_x(dphidy1);
			call alloc_x(dphidz1);
		endif
	endif
    allocate(sx(xsize(2),xsize(3)),vx(xsize(2),xsize(3)))
    !inflow/ouflow 2d arrays
    allocate(bxx1(xsize(2),xsize(3)),bxy1(xsize(2),xsize(3)))
    allocate(bxz1(xsize(2),xsize(3)),bxxn(xsize(2),xsize(3)))
    allocate(bxyn(xsize(2),xsize(3)),bxzn(xsize(2),xsize(3)))
    allocate(bxo(xsize(2),xsize(3)),byo(xsize(2),xsize(3)))
    allocate(bzo(xsize(2),xsize(3)))
    allocate(byx1(xsize(1),xsize(3)),byy1(xsize(1),xsize(3)))
    allocate(byz1(xsize(1),xsize(3)),byxn(xsize(1),xsize(3)))
    allocate(byyn(xsize(1),xsize(3)),byzn(xsize(1),xsize(3)))   
    allocate(bzx1(xsize(1),xsize(2)),bzy1(xsize(1),xsize(2)))
    allocate(bzz1(xsize(1),xsize(2)),bzxn(xsize(1),xsize(2)))
    allocate(bzyn(xsize(1),xsize(2)),bzzn(xsize(1),xsize(2)))
    !pre_correc 2d array
    allocate(dpdyx1(xsize(2),xsize(3)),dpdyxn(xsize(2),xsize(3)))
    allocate(dpdzx1(xsize(2),xsize(3)),dpdzxn(xsize(2),xsize(3)))
    allocate(dpdxy1(xsize(1),xsize(3)),dpdxyn(xsize(1),xsize(3)))
    allocate(dpdzy1(xsize(1),xsize(3)),dpdzyn(xsize(1),xsize(3)))
    allocate(dpdxz1(xsize(1),xsize(2)),dpdxzn(xsize(1),xsize(2)))
    allocate(dpdyz1(xsize(1),xsize(2)),dpdyzn(xsize(1),xsize(2)))

!arrays for statistic collection!pay attention to the size!
    allocate (umean(xstS(1):xenS(1),xstS(2):xenS(2),xstS(3):xenS(3)))
    allocate (vmean(xstS(1):xenS(1),xstS(2):xenS(2),xstS(3):xenS(3)))
    allocate (wmean(xstS(1):xenS(1),xstS(2):xenS(2),xstS(3):xenS(3)))
    allocate (uumean(xstS(1):xenS(1),xstS(2):xenS(2),xstS(3):xenS(3)))
    allocate (vvmean(xstS(1):xenS(1),xstS(2):xenS(2),xstS(3):xenS(3)))
    allocate (wwmean(xstS(1):xenS(1),xstS(2):xenS(2),xstS(3):xenS(3)))
    allocate (uvmean(xstS(1):xenS(1),xstS(2):xenS(2),xstS(3):xenS(3)))
    allocate (uwmean(xstS(1):xenS(1),xstS(2):xenS(2),xstS(3):xenS(3)))
    allocate (vwmean(xstS(1):xenS(1),xstS(2):xenS(2),xstS(3):xenS(3)))
    allocate (tmean(xstS(1):xenS(1),xstS(2):xenS(2),xstS(3):xenS(3)))    
    if (iscalar==1) then
       allocate (phimean(xstS(1):xenS(1),xstS(2):xenS(2),xstS(3):xenS(3)))
       allocate (phiphimean(xstS(1):xenS(1),xstS(2):xenS(2),xstS(3):xenS(3)))
    else
       allocate (phimean(1,1,1))
       allocate (phiphimean(1,1,1))
    endif

!arrays for visualization!pay attention to the size!
    allocate (uvisu(xstV(1):xenV(1),xstV(2):xenV(2),xstV(3):xenV(3)))

!Y PENCILS
    call alloc_y(ux2);call alloc_y(uy2);call alloc_y(uz2)
    call alloc_y(ta2);call alloc_y(tb2);call alloc_y(tc2)
    call alloc_y(td2);call alloc_y(te2);call alloc_y(tf2)
    call alloc_y(tg2);call alloc_y(th2);call alloc_y(ti2)
    call alloc_y(tj2)
    call alloc_y(di2);call alloc_y(phi2)
    if (iles > 0) then
		! allocation for explicit LES 
		call alloc_y(duxdx2);  call alloc_y(duydx2); call alloc_y(duzdx2);
		call alloc_y(duxdy2);  call alloc_y(duydy2); call alloc_y(duzdy2);
		call alloc_y(duydz2);  call alloc_y(duxdz2); 
		call alloc_y(div_tau_x2);  call alloc_y(div_tau_y2); call alloc_y(div_tau_z2);
		call alloc_y(xnu_sgs2);  call alloc_y(les_a2); call alloc_y(les_b2); 
		if (iscalar.eq.1) then 
			call alloc_y(tau_phi_y2); 
		endif
		if ((iles.eq.4).or.(iles.eq.5)) then
			call alloc_y(k_sgs2);  call alloc_y(e_svm_x2); 
			call alloc_y(e_svm_y2);  call alloc_y(e_svm_z2);
			call alloc_y(dphidx2);  call alloc_y(dphidy2);
			call alloc_y(dphidz2);
		endif		
	endif
    allocate(sy(ysize(1),ysize(3)),vy(ysize(1),ysize(3)))
!Z PENCILS
    call alloc_z(ux3);call alloc_z(uy3);call alloc_z(uz3)
    call alloc_z(ta3);call alloc_z(tb3);call alloc_z(tc3)
    call alloc_z(td3);call alloc_z(te3);call alloc_z(tf3)
    call alloc_z(tg3);call alloc_z(th3);call alloc_z(ti3)
    call alloc_z(di3);call alloc_z(phi3)
    if (iles > 0) then
		! allocation for explicit LES 
		call alloc_z(duxdx3);  call alloc_z(duydx3); call alloc_z(duzdx3);
		call alloc_z(duxdy3);  call alloc_z(duydy3); call alloc_z(duzdy3);
		call alloc_z(duxdz3);  call alloc_z(duydz3); call alloc_z(duzdz3);
		call alloc_z(div_tau_x3);  call alloc_z(div_tau_y3); call alloc_z(div_tau_z3);
		call alloc_z(xnu_sgs3);  call alloc_z(les_a3); 
		if (iscalar.eq.1) then 
			call alloc_z(tau_phi_z3); 
		endif
		if ((iles.eq.4).or.(iles.eq.5)) then
			call alloc_z(k_sgs3);  call alloc_z(e_svm_x3); 
			call alloc_z(e_svm_y3);  call alloc_z(e_svm_z3);
			call alloc_z(dphidx3);  call alloc_z(dphidy3);
			call alloc_z(dphidz3);
		endif		
	endif
    allocate(sz(zsize(1),zsize(2)),vz(zsize(1),zsize(2)))

 ! if all periodic
 !   allocate (pp3(ph%zst(1):ph%zen(1),ph%zst(2):ph%zen(2),ph%zst(3):ph%zen(3)))
 !   allocate (dv3(ph%zst(1):ph%zen(1),ph%zst(2):ph%zen(2),ph%zst(3):ph%zen(3)))
 !   allocate (po3(ph%zst(1):ph%zen(1),ph%zst(2):ph%zen(2),ph%zst(3):ph%zen(3)))
    call alloc_z(pp3,ph,.true.)
    call alloc_z(dv3,ph,.true.)
    call alloc_z(po3,ph,.true.)

    return
  end subroutine init_variables

end module var

