

!********************************************************************
!
subroutine xnu_smag(duxdx3,duydx3,duzdx3,&
     duxdy3,duydy3,duzdy3,&
     duxdz3,duydz3,duzdz3,&
     xnu_sgs3)
! 
!********************************************************************
USE param
USE variables
USE decomp_2d

!Z pencils 

implicit none

real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: duxdx3,duydx3,duzdx3
real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: duxdy3,duydy3,duzdy3
real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: duxdz3,duydz3,duzdz3
real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: xnu_sgs3


integer :: ijk,nvect1,nvect2,nvect3,i,j,k
real(mytype) :: x,y,z

nvect3=zsize(1)*zsize(2)*zsize(3)
!calculate nu_SGS=(delta*C_S)^2 * sqrt(2*Sij*Sij)
    do ijk=1,nvect3
       xnu_sgs3(ijk,1,1) = 2.0 * duxdx3(ijk,1,1) * duxdx3(ijk,1,1)
       xnu_sgs3(ijk,1,1) = 2.0 * duydy3(ijk,1,1) * duydy3(ijk,1,1) + xnu_sgs3(ijk,1,1)
       xnu_sgs3(ijk,1,1) = 2.0 * duzdz3(ijk,1,1) * duzdz3(ijk,1,1) + xnu_sgs3(ijk,1,1)
       xnu_sgs3(ijk,1,1) = 4.0 * ( 0.5*(duxdy3(ijk,1,1) + duydx3(ijk,1,1) ) ) ** 2.0 + xnu_sgs3(ijk,1,1)
       xnu_sgs3(ijk,1,1) = 4.0 * ( 0.5*(duxdz3(ijk,1,1) + duzdx3(ijk,1,1) ) ) ** 2.0 + xnu_sgs3(ijk,1,1)
       xnu_sgs3(ijk,1,1) = 4.0 * ( 0.5*(duydz3(ijk,1,1) + duzdy3(ijk,1,1) ) ) ** 2.0 + xnu_sgs3(ijk,1,1)
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
		do ijk=1,nvect3
			xnu_sgs3(ijk,1,1) = (delta_bar(1) * Csmag)**2.0 * sqrt(xnu_sgs3(ijk,1,1))
		enddo
    endif
!~     if (nrank.eq.1.and.itime.gt.5) then
!~     print *,duxdx3(5,5,5),duxdy3(5,5,5),duxdz3(5,5,5)
!~     print *,duydx3(5,5,5),duydy3(5,5,5),duydz3(5,5,5)
!~     print *,duzdx3(5,5,5),duzdy3(5,5,5),duzdz3(5,5,5)
!~     print *,delta_bar(1),Csmag
!~     print *,xnu_sgs3(5,5,5)
!~     print *,dx,dy,dz
!~     stop
!~     endif
	return
end subroutine xnu_smag



!********************************************************************
!
subroutine xnu_sigma(duxdx3,duydx3,duzdx3,&
     duxdy3,duydy3,duzdy3,&
     duxdz3,duydz3,duzdz3,&
     xnu_sgs3)
! 
!********************************************************************
USE param
USE variables
USE decomp_2d
USE mymath

!Z pencils 

implicit none

real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: duxdx3,duydx3,duzdx3
real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: duxdy3,duydy3,duzdy3
real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: duxdz3,duydz3,duzdz3
real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: xnu_sgs3
real(mytype) :: G11,G12,G13,G22,G23,G33 !elements of the G symmetric matrix

integer :: ijk,nvect1,nvect2,nvect3,i,j,k
real(mytype) :: x,x1,x2,x3

do j=1,zsize(2)
	do i=1,zsize(1)
		do k=1,zsize(3) 
			G11 = duxdx3(i,j,k) * duxdx3(i,j,k)  + duydx3(i,j,k) * duydx3(i,j,k)&
					+ duzdx3(i,j,k) * duzdx3(i,j,k)
			G22 = duxdy3(i,j,k) * duxdy3(i,j,k)  + duydy3(i,j,k) * duydy3(i,j,k)&
					+ duzdy3(i,j,k) * duzdy3(i,j,k)
			G33 = duxdz3(i,j,k) * duxdz3(i,j,k)  + duydz3(i,j,k) * duydz3(i,j,k)&
					+ duzdz3(i,j,k) * duzdz3(i,j,k)
			G12 = duxdx3(i,j,k) * duxdy3(i,j,k)  + duydx3(i,j,k) * duydy3(i,j,k)&
					+ duzdx3(i,j,k) * duzdy3(i,j,k)
			G13 = duxdx3(i,j,k) * duxdz3(i,j,k)  + duydx3(i,j,k) * duydz3(i,j,k)&
					+ duzdx3(i,j,k) * duzdz3(i,j,k)	
			G23 = duxdy3(i,j,k) * duxdz3(i,j,k)  + duydy3(i,j,k) * duydz3(i,j,k)&
					+ duzdy3(i,j,k) * duzdz3(i,j,k)
			call eig_sym_smith(G11,G12,G13,G22,G23,G33,x1,x2,x3)
			x1 = sqrt(abs(x1))
			x2 = sqrt(abs(x2))
			x3 = sqrt(abs(x3))
			call sort3(x1,x2,x3)
			xnu_sgs3(i,j,k) = x3*(x1-x2)*(x2-x3)/x1/x1
!~ 			if (j.eq.33) then
!~ 			print *,xnu_sgs3(i,j,k),x1,x2,x3,G11,G12,G13,G22,G23,G33
!~ 			print *,duxdx3(i,j,k) , duxdy3(i,j,k) , duxdz3(i,j,k)
!~ 			print *,duydx3(i,j,k) , duydy3(i,j,k) , duydz3(i,j,k)
!~ 			print *,duzdx3(i,j,k) , duzdy3(i,j,k) , duzdz3(i,j,k)
!~ 			endif
			if (istret>0) then
				xnu_sgs3(i,j,k) = (delta_bar(zstart(2)-1+j) * Csigma)**2.0 * xnu_sgs3(i,j,k)
			else
				xnu_sgs3(i,j,k) = (delta_bar(1) * Csigma)**2.0 * xnu_sgs3(i,j,k)
			endif
		enddo
	enddo
enddo   
return
end subroutine xnu_sigma

subroutine q_integral(r,delta,cos_phi,q)
	USE param
	USE variables
	USE decomp_2d
	USE mymath
	
	implicit none
	
	real(mytype),intent(IN) :: r,delta,cos_phi
	real(mytype) :: G1,C1,d,f,g,s,numerator,denominator,sigma,delta_eff
	real(mytype),intent(OUT) :: q
	
	delta_eff = min(delta,2.*dy)
	sigma = 1.0-cos_phi**2.0
	G1 = pi*(2.0 - sigma/3.0 - sigma**2/12.0 - &
			25.0*sigma**3.0/648.0-175.0*sigma**4.0/7776.0)
    C1 = 1.432! pi / 2**(2./3.) / sqrt(3.0) / gamma(4./3) ** 2.0
    d = r/delta_eff
    f = 2.710301 !3.0*pi**(7./3)/16
    g = C1*G1
    s = 2.0-sigma
    numerator = f*s*g*(d**2.0)
    denominator = g + f*s*d**(4./3.0)
    q = numerator/denominator
!~     print *,f,s,g,d
    
    return
    
end subroutine q_integral

!********************************************************************
!
subroutine k_sgs_circular(ux,uy,uz,k_sgs,ix,iy,iz,evec_y,pencil)
! 
!********************************************************************
USE param
USE variables
USE decomp_2d
USE mymath

!Z pencils 

implicit none

integer :: ix,iy,iz,pencil
real(mytype),dimension(ix,iy,iz) :: ux,uy,uz,k_sgs,evec_y
real(mytype) :: r,f2,d,factor,delta
real(mytype) :: kc,q,K0eps23
integer :: i,j,k,ip,im,kp,km

r = sqrt(dx*dz)
if (pencil.eq.1) then
	d = dx
	k_sgs = 0.0 !initialize k_sgs, first call
endif
if (pencil.eq.3) then
	d = dz
endif
factor = (r/d) ** (2.0/3.0)
do i = 1,ix
	if (pencil.eq.1) then
		ip = i+1
		im = i-1
		if (ip > ix) then
			ip = 1
		endif
		if (im < 1) then
			im = ix
		endif	
	endif
	do j = 1,iy
		if (pencil.eq.1) then
			if (j+xstart(2) > xend(2)) then
				delta = yp(j+xstart(2)-1)-yp(j+xstart(2)-2)
			else
				delta = yp(j+xstart(2))-yp(j+xstart(2)-1)
			endif
		endif
		if (pencil.eq.3) then
			if (j+zstart(2) > zend(2)) then
				delta = yp(j+zstart(2)-1)-yp(j+zstart(2)-2)
			else
				delta = yp(j+zstart(2))-yp(j+zstart(2)-1)
			endif
		endif
		kc = 3.1415 / delta
		do k = 1,iz
				if (pencil.eq.3) then
					kp = k+1
					km = k-1
					if (kp > iz) then
						kp = 1
					endif
					if (km < 1) then
						km = iz
					endif	
				endif
			if (pencil.eq.1) then
				f2 =      factor*(ux(i,j,k) - ux(im,j,k)) ** 2.0
				f2 = f2 + factor*(ux(i,j,k) - ux(ip,j,k)) ** 2.0
				f2 = f2 + factor*(uy(i,j,k) - uy(im,j,k)) ** 2.0
				f2 = f2 + factor*(uy(i,j,k) - uy(ip,j,k)) ** 2.0
				f2 = f2 + factor*(uz(i,j,k) - uz(im,j,k)) ** 2.0
				f2 = f2 + factor*(uz(i,j,k) - uz(ip,j,k)) ** 2.0
			endif
			if (pencil.eq.3) then
				f2 =      factor*(ux(i,j,k) - ux(i,j,km)) ** 2.0
				f2 = f2 + factor*(ux(i,j,k) - ux(i,j,kp)) ** 2.0
				f2 = f2 + factor*(uy(i,j,k) - uy(i,j,km)) ** 2.0
				f2 = f2 + factor*(uy(i,j,k) - uy(i,j,kp)) ** 2.0
				f2 = f2 + factor*(uz(i,j,k) - uz(i,j,km)) ** 2.0
				f2 = f2 + factor*(uz(i,j,k) - uz(i,j,kp)) ** 2.0
			endif
			f2 = f2 / 4.0
			call q_integral(r, delta, evec_y(i,j,k), q )
			K0eps23 = 3.1415 * f2/(2.0 * delta**(2.0/3.0) * q)
			k_sgs(i,j,k) = k_sgs(i,j,k) + 3.0/2.0 * K0eps23 * kc ** (-2./3.0)
		enddo
	enddo
enddo

return

end subroutine k_sgs_circular


!********************************************************************
!
subroutine e_svm_v(duxdx3,duydx3,duzdx3,duxdy3,duydy3,duzdy3,&
		 duxdz3,duydz3,duzdz3,e_svm_x3,e_svm_y3,e_svm_z3)
! 
!********************************************************************
! calculate SGS vortex orientation according to resolved vorticity 
USE param
USE variables
USE decomp_2d
USE mymath

!Z pencils 

implicit none

real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: duxdx3,duydx3,duzdx3,&
		duxdy3,duydy3,duzdy3,&
		duxdz3,duydz3,duzdz3
real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: e_svm_x3,e_svm_y3,e_svm_z3
real(mytype) :: f
integer :: i,j,k


e_svm_x3(:,:,:) = duzdy3(:,:,:) - duydz3(:,:,:)
e_svm_y3(:,:,:) = duxdz3(:,:,:) - duzdx3(:,:,:)
e_svm_z3(:,:,:) = duydx3(:,:,:) - duxdy3(:,:,:)

!normalization
do i=1, zsize(1)
	do j=1, zsize(2)
		do k=1, zsize(3)
			f =     e_svm_x3(i,j,k)**2.0
			f = f + e_svm_y3(i,j,k)**2.0
			f = f + e_svm_z3(i,j,k)**2.0
			f = sqrt(f)
			e_svm_x3(i,j,k) = e_svm_x3(i,j,k) / f
			e_svm_y3(i,j,k) = e_svm_y3(i,j,k) / f
			e_svm_z3(i,j,k) = e_svm_z3(i,j,k) / f
		enddo
	enddo
enddo

return



end subroutine e_svm_v

!********************************************************************
!
subroutine e_svm_s(duxdx3,duydx3,duzdx3,duxdy3,duydy3,duzdy3,&
		 duxdz3,duydz3,duzdz3,e_svm_x3,e_svm_y3,e_svm_z3)
! 
!********************************************************************
! calculate SGS vortex orientation according to resolved vorticity 
USE param
USE variables
USE decomp_2d
USE mymath

!Z pencils 

implicit none

real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: duxdx3,duydx3,duzdx3,&
		duxdy3,duydy3,duzdy3,&
		duxdz3,duydz3,duzdz3
real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: e_svm_x3,e_svm_y3,e_svm_z3
real(mytype) :: f,x1,x2,x3,S11,S12,S13,S22,S23,S33
real(mytype),dimension(3,3) :: v
integer :: i,j,k

do i=1, zsize(1)
	do j=1, zsize(2)
		do k=1, zsize(3)
			!constructing S_ij
			S11 = duxdx3(i,j,k)
			S22 = duydy3(i,j,k)
			S33 = duzdz3(i,j,k)
			S12 = 0.5*(duxdy3(i,j,k) + duydx3(i,j,k))
			S13 = 0.5*(duxdz3(i,j,k) + duzdx3(i,j,k))
			S23 = 0.5*(duydz3(i,j,k) + duzdy3(i,j,k))
			!calculate eigenvalues/eigenvectors
			call eig_sym_smith(S11,S12,S13,S22,S23,S33,x1,x2,x3)
			call eigvec(S11,S12,S13,S22,S23,S33,x1,x2,x3,v)
			e_svm_x3(i,j,k) = v(1,1)
			e_svm_y3(i,j,k) = v(2,1)
			e_svm_z3(i,j,k) = v(3,1)
			!normalization
			f =     e_svm_x3(i,j,k)**2.0
			f = f + e_svm_y3(i,j,k)**2.0
			f = f + e_svm_z3(i,j,k)**2.0
			f = sqrt(f)
			e_svm_x3(i,j,k) = e_svm_x3(i,j,k) / f
			e_svm_y3(i,j,k) = e_svm_y3(i,j,k) / f
			e_svm_z3(i,j,k) = e_svm_z3(i,j,k) / f
		enddo
	enddo
enddo

return



end subroutine e_svm_s



