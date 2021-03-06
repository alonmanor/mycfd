module mymath
!contains my added mathematical subroutines
contains

SUBROUTINE SWAP(A, B)
		USE param
		USE variables
		USE decomp_2d
       IMPLICIT NONE
       REAL(mytype), INTENT (INOUT)      :: A, B
       REAL(mytype)                      :: TEMP
               TEMP = A ; A = B ; B = TEMP
       END SUBROUTINE SWAP
		
subroutine sort3(x1,x2,x3)
		USE param
		USE variables
		USE decomp_2d
       IMPLICIT NONE
       REAL(mytype), INTENT (INOUT) :: x1,x2,x3
        if (x1<x2) then
			call swap(x1,x2)
		endif
		if (x1<x3) then
			call swap(x1,x3)
		endif
		if (x2<x3) then
			call swap(x2,x3)
		endif
       END SUBROUTINE sort3
 
 subroutine det(a11,a12,a13,a21,a22,a23,a31,a32,a33,d)
! return a 3x3 determinant in d
		USE decomp_2d
		real(mytype), intent(in) :: a11,a12,a13,a21,a22,a23,a31,a32,a33
		real(mytype), intent(out) :: d
		
		d = a11*(a22*a33-a23*a32) - a12*(a21*a33-a23*a31) + a13*(a21*a32-a22*a31)
 end subroutine det
 
 subroutine eigvec(G11,G12,G13,G22,G23,G33,x1,x2,x3,v)
 ! given the symmetric 3X3 matrix Gij and it's eigenvalues xi
 ! return the eigenvectors in the columns of v
		USE decomp_2d
       IMPLICIT NONE
       REAL(mytype), INTENT (IN) :: G11,G12,G13,G22,G23,G33
       REAL(mytype), INTENT (IN) :: x1,x2,x3
       REAL(mytype), dimension(3,3), INTENT (OUT) :: v
       REAL(mytype) ::d,d1,d2,d3
       REAL(mytype), dimension(3):: x
       real(mytype),parameter::zero=0.0
       integer i,j
       REAL(mytype) ::B11,B22,B33
       x(1) = x1;x(2)=x2;x(3)=x3
       do i=1,3
			B11 = G11-x(i)
			B22 = G22-x(i)
			B33 = G33-x(i)
			v(3,i)=1.0
			v(2,i)=(G13*G12/B11-G23)/(B22-G12*G12/B11)
			v(1,i)=(-G12*v(2,i)-G13)/B11
			v(:,i) = v(:,i) / sqrt(v(1,i)**2.0 + v(2,i)**2.0 + v(3,i)**2.0)
		enddo
		return
end subroutine eigvec
       
		
  
  
 subroutine eig_sym_smith(G11,G12,G13,G22,G23,G33,x1,x2,x3)
 ! for the symmetric 3X3 matrix Gij, return 
 ! the 3 eigenvalues, such that x1>x2>x3
 !   based on Smith, 1961
		USE param
		USE variables
		USE decomp_2d
       IMPLICIT NONE
       REAL(mytype), INTENT (IN) :: G11,G12,G13,G22,G23,G33
       REAL(mytype), INTENT (OUT) :: x1,x2,x3
       REAL(mytype) ::Q11,Q22,Q33,q,m,eta,p
       real(mytype), parameter::three=3.0

		m = (G11+G22+G33)/3.0
		Q11=G11-m
		Q22=G22-m
		Q33=G33-m
		q = 0.5*(Q11*(Q22*Q33-G23**2.0) - G12*(G12*Q33-G23*G13) + G13*(G12*G23-Q22*G13))
		
		p = 1.0/6.0*(Q11**2.0+Q22**2.0+Q33**2.0+2*G12**2.0+2*G13**2.0+2*G23**2.0)
		eta=1.0/3.0*atan(sqrt(p**3.0-q**2.0)/q)
!~ 		print *,'m=',m
!~ 		print *,'q=',q
!~ 		print *,'p=',p
!~ 		print *,'eta=',eta
		x1 = m+2*sqrt(p)*cos(eta)
		x2 = m-sqrt(p)*(cos(eta) + sqrt(three)*sin(eta))
		x3 = m-sqrt(p)*(cos(eta) - sqrt(three)*sin(eta))
        call sort3(x1,x2,x3)
        return
end subroutine eig_sym_smith


subroutine nanhunt(A,nx,ny,nz)
	USE decomp_2d
	real(mytype),dimension(:,:,:)::A
	integer nx,ny,nz
	integer i,j,k
	do i=1,nx
	do j=1,ny
	do k=1,nz
	if (A(i,j,k).ne.A(i,j,k)) then
	print *,'proc. ',nrank,'Nan at ',i,j,k
	stop
	endif
	enddo
	enddo
	enddo
end subroutine nanhunt

! show print the element with the absolute index i,j,k. Xpencils
subroutine showval1(a, i, j, k)
	USE decomp_2d
	real(mytype),dimension(xstart(1):xend(1),xstart(2):xend(2),&
						xstart(3):xend(3)) :: a
	integer i,j,k
	if ( (i.ge.xstart(1)).and.(i.le.xend(1)) ) then
		if ( (j.ge.xstart(2)).and.(j.le.xend(2)) ) then
			if ( (k.ge.xstart(3)).and.(k.le.xend(3)) ) then
				print *,'Xpencil- nrank=',nrank,'relative ijk',i+1-xstart(1),&
				j+1-xstart(2),k+1-xstart(3),'abs ijk',i,j,k,a(i,j,k)
			endif
		endif
	endif
end subroutine showval1

! show print the element with the absolute index i,j,k. Ypencils
subroutine showval2(a, i, j, k)
	USE decomp_2d
	real(mytype),dimension(ystart(1):yend(1),ystart(2):yend(2),&
						ystart(3):yend(3)) :: a
	integer i,j,k
	if ( (i.ge.ystart(1)).and.(i.le.yend(1)) ) then
		if ( (j.ge.ystart(2)).and.(j.le.yend(2)) ) then
			if ( (k.ge.ystart(3)).and.(k.le.yend(3)) ) then
				print *,'Ypencil- nrank=',nrank,'relative ijk',i+1-ystart(1),&
				j+1-ystart(2),k+1-ystart(3),'abs ijk',i,j,k,a(i,j,k)
			endif
		endif
	endif
end subroutine showval2

! show print the element with the absolute index i,j,k. Zpencils
subroutine showval3(a, i, j, k)
	USE decomp_2d
	real(mytype),dimension(zstart(1):zend(1),zstart(2):zend(2),&
						zstart(3):zend(3)) :: a
	integer i,j,k
	if ( (i.ge.zstart(1)).and.(i.le.zend(1)) ) then
		if ( (j.ge.zstart(2)).and.(j.le.zend(2)) ) then
			if ( (k.ge.zstart(3)).and.(k.le.zend(3)) ) then
				print *,'Zpencil- nrank=',nrank,'relative ijk',i+1-ystart(1),&
				j+1-ystart(2),k+1-ystart(3),'abs ijk',i,j,k,a(i,j,k)
			endif
		endif
	endif
end subroutine showval3

end module mymath
