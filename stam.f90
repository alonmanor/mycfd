program stam

implicit none

!integer(8),dimension(3:8,2:5,2:6) :: a
real(8),dimension(3,4,5) :: a,b
integer i,j,k,ijk

do i=1,3
do j=1,4
do k=1,5
a(i,j,k) = i*j*k
enddo
enddo
enddo
!~ do ijk=1,3*4*5
!~ write(*,*) a(ijk,1,1)
!~ enddo
!~ b(:,:,:) = a(:,:,:)*0.5*a(:,:,:)*2.0
!~ do ijk=1,3*4*5
!~ write(*,*) a(ijk,1,1),b(ijk,1,1)
!~ enddo
b = a(:,:,:)**2.0
do i=1,3
do j=1,4
do k=1,5
print *,i,j,k,b(i,j,k)
enddo
enddo
enddo
end
