program stam

implicit none

!integer(8),dimension(3:8,2:5,2:6) :: a
integer(8),dimension(3:8,2:5,2:6) :: a
integer i,j,k,ijk

do i=3,8
do j=2,5
do k=2,6
a(i,j,k) = i*j*k
enddo
enddo
enddo
do ijk=1,4*4*5
write(*,*) a(ijk,1,1)
enddo

end
