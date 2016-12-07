
USE param
USE variables
USE decomp_2d

!********************************************************************
!
subroutine show_single_value1(array,i,j,k)
! 
!********************************************************************
integer i,j,k
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: array
if 
