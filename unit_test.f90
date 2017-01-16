

 subroutine test_EV_smith()
    USE param
    USE variables
    USE decomp_2d
    use mymath
   IMPLICIT NONE
   REAL(mytype):: G11,G12,G13,G22,G23,G33
   REAL(mytype) :: x1,x2,x3,y1,y2,y3
   integer i
   real(mytype),dimension(3,3)::v
   
   G11=4.0; G22=5.0;G33=-2
   G23=6.0; G12=-10;G13=0.2
    call eig_sym_smith(G11,G12,G13,G22,G23,G33,x1,x2,x3)
    print *,G11,G12,G13
    print *,G12,G22,G23
    print *,G13,G23,G33
    print *,' '
    print *,x1,x2,x3
    y1 = abs(x1-15.573253857)
    y2 = abs(x2-6.0983473713811520E-002)
    y3 = abs(x3+8.6342373307883)
    
!~     print *,y1,y2,y3
    call eigvec(G11,G12,G13,G22,G23,G33,x1,x2,x3,v)
    print *,x1,(v(i,1),i=1,3)
    print *,x2,(v(i,2),i=1,3)
    print *,x3,(v(i,3),i=1,3)
    if ((y1 < 1.e-6).and.(y2 < 1.e-6).and.(y3<1.e-6)) then
    print *,'pass'
    endif
    
    
  end subroutine


program unit_test
call test_EV_smith()
end program
