! From COMMAND line:
! gfortran-mp-5 -c lu_solver.f90
! gfortran-mp-5 -o blasius blasius.f90 lu_solver.o 
! ./blasius
! gnuplot f.gpl

program BLASIUS


use LU

implicit none

integer :: i,j,k,pointsETA,iterations
integer :: RC, d
real*8 :: stopETA, h ! real*8 double precision, 64 bit number (not 32 which is single)

parameter(pointsETA=31,stopETA=5.0,iterations=20)

real*8, dimension(pointsETA) :: eta,f1,f2,f3
real*8, dimension(pointsETA-1) :: halfF1,halfF2,halfF3
real*8, dimension(3*pointsETA) :: b ! breaking into a 3x3 system @ each point
integer, dimension(3*pointsETA) :: INDX
real*8, dimension(3*pointsETA,3*pointsETA) :: A


h = stopETA/(pointsETA-1.0) ! step size

do i=1,pointsETA
  eta(i) = (i-1)*h !(i-1) starts at zero
end do

! these are our initial guesses that satisfy the bcs
f1 = eta  ! this is f(eta)
f2 =1.0-exp(-eta)  ! this is f'(eta)
f3 =(1.0-tanh(eta -3.0))/6.0 ! this is f''(eta)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do j=1,iterations
do i=1,pointsETA-1 ! evaluate our guesses/solutions at their centered values
  halfF1(i) = ( f1(i+1)+f1(i) ) / 2.0
  halfF2(i) = ( f2(i+1)+f2(i) ) / 2.0
  halfF3(i) = ( f3(i+1)+f3(i) ) / 2.0
end do

A=0.0
b=0.0
do i=1, pointsETA-1
  k=(i-1)*3
  A( k+1, (/k+1,k+2,k+3,k+4,k+5,k+6/) ) = (/-1/h, -1/2.0d0, 0.0d0, 1/h, -1/2.0d0, 0.0d0/)

  A( k+2, (/k+1,k+2,k+3,k+4,k+5,k+6/) ) = (/0.0d0, -1/h, -1/2.0d0, 0.0d0,1/h, -1/2.0d0/)

  A( k+3,   (/k+1,k+2,k+3,k+4,k+5,k+6/) ) = (/halfF3(i)/4.0, 0.0d0, halfF1(i)/4.0-1/h, &
                                                halfF3(i)/4.0, 0.0d0, halfF1(i)/4.0+1/h/)

  b( (/k+1,k+2,k+3/) ) = (/ -(f1(i+1)-f1(i))/h +halfF2(i), -(f2(i+1)-f2(i))/h +halfF3(i), &
                      -( f3(i+1)-f3(i) )/h - halfF1(i)*halfF3(i)/2.0 /)
end do

! do the last 3 rows to assign bcs
k=pointsETA*3-3
A(k+1,1)=1 ! df1(1)=0
A(k+2,2)=1 ! df2(1)=0
A(k+3,pointsETA*3-1)=1 ! df2(pointsEta)=0


! take A and b and put solution into b
call LUDCMP(A, pointsETA*3,INDX,d,RC)
if (RC.eq.0) then
  call LUBKSB(A, pointsETA*3,INDX,b)
end if

do i=1,pointsETA
  k=(i-1)*3
  f1(i)=f1(i)+b(k+1)
  f2(i)=f2(i)+b(k+2)
  f3(i)=f3(i)+b(k+3)
end do

end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! clean up round off error, to prevent negative values when plotting
f1(1)=0
f2(1)=0

102 FORMAT(f12.6,1x,f12.6,1x,f12.6,1x,f12.6)

! output with each value having its own line
open(1,file = 'f.dat', status = 'replace')
!write(1,*)'#    eta                  f1                    f2                    f3'
write(1,*)'#    eta          f1          f2            f3'
do i = 1, pointsETA
  write(1,102) eta(i),f1(i),f2(i),f3(i)
end do


end program







