!---------------------------------------------------------------------------
! This is a fortran90 program to evaluate the 2 dimension laplace equation:
!   Fxx + Fyy = 0
!   {F(Xi+1,Yj)-2*(Fxi,Yj)+F(Xi-1,Yj)}/dX^2 + 
!   {F(Xi,Yj+1)-2*F(Xi,Yj)+F(Xi,Yj-1)}/dY^2 = 0
!   for 1/2 >= {x,y} >= 0
!
! The boundary conditions are:
!   f(x=0,y) = 0
!   f(x,y=0) = 0
!   f(x=1/2,y) = 200y
!   f(x,y=1/2) = 200x
! In this example the boundary condition represents temperature in degC.
! 
! Mapping:
!   i,j goes from 1 to N-1
!   i,j at the edges (0 and N) are constrained by the boundary conditions
!
! 2017-03-15
!      - working now, 171 iterations needed
!
! 2017-03-06
!      - Initial try.
!      - jacobi method:
!        1. fill f square with boundary conditions,
!        2. iterate through each interior square, where:
!           U(i,j) = (1/4)*{(U(i-1,j)+U(i+1,j)+U(i,j-1)+U(i,j+1)}
!        3. how many iterations does it take to converge?
!        4. compare error.
!      - based on professor xxx's code,
!      - written by xxx.
!---------------------------------------------------------------------------

program part3_jacobi_10
  implicit none

  integer, parameter :: N=10,NM1=N-1,M=NM1*NM1

  real ALX,ALY,dx,dy
  real tolerance, error2, error

  real, dimension (0:N,0:N)  :: F,Fold,Fexact
  real, dimension (0:N)      :: X,Y
  
  integer i,j,k,iter
  
  open (2, file = 'part3_jacobi_10.out')

  tolerance = 1.0e-4
  error=0; error2=0

  ALX = 0.5 ! boundary
  ALY = 0.5 
  dx = ALX/N ! discretize
  dy = ALY/N
  do i=0,N,1
     X(i) = i*dx
     Y(i) = i*dy
  enddo

  ! setup Fexact matrix
  do i=0,N,1
     do j=0,N,1
        Fexact(i,j) = 400*X(i)*Y(j) ! exact solution
     enddo                         
  enddo

  ! setup F matrix, boundary conditions
  F=0
  do i=0,N,1
     F(i,0) = 0
     F(0,i) = 0
     F(N,i) = 200*Y(i)
     F(i,N) = 200*X(i)
  enddo

  ! iterate using Jacobi method
  Fold = F
  do iter=1,200,1
     do i=1,NM1,1
        do j=1,N-1,1
           F(i,j) = 0.25*( &
                          Fold(i-1,j) + Fold(i+1,j) + &
                          Fold(i,j-1) + Fold(i,j+1)   &
                         )
        enddo
     enddo

     ! calculate error
     error2 = 0 
     do i=1,NM1,1
        do j=1,NM1,i
           error2 = error2 + (F(i,j) - Fold(i,j))**2 
        enddo
     enddo
     error = sqrt(error2)/NM1 ! this calculates the L2 norm error / #elements
     write(2,*), iter, 'l2 norm error / #elements:',  error
     
     Fold = F
     if (error.le.tolerance) goto 1000
  enddo
1000 continue

  ! write F and Fexact to output
  write(2,*), 'F:'
  do i=0,N,1
     write(2,*) (F(i,j),j=0,N)
  enddo
  write(2,*), 'Fexact:'
  do i=0,N,1
     write(2,*) (Fexact(i,j),j=0,N)
  enddo
   
end program
