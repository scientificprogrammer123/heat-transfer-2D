# heat-transfer-2D

This is a fortran90 program to evaluate 2-D heat transfer using Laplace equation, 
via the Jacobi method.
This program computes how the temperature at a boundary propogates to another boundary 
with a different temperature.

The Laplace equation in 2-D is of the form:
  div div f(x,y) = 0

Expanding it out gives:
  partial^2 f / partial^2 x + partial^2 f / partial^2 y = 0
  
The Jacobi method works by setting up the boundary conditions, then iteratively 
computing the values between the boundaries until the values converge.

The boundary conditions for this problem are:
  f(x=0,   y    ) = 0
  f(x,     y=0  ) = 0
  f(x=0.5, y    ) = 200y
  f(x,     y=0.5) = 200x
  
In each iteration, the value at a particular location is computed by taking
the value of the 4 point around the the point of interest at the previous 
iteration, and averaging them:
  do i=1,NM1,1
    do j=1,N-1,1
      F(i,j) = 0.25*( &
                      Fold(i-1,j) + Fold(i+1,j) + &
                      Fold(i,j-1) + Fold(i,j+1)   &
                    )
    enddo
  enddo

The exact solution of this equation is:
  Fexact(i,j) = 400 * X(i) * Y(j)
  
At the end, the program outputs the number of iterations taken, and the value
of the numerically computed solution and the exact solution.

To run the program, go to a termina and type: gfortran heat_transfer_2D.f90
This program has been verified on a computer running Scientific Linux 6.9.
