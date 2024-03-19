Program Driver_LinAl

  use LinAl

  implicit none
  
  integer :: i,j,k
  real, parameter :: tol = 10.d-5
  integer, parameter :: ma = 10, nb = 1
  real, dimension(ma, ma) :: A, eye, As
  real, dimension(ma, nb) :: B, X, Bs
  real :: norm
  character :: jors

  print *, "--------------------------------------------------------------"
  print *, " "
  print *, "Question 7: Iterative Methods for solving systems of equations"
  print *, " "

  A = 1.0
  call ident(eye, ma)
  A = A - eye

 print *, "Please enter a value for A(i, i): "
 read *, norm

    
  do i = 1, ma
    B(i, 1) = i
    A(i, i) = norm
  end do
  X = 0.0

  As = A
  Bs = B

 print *, "Would you like the code to perform Gauss-Jordan (J), Gauss-Seidel (S), or Conjugate Gradient (C)?"
 read *, jors

  if (jors .eq. "J" .or. jors .eq.  "j") then

    print *, " "
    print *, "--------------------------------------------------------------"
    print *, " "

    call GJ(A, B, X, ma, nb, tol, .false., 0)

    print *, "This is the solution matrix X"
    call printmat(X, ma, nb)
    print *, " "

    B = B - matmul(As, X)
    call twonorm(B(:, 1), norm)
    print *, "ERROR :", norm
    
  else if (jors .eq. "S" .or. jors .eq. "s") then

    A = As
    B = Bs
    X = 0.0 

    call GS(A, B, X, ma, nb, tol, .false., 0)

    print *, "This is the solution matrix X"
    call printmat(X, ma, nb)
    print *, " "

    B = B - matmul(As, X)
    call twonorm(B(:, 1), norm)
    print *, "ERROR :", norm

  else if (jors .eq. "C" .or. jors .eq. "c") then
    A = As
    B = Bs
    X = 0.0 
    call CG(A, B, X, ma, nb, tol)

    print *, "This is the solution matrix X"
    call printmat(X, ma, nb)
    print *, " "

    B = B - matmul(As, X)
    call twonorm(B(:, 1), norm)
    print *, "ERROR :", norm
   
   else 
     print *, "That was not one of the listed options please restart the code and try again "
   end if


End Program Driver_LinAl
