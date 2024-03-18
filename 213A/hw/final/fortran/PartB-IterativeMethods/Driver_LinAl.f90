Program Driver_LinAl

  use LinAl

  implicit none
  
  integer :: i,j,k
  real, parameter :: tol = 10.d-5
  integer, parameter :: ma = 3, nb = 1
  real, dimension(ma, ma) :: A, eye, As
  real, dimension(ma, nb) :: B, X
  real :: norm
  character :: jors

  print *, "--------------------------------------------------------------"
  print *, " "
  print *, "Question 7: Iterative Methods for solving systems of equations"
  print *, " "

  A = 1.0
  call ident(eye, ma)
  A = A - eye

  do i = 1, ma
    print *, "Please enter a value for A("//trim(str(i))//trim(str(i))//"): "
    read *, A(i, i)
    B(i, 1) = i
  end do
  X = 0.0

  As = A

  print *, "Would you like the code to perform Gauss-Jordan (J) or Gauss-Seidel (S) ?"
  read *, jors

  if (jors .eq. "J" .or. jors .eq.  "j") then

    call GJ(A, B, X, ma, nb, tol)

    print *, "This is the solution matrix X"
    call printmat(X, ma, nb)
    print *, " "

    B = B - matmul(As, X)
    call twonorm(B(:, 1), norm)
    print *, "ERROR :", norm
    
  else if (jors .eq. "S" .or. jors .eq. "s") then

    call GS(A, B, X, ma, nb, tol)

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
