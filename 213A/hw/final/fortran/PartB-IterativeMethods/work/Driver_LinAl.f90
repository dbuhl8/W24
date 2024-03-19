Program Driver_LinAl

  use LinAl

  implicit none
  
  integer :: i,j,k
  real, parameter :: tol = 10.d-5
  integer, parameter :: ma = 10, nb = 1, nd = 5, mc = 100
  real, dimension(ma, ma) :: A, eye, As
  real, dimension(ma, nb) :: B, X, Bs
  real, dimension(mc, mc) :: C, eye2, Cs
  real, dimension(mc, nb) :: B2, X2, B2s
  real, dimension(nd) :: d
  integer, dimension(nd) ::di
  real :: norm
  character :: jors

  print *, "--------------------------------------------------------------"
  print *, " "
  print *, "Question 7: Iterative Methods for solving systems of equations"
  print *, " "

  A = 1.0
  call ident(eye, ma)
  A = A - eye
  d = (/ 2., 5., 10., 100., 1000. /)
  di = (/ 2, 5, 10, 100, 1000 /)

! print *, "Please enter a value for A(i, i): "
! read *, norm

do k = 1, nd 
    
  do i = 1, ma
    B(i, 1) = i
    A(i, i) = d(k)
  end do
  X = 0.0

  As = A
  Bs = B

! print *, "Would you like the code to perform Gauss-Jordan (J), Gauss-Seidel (S), or Conjugate Gradient (C)?"
! read *, jors

!  if (jors .eq. "J" .or. jors .eq.  "j") then

    print *, " "
    print *, "--------------------------------------------------------------"
    print *, " "
    print *, "Results for D: ", d(k)
    print *, " "

    call GJ(A, B, X, ma, nb, tol, .true., di(k))

    print *, "This is the solution matrix X"
    call printmat(X, ma, nb)
    print *, " "

    B = B - matmul(As, X)
    call twonorm(B(:, 1), norm)
    print *, "ERROR :", norm
    
!  else if (jors .eq. "S" .or. jors .eq. "s") then

    A = As
    B = Bs
    X = 0.0 

    call GS(A, B, X, ma, nb, tol, .true., di(k))

    print *, "This is the solution matrix X"
    call printmat(X, ma, nb)
    print *, " "

    B = B - matmul(As, X)
    call twonorm(B(:, 1), norm)
    print *, "ERROR :", norm

!  else if (jors .eq. "C" .or. jors .eq. "c") then
    A = As
    B = Bs
    X = 0.0 
    call CG(A, B, X, ma, nb, tol, .false.)

    print *, "This is the solution matrix X"
    call printmat(X, ma, nb)
    print *, " "

    B = B - matmul(As, X)
    call twonorm(B(:, 1), norm)
    print *, "ERROR :", norm
   
!   else 
!     print *, "That was not one of the listed options please restart the code and try again "
!   end if

end do

  do i = 1, ma
    B(i, 1) = i
    A(i, i) = i
  end do
  X = 0.0

  As = A
  Bs = B
  C = 1.0
  call ident(eye2, mc)
  C = C - eye2

  do i = 1, mc
    B2(i, 1) = i
    C(i, i) = i
  end do
  Cs = C
  B2s = B2
  X2 = 0.0

    print *, " "
    print *, "--------------------------------------------------------------"
    print *, " "
    print *, "Results for 10 x 10 matrices with A(i, i) = i "
    print *, " "

    call GJ(A, B, X, ma, nb, tol, .false., 0)

    print *, "This is the solution matrix X"
    call printmat(X, ma, nb)
    print *, " "

    B = B - matmul(As, X)
    call twonorm(B(:, 1), norm)
    print *, "ERROR :", norm
    
!  else if (jors .eq. "S" .or. jors .eq. "s") then

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

!  else if (jors .eq. "C" .or. jors .eq. "c") then
    A = As
    B = Bs
    X = 0.0 
    call CG(A, B, X, ma, nb, tol, .false.)

    print *, "This is the solution matrix X"
    call printmat(X, ma, nb)
    print *, " "

    B = B - matmul(As, X)
    call twonorm(B(:, 1), norm)
    print *, "ERROR :", norm

!   A = As
!   B = Bs
!   X = 0.0 
!   call CG(A, B, X, ma, nb, tol, .true.)

!   print *, "This is the solution matrix X"
!   call printmat(X, ma, nb)
!   print *, " "

!   B = B - matmul(As, X)
!   call twonorm(B(:, 1), norm)
!   print *, "ERROR :", norm

    print *, " "
    print *, "--------------------------------------------------------------"
    print *, " "
    print *, "Results for 100 x 100 matrices with A(i, i) = i "
    print *, " "

  
    call GJ(C, B2, X2, mc, nb, tol, .false., 0)

!   print *, "This is the solution matrix X"
!   call printmat(X2, mc, nb)
!   print *, " "

    B2 = B2 - matmul(Cs, X2)
    call twonorm(B2(:, 1), norm)
    print *, "ERROR :", norm
    
    C = Cs
    B2 = B2s
    X2 = 0.0 

    call GS(C, B2, X2, mc, nb, tol, .false., 0)

!   print *, "This is the solution matrix X"
!   call printmat(X2, mc, nb)
!   print *, " "

    B2 = B2 - matmul(Cs, X2)
    call twonorm(B2(:, 1), norm)
    print *, "ERROR :", norm

    C = Cs
    B2 = B2s
    X = 0.0 
    call CG(C, B2, X2, mc, nb, tol, .false.)

!   print *, "This is the solution matrix X"
!   call printmat(X2, mc, nb)
!   print *, " "

!   B2 = B2 - matmul(Cs, X2)
!   call twonorm(B2(:, 1), norm)
!   print *, "ERROR :", norm

!   C = Cs
!   B2 = B2s
!   X = 0.0 
!   call CG(C, B2, X2, mc, nb, tol, .true.)

!   print *, "This is the solution matrix X"
!   call printmat(X2, mc, nb)
!   print *, " "

!   B2 = B2 - matmul(Cs, X2)
!   call twonorm(B2(:, 1), norm)
!   print *, "ERROR :", norm

End Program Driver_LinAl
