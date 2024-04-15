Program Driver_LinAl

  use LinAl

  implicit none
  
  integer :: i,j
  real, parameter :: pi = 4.*atan(1.0), tol = 10.0**(-14)
  integer, parameter :: ma1 = 4, ma2 = 3
  real, dimension(4, 4) :: A1, A4, As, Eigmat, D4
  real, dimension(3, 3) :: A2, A3, D, A8
  real, dimension(3, 1) :: B8, X8
  integer, dimension(3) :: P
  real, dimension(4, 1) :: EA4, eigvec
  real :: norm
  logical :: bool

  A8(1, :) = (/1.0, 4.0, 9.0/)
  A8(2, :) = (/1.0, 2.0, 3.0/)
  A8(3, :) = (/1.0, 1.0, 1.0/)

  B8(1, 1) = 0.0
  B8(2, 1) = 0.0
  B8(3, 1) = 1.0

  P = 0

  A1(1, :) = (/5.0, 4.0, 1.0, 1.0/)
  A1(2, :) = (/4.0, 5.0, 1.0, 1.0/)
  A1(3, :) = (/1.0, 1.0, 4.0, 2.0/)
  A1(4, :) = (/1.0, 1.0, 2.0, 4.0/)

  A2(1, :) = (/3.0, 1.0, 0.0/)
  A2(2, :) = (/1.0, 2.0, 1.0/)
  A2(3, :) = (/0.0, 1.0, 1.0/)

  A3 = A2

  A4(1, :) = (/2.0, 1.0, 3.0, 4.0/)
  A4(2, :) = (/1.0, -3.0, 1.0, 5.0/)
  A4(3, :) = (/3.0, 1.0, 6.0, -2.0/)
  A4(4, :) = (/4.0, 5.0, -2.0, -1.0/)

  EA4(:, 1) = (/-8.0286, 7.9329, 5.6689, -1.5732/)

  print *, " "
  print *, "--------------------------------------------------------------"
  print *, " "
  print *, " Random quesiton from 213B"
  print *, " "

    call LU(A8, 3, bool, P)
    call LUsolve(A8, 3, B8, X8, 1, P)

    call printmat(X8, 3, 1)

  print *, " "
  print *, "--------------------------------------------------------------"
! print *, " "
! print *, "Question 1: Hessenburg Form (Tridiagonalization)"
! print *, " "

! print *, "A before tridiagonalization"
! call printmat(A1, ma1, ma1)
! print *, " "

! call tridiagonal(A1, ma1)

! print *, "A after Tridiagonalization"
! call printmat(A1, ma1, ma1)
! print *, " "

! do i = 1, ma1
!   call twonorm(A1(:, i), norm)
!   print *, "Two norm of column "//trim(str(i))//":", norm
! end do

! print *, " "
! print *, "--------------------------------------------------------------"
! print *, " "
! print *, "Question 2: QR Algorithm with and without Shifts"
! print *, " "

! print *, "Matrix A before QR w/o shift"
! call printmat(A2, ma2, ma2)  
! print *, " "

! call eigQR(A2, ma2, .false., 10.d-14)

! call diag(A2, ma2, D)

! print *, "Matrix D after QR w/o shift"
! call printmat(D, ma2, ma2)  
! print *, " "

! print *, "Matrix A before QR w shift"
! call printmat(A3, ma2, ma2)  
! print *, " "

! call eigQR(A3, ma2, .true., 10.d-14)

! call diag(A3, ma2, D)

! print *, "Matrix D after QR w shift"
! call printmat(D, ma2, ma2)  
! print *, " "

! print *, " "
! print *, "--------------------------------------------------------------"
! print *, " "
! print *, "Question 3: Inverse Iteration Method"
! print *, " "

! print *, "Matrix A before Inverse Iteration"
! call printmat(A4, ma1, ma1)  
! print *, " "

! As = A4
! print *, "Calculating Eigenvalues of A with QR w shift with high accuracy"

! call eigQR(A4, ma1, .false., 10.d-14)
! call diag(A4, ma1, D4)
!   
! call printmat(D4, ma1, ma1)

! A4 = As

! print *, "Eigenvalues of A" 
! call printmat(EA4, ma1, 1)  
! print *, " "

! do i = 1, ma1

!   A4 = As

!   call inviter(A4, D4(i, i), eigvec, ma1, 10.d-14)

!   Eigmat(:, i) = eigvec(:, 1)
!   
! end do
! 
! print *, " "
! print *, "Eigenvectors (in columns)"
! call printmat(eigmat, ma1, ma1)
! print *, " "
! do i = 1, ma1
!   call twonorm(matmul(As, eigmat(:, i)) - D4(i, i)*eigmat(:, i), norm)
!   print *, "Error in eigenvector "//trim(str(i))//" calculation norm: ", norm
! end do

! print *, "--------------------------------------------------------------"


End Program Driver_LinAl
