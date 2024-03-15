Program Driver_LinAl

  use LinAl

  implicit none
  
  integer :: i,j
  real, parameter :: pi = 4.*atan(1.0), tol = 10.0**(-14)
  integer, parameter :: ma1 = 4, ma2 = 3
  real, dimension(4, 4) :: A1
  real, dimension(3, 3) :: A2, A3, D
  real, dimension(4, 4) :: A4
  real :: norm


  A1(1, :) = (/5.0, 4.0, 1.0, 1.0/)
  A1(2, :) = (/4.0, 5.0, 1.0, 1.0/)
  A1(3, :) = (/1.0, 1.0, 4.0, 2.0/)
  A1(4, :) = (/1.0, 1.0, 2.0, 4.0/)

  A2(1, :) = (/3.0, 1.0, 0.0/)
  A2(2, :) = (/1.0, 2.0, 1.0/)
  A2(3, :) = (/0.0, 1.0, 1.0/)

  A3 = A2


  print *, " "
  print *, "--------------------------------------------------------------"
  print *, " "
  print *, "Question 1: Hessenburg Form (Tridiagonalization)"
  print *, " "

  print *, "A before tridiagonalization"
  call printmat(A1, ma1, ma1)
  print *, " "

  call tridiagonal(A1, ma1)

  print *, "A after Tridiagonalization"
  call printmat(A1, ma1, ma1)
  print *, " "

  do i = 1, ma1
    call twonorm(A1(:, i), norm)
    print *, "Two norm of column "//trim(str(i))//":", norm
  end do

  print *, " "
  print *, "--------------------------------------------------------------"
  print *, " "
  print *, "Question 2: QR Algorithm with and without Shifts"
  print *, " "

  print *, "Matrix A before QR w/o shift"
  call printmat(A2, ma2, ma2)  
  print *, " "

  call eigQR(A2, ma2, .false., 10.d-10)

  call diag(A2, ma2, D)

  print *, "Matrix D after QR w/o shift"
  call printmat(D, ma2, ma2)  
  print *, " "

  print *, "Matrix A before QR w shift"
  call printmat(A3, ma2, ma2)  
  print *, " "

  call eigQR(A3, ma2, .true., 10.d-10)

  call diag(A3, ma2, D)

  print *, "Matrix D after QR w shift"
  call printmat(D, ma2, ma2)  
  print *, " "


!    print *, "Matrix R, after QR"
!    call printmat(R, na, na)
!    print *, "Matrix Q, after QR"
!    call printmat(transpose(Q), ma, na)
!    print *, "Matrix QTQ, after QR"
!    call printmat(matmul(transpose(Q), Q), na, na)
!   call ident(eye, na)
!   print *, "Matrix I - QTQ"
!   call printmat(eye - matmul(transpose(Q), Q), na, na)
!   call frobnorm(eye - matmul(transpose(Q), Q), norm)
!   print *, "Frob Norm of I - QTQ: ", norm
!   print *, "Matrix A - QR"
!   call printmat(As - matmul(Q, R), ma, na)
!   call frobnorm(As - matmul(Q, R), norm)
!   print *, "Frob Norm of A - QR: ", norm

!   B1 = matmul(transpose(Q), B1)

!   call backsub(R, B1, X1, na, mb) 

!   E1 = matmul(As, X1) - Bs
!   print *, "Matrix X"
!   call printmat(X1, matb, nb)
!    print *, "Matrix E"
!    call printmat(E1, mb, nb)
!   
!   call frobnorm(E1, norm)
!   print *, "Frobenious norm of the error is: ", norm

!   regression(:, 2) = 0.0
!   do i = 1, na
!       regression(:, 2) = regression(:, 2) + (regression(:, 1)**(i-1))*X1(i, 1)
!   end do

!   open(15, file="qr.dat")
!   do i = 1, nsteps
!       write(15, "(2F10.3)") regression(i, :)
!   end do
!   close(15)


!else 
!  print *, "The Matrix is Singular"
!end if

!print *, " "
!print *, "--------------------------------------------------------------"

! deallocate(A1, B1, X1, E1, As, Bs, ATA, ATB, R, Q, eye)

End Program Driver_LinAl
