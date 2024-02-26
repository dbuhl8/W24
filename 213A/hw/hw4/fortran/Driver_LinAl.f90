Program Driver_LinAl

  use LinAl

  implicit none
  
  character(len=100) :: myFileName
  integer :: i,j
  real, parameter :: pi = 4.*atan(1.0), tol = 10.0**(-14)
  integer, parameter :: ma = 20, mb = 20, na = 20, nb = 1
  logical :: isSingular
  real, allocatable :: A1(:, :), B1(:, :), X1(:, :), As(:, :), Bs(:, :), E1(:, :)
  real, allocatable :: R(:, :), Q(:, :)
  real, dimension(ma) :: rvec
  real, dimension(20, 2) :: data, data2
  real, dimension(1, 2) :: datapivot
  real :: norm

  isSingular = .false.

  
  myFileName = 'atkinson.dat'

!  open(10,file=myFileName)
!  read(10,*) msize,nsize
!  close(10)

!  ma = msize
  allocate(A1(ma, na), B1(mb, nb), X1(mb, nb), E1(mb, nb), As(ma, ma), Bs(mb, nb))
  allocate(R(na, na), Q(ma, na))
  A1 = 0.0
  As = 0.0
  B1 = 0.0
  Bs = 0.0
  X1 = 0.0
  E1 = 0.0
  R = 0.0
  Q = 0.0
  rvec = 0.0
  data = 0.0
  datapivot = 0.0

  open(10, file=myFileName)
    do i = 1, ma
        read(10, *)  data(i, :)
    end do
  close(10)
  
  !pivoting data to improve cholesky
  !datapivot(1, :) = data(1, :)
  do i = 1, 20
    data2(i, :) = data(21 - i, :)
  end do
 ! data(20, :) = datapivot(1, :)


  !constructing Vandermonde matrix A1
  do j = 1, na
    As(:, j) = data2(:, 1)**(20 - j)
    Bs(j, 1) = data2(j, 2)
  end do

  A1 = As
  B1 = Bs


  print *, " "
  print *, "--------------------------------------------------------------"
  print *, " "
  print *, "Question 1: Cholesky Factorization"
  print *, " "

  print *, "Vandermonde Matrix for from Atkinson.dat"
  call printmat(As, ma, na)
  print *, " "

  call cholesky(A1, ma, isSingular, tol)

  call printmat(A1, ma, ma)
  print *, " "
    
  if (.not. isSingular) then  
      call LLsolve(A1, B1, X1, ma, nb)
    
      print *, "Matrix X1"
      call printmat(X1, mb, nb)
    
      E1 = matmul(As, X1) - Bs
    
      print *, "Matrix E1"
      call printmat(E1, mb, nb)

      call frobnorm(E1, norm)
      print *, "Frobenious norm of the error is: ", norm
  else 
    print *, "The matrix is singular"
  end if

  print *, " "
  print *, "--------------------------------------------------------------"
  print *, " "
  print *, "Question 2: QR Solution of the Least-Squares Problem"
  print *, " "

  A1 = As
  B1 = Bs
  X1 = 0.0
  E1 = 0.0

  print *, "Matrix A before QR"
  call printmat(A1, ma, na)  
  print *, "does this even get printed?"

  call householderQR(A1, rvec, ma, na, isSingular, tol)
  if (.not. isSingular) then
    print *, "Matrix A1, after QR"
    call printmat(A1, ma, na)

    call formR(A1, R, rvec, na)
    call formQstar(A1, Q, ma, na)

    print *, "Matrix R, after QR"
    call printmat(R, na, na)
    print *, "Matrix Q, after QR"
    call printmat(transpose(Q), ma, na)

    B1 = matmul(Q, B1)

    call backsub(R, B1, X1, ma, mb) 

    E1 = matmul(As, X1) - Bs
    print *, "Matrix X"
    call printmat(X1, mb, nb)
    print *, "Matrix E"
    call printmat(E1, mb, nb)
    
    call frobnorm(E1, norm)
      print *, "Frobenious norm of the error is: ", norm
 else 
   print *, "The Matrix is Singular"
 end if

 print *, " "
 print *, "--------------------------------------------------------------"
! print *, " "
! print *, "Question 5: Application Problem"
! print *, " "
! 
! !Aa = reshape((/ -3., 2., 1., 1., 2., 1., pi, exp(1.0), 1. /), shape(Aa))
! Aa = reshape((/ -3., 1., pi, 2., 2., exp(1.), 1., 1., 1. /), shape(Aa))
! Ba = reshape((/ 5., 3., -sqrt(2.0) /), shape(Ba))
! call printmat(Ba, 3, 1)
! Xa = 0.0
! print *, "In order to find the equation for the plane, we simply need to solve the Linear System, AX = B"
! print *, "Here is matrix A"
! call printmat(Aa, 3, 3)

! call LU(Aa, 3, isSingular, Pa)
! if (.not. isSingular) then
!   call LUsolve(Aa, 3, Ba, Xa, 1, Pa)
!   print *, "The solution is:"
!   call printmat(Xa, 3, 1)

!   Aa = reshape((/ -3., 1., pi, 2., 2., exp(1.), 1., 1., 1. /), shape(Aa))
! Ba = reshape((/ 5., 3., -sqrt(2.0) /), shape(Ba))

!   Xa = Ba - matmul(Aa, Xa)
!   print *, "The error is:"
!   call printmat(Xa, 3, 1)

!   call twonorm(Xa(:, 1), 3, norm)
!   print *, norm
!   
!   print *, "The permutation vector is"
!   print *, Pa(:)
! else
!   print *, "The matrix A is singular!"
! end if
!   

  deallocate(A1, As, B1, Bs, X1, E1)

End Program Driver_LinAl
