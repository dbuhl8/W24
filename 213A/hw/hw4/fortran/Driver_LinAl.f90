Program Driver_LinAl

  use LinAl

  implicit none
  
  character(len=100) :: myFileName
  integer :: i,j
  real, parameter :: pi = 4.*atan(1.0), tol = 10.0**(-14)
  integer, parameter :: ma = 20, mb = 20, na = 10, nb = 1, mata = 10, matb = 10, nsteps = 1000
  logical :: isSingular
  real, allocatable :: A1(:, :), B1(:, :), X1(:, :), As(:, :), Bs(:, :), E1(:, :), ATA(:, :), ATB(:, :)
  real, allocatable :: R(:, :), Q(:, :)
  real, allocatable :: eye(:, :)
  real, dimension(ma) :: rvec
  real, dimension(20, 2) :: data, data2
  real, dimension(1, 2) :: datapivot
  real, dimension(nsteps, 2) :: regression
  real :: norm, xstart, xstop, dx

  isSingular = .false.

  
  myFileName = 'atkinson.dat'

!  open(10,file=myFileName)
!  read(10,*) msize,nsize
!  close(10)

!  ma = msize
  allocate(A1(ma, na), B1(mb, nb), X1(na, nb), E1(mb, nb), As(ma, na), Bs(mb, nb), ATA(na, na), ATB(na, nb))
  allocate(R(na, na), Q(ma, na), eye(na, na))
  !note that ma = mb, na = mata = matb
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

  !constructing Vandermonde matrix A1
  do j = 1, na
    As(:, j) = data(:, 1)**(j-1)
    !Bs(j, 1) = data(j, 2)
  end do
  Bs(:, 1) = data(:, 2)

  A1 = As
  ATA = matmul(transpose(As), As)
  B1 = Bs
  ATB = matmul(transpose(As), Bs)

  print *, "Vandermonde Matrix from Atkinson.dat"
  call printmat(As, ma, na)
  print *, " "
    
  print *, "Y val Matrix from Atkinson.dat"
  call printmat(Bs, mb, nb)
  print *, " "



  print *, " "
  print *, "--------------------------------------------------------------"
  print *, " "
  print *, "Question 1: Cholesky Factorization"
  print *, " "

!  print *, "Vandermonde Matrix for from Atkinson.dat"
!  call printmat(As, ma, na)
!  print *, " "
    
  print *, "Vandermonde Matrix (Normalized) from Atkinson.dat"
  call printmat(ATA, mata, na)
  print *, " "

! print *, " ATB before cholesky factorization"
! call printmat(ATB, matb, nb)
! print *, " "

  call cholesky(ATA, mata, isSingular, tol)

! print *, "LL* after cholesky factorization"
! call printmat(matmul(ATA, transpose(ATA)), mata, na)
! print *, " "

! print *, "ATA after cholesky factorization"
! call printmat(ATA, mata, na)
! print *, " "

  if (.not. isSingular) then  
      call LLsolve(ATA, ATB, X1, mata, matb)
    
      print *, "Matrix X1"
      call printmat(X1, mb, nb)
    
      E1 = matmul(As, X1) - Bs
    
!     print *, "Matrix E1"
!     call printmat(E1, mb, nb)

      call frobnorm(E1, norm)
      print *, "Frobenious norm of the error is: ", norm

    !generating regression
    xstart = data(1, 1)
    xstop = data(ma, 1)
    dx = (xstop - xstart)/nsteps

    regression = 0.0
    
    do i = 1, nsteps
        regression(i, 1) = xstart + dx*(i-1)
    end do
    do i = 1, na
        regression(:, 2) = regression(:, 2) + (regression(:, 1)**(i-1))*X1(i, 1)
    end do
    
    !writing data and regression to a file for gnuplot
    open(15, file="plot.dat")
    do i = 1, ma
        write(15, "(2F10.3)") data(i, :)
    end do
    close(15)
    open(15, file="regression.dat")
    do i = 1, nsteps
        write(15, "(2F10.3)") regression(i, :)
    end do
    close(15)
   
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

!  print *, "Matrix A before QR"
!  call printmat(A1, ma, na)  

  call householderQR(A1, rvec, ma, na, isSingular, tol)
  if (.not. isSingular) then
!    print *, "Matrix A1, after QR"
!    call printmat(A1, ma, na)

    call formR(A1, R, rvec, na)
    call formQstar(A1, Q, ma, na)

!    print *, "Matrix R, after QR"
!    call printmat(R, na, na)
!    print *, "Matrix Q, after QR"
!    call printmat(transpose(Q), ma, na)
!    print *, "Matrix QTQ, after QR"
!    call printmat(matmul(transpose(Q), Q), na, na)
    call ident(eye, na)
    print *, "Matrix I - QTQ"
    call printmat(eye - matmul(transpose(Q), Q), na, na)
    call frobnorm(eye - matmul(transpose(Q), Q), norm)
    print *, "Frob Norm of I - QTQ: ", norm
    print *, "Matrix A - QR"
    call printmat(As - matmul(Q, R), ma, na)
    call frobnorm(As - matmul(Q, R), norm)
    print *, "Frob Norm of A - QR: ", norm





    B1 = matmul(transpose(Q), B1)

    call backsub(R, B1, X1, na, mb) 

    E1 = matmul(As, X1) - Bs
    print *, "Matrix X"
    call printmat(X1, mb, nb)
!    print *, "Matrix E"
!    call printmat(E1, mb, nb)
    
    call frobnorm(E1, norm)
      print *, "Frobenious norm of the error is: ", norm
 else 
   print *, "The Matrix is Singular"
 end if

 print *, " "
 print *, "--------------------------------------------------------------"

  deallocate(A1, B1, X1, E1, As, Bs, ATA, ATB, R, Q, eye)

End Program Driver_LinAl
