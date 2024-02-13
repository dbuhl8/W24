Program Driver_LinAl

  use LinAl

  implicit none
  
  character(len=100) :: myFileName
  integer :: i,j
  real :: traceA=0., norm=0.
  real, allocatable :: A1(:, :), A2(:,:), As(:,:), B1(:,:), B2(:,:), Bs(:,:), X1(:,:), X2(:,:), E1(:,:), E2(:,:), E3(:, :)
  real, allocatable :: L(:, :), U(:, :)
  integer, allocatable :: P1(:)
  integer :: ma, mb, nb
  logical :: isSingular

  
  myFileName = 'Amat.dat'

  open(10,file=myFileName)
  read(10,*) msize,nsize
  close(10)

  allocate(A1(msize,nsize), A2(msize,nsize), As(msize, nsize), mat(msize, nsize), P1(msize), L(msize, nsize), U(msize, nsize),&
    & E3(msize, nsize))
  ma = msize
  ! Always initialize with zeros
  mat = 0.0

  call readMat(myFileName)
  A1 = mat
  A2 = A1
  As = A1
  L = 0.
  U = 0.
  deallocate(mat)

  myFileName = 'Bmat.dat'

  open(10,file=myFileName)
  read(10,*) msize,nsize
  close(10)

  allocate(B1(msize,nsize), B2(msize, nsize), Bs(msize, nsize), mat(msize, nsize), X1(msize, nsize), X2(msize, nsize))
  allocate(E1(msize, nsize), E2(msize, nsize))
  mb = msize
  nb = nsize
  mat = 0.0

  call readMat(myFileName)
  B1 = mat
  B2 = B1
  Bs = B1
  E1 = 0.0
  E2 = 0.0
  X1 = 0.0
  X2 = 0.0
  E3 = 0.0

  print *, "Question 2: Basic Fortran Routines"
  call printmat(A1, ma, nsize)

  call trace(A1, ma, traceA)
  print *, "Trace of the matrix is ", traceA

  do i = 1, ma
    call twonorm(A1(:, i), ma, norm)
    print *, "Norm of "//trim(str(i))//"th column is ", norm
  end do

  do i = 1, msize
     write(*,*) (A1(i,j) , j = 1, ma)
  end do
  
  print *, "Question 3: Gaussian Elimination"
  print *, "Matrix A1, before GE"
  call printmat(A1, ma, ma)
  print *, "Matrix B1, before GE"
  call printmat(B1, mb, nb)
  call GE(A1, B1, ma, nb, isSingular)
  print *, "Matrix A1, after GE"
  call printmat(A1, ma, ma)
  print *, "Matrix B1, after GE"
  call printmat(B1, mb, nb)
  call backsub(A1, B1, X1, ma, nb)
  print *, "Matrix X1"
  call printmat(X1, mb, nb)
  E1 = matmul(As, X1) - Bs
  print *, "Matrix E1"
  call printmat(E1, mb, nb)
  do i = 1, msize
     write(*,*) (E1(i,j) , j = 1, nb)
  end do
  do i = 1, nb
    call twonorm(E1(:, i), mb, norm)
    print *, "Norm of "//trim(str(i))//"th column is ", norm
  end do

  print *, "Question 4: LU Method for Solving Systems of Equations"
  print *, "Matrix A2 before LU"
  call printmat(A2, ma, ma)
  call LU(A2, ma, isSingular, P1)
  print *, "Matrix A2, after LU"
  call printmat(A2, ma, ma)
  do j = 1, ma
    do i = 1, ma
        if (i .eq. j) then
            L(i, j) = 1.
            U(i, j) = A2(i, j)
        else if (i .gt. j) then
            L(i, j) = A2(i, j)
        else 
            U(i, j) = A2(i, j)
        end if
    end do
  end do
  E3 = matmul(L, U)
  print *, "DEBUGGING: This is the error in the LU decomp"
  call printmat(E3, ma, ma)
  print *, "L after LU decomposition"
  call printmat(L, ma, ma)
  print *, "U after LU Decomposition"
  call printmat(U, ma, ma)
  print *, "A1 from GE (should be same as U)"
  call printmat(A1, ma, ma)
  call LUsolve(A2, ma, B2, X2, nb, P1)
  E2 = matmul(As, X2) - Bs
  print *, "Matrix X2"
  call printmat(X2, mb, nb)
  print *, "Matrix E2"
  call printmat(E2, mb, nb)
  do i = 1, nb
    call twonorm(E2(:, i), mb, norm)
    print *, "Norm of "//trim(str(i))//"th column is ", norm
  end do

  deallocate(mat, A1, A2, As, B1, B2, Bs, X1, X2, E1, E2, P1, L, U)

End Program Driver_LinAl
