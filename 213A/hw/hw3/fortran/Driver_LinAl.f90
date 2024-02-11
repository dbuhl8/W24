Program Driver_LinAl

  use LinAl, only: mat, msize, nsize, readMat, str, trace, printmat, twonorm

  implicit none
  
  character(len=100) :: myFileName
  integer :: i,j
  real :: traceA=0., norm=0.

  
  myFileName = 'Amat.dat'

  open(10,file=myFileName)
  read(10,*) msize,nsize
  close(10)

  allocate(mat(msize,nsize))
  ! Always initialize with zeros
  mat = 0.0

  
  
  call readMat(myFileName)

  call printmat(mat, msize, nsize)

  call trace(mat, msize, traceA)
  print *, "Trace of the matrix is ", traceA

  do i = 1, nsize
    call twonorm(mat(:, i), msize, norm)
    print *, "Norm of "//trim(str(i))//"th column is ", norm
  end do

  do i = 1, msize
     write(*,*) (mat(i,j) , j = 1, nsize )
  end do
  

  deallocate(mat)





End Program Driver_LinAl
