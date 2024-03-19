Program Driver_LinAl

  use LinAl

  implicit none
  
  integer :: i,j,k
  real, parameter :: tol = 10.0**(-14)
  integer, parameter :: ma = 1920, na = 1279, ml = 3, nb = 1
  integer, parameter :: LDA=ma, LDU=ma, LDV=na
  real, dimension(ma, na) :: A, As, At
  real, dimension(ml, ml) :: L, eye, Ls
  real, dimension(ml, nb) :: B, X
  character :: jobu='A', jobvt='A'
  integer :: lwork, info
  real, dimension(min(ma, na)) :: S
  real, dimension(ma, na) :: Sigma, Sigma10, Sigma20, Sigma40, Sigma80, Sigma160
  real, dimension(ma, na) :: Sigma320, Sigma640
  real, dimension(LDU, ma) :: U
  real, dimension(LDV, na) :: VT
  real, allocatable :: work(:)
  real :: norm
  character :: jors

  open(10, file="dog_bw_data.dat")
  do i = 1, na
    read (10, *) A(:, i)
  end do
  close(10) 

  print *, " "
  print *, "--------------------------------------------------------------"
  print *, " "
  print *, "Question 1: Reading in dog_bw_data.dat into A"
  print *, " "
  print *, " There was a print statement for A here before. It was ridiculously large so it is omitted"
  print *, " "
  print *, "--------------------------------------------------------------"
  print *, " "
  print *, "Question 2: Performing SVD on A"
  print *, " "

  As = A

  ! Finding the optimal sizae for work
  info = 0
  lwork = -1
  allocate(work(1))
  call dgesvd(jobu, jobvt, ma, na, A, LDA, S, U, LDU, VT, LDV, work, lwork, info)
  lwork = work(1)
  deallocate(work)
  allocate(work(lwork))

  call dgesvd(jobu, jobvt, ma, na, A, LDA, S, U, LDU, VT, LDV, work, lwork, info)

  !constructing sigma
  Sigma = 0.0
  do i = 1, min(ma, na)
    Sigma(i, i) = S(i)
  end do


  print *, " "
  print *, "--------------------------------------------------------------"
  print *, " "
  print *, "Question 3-4: Reduced Rank Reconstructions"
  print *, " "

  call reducedrank(Sigma, Sigma10, 10, ma, na)
  call reducedrank(Sigma, Sigma20, 20, ma, na)
  call reducedrank(Sigma, Sigma40, 40, ma, na)
  call reducedrank(Sigma, Sigma80, 80, ma, na)
  call reducedrank(Sigma, Sigma160, 160, ma, na)
  call reducedrank(Sigma, Sigma320, 320, ma, na)
  call reducedrank(Sigma, Sigma640, 640, ma, na)

  Sigma10 = matmul(U, matmul(Sigma10, VT))
  Sigma20 = matmul(U, matmul(Sigma20, VT))
  Sigma40 = matmul(U, matmul(Sigma40, VT))
  Sigma80 = matmul(U, matmul(Sigma80, VT))
  Sigma160 = matmul(U, matmul(Sigma160, VT))
  Sigma320 = matmul(U, matmul(Sigma320, VT))
  Sigma640 = matmul(U, matmul(Sigma640, VT))
  Sigma = matmul(U, matmul(Sigma, VT))

  open(11, file="Image_appn_100010.dat")
    do i = 1, na
        write(11, *) Sigma10(:, i)
    end do
  close(11)
  open(12, file="Image_appn_100020.dat")
    do i = 1, na
        write(12, *) Sigma20(:, i)
    end do
  close(12)
  open(13, file="Image_appn_100040.dat")
    do i = 1, na
        write(13, *) Sigma40(:, i)
    end do
  close(13)
  open(14, file="Image_appn_100080.dat")
    do i = 1, na
        write(14, *) Sigma80(:, i)
    end do
  close(14)
  open(15, file="Image_appn_100160.dat")
    do i = 1, na
        write(15, *) Sigma160(:, i)
    end do
  close(15)
  open(16, file="Image_appn_100320.dat")
    do i = 1, na
        write(16, *) Sigma320(:, i)
    end do
  close(16)
  open(17, file="Image_appn_100640.dat")
    do i = 1, na
        write(17, *) Sigma640(:, i)
    end do
  close(17)
  open(18, file="Image_appn_101279.dat")
    do i = 1, na
        write(18, *) Sigma(:, i)
    end do
  close(18)

  print *, " "
  print *, "--------------------------------------------------------------"
  print *, " "
  print *, "Question 6: Comparing Reduced Rank Approximation Errors"
  print *, " "

  print *, "Here are the first 10 singular values: "
  print *, S(1:10)
  print *, " "

  print *, "Here are the 'rest singular values' for each k (in order of ascending k): "
  print "(8F16.8)", S(10), S(20), S(40), S(80), S(160), S(320), S(640), S(1279)

  At = As - Sigma10
  call frobnorm(At, norm)
  norm = norm/(ma*na)
  print *, "Here is the Frob Norm of A - U(Sigma10)VT. ERROR: ", norm
  print *, " "
  At = As - Sigma20
  call frobnorm(At, norm)
  norm = norm/(ma*na)
  print *, "Here is the Frob Norm of A - U(Sigma20)VT. ERROR: ", norm
  print *, " "
  At = As - Sigma40
  call frobnorm(At, norm)
  norm = norm/(ma*na)
  print *, "Here is the Frob Norm of A - U(Sigma40)VT. ERROR: ", norm
  print *, " "
  At = As - Sigma80
  call frobnorm(At, norm)
  norm = norm/(ma*na)
  print *, "Here is the Frob Norm of A - U(Sigma80)VT. ERROR: ", norm
  print *, " "
  At = As - Sigma160
  call frobnorm(At, norm)
  norm = norm/(ma*na)
  print *, "Here is the Frob Norm of A - U(Sigma160)VT. ERROR: ", norm
  print *, " "
  At = As - Sigma320
  call frobnorm(At, norm)
  norm = norm/(ma*na)
  print *, "Here is the Frob Norm of A - U(Sigma320)VT. ERROR: ", norm
  print *, " "
  At = As - Sigma640
  call frobnorm(At, norm)
  norm = norm/(ma*na)
  print *, "Here is the Frob Norm of A - U(Sigma640)VT. ERROR: ", norm
  print *, " "
  At = As - Sigma
  call frobnorm(At, norm)
  norm = norm/(ma*na)
  print *, "Here is the Frob Norm of A - U(Sigma1279)VT. ERROR: ", norm
  print *, " "

  print *, " "
  print *, "--------------------------------------------------------------"

End Program Driver_LinAl
