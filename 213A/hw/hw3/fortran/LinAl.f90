module LinAl 

implicit none 
integer, save :: msize, nsize 
real, dimension(:,:), allocatable, save :: mat
contains

  !********************************************************

  subroutine readMat(filename)

    implicit none
    character(len=*) :: filename

    integer :: i,j

    ! Reads a file containing the matrix A 
    ! Sample file:
    !
    ! 4 4 
    ! 2.0 1.0 1.0 0.0
    ! 4.0 3.0 3.0 1.0
    ! 8.0 7.0 9.0 5.0
    ! 6.0 7.0 9.0 8.0
    !
    ! Note that the first 2 numbers in the first line are the matrix dimensions, i.e., 4x4,
    ! then the next msize lines are the matrix entries. This matrix is found in Eq. 2.18 of the lecture note.
    ! Note that entries must be separated by a tab.


    open(10,file=filename)

    ! Read the matrix dimensions
    read(10,*) i,j

    ! Read matrix
    do i=1,msize
       read(10,*) ( mat(i,j), j=1,nsize )
    enddo

    close(10)
    
  end subroutine readMat

  subroutine trace(A, m, tr)
    
    implicit none
    
    real, intent(in) :: A(:, :)
    integer, intent(in) :: m
    real, intent(out) :: tr 
    integer :: i

    do i = 1, m
        tr = tr + A(i, i)
    end do

  end subroutine trace

  subroutine twonorm(v, m, nrm)
    
    implicit none

    real, intent(in) :: v(:)
    integer, intent (in) :: m
    real, intent(out) :: nrm 
    integer :: i

    do i = 1, m
        nrm = nrm + v(i)**2
    end do

    nrm = sqrt(nrm)
    !nrm = sqrt(sum(v**2))

  end subroutine twonorm

  subroutine printmat(A, m, n)

    real, intent(in) :: A(:, :)
    integer, intent(in) :: m, n
    integer :: i, j

    print "(A, I3, A, I3)", "This is an ", m," by ", n," matrix."
    do i = 1, m
        print "("//trim(str(n))//"F8.2)", A(i, :)
    end do

  end subroutine printmat

character(len=20) function str(k)
!   "Convert an integer to string."
    integer, intent(in) :: k
    write (str, *) k
    str = adjustl(str)
end function str

  subroutine GE(A, B, ma, bool)
        
    implicit none
    
    integer, intent(in) :: ma
    logical :: bool
    real :: A(:, :), B(:, :)
    real, allocatable :: rowapivot(:)
    !real, dimension(1, mb) :: rowbpivot
    real, allocatable :: rowlpivot(:)
    real, dimension(ma) :: rowppivot
    real, dimension(ma, ma) :: L
    integer, dimension(ma, ma) :: P
    integer :: i, j, k
    bool = .false. 
    !initializing our matrices. 
    L = 0.0
    P = 0
    do i = 1, ma
        L(i, i) = 1.0 
        P(i, i) = 1
    end do

    !begin GE 
    do i = 1, ma-1
        !pivoting step
        k = i
        do j = i+1, ma
            if (abs(A(j, i)) .gt. abs(A(k, i))) then
                k = j
            end if
        end do
        allocate(rowlpivot(1:(k-1)), rowapivot(k:ma))
        rowapivot = A(i, k:ma)
        !rowbpivot = B(i, k:m)
        rowlpivot = L(i, 1:(k-1))
        rowppivot = P(i, :)
        A(i, k:ma) = A(k, k:ma)
        !B(i, k:m) = B(k, k:m)
        L(i, 1:(k-1)) = L(k, 1:(k-1))
        P(i, :) = P(k, :)
        A(k, k:ma) = rowapivot
        !B(k, k:m) = rowbpivot
        L(k, 1:(k-1)) = rowlpivot
        P(k, :) = rowppivot
        deallocate(rowlpivot, rowapivot)
        if (A(i, i) .eq. 0) then
            bool = .true.
            return
        end if
        !L U step
        L(i+1:ma, i) = A(i+1:ma,i)/A(i, i)
        do  j = i+1, ma
            A(j, i:ma) = A(j, i:ma) - L(j, i)*A(i, i:ma)
        end do
    end do 
    B = matmul(L, matmul(P, B))
    
  end subroutine GE

  subroutine backsub(U, B, X, mu, mb)
    implicit none
    
    integer, intent(in) :: mu, mb
    real, intent(in) :: B(:, :), U(:, :)
    real :: X(:, :)

    integer :: i, j, k

    X(mu, :) = B(mu, :)/U(mu,mu)
    do j = 1, mb
        do i = 1, mu-1
            k = mu - i
            X(k, j) = (B(k, j) - dot_product(U(k, k+1:mu), X(k+1:mu, j)))/U(k, k)
        end do
    end do

  end subroutine backsub

  subroutine LU(A, ma, bool, P)
    
    implicit none
    
    integer, intent(in) :: ma
    logical :: bool
    real :: A(:, :)
    real, allocatable :: rowapivot(:)
    real, allocatable :: rowlpivot(:)
    real, dimension(ma) :: rowppivot
    real, dimension(ma, ma) :: L
    integer :: P(:, :)
    integer :: i, j, k
    bool = .false. 
    !initializing our matrices. 
    L = 0.0
    P = 0
    do i = 1, ma
        L(i, i) = 1.0 
        P(i, i) = 1
    end do

    !begin GE 
    do i = 1, ma-1
        !pivoting step
        k = i
        do j = i+1, ma
            if (abs(A(j, i)) .gt. abs(A(k, i))) then
                k = j
            end if
        end do
        allocate(rowlpivot(1:(k-1)), rowapivot(k:ma))
        rowapivot = A(i, k:ma)
        rowlpivot = L(i, 1:(k-1))
        rowppivot = P(i, :)
        A(i, k:ma) = A(k, k:ma)
        L(i, 1:(k-1)) = L(k, 1:(k-1))
        P(i, :) = P(k, :)
        A(k, k:ma) = rowapivot
        L(k, 1:(k-1)) = rowlpivot
        P(k, :) = rowppivot
        deallocate(rowlpivot, rowapivot)
        if (A(i, i) .eq. 0) then
            bool = .true.
            return
        end if
        !L U step
        L(i+1:ma, i) = A(i+1:ma,i)/A(i, i)
        do  j = i+1, ma
            A(j, i:ma) = A(j, i:ma) - L(j, i)*A(i, i:ma)
        end do
        A = A + L !this gives A as L + U where the only overlap is along the diagonals
    end do 

  end subroutine LU

  subroutine LUsolve(LU, ma, B, X, mb, P)

    implicit none

    integer, intent(in) :: ma, mb, P(:, :)
    real, intent(in) :: LU(:, :), X(:,:)
    real :: B(:, :)
    real, dimension(ma, ma) :: L, U
    real, dimension(ma, mb) :: Y
    integer :: i , j, k

    L = 0.
    U = 0.

    do j = 1, ma
        do i = 1, ma
            if(i .eq. j) then   
                L(i, j) = 1.
                U(i, j) = LU(i, j) - 1
            else if(i .gt. j) then
                L(i, j) = LU(i, j)
            else
                U(i, j) = LU(i, j)
            end if
        end do
    end do
    B = matmul(P, B)
    !forward sub for Ly = b
    call forwardsub(L, B, Y, ma, mb)

    !Backward sub for Ux = y
    call backsub(U, Y, X, ma, mb)

  end subroutine LUsolve
    
  subroutine forwardsub(L, B, Y, ma, mb) 

    implicit none

    integer, intent(in) :: ma, mb
    real :: L(:, :), B(:, :), Y(:, :)
    integer :: i, j
    Y = 0. 
     
    Y(1, :) = B(1, :)
    do j = 1, mb
        do i = 2, ma
            ! THIS NEEDS TO BE A DOT PRODUCT INSTEAD OF COORESP SCALAR MULT
            Y(i, j) = (B(i, j) - dot_product(L(i, 1:i-1), Y(1:i-1, j)))/L(i, i)
        end do
    end do

  end subroutine forwardsub

end module LinAl
