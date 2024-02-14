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

    nrm = sqrt(sum(v**2))

  end subroutine twonorm

  subroutine printmat(A, m, n)

    real, intent(in) :: A(:, :)
    integer, intent(in) :: m, n
    integer :: i, j

    print "(A, I3, A, I3)", "This is ", m," by ", n," matrix."
    do i = 1, m
        print "("//trim(str(n))//"F10.3)", A(i, :)
    end do

  end subroutine printmat

character(len=20) function str(k)
!   "Convert an integer to string."
    integer, intent(in) :: k
    write (str, *) k
    str = adjustl(str)
end function str

  subroutine GE(A, B, ma, nb, bool)
        
    implicit none
    
    integer, intent(in) :: ma, nb
    logical :: bool
    real :: A(:, :), B(:, :)
    real, allocatable :: rowapivot(:), rowbpivot(:)
    real, dimension(ma, ma) :: L
    integer :: i, j, k
    bool = .false. 
    
    !initializing our matrices. 
    L = 0.0
    do i = 1, ma
        L(i, i) = 1.0 
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
        allocate(rowapivot(k:ma), rowbpivot(nb))
        rowapivot = A(i, i:ma)
        rowbpivot = B(i, :)
        A(i, i:ma) = A(k, i:ma)
        B(i, :) = B(k, :)
        A(k, i:ma) = rowapivot
        B(k, :) = rowbpivot
        deallocate(rowapivot, rowbpivot)
        if (A(i, i) .eq. 0) then
            bool = .true.
            return
        end if
        !L Step
        L(i+1:ma, i) = A(i+1:ma,i)/A(i, i)
    
        A(i+1:ma, i:ma) = A(i+1:ma, i:ma) - matmul(L(i+1:ma, i:i), A(i:i, i:ma))
        !do  j = i+1, ma
            ! this can be optimized to not require a loop 
            ! A(i+1:ma, i:ma) = A(i+1:ma, i:ma) - matmul(L(i+1:ma, i), A(i, i:ma))
        !    A(j, i:ma) = A(j, i:ma) - L(j, i)*A(i, i:ma)
        !end do
        do j = i+1, ma
            B(j, :) = B(j, :) - L(j, i)*B(i, :)
        end do
    end do 
    
  end subroutine GE

  subroutine backsub(U, B, X, mu, mb)
    implicit none
    
    integer, intent(in) :: mu, mb
    real, intent(in) :: B(:, :), U(:, :)
    real :: X(:, :)

    integer :: i, j, k

    X(mu, :) = B(mu, :)/U(mu,mu)
    do i = 1, mu-1
        k = mu - i
        X(k, :) = (B(k, :) - matmul(U(k, k+1:mu), X(k+1:mu, :)))/U(k, k)
    end do

  end subroutine backsub

  subroutine LU(A, ma, bool, P)
    
    implicit none
    
    integer, intent(in) :: ma
    logical :: bool
    real :: A(:, :)
    real, allocatable :: rowapivot(:)
    real, allocatable :: rowlpivot(:)
    real :: rowppivot
    real, dimension(ma, ma) :: L
    integer :: P(:)
    integer :: i, j, k
    bool = .false. 
    !initializing our matrices. 
    L = 0.0
    P = 0
    do i = 1, ma
        P(i) = i
    end do
    k = 0
    !begin GE 
    do i = 1, ma-1
        !pivoting step
        k = i
        do j = i+1, ma
            if (abs(A(j, i)) .gt. abs(A(k, i))) then !set pivot row to k
                k = j
            end if
        end do
        if (k .ne. i) then !pivot A and P
            allocate(rowapivot(k:ma))
            rowapivot = A(i, i:ma)
            rowppivot = P(i)
            A(i, i:ma) = A(k, i:ma)
            P(i) = P(k)
            A(k, i:ma) = rowapivot
            P(k) = rowppivot
            deallocate(rowapivot)
            if (i .ne. 1) then !pivot l
                allocate(rowlpivot(1:i-1))
                rowlpivot = L(i, 1:i-1)
                L(i, 1:i-1) = L(k, 1:i-1)
                L(k, 1:i-1) = rowlpivot
                deallocate(rowlpivot)
            end if
        end if
        if (A(i, i) .eq. 0) then !return singular
            bool = .true.
            return
        end if
        !L U step
        L(i+1:ma, i) = A(i+1:ma,i)/A(i, i)
        do  j = i+1, ma
            A(j, i:ma) = A(j, i:ma) - L(j, i)*A(i, i:ma)
        end do
    end do 
!    print *, "Printing L before adding L to A"
!    call printmat(L, ma, ma)
    A = A + L !this gives A as L + U where the only overlap is along the diagonals

  end subroutine LU

  subroutine LUsolve(LU, ma, B, X, mb, P)

    implicit none

    integer, intent(in) :: ma, mb, P(:)
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
                U(i, j) = LU(i, j)
            else if(i .gt. j) then
                L(i, j) = LU(i, j)
            else
                U(i, j) = LU(i, j)
            end if
        end do
    end do
    ! need to permute B
    do i = 1, ma
        Y(i, :) = B(P(i), :)
    end do
    B = Y
    Y = 0.
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
    do i = 2, ma
        Y(i, :) = B(i, :) - matmul(L(i, 1:i-1), Y(1:i-1, :))
    end do

  end subroutine forwardsub

end module LinAl
