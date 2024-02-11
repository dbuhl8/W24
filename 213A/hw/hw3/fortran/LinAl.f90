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

  end subroutine twonorm

  subroutine printmat(A, m, n)

    real, intent(in) :: A(:, :)
    integer, intent(in) :: m, n
    integer :: i, j

    print "(A, I3, A, I3)", "A is an ", m," by ", n," matrix."
    do i = 1, m
        print "("//trim(str(n))//"F6.2)", A(i, :)
    end do

  end subroutine printmat

character(len=20) function str(k)
!   "Convert an integer to string."
    integer, intent(in) :: k
    write (str, *) k
    str = adjustl(str)
end function str

  subroutine GE(A, B, ma, mb, nb, bool)

  end subroutine GE

  subroutine UBsolver(U, B, X, mu, mb, n)

  end subroutine UBsolver

  subroutine LUdecomp(A, m, bool, s)

  end subroutine LUdecomp

  subroutine LUsolver(LU, m, B, mb, s)

  end subroutine LUsolver

end module LinAl
