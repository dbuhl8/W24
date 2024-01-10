

program betteraddmats

    implicit none

    integer, parameter :: dimmat = 10000
    integer :: i, j
    real, dimension(dimmat, dimmat) :: a, b, c
    a = 0.
    b = 0.
    c = 0.
    ! This creates the matrices.
    a(1,2) = 2.0
    do i=2,dimmat-1
        a(i,i+1) = 2.0
        b(i,i-1) = 1.0
    enddo
    b(dimmat,dimmat-1) = 1.0

    c = a + b

!    do i = 1, dimmat
!        write(*,*) (c(i,j), j=1,dimmat)
!    enddo

end program betteraddmats
