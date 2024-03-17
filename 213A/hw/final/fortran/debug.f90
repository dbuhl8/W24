


program debug

    use LinAl

    implicit none

    integer, parameter :: ma = 4
    real, parameter :: tol=10.d-14

    real, dimension(ma, ma) :: A, R, Q, eye, B
    real, dimension(ma) :: rvec

    logical :: isSingular

    R = 0.0
    Q = 0.0
    eye = 0.0
    
    A(1, :) = (/2.0, 0.0, 0.0, 0.0/)
    A(2, :) = (/0.0, 1.0, 0.0, 0.0/)
    A(3, :) = (/0.0, 0.0, 3.0, 0.0/)
    A(4, :) = (/0.0, 0.0, 0.0, 0.0/)

    call ident(eye, ma)

    call householderQR(A, rvec, ma, ma, isSingular, tol)
    call formR(A, R, rvec, ma)
    call formQstar(A, Q, ma, ma)

    B = matmul(R, Q)
    call printmat(B, ma, ma)
    
end program debug
