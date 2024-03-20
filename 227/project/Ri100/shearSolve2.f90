

program shearSolve

    use LinAl

    implicit none

    integer, parameter :: M = 5, N = 5, nk = 40
    integer, parameter :: LDA = 5*(2*N+1)*(2*M+1)
    real, dimension(nk, 1) :: kz
    real :: ky = 1, kx = 0.5
    real :: kzmin = 0.01, kzmax = 1.5
    integer :: i, j, k, indu, indv, indw, indt, indp, indm, l
    real :: Pe, Re, Ri, f, dispterm
    real :: v_1 = 0.5, v0 = -0.915, v1 = -0.5, alpha=0.1
    character :: jobvl = "N", jobvr = "V"
!    real, dimension(LDA) :: ALPHAR, ALPHAI, BETA
!    real, dimension(LDA, LDA) :: A, B, D, V, VL
!    real, dimension(LDA, LDA) :: VR
!    real, dimension(nk, nk) :: lambda
!    real, allocatable :: work(:)
    complex, dimension(LDA) :: ALFA, BETA
    complex, allocatable :: A(:, :), B(:, :), V(:, :), VL(:, :) !dimension(LDA, LDA) :: A, B, V, VL
    real, allocatable :: D(:, :) !dimension(LDA, LDA) :: D !Note that we only take the real part of the eigenvalue when we find the max
    complex, allocatable :: VR(:, :) !dimension(LDA, LDA) :: VR
    real, dimension(nk, nk) :: lambda !Note that we only store real parts of the eigenvalues. 
    real, dimension(8*LDA) :: RWORK
    complex, allocatable :: work(:)
    integer :: lwork, info

    Pe = 1000.0
    Re = 1000.0
    Ri = 100.0

    f = 0.0

    allocate(A(LDA, LDA), B(LDA, LDA), D(LDA, LDA), V(LDA, LDA), VL(LDA, LDA), VR(LDA, LDA))

    do j = 1, nk
        
        kz(j, 1) = kzmin + (j-1)*(kzmax - kzmin)/(nk-1)
        A = 0.0
        B = 0.0
    
        do k = -M, M

            do i = -N, N

                indu = (i + N + 1) + (2*N+1)*(M+k)
                indv = (2*N + 1)*(2*M + 1) + indu
                indw = (2*N + 1)*(2*M + 1) + indv
                indt = (2*N + 1)*(2*M + 1) + indw
                indp = (2*N + 1)*(2*M + 1) + indt

                dispterm = ((f+i)**2)*kx**2 + ((f+k)**2)*ky**2 + kz(j, 1)**2
                ! check if i > -N
                if (i > -N) then
                    !w equation
                    A(indw, indw-1) = -(f+k)*ky*alpha*v0*(0.0, 1.0)
                    !t equation
                    A(indt, indt-1) = -(f+k)*ky*alpha*v0*(0.0, 1.0)
                    !v equation
                    A(indv, indu-1) = -kx*alpha*v0*(0.0, 1.0)
                    A(indv, indv-1) = -(f+k)*ky*alpha*v0*(0.0, 1.0)
                    !u equation
                    A(indu, indu-1) = -(f+k)*ky*alpha*v0*(0.0, 1.0)

                    ! check if k > -M (i > -N)
                    if (k > -M) then
                        A(indw, indw-2*N-2) = -(i+f)*kx*alpha*(0.0, 1.0) - (f+k)*ky*alpha*v1*(0.0, 1.0)
                        A(indt, indt-2*N-2) = -(i+f)*kx*alpha*(0.0, 1.0) - (f+k)*ky*alpha*v1*(0.0, 1.0)
                        A(indv, indv-2*N-2) = -(i+f)*kx*alpha*(0.0, 1.0) - (f+k)*ky*alpha*v1*(0.0, 1.0) + (ky*alpha/2.0)*(0.0, 1.0)
                        A(indv, indu-2*N-2) =  (kx*alpha/2.0)*(0.0, 1.0)
                        A(indu, indu-2*N-2) = -(kx*alpha)*(0.0, 1.0) - (f+k)*ky*alpha*v1*(0.0, 1.0) - (f+i)*kx*alpha*(0.0, 1.0)
                        A(indu, indv-2*N-2) = -ky*alpha*(0.0, 1.0)
                    end if
                    ! check if k < M (i > -N)
                    if (k < M) then
                        A(indw, indw+2*N) = -(i+f)*kx*alpha*(0.0, 1.0) - (f+k)*ky*alpha*v_1*(0.0, 1.0) 
                        A(indt, indt+2*N) = -(i+f)*kx*alpha*(0.0, 1.0) - (f+k)*ky*alpha*v_1*(0.0, 1.0)
                        A(indv, indv+2*N) = -(i+f)*kx*alpha*(0.0, 1.0) - (f+k)*ky*alpha*v_1*(0.0, 1.0) + (ky*alpha/2.0)*(0.0, 1.0)
                        A(indv, indu+2*N) = -(kx*alpha/2.0)*(0.0, 1.0)
                        A(indu, indu+2*N) = -(kx*alpha)*(0.0, 1.0) - (f+k)*ky*alpha*v_1*(0.0, 1.0) - (f+i)*kx*alpha*(0.0, 1.0)
                        A(indu, indv+2*N) =  ky*alpha*(0.0, 1.0)
                    end if

                end if

                ! this logic section is not needed yet for this problem
!               ! check if i < N
!               if(i < N) then
!                   ! check if k > -M (i < N)
!                   if (k > -M) then
!   
!                   end if
!                   ! check if k < M (i < N)
!                   if (k < M) then
!   
!                   end if
!   
!               end if 
 
                if (k > -M) then
                    A(indw, indw-(2*N+1)) = -((f+i)*kx/2.0)*(1.0, 0.0)
                    A(indt, indt-(2*N+1)) = -((f+i)*kx/2.0)*(1.0, 0.0)
                    A(indv, indv-(2*N+1)) = -((f+i)*kx/2.0)*(1.0, 0.0)
                    A(indu, indu-(2*N+1)) = -((f+i)*kx/2.0)*(1.0, 0.0)
                    A(indu, indv-(2*N+1)) = -(ky/2.0)*(1.0, 0.0)
    
                end if
                ! check if k < M
                if (k < M) then
                    A(indw, indw+(2*N+1)) = ((f+i)*kx/2.0)*(1.0, 0.0)
                    A(indt, indt+(2*N+1)) = ((f+i)*kx/2.0)*(1.0, 0.0)
                    A(indv, indv+(2*N+1)) = ((f+i)*kx/2.0)*(1.0, 0.0)
                    A(indu, indu+(2*N+1)) = ((f+i)*kx/2.0)*(1.0, 0.0)
                    A(indu, indv+(2*N+1)) = -(ky/2.0)*(1.0, 0.0)
    
                end if
               

                ! Dissipitive Terms
                A(indu,indu)= -dispterm/Re;
                A(indw,indw)= -dispterm/Re;
                A(indv,indv)= -dispterm/Re;
                A(indt,indt)= -dispterm/Pe;
                
                ! Stratification Terms here
                A(indt,indw)=-1.0;
                A(indw,indt)= Ri;

                ! Pressure Terms
                A(indu,indp)=-kx*(f+i);  ! avoid complex matrix by defining variable ip
                A(indv,indp)=-ky*(f+k);
                A(indw,indp)=-kz(j, 1);

                ! Eigenvalue Matrix B, 1 on diagonal
                B(indu,indu)=1.0;
                B(indv,indv)=1.0;
                B(indw,indw)=1.0;
                B(indt,indt)=1.0;
                B(indp,indp)=0.0;
                
                ! Continuity Equation
                A(indp,indu) = kx*(f+i);
                A(indp,indv) = ky*(f+k);
                A(indp,indw) = kz(j, 1);

            end do

            info = 0
            VL = 0.0
            VR = 0.0
            ALFA = 0.0
            BETA = 0.0
            lwork = -1
            allocate(work(1))
            work = 0.0
            call zggev(jobvl, jobvr, LDA, A, LDA, B, LDA, ALFA, BETA, VL, LDA, VR, LDA, WORK, LWORK, RWORK, info)
            lwork = work(1)
            deallocate(work)  
            allocate(work(lwork))
            call zggev(jobvl, jobvr, LDA, A, LDA, B, LDA, ALFA, BETA, VL, LDA, VR, LDA, WORK, LWORK, RWORK, info)
            deallocate(work)
            do l = 1, LDA
                if(abs(BETA(l)) > 10.d-6) then
                    D(l, l) = real(ALFA(l)/BETA(l))
                end if
            end do
            indm = 1 
            do l = 1, LDA
                if(D(l, l) .le. 0) then
                    D(l, l) = 0.;
                end if
                if(D(l, l) .ge. 100.0) then
                    D(l, l) = 0.;
                end if
                if(D(l, l) > D(indm, indm)) then
                    indm = l
                end if
            end do
            lambda(j, j) = D(indm, indm)
        end do
    end do

    

    call printmat(lambda, nk, nk)
!    call printmat(kx, nk, 1)
    call printmat(kz, nk, 1)

    open(10, file="kz.dat")
        do i = 1, nk
           write (10, *) kz(i, 1) 
        end do
    close(10)
!   open(11, file="kx.dat")
!       do i = 1, nk
!          write (11, *) kx(i, 1) 
!       end do
!   close(11)
    open(12, file="lambda.dat")
        do i = 1, nk
           write (12, *) lambda(i, :)
        end do
    close(12)
    open(13, file="nk.dat")
        write (13, *) nk
    close(13)

    deallocate(A, B, V, VL, D, VR)

  
end program shearSolve




