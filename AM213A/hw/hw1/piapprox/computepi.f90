

program computepi

    implicit none
    
    integer :: numiter
    real :: threshold, difference

    threshold = 10.d0**(-4)
    call piapprox(threshold, numiter, difference)
    print "(A, e8.2, A, I5, A, e10.5)", 'For the threshold ', threshold, &
                                        ' the sum contained ', numiter, &
         ' iterations, and returned an absolute difference of ', difference   

    threshold = 10.d0**(-8)
    call piapprox(threshold, numiter, difference)
    print "(A, e8.2, A, I5, A, e10.5)", 'For the threshold ', threshold, &
                                        ' the sum contained ', numiter, &
        ' iterations, and returned an absolute difference of ', difference   


    threshold = 10.d0**(-12)
    call piapprox(threshold, numiter, difference)
    print "(A, e8.2, A, I5, A, e10.5)", 'For the threshold ', threshold, &
                                        ' the sum contained ', numiter, &
        ' iterations, and returned an absolute difference of ', difference   


    threshold = 10.d0**(-16)
    call piapprox(threshold, numiter, difference)
    print "(A, e8.2, A, I5, A, e10.5)", 'For the threshold ', threshold, &
                                        ' the sum contained ', numiter, &
        ' iterations, and returned an absolute difference of ', difference   


    contains

        subroutine piapprox(thres, n, diff)

            implicit none

            integer :: n
            real :: thres
            real :: diff
            real :: pisum, pitrue = acos(-1.d0)

            n = 0
            diff = thres + 1.
            pisum = 0.

            do while (thres .lt. diff)
                pisum = pisum + sumterm(n)
                n = n + 1

                diff = abs(pisum - pitrue)
            end do 

        end subroutine piapprox

        function sumterm(i) result(term)

            implicit none
            
            integer :: i
            real :: term
    
            term = (16.d0**(-i)) * ((4.d0/(8.d0*i + 1.d0)) - (2.d0/(8.d0*i + 4.d0)) - (1.d0/(8.d0*i + 5.d0)) - (1.d0/(8.d0*i + 6.d0)))

        end function sumterm

end program computepi



