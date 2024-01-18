program solbybisection

    implicit none
    real, parameter :: xmin=-4.0, xmax=4.0
    integer, parameter :: iter = 10
    real :: sol,err,fcosx
    external fcosx
    
    call bisect(xmin,xmax,fcosx,sol,iter,err)
    write(*,*) sol,err

end program solbybisection
