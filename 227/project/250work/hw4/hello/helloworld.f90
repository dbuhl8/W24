


program helloworld

include 'mpif.h'

integer ierr, myid, numprocs

call MPI_INIT(ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)


print *,"Hello from processor #",myid," out of ",numprocs

call MPI_FINALIZE(ierr)
stop
end program helloworld
