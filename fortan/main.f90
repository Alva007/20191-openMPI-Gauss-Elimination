program main
    use mpi
    use mpi_utils, only: mpi_matinv, mpi_rank, mpi_nproc, mpi_ierr

    implicit none

    integer, parameter :: n = 1000 ! size of the matrix
    integer :: i
    double precision :: t0, tf
    double precision, dimension(n * n) :: array
    double precision, dimension(n, n) :: A, A_inv

    ! Initialize MPI (needed in all MPI codes)
    call MPI_INIT(mpi_ierr)

    ! Get the rank and the number of process (needed for mpi_matinv)
    call MPI_COMM_RANK(MPI_COMM_WORLD,mpi_rank,mpi_ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,mpi_nproc,mpi_ierr)

    ! Fill the matrix to inverse with random numbers
    do i = 1, n * n
        array(i) = rand(0) * 1D0
    enddo
    A = transpose(reshape(array, shape(A)))
  
    ! Initialize the inversed matrix
    A_inv(:, :) = 0D0
    
    call cpu_time(t0) ! store start inversion time

    ! Do the inversion
    call mpi_matinv(A(:, :),A_inv(:, :), n)

    call cpu_time(tf) ! store end inversion time

    ! Check if the result is correct
    if (mpi_rank == 0) print *, 'result =', sum(matmul(A(:,:), A_inv(:,:))) / n  ! ideally = 1
    if (mpi_rank == 0) print *, 'elapsed time (s) = ', tf - t0

    ! Wait for all the process (avoid program ending before the display of the results)
    call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)

    ! Put this at the end of all MPI codes
    call MPI_FINALIZE(mpi_ierr)
end program main

