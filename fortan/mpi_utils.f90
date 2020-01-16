module mpi_utils
    ! """
    ! Contains the parameters needed for MPI and some useful MPI functions.
    ! """
    use mpi

    implicit none

    integer, parameter :: &
        mpi_tag1 = 100  ! tag

    integer :: &
        mpi_nproc,  & ! number of processes in the run
        mpi_rank,   & ! rank of the process
        mpi_ierr      ! run code

    integer, dimension(MPI_STATUS_SIZE) :: &
        mpi_status1  ! status

    save

    contains
        subroutine mpi_fair_share(n0, n, displ)
            ! """
            ! Redistributes equitably the size and the indexes of an array to all tasks in a group.
            ! Useful for the paralellization of a loop over an array and for all *v mpi functions, like MPI_SCATTERV.
            !
            ! input:
            !   n0: original size
            !
            ! outputs:
            !   n: new size for a process can be used to build the sendcounts array for a *v mpi function
            !   displ: displacement relative to n0 for a process can be used to build the displs array for a *v mpi
            !          function
            !
            ! notes:
            !   mpi_rank and mpi_nproc must have been initialized
            ! """
            implicit none

            integer, intent(in) :: n0
            integer, intent(out) :: displ, n

            if (mpi_nproc > n0 .and. mpi_rank == 0) then
                print '("WARNING: number of processors (",I6,") > number of elements (",I6,")")',mpi_nproc,n0
                print '(1X,"> prepare for unforeseen consequences...")'
            endif

            if (mod(n0, mpi_nproc) == 0) then
                n = n0  /mpi_nproc
                displ = n * mpi_rank
            else
                if (mpi_rank < mod(n0, mpi_nproc)) then
                    n = floor(real(n0 / mpi_nproc)) + 1
                    displ = n * mpi_rank
                else
                    n = floor(real(n0 / mpi_nproc))
                    if (mpi_rank == mod(n0, mpi_nproc)) then
                        displ = (n + 1) * mpi_rank
                    else
                        displ = (n + 1) * mpi_rank + mod(n0, mpi_nproc) - mpi_rank
                    endif
                endif
            endif

            return

        end subroutine mpi_fair_share


        subroutine mpi_fair_sharev(recvcount, sendcount, sendcounts, displ, displs)
            ! """
            ! Redistributes equitably the size and the indexes of an array for all tasks in a group, then build the
            ! input arrays needed for any *v MPI function.
            ! Useful for the paralellization of a loop over an array and for all *v mpi functions, like MPI_SCATTERV.
            !
            ! input:
            !   recvcount: number of elements in receive buffer (integer)
            !
            ! outputs:
            !   sendcount: number of elements in send buffer (integer)
            !   sendcounts: integer array (of length group size) specifying the number of elements to send to each
            !               processor
            !   displ: displacement relative to sendbuf from which to take the outgoing data to the current process
            !   displs: integer array (of length group size). Entry i specifies the displacement relative to sendbuf
            !           from which to take the outgoing data to process i
            !
            ! notes:
            !   mpi_rank and mpi_nproc must have been initialized
            ! """
            implicit none

            integer, intent(in) :: recvcount
            integer, intent(out) :: sendcount, displ
            integer, intent(out), dimension(mpi_nproc) :: sendcounts, displs

            call mpi_fair_share(recvcount, sendcount, displ)
            call MPI_ALLGATHER(sendcount, 1, MPI_INT, sendcounts, 1, MPI_INT, MPI_COMM_WORLD, mpi_ierr)
            call MPI_ALLGATHER(displ, 1, MPI_INT, displs, 1, MPI_INT, MPI_COMM_WORLD, mpi_ierr)

            return

        end subroutine mpi_fair_sharev


        subroutine mpi_matinv(A0, A0_inv, n0)
            ! """
            ! Inverse a 2D matrix of dimension (n0, n0) using MPI.
            !
            ! inputs:
            !   A0(n0, n0): matrix to inverse
            !   n0: dimension of the matrix
            !
            ! output:
            !   A0_inv(n0, n0): inverse of matrix A0
            !
            ! notes:
            !   Method: Parallel gaussian ellimination
            !   Based on: CSE 633 Parallel Algorithms (Spring 2014) by Aravindhan Thanigachalam (athaniga@buffalo.edu)
            !   Dependent of subroutine mpi_fair_sharv, which redistributes equitably indexes to all tasks in a group.
            ! """
            implicit none

            integer, intent(in) :: n0
            double precision, dimension(n0, n0), intent(in) :: A0
            double precision, dimension(n0, n0), intent(inout) :: A0_inv

            integer :: i, i_para, j, j_min, j_max, jz, jz_para, ncol, r, root, size_A
            integer, dimension(mpi_nproc) :: ncols, sizes, displs_j, displs_j0  ! mpi_fair_sharev utils
            double precision :: scale
            double precision, dimension(n0) :: factor
            double precision, dimension(:), allocatable :: row
            double precision, dimension(n0, n0) :: Id0  ! identity matrix
            double precision, dimension(:, :), allocatable :: A, A_inv

            factor(:) = 0D0

            ! Build identity matrix
            Id0(:, :) = 0D0
            do i = 1, n0
                Id0(i, i) = 1D0
            enddo

            ! Column wise distribution
            call mpi_fair_sharev(n0, ncol, ncols, j_min, displs_j0)
            j_min = j_min + 1 ! fortran index format
            j_max = j_min + ncol - 1

            allocate(A(n0, ncol), A_inv(n0, ncol), row(ncol))
            A(:, :) = A0(:, j_min:j_max)
            A_inv(:, :) = Id0(:, j_min:j_max)
            row(:) = 0D0

            ! Gaussian elimination phase
            i_para = 1
            root = 0

            do i=1,n0
                ! Step 1: swap row i with the nearest subsequent row r such that after swapping A(i,i) /= 0
                ! Step 1a: find the rank (root) of the process where the column i is
                if (i > displs_j0(root + 1) + ncols(root + 1)) then
                    root = root + 1
                endif

                ! Step 1b: find the nearest subsequent row r where A(r,i) /= 0
                if (mpi_rank == root) then
                    r = i
                    i_para = i - sum(ncols(1:root))
                    do while (abs(A(r, i_para)) < tiny(0.))
                        r = r + 1
                        if (r > n0) then
                            if (mpi_rank == 0) print '("mpi_matinv: Inverse does not exist")'
                            return
                        endif
                    enddo
                endif

                call MPI_BCAST(r, 1, MPI_INTEGER, root, MPI_COMM_WORLD, mpi_ierr)

                ! Step 1c: swap row i with row r
                if (r > i) then
                    row(:) = A(i, :)
                    A(i, :) = A(r, :)
                    A(r, :) = row(:)

                    row(:) = A_inv(i, :)
                    A_inv(i, :) = A_inv(r, :)
                    A_inv(r, :) = row(:)
                endif

                ! Step 2: divide row i by scale = A(i,i)
                scale = A(i, i_para)
                call MPI_BCAST(scale, 1, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, mpi_ierr)

                do j = 1, ncol
                    A(i, j) = A(i, j) / scale
                    A_inv(i, j) = A_inv(i, j) / scale
                enddo

                if (mpi_rank == root) then
                    factor(:) = A(:, i_para)
                endif

                call MPI_BCAST(factor(:), n0, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, mpi_ierr)

                ! Step 3: change all rows below row i to take the previous modifications into account
                if (i < n0) then
                    do r = i + 1, n0
                        do j=1,ncol
                            A(r, j) = A(r, j) - factor(r) * A(i, j)
                            A_inv(r, j) = A_inv(r, j) - factor(r) * A_inv(i, j)
                        enddo
                    enddo
                endif
            enddo

            ! Back substitution phase
            root = mpi_nproc - 1
            jz_para = 1

            do jz = n0, 2, -1
                ! Step 1: find the rank of the process where column jz (zeroing column) is
                if (jz < displs_j0(root + 1) + 1) then
                    root = root - 1
                endif

                if (mpi_rank == root) then
                    jz_para = jz - sum(ncols(1:root))
                    factor(:) = A(:, jz_para)
                endif
                call MPI_BCAST(factor(:), n0, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, mpi_ierr)

                ! Step 2: transform matrix A_inv into A^-1 (doing the same transformation on A gives the Id matrix)
                do i = jz - 1, 1, -1
                    do j = 1, ncol
                        A_inv(i, j) = A_inv(i, j) - factor(i) * A_inv(jz, j)
                    enddo
                enddo
            enddo

            ! Gather inversed matrix
            size_A = size(A(:, :))
            call MPI_ALLGATHER(size_A, 1, MPI_INT, sizes, 1, MPI_INT, MPI_COMM_WORLD, mpi_ierr)

            j = sum(sizes(1:mpi_rank))
            call MPI_ALLGATHER(j, 1, MPI_INT, displs_j, 1, MPI_INT, MPI_COMM_WORLD, mpi_ierr)

            call MPI_ALLGATHERV(A_inv(:, :), size_A, MPI_DOUBLE_PRECISION, &
                                A0_inv(:, :), sizes(:), displs_j(:), MPI_DOUBLE_PRECISION, &
                                MPI_COMM_WORLD, mpi_ierr)

            deallocate(A, A_inv, row)

        end subroutine mpi_matinv

end module mpi_utils


