! Module to implement LU decomposition
! Not to be used. Use `stdlib_linalg` -> `solve_lu`
module lufac
    implicit none
    real, allocatable :: a(:, :), l(:, :), u(:, :), b(:), y(:), x(:)
    ! External Input: N, A, b
    contains

    subroutine free_arrays()
        if (allocated(a)) deallocate(a)
        if (allocated(l)) deallocate(l)
        if (allocated(u)) deallocate(u)
        if (allocated(b)) deallocate(b)
        if (allocated(y)) deallocate(y)
        if (allocated(x)) deallocate(x)
    end subroutine

    subroutine allocate_arrays(n)
        implicit none
        integer, intent(in) :: n

        if (.not. allocated(a)) allocate(a(n, n))
        if (.not. allocated(l)) allocate(l(n, n))
        if (.not. allocated(u)) allocate(u(n, n))
        if (.not. allocated(b)) allocate(b(n))
        if (.not. allocated(y)) allocate(y(n))
        if (.not. allocated(x)) allocate(x(n))
    end subroutine allocate_arrays

    ! Uses Dolittle's Method
    subroutine solve_linear_system(a_in, b_in, n)
        integer :: n, i, j, k
        real :: sum
        real, intent(in) :: a_in(:, :)
        real, intent(in) :: b_in(:)
        call allocate_arrays(n)

        a = a_in
        b = b_in

        ! Initialize L and U
        l = 0.0d0
        u = 0.0d0
        do i = 1, n
            l(i, i) = 1.0d0
        end do

        ! First row of U
        u(1, 1:n) = a(1, 1:n)

        ! First column of L
        do i = 2, n
            l(i, 1) = a(i, 1) / u(1, 1)
        end do

        ! LU decomposition
        do i = 2, n
            do j = i, n
                sum = 0.0d0
                do k = 1, i - 1
                    sum = sum + l(i, k) * u(k, j)
                end do
                u(i, j) = a(i, j) - sum
            end do
            do k = i + 1, n
                sum = 0.0d0
                do j = 1, i - 1
                    sum = sum + l(k, j) * u(j, i)
                end do
                l(k, i) = (a(k, i) - sum) / u(i, i)
            end do
        end do

        ! Forward substitution: solve Ly = b
        y(1) = b(1)
        do i = 2, n
            sum = 0.0d0
            do j = 1, i - 1
                sum = sum + l(i, j) * y(j)
            end do
            y(i) = b(i) - sum
        end do

        ! Backward substitution: solve Ux = y
        x(n) = y(n) / u(n, n)
        do i = n - 1, 1, -1
            sum = 0.0d0
            do j = i + 1, n
                sum = sum + u(i, j) * x(j)
            end do
            x(i) = (y(i) - sum) / u(i, i)
        end do
    end subroutine solve_linear_system

    subroutine print_arrays(n)
        integer :: i, n
        print *, 'Lower triangular matrix L:'
        do i = 1, n
            print '(100F10.4)', l(i, 1:n)
        end do

        print *, 'Upper triangular matrix U:'
        do i = 1, n
            print '(100F10.4)', u(i, 1:n)
        end do

        print *, 'Solution vector x:'
        do i = 1, n
            print *, x(i)
        end do
    end subroutine print_arrays

end module lufac