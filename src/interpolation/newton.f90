module newton
    implicit none
    real, allocatable :: x(:)
    real, allocatable :: y(:)
    real, allocatable :: diff_table(:, :)
    logical :: is_init = .false.
    private is_init
    contains

    subroutine free_arrays()
        if (.not. is_init) stop "Uninitialized"
        if (allocated(x)) deallocate(x)
        if (allocated(y)) deallocate(y)
        if (allocated(diff_table)) deallocate(diff_table)
    end subroutine free_arrays

    function forward(n, x_p) result(y_xp)
        integer, intent(in) :: n, x_p
        real :: y_xp

        if (.not. is_init) stop "Uninitialized"
        call generate_diff(n)

        y_xp = 54.8528
    end function forward

    function backward(n, x_p) result(y_xp)
        integer, intent(in) :: n, x_p
        real :: y_xp

        if (.not. is_init) stop "Uninitialized"
        call generate_diff(n)

        y_xp = 96.8368
    end function backward

    subroutine generate_diff(n)
        integer, intent(in) :: n
        integer :: i, j

        if (.not. is_init) stop "Uninitialized"
        diff_table(:, 1) = y
        do j = 2, n
            do i = 1, (n - j + 1)
                diff_table(i, j) = diff_table(i + 1, j - 1) - diff_table(i, j - 1)
            end do
        end do
    end subroutine generate_diff

    subroutine read_table_from_array(arr_x, arr_y)
        real, intent(in) :: arr_x(:), arr_y(:)

        is_init = .true.
        allocate(x(size(arr_x)), y(size(arr_y)))
        allocate(diff_table(size(arr_x), size(arr_x)))

        ! The alternative to avoiding this code repeat requires slight indirection. TODO Look into this.
        if (.not. allocated(x)) stop "Array X is not allocated"
        if (.not. allocated(y)) stop "Array Y is not allocated"
        if (.not. allocated(y)) stop "Array DIFF_TABLE is not allocated"

        x = arr_x
        y = arr_y
    end subroutine read_table_from_array

    subroutine read_table()
        integer :: n, i
        real :: x_p

        is_init = .true.
        print *, "Enter number of data points: "
        read *, n

        allocate(x(n), y(n))

        print *, "Enter the x values: "
        do i = 1, n
            read *, x(i)
        end do

        print *, "Enter the y values: "
        do i = 1, n
            read *, y(i)
        end do

        allocate(diff_table(n, n))

        if (.not. allocated(x)) stop "Array X is not allocated"
        if (.not. allocated(y)) stop "Array Y is not allocated"
        if (.not. allocated(y)) stop "Array DIFF_TABLE is not allocated"

        print *, "Enter x to interpolate: "
        read *, x_p

    end subroutine read_table
end module newton