module newton
    implicit none
    real, allocatable :: x(:)
    real, allocatable :: y(:)
    real, allocatable :: diff_table(:, :)
    logical :: is_init = .false.
    private :: is_init

    interface
        module subroutine generate_div_diff(n)
            integer, intent(in) :: n
        end subroutine generate_div_diff

        module function divdiff(n, x_p) result(y_xp)
            integer, intent(in) :: n
            real, intent(in) :: x_p
            real :: y_xp
        end function divdiff
    end interface

    contains

    ! Free the allocated arrays
    subroutine free_arrays()
        if (.not. is_init) stop "Uninitialized"
        if (allocated(x)) deallocate(x)
        if (allocated(y)) deallocate(y)
        if (allocated(diff_table)) deallocate(diff_table)
    end subroutine free_arrays

    ! Newton's Forward Interpolation
    ! @params n: size of x array
    ! @params x_p : Point to interpolate
    ! @result y_xp: The value of x_p post interpolation
    ! Picked when |x_p - x(0)| < eps
    function forward(n, x_p) result(y_xp)
        integer, intent(in) :: n
        real :: x_p, y_xp, h, p, fact, top
        integer :: i
        ! h: step length

        if (.not. is_init) stop "Uninitialized"
        call generate_diff(n)

        h = x(2) - x(1)
        p = (x_p - x(1)) / h

        y_xp = diff_table(1, 1) ! First iteration
        fact = 1.0d0
        top = 1.0d0

        ! 1 to n - 1 => n - 1 iterations.
        do i = 1, n - 1
            top = top * (p - i + 1)
            fact = fact * i
            y_xp = y_xp + (top / fact) * diff_table(1, i + 1)
        end do
        call free_arrays()
    end function forward

    ! Newton's Backward Interpolation
    ! @params n: size of x array
    ! @params x_p : Point to interpolate
    ! @result y_xp: The value of x_p post interpolation
    ! picked when |x_p - x(n)| < eps
    function backward(n, x_p) result(y_xp)
        integer, intent(in) :: n
        real :: x_p, y_xp, h, p, fact, top
        integer :: i
        ! h: step length

        if (.not. is_init) stop "Uninitialized"
        call generate_diff(n)

        h = x(2) - x(1)
        p = (x_p - x(n)) / h

        y_xp = diff_table(n, 1) ! First iteration
        fact = 1.0d0
        top = 1.0d0

        ! 1 to n - 1 => n - 1 iterations.
        do i = 1, n - 1
            top = top * (p + i - 1)
            fact = fact * i
            y_xp = y_xp + (top / fact) * diff_table(n - i, i + 1)
        end do
        call free_arrays()
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

submodule (newton) divided_difference
implicit none
contains
subroutine generate_div_diff(n)
    integer, intent(in) :: n
    integer :: i, j

    diff_table(:, 1) = y
    do j = 2, n
        do i = 1, n - j + 1
            ! In the code below, the j - 1 in x(i + j - 1) defines the "gap" between the two variables in the denominator
            diff_table(i, j) = (diff_table(i + 1, j - 1) - diff_table(i, j - 1)) / (x(i + j - 1) - x(i))
        end do
    end do
end subroutine generate_div_diff

! Newton's Divided Difference Interpolation
! @params n: size of x array
! @params x_p : Point to interpolate
! @result y_xp: The value of x_p post interpolation
! picked when unequal intervals
function divdiff(n, x_p) result(y_xp)
    integer, intent(in) :: n
    real, intent(in) :: x_p 
    real :: y_xp, side
    integer :: i

    if (.not. is_init) stop "Uninitialized"

    call generate_div_diff(n)

    y_xp = diff_table(1, 1) ! First iteration
    side = 1.0  

    do i = 2, n
        side = side * (x_p - x(i - 1))
        y_xp = y_xp + diff_table(1, i) * side 
    end do
    call free_arrays()
end function divdiff
end submodule divided_difference