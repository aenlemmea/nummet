module interpolation_test
    use newton
    use lagrange
    use testdrive, only : new_unittest, unittest_type, error_type, check
    implicit none
    private

    public :: collect_suite_interpolation
    contains

    subroutine collect_suite_interpolation(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [&
        new_unittest("Read Test", test_valid_read), &
        new_unittest("Diff Table Gen", test_valid_gen_diff_table), &
        new_unittest("Backward y_xp", test_valid_backward), &
        new_unittest("Forward y_xp", test_valid_forward), &
        new_unittest("Divided Diff Table Gen", test_valid_gen_div_diff_table), &
        new_unittest("Divided Diff y_xp", test_valid_divdiff), &
        new_unittest("Lagrange y_xp", test_valid_lagrange) &
        !new_unittest("invalid", test_invalid_read, should_fail=.true.) &
        ]
    end subroutine collect_suite_interpolation

    subroutine setup(x_sz)
        real, parameter :: x_in(5) = [1891, 1901, 1911, 1921, 1931]
        real, parameter :: y_in(5) = [46, 66, 81, 93, 101]
        integer, intent(out), optional :: x_sz

        x_sz = size(x_in)

        call read_table_from_array(x_in, y_in)
    end subroutine setup
    
    subroutine test_valid_lagrange(error) 
        type(error_type), allocatable, intent(out) :: error
        integer :: x_sz
        real :: expected 
        real, parameter :: tolerance = 1e-6

        call setup(x_sz)
        expected = lagrange_interpolate(x_sz, 1896.0)
        call check(error, abs(expected - 56.8671875) < tolerance)
        call free_arrays()
        if (allocated(error)) return
    end subroutine test_valid_lagrange

    subroutine test_valid_divdiff(error) 
        type(error_type), allocatable, intent(out) :: error
        integer :: x_sz
        real :: expected
        real, parameter :: tolerance = 1e-6

        call setup(x_sz)
        expected = divdiff(x_sz, 1896.0)
        call check(error, abs(expected - 56.8671874875) < tolerance)
        call free_arrays()
        if (allocated(error)) return
    end subroutine test_valid_divdiff

    subroutine test_valid_gen_div_diff_table(error)
        type(error_type), allocatable, intent(out) :: error
        integer :: x_sz
        real, parameter :: tolerance = 1e-6

        call setup(x_sz)
        call generate_div_diff(x_sz)

        ! Strangely the values of this diff_table are similar to
        ! the one in forward and backward tables,
        ! except they are divided by 10 for the first column
        call check(error, abs(diff_table(1, 2) - 2.0000) < tolerance)
        call check(error, abs(diff_table(2, 2) - 1.5000) < tolerance)
        call check(error, abs(diff_table(3, 2) - 1.2000) < tolerance)
        call check(error, abs(diff_table(4, 2) - 0.8000) < tolerance)

        call check(error, abs(diff_table(1, 3) - (-0.0250)) < tolerance)
        call check(error, abs(diff_table(2, 3) - (-0.0150)) < tolerance)
        call check(error, abs(diff_table(3, 3) - (-0.0200)) < tolerance)

        call check(error, abs(diff_table(1, 4) - 0.0003333) < tolerance)
        call check(error, abs(diff_table(2, 4) - (-0.000166)) < tolerance)

        call check(error, abs(diff_table(1, 5) - (-0.0000125)) < tolerance)

        call free_arrays()
        if (allocated(error)) return
    end subroutine test_valid_gen_div_diff_table


    subroutine test_valid_forward(error)
        type(error_type), allocatable, intent(out) :: error
        integer :: x_sz
        real, parameter :: tolerance = 1e-6
        real :: expected

        call setup(x_sz)
        expected = forward(x_sz, 1895.0)
        call check(error, abs(expected - 54.8528) < tolerance)
        call free_arrays()
        if (allocated(error)) return
    end subroutine test_valid_forward


    subroutine test_valid_backward(error)
        type(error_type), allocatable, intent(out) :: error
        integer :: x_sz
        real, parameter :: tolerance = 1e-6
        real :: expected

        call setup(x_sz)
        expected = backward(x_sz, 1925.0)
        call check(error, abs(expected - 96.8368) < tolerance)
        call free_arrays()
    end subroutine test_valid_backward

    ! | A(1,1)   A(2,1)   A(3,1)   A(4,1)   A(5,1)  |   -> First column = y
    ! | A(1,2)   A(2,2)   A(3,2)   A(4,2)           |   -> Second column
    ! | A(1,3)   A(2,3)   A(3,3)                    |   -> Third column
    ! | A(1,4)   A(2,4)                             |   -> Fourth column
    ! | A(1,5)                                      |   -> Fifth column
    subroutine test_valid_gen_diff_table(error)
        type(error_type), allocatable, intent(out) :: error
        real, parameter :: tolerance = 1e-6
        integer :: x_sz

        call setup(x_sz)
        call generate_diff(x_sz)

        call check(error, abs(diff_table(1, 2) - 20.0000) < tolerance)
        call check(error, abs(diff_table(2, 2) - 15.0000) < tolerance)
        call check(error, abs(diff_table(3, 2) - 12.0000) < tolerance)
        call check(error, abs(diff_table(4, 2) - 8.0000) < tolerance)

        call check(error, abs(diff_table(1, 3) - (-5.000)) < tolerance)
        call check(error, abs(diff_table(2, 3) - (-3.0000)) < tolerance)
        call check(error, abs(diff_table(3, 3) - (-4.0000)) < tolerance)

        call check(error, abs(diff_table(1, 4) - 2.0000) < tolerance)
        call check(error, abs(diff_table(2, 4) - (-1.0000)) < tolerance)

        call check(error, abs(diff_table(1, 5) - (-3.0000)) < tolerance)

        call free_arrays()

        if (allocated(error)) return
    end subroutine test_valid_gen_diff_table


    subroutine test_valid_read(error)
        type(error_type), allocatable, intent(out) :: error
        real, parameter :: tolerance = 1e-6

        real :: x_in(3) = [1, 2, 3]
        real :: y_in(3) = [10, 20, 30]

        call read_table_from_array(x_in, y_in)

        call check(error, abs(x(1) - 1) < tolerance)
        call check(error, abs(x(2) - 2) < tolerance)
        call check(error, abs(x(3) - 3) < tolerance)

        if (allocated(error)) return
        call check(error, abs(y(1) - 10) < tolerance)
        call check(error, abs(y(2) - 20) < tolerance)
        call check(error, abs(y(3) - 30) < tolerance)
        if (allocated(error)) return

        call free_arrays()
    end subroutine test_valid_read

    subroutine test_invalid_read(error)
        type(error_type), allocatable, intent(out) :: error
        ! ...
    end subroutine test_invalid_read
end module interpolation_test