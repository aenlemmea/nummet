module sysqn_test
    use lufac
    use testdrive, only : new_unittest, unittest_type, error_type, check
    implicit none
    private

    public :: collect_suite_sysqn
    contains

    subroutine collect_suite_sysqn(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [&
        new_unittest("valid", test_valid) &
        !   new_unittest("invalid", test_invalid, should_fail=.true.) &
        ]

    end subroutine collect_suite_sysqn

    subroutine test_valid(error)
        type(error_type), allocatable, intent(out) :: error
        integer :: n
        real :: a_in(3, 3), b_in(3), l_expect(3, 3), u_expect(3, 3), x_expect(3)
        real, parameter :: tolerance = 1.0e-4
        n = 3

        l_expect = reshape([1.0, 0.0, 0.0, &
        0.5, 1.0, 0.0, &
        0.5, 0.3333, 1.0], &
        [n, n], order = [2, 1])

        u_expect = reshape([2.0, 4.0, -6.0, &
        0.0, 3.0, 6.0, &
        0.0, 0.0, 3.0], &
        [n, n], order = [2, 1])

        x_expect = [ -3.0, 2.0, 1.0]

        a_in = reshape([2.0, 4.0, -6.0, &
        1.0, 5.0, 3.0, &
        1.0, 3.0, 2.0], [n, n], order = [2, 1])

        b_in = [ -4.0, 10.0, 5.0]

        call solve_linear_system(a_in, b_in, n)
        call check(error, all(abs(l_expect - l) < tolerance))
        if (allocated(error)) return
        call check(error, all(abs(u_expect - u) < tolerance))
        if (allocated(error)) return
        call check(error, all(abs(x_expect - x) < tolerance))
        if (allocated(error)) return

        call free_arrays()
    end subroutine test_valid

    !   subroutine test_invalid(error)
    !     type(error_type), allocatable, intent(out) :: error
    !     ! ...
    !   end subroutine test_invalid

end module sysqn_test