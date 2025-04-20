module algebqn_test
    use bisection
    use testdrive, only : new_unittest, unittest_type, error_type, check
    implicit none
    private

    public :: collect_suite_algebqn
    contains
    subroutine collect_suite_algebqn(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [&
        new_unittest("Bisection: init_data", test_valid_init), &
        new_unittest("Bisection: Iter 0 bruteforce", test_valid_iter0_interval_bruteforce), &
        new_unittest("Bisection: Next interval", test_valid_interval_find), &
        new_unittest("Bisection: root", test_valid_bisection), &
        new_unittest("Bisection: Iter 0 bruteforce > Invalid", test_invalid_iter0_interval_bruteforce, should_fail = .true.) &
        ]

    end subroutine collect_suite_algebqn

    subroutine setup(test_prob)
        type(bisection_solver), intent(out) :: test_prob
        procedure(eqn_interface), pointer :: f_ptr => null()
        f_ptr => test_f
        call test_prob%init(f_ptr, -1.0, -1.0)
    end subroutine setup

    function test_f(x) result(out)
        real, intent(in) :: x
        real :: out
        out = (x * x * x) - x - 1
    end function test_f

    function test_g(x) result(out)
        real, intent(in) :: x
        real :: out
        out = (2 * (x * x * x)) - (2 * x) - 5
    end function test_g

    subroutine setup_sec(test_prob)
        type(bisection_solver), intent(out) :: test_prob
        procedure(eqn_interface), pointer :: g_ptr => null()
        g_ptr => test_g
        call test_prob%init(g_ptr, 1.0, 2.0, .false.)
    end subroutine setup_sec

    subroutine test_valid_init(error)
        type(error_type), allocatable, intent(out) :: error
        type(bisection_solver) :: test_prob
        call setup(test_prob)

        call check(error, test_prob%a, -1.0)
        call check(error, test_prob%b, -1.0)
        call check(error, test_prob%iter_zero, .true.)
        call check(error, associated(test_prob%f), .true.)
        if (allocated(error)) return
    end subroutine test_valid_init

    subroutine test_valid_iter0_interval_bruteforce(error)
        type(error_type), allocatable, intent(out) :: error
        type(bisection_solver) :: test_prob
        real, parameter :: tolerance = 1e-5

        call setup(test_prob)
        call test_prob%set_first_interval(1, 100)
        call check(error, abs(test_prob%a - 1.0) < tolerance)
        if (allocated(error)) return
        call check(error, abs(test_prob%b - 2.0) < tolerance)
        if (allocated(error)) return
    end subroutine test_valid_iter0_interval_bruteforce

    subroutine test_valid_interval_find(error)
        type(error_type), allocatable, intent(out) :: error
        type(bisection_solver) :: test_prob
        real, parameter :: tolerance = 1e-5

        call setup(test_prob)
        call test_prob%set_first_interval(1, 5)
        call test_prob%next_set()
        call check(error, abs(test_prob%a - 1.0) < tolerance)
        if (allocated(error)) return
        call check(error, abs(test_prob%b - 1.5) < tolerance)
        if (allocated(error)) return
    end subroutine test_valid_interval_find

    subroutine test_invalid_iter0_interval_bruteforce(error)
        type(error_type), allocatable, intent(out) :: error
        type(bisection_solver) :: test_prob

        call setup(test_prob)
        call test_prob%set_first_interval(1, -100)
        call check(error, test_prob%iter_zero, .false.)
        if (allocated(error)) return
    end subroutine test_invalid_iter0_interval_bruteforce

    subroutine test_valid_bisection(error)
        type(error_type), allocatable, intent(out) :: error
        type(bisection_solver) :: test_prob
        real, parameter :: tolerance = 1.0e-5

        call setup(test_prob)
        call test_prob%set_first_interval(1, 5)
        call check(error, abs(test_prob%bisect(20) - 1.32471) < tolerance)
        if (allocated(error)) return

        call setup_sec(test_prob)
        call check(error, abs(test_prob%bisect(25) - 1.60059) < tolerance)
        if (allocated(error)) return
    end subroutine test_valid_bisection



    ! subroutine test_invalid(error)
    !     type(error_type), allocatable, intent(out) :: error
    !     ! ...
    ! end subroutine test_invalid

end module algebqn_test