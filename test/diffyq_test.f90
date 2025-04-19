module diffyq_test
    use diffyq
    use rk_two
    use rk_fourth
    use testdrive, only : new_unittest, unittest_type, error_type, check
    implicit none
    private

    public :: collect_suite_diffyq
    contains
    subroutine collect_suite_diffyq(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [&
        new_unittest("RK: init_data", test_valid_init), &
        new_unittest("RK_Two: y_req", test_valid_two), &
        new_unittest("RK_Two: y_req(heun)", test_valid_two_heun), &
        new_unittest("RK_Fourth: y_req", test_valid_fourth), &
        new_unittest("Euler: y_req", test_valid_euler) &
        ! new_unittest("invalid", test_invalid, should_fail=.true.) &
        ]

    end subroutine collect_suite_diffyq

    subroutine setup(test_prob)
        type(diffyq_prob), intent(out) :: test_prob
        procedure(eqn_interface), pointer :: f_ptr => null()

        f_ptr => test_f
        test_prob = init_prob(0.0, 1.0, 0.1, 0.2, f_ptr)
    end subroutine setup

    function test_f(t, y) result(out)
        real, intent(in) :: t, y
        real :: out
        out = (t - y) / 2.0
    end function test_f

    subroutine test_valid_init(error)
        type(error_type), allocatable, intent(out) :: error
        type(diffyq_prob) :: test_prob
        call setup(test_prob)

        call check(error, test_prob%t, 0.0)
        call check(error, test_prob%y, 1.0)
        call check(error, test_prob%h, 0.1)
        call check(error, test_prob%x_p, 0.2)
        call check(error, associated(test_prob%f), .true.)
        if (allocated(error)) return
    end subroutine test_valid_init

    subroutine test_valid_euler(error)
        type(error_type), intent(out), allocatable :: error
        
        ! call setup(test_prob)
        ! call check(error, abs(solver%solve - ))
    end subroutine test_valid_euler

    subroutine test_valid_two(error)
        type(error_type), allocatable, intent(out) :: error
        class(solver_base), allocatable :: solver
        type(diffyq_prob) :: test_prob
        real, parameter :: tolerance = 1e-5

        call setup(test_prob)
        allocate(rk_two_solve :: solver)
        if (allocated(solver)) then
            call check(error, abs(solver%solve(test_prob) - 0.91463) < tolerance)
            if (allocated(error)) return
        end if
    end subroutine test_valid_two

    subroutine test_valid_two_heun(error)
        type(error_type), allocatable, intent(out) :: error
        class(solver_base), allocatable :: solver
        type(diffyq_prob) :: test_prob
        real, parameter :: tolerance = 1e-5

        call setup(test_prob)
        allocate(rk_two_solve :: solver)
        if (allocated(solver)) then
            call check(error, abs(solver%solve(test_prob) - 0.91463) < tolerance)
            if (allocated(error)) return
        end if
    end subroutine test_valid_two_heun

    subroutine test_valid_fourth(error)
        type(error_type), allocatable, intent(out) :: error
        class(solver_base), allocatable :: solver
        type(diffyq_prob) :: test_prob
        real, parameter :: tolerance = 1e-5

        call setup(test_prob)
        allocate(rk_fourth_solve :: solver)
        if (allocated(solver)) then
            call check(error, abs(solver%solve(test_prob) - 0.91451) < tolerance)
            if (allocated(error)) return
        end if
    end subroutine test_valid_fourth

    ! subroutine test_invalid(error)
    !     type(error_type), allocatable, intent(out) :: error
    !     ! ...
    ! end subroutine test_invalid

end module diffyq_test