module bisection
    implicit none
    type :: algeb_prob
        logical :: iter_zero = .true.
        procedure(eqn_interface), pointer, nopass :: f => null()
        real :: a, b, c
        contains
        procedure :: print => print_prob
        procedure :: set_first_interval
        procedure :: next_set => next_interval_set
        procedure :: bisect => do_bisection
    end type algeb_prob

    abstract interface
        function eqn_interface(x) result(res)
            real, intent(in) :: x
            real :: res
        end function eqn_interface
    end interface

    contains

    ! Performs bisection on the current [a, b] interval
    ! @params iter: Max number of iterations (default: 15)
    ! @params tolerance: Convergence threshold (default: 1e-6)
    ! @result root: Approximated root of the equation
    function do_bisection(this, iter, tolerance) result(root)
        class(algeb_prob), intent(inout) :: this
        integer, optional, intent(in) :: iter
        real, optional, intent(in) :: tolerance
        integer :: trial, i
        real :: tol, root

        if (this%f(this%a) * this%f(this%b) > 0.0) then
            print *, "Warning: f(a) and f(b) have the same sign. Bisection might fail."
        end if

        if (present(iter)) then
            trial = iter
        else
            trial = 15
        end if

        if (present(tolerance)) then
            tol = tolerance
        else
            tol = 1.0e-6
        end if

        bisect: do i = 1, trial
            call this%next_set()
            ! call this%print()
            if (abs(this%f(this%c)) < tol) then
                root = this%c
                print *, " f(c) < tol @ iter: ", i
                return
            end if

            if (abs(this%b - this%a) < tol) then
                root = this%c
                print *, " b - a < tol @ iter: ", i
                return
            end if
        end do bisect
    end function do_bisection

    ! This method is not really to be used. It only illustrates the idea
    ! that the first interval is to be bruteforced
    subroutine set_first_interval(this, from, to)
        class(algeb_prob), intent(inout) :: this
        integer, intent(in) :: from, to
        integer :: i
        if (.not. this%iter_zero) then
            call next_interval_set(this)
            return
        end if
        ! Searching iter 0 for a,b inits currently using symmteric and monotonic bounds
        ! Basically search in a predefined range for now
        ! TODO double loops for bruteforcing.
        first: do i = from, to
            if (this%f(real(i - 1)) * this%f(real(i)) < 0) then
                this%a = i - 1
                this%b = i
                this%iter_zero = .false.
                exit first
            end if
        end do first

        if (this%iter_zero) then
            second: do i = from, to
                if (this%f(real(i - 1)) * this%f(real(-(i - 1))) < 0) then
                    this%a = i - 1
                    this%b = -(i - 1)
                    this%iter_zero = .false.
                    exit second
                end if
            end do second
        end if

    end subroutine set_first_interval

    ! Actual method to change poisition of a, b post iter 0
    subroutine next_interval_set(this)
        class(algeb_prob), intent(inout) :: this

        if (this%iter_zero) stop "Please use set_first_interval or set iter_zero false by specifying first interval."

        this%c = (this%a + this%b) / 2
        if (this%f(this%a) * this%f(this%c) < 0) then
            this%b = this%c
        else
            this%a = this%c
        end if
    end subroutine next_interval_set

    subroutine print_prob(this)
        class(algeb_prob), intent(in) :: this
        print *, "Problem:[a= ", this%a, " b= ", this%b, &
        " c= ", this%c, " iter_zero= ", this%iter_zero, "]"
    end subroutine print_prob

    function init_prob(eqn, a, b, iter_zero) result(aprob)
        type(algeb_prob) :: aprob
        procedure(eqn_interface), pointer :: eqn
        real, intent(in) :: a, b
        logical, optional, intent(in) :: iter_zero

        aprob%a = a
        aprob%b = b
        if (present(iter_zero)) aprob%iter_zero = iter_zero
        aprob%f => eqn
    end function init_prob
end module bisection