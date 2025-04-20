!! Highly coupled module for the bisection method. Totally not sarcastic :sadface: Yes I am aware of it.
! The problem with decoupling this is that next_set is very tightly used in the other procedures.
! It is not worth to decouple for extraction of a single type.
module bisection
    use algebdef
    implicit none
    type, extends(algeb_prob) :: bisection_solver
        contains
        procedure :: set_first_interval
        procedure :: next_set => next_interval_set
        procedure :: bisect => do_bisection
    end type bisection_solver

    contains

    ! Performs bisection on the current [a, b] interval
    ! @params iter: Max number of iterations (default: 15)
    ! @params tolerance: Convergence threshold (default: 1e-6)
    ! @result root: Approximated root of the equation
    function do_bisection(this, iter, tolerance) result(root)
        class(bisection_solver), intent(inout) :: this
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
            tol = 1.0e-5 !!  NOTE: This affects testing. This indirectly 
            !!        controls iterations done. Hence adjust accordingly.
        end if

        bisect: do i = 1, trial
            call this%next_set()
            ! call this%print()
            if (abs(this%f(this%c)) < tol .or. abs(this%b - this%a) < tol) then
                root = this%c
                print *, " f(c) < tol @ iter: ", i
                return
            end if

        end do bisect
    end function do_bisection

    ! This method is not really to be used. It only illustrates the idea
    ! that the first interval is to be bruteforced
    subroutine set_first_interval(this, from, to)
        class(bisection_solver), intent(inout) :: this
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
        class(bisection_solver), intent(inout) :: this

        if (this%iter_zero) stop "Please use set_first_interval or set iter_zero false by specifying first interval."

        this%c = (this%a + this%b) / 2
        if (this%f(this%a) * this%f(this%c) < 0) then
            this%b = this%c
        else
            this%a = this%c
        end if
    end subroutine next_interval_set
end module bisection