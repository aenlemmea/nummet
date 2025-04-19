module range_kutta ! For first order derivatives
    implicit none

    type, public :: rk_prob
        logical :: is_first_iter = .true.
        real :: t, y, h, x_p, y_req = -1.0 ! Default "shadow" value for error handling
        procedure(eqn_interface), pointer, nopass :: f => null()
        contains
        procedure :: print => print_problem
        procedure :: next => next_iter
    end type rk_prob

    type, abstract, public :: range_kutta_base
        contains
        procedure(solve_interface), deferred :: solve
    end type range_kutta_base

    abstract interface ! Mostly for type safety
        function solve_interface(self, prob) result(y_req)
            import :: range_kutta_base, rk_prob
            class(range_kutta_base), intent(inout) :: self
            type(rk_prob), intent(inout) :: prob
            real :: y_req
        end function solve_interface

        function eqn_interface(t, y)
            real :: eqn_interface
            real, intent(in) :: t, y
        end function eqn_interface
    end interface
    contains

    subroutine next_iter(this)
        class(rk_prob), intent(inout) :: this
        if (this%is_first_iter) stop "Trying to go to next iteration while first is not complete."
        this%y = this%y_req
        this%t = this%t + this%h
    end subroutine next_iter

    function init_prob(t, y, h, x_p, f) result(aprob)
        real, intent(in) :: t, y, h, x_p
        procedure(eqn_interface), pointer :: f
        type(rk_prob) :: aprob

        aprob%t = t
        aprob%y = y
        aprob%h = h
        aprob%x_p = x_p
        aprob%f => f
    end function init_prob

    subroutine print_problem(this)
        class(rk_prob), intent(in) :: this
        print *, "Problem:[eqn= ", this%f(this%t, this%y), " t= ", this%t, " y= ", this%y, " h= ", this%h, " x_p= ", this%x_p, "]"
        ! if (associated(this%f)) then
        !     print *, "Found an associated f"
        ! else
        !     print *, "f is not associated. Please init."
        ! end if
    end subroutine print_problem
end module range_kutta