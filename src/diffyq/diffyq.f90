module diffyq ! For first order derivatives
    implicit none

    type, public :: diffyq_prob
        logical :: is_first_iter = .true.
        real :: t, y, h, x_p, y_req = -1.0 ! Default "shadow" value for error handling
        procedure(eqn_interface), pointer, nopass :: f => null()
        contains
        procedure :: print => print_problem
        procedure :: next => next_iter
    end type diffyq_prob

    ! Base solver type to offer a single point of entry "solve"
    type, abstract, public :: solver_base
        contains
        procedure(solve_interface), deferred :: solve
    end type solver_base

    abstract interface ! Mostly for type safety
        function solve_interface(self, prob) result(y_req)
            import :: solver_base, diffyq_prob
            class(solver_base), intent(inout) :: self
            type(diffyq_prob), intent(inout) :: prob
            real :: y_req
        end function solve_interface

        function eqn_interface(t, y)
            real :: eqn_interface
            real, intent(in) :: t, y
        end function eqn_interface
    end interface
    contains

    subroutine next_iter(this)
        class(diffyq_prob), intent(inout) :: this
        if (this%is_first_iter) stop "Trying to go to next iteration while first is not complete."
        this%y = this%y_req
        this%t = this%t + this%h
    end subroutine next_iter

    function init_prob(t, y, h, x_p, f) result(aprob)
        real, intent(in) :: t, y, h, x_p
        procedure(eqn_interface), pointer :: f
        type(diffyq_prob) :: aprob

        aprob%t = t
        aprob%y = y
        aprob%h = h
        aprob%x_p = x_p
        aprob%f => f
    end function init_prob

    subroutine print_problem(this)
        class(diffyq_prob), intent(in) :: this
        print *, "Problem:[eqn= ", this%f(this%t, this%y), " t= ", this%t, " y= ", this%y, " h= ", this%h, " x_p= ", this%x_p, "]"
        ! if (associated(this%f)) then
        !     print *, "Found an associated f"
        ! else
        !     print *, "f is not associated. Please init."
        ! end if
    end subroutine print_problem
end module diffyq