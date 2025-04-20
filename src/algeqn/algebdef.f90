module algebdef
    implicit none
    type :: algeb_prob
        logical :: iter_zero = .true.
        procedure(eqn_interface), pointer, nopass :: f => null()
        real :: a, b, c
        contains
        procedure :: init => init_prob
        procedure :: print => print_prob

    end type algeb_prob

    abstract interface
        function eqn_interface(x) result(res)
            real, intent(in) :: x
            real :: res
        end function eqn_interface
    end interface

    contains

    subroutine print_prob(this)
        class(algeb_prob), intent(in) :: this
        print *, "Problem:[a= ", this%a, " b= ", this%b, &
        " c= ", this%c, " iter_zero= ", this%iter_zero, "]"
    end subroutine print_prob

    subroutine init_prob(this, eqn, a, b, iter_zero)
        class(algeb_prob), intent(inout) :: this
        procedure(eqn_interface), pointer :: eqn
        real, intent(in) :: a, b
        logical, optional, intent(in) :: iter_zero

        this%a = a
        this%b = b
        if (present(iter_zero)) this%iter_zero = iter_zero
        this%f => eqn
    end subroutine init_prob

end module algebdef