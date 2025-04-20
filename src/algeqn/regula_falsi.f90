! Converges faster than bisection
module regula_falsi
    use algebdef
    implicit none

    type, extends(algeb_prob) :: regula_falsi_solver
    contains 
        procedure :: find_first
        procedure :: next_set => next_interval_set
        procedure :: regula_falsi => do_regula_falsi
    end type regula_falsi_solver
contains

    subroutine find_first(this)
        class(regula_falsi_solver), intent(inout) :: this
        
    end subroutine find_first

    subroutine next_interval_set(this)
        class(regula_falsi_solver), intent(inout) :: this
        real :: a, b, fa, fb, fc
        
        if (this%iter_zero) error stop "First interval is not set. Either override in algeb_prob or call set_first_interval."

        a = this%a
        b = this%b
        fa = this%f(a)
        fb = this%f(b)

        if (abs(fb - fa) < 1.0e-12) error stop "Division by zero in Regula Falsi."
    
        ! Unlike in bisection, here we don't do signs 
        ! We compute the secant and c is where the secant cuts the x axis
        this%c = (a * fb - b * fa) / (fb - fa) 
        fc = this%f(this%c)
    
        if (fa * fc < 0.0) then
            this%b = this%c
        else
            this%a = this%c
        end if
    end subroutine next_interval_set

    function do_regula_falsi(this, iter, tolerance) result(root)
        class(regula_falsi_solver), intent(inout) :: this
        integer, optional, intent(in) :: iter
        real, optional, intent(in) :: tolerance
        integer :: i, trial
        real :: tol, root

        if (present(iter)) then
            trial = iter
        else
            trial = 15
        end if

        if (present(tolerance)) then
            tol = tolerance
        else
            tol = 1.0e-5    !!  NOTE: This affects testing. This indirectly 
                            !!        controls iterations done. Hence adjust accordingly.
        end if

        do i = 1, trial
            call this%next_set()
            if (abs(this%f(this%c)) < tol) then
                root = this%c
                ! print *, "root: ", root, "@ iter: ", i
                return
            end if
        end do
    end function do_regula_falsi
end module regula_falsi