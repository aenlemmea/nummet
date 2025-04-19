! Module for Euler's Method
! y_req = y + hf(t,y)
module euler
    use diffyq
    implicit none

    type, extends(solver_base) :: euler_solve
        contains
        procedure :: solve => solve_euler
    end type euler_solve

    contains

    function solve_euler(self, prob) result(y_req)
        class(euler_solve), intent(inout) :: self
        type(diffyq_prob), intent(inout) :: prob
        real :: y_req
        integer :: steps, i

        steps = int((prob%x_p - prob%t) / prob%h + 0.5)

        do i = 1, steps
            prob%y_req = solve_euler_once(self, prob) 
            call prob%next()
        end do
        y_req = prob%y_req
    end function solve_euler

    function solve_euler_once(self, prob) result(y_req)
        class(euler_solve), intent(inout) :: self
        type(diffyq_prob), intent(inout) :: prob
        real :: y_req

        prob%is_first_iter = .false.
        y_req = prob%y + prob%h * prob%f(prob%t, prob%y)
    end function solve_euler_once
end module euler