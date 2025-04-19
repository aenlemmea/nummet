! Module for the Range Kutta 4th Order Method
! k1 = h * f(t, y)
! k2 = h * f(t + h/2, y + k1/2)
! k3 = h * f(t + h/2, y + k2/2)
! k4 = h * f(t + h, y + k3)

! y_req = y + (1/6)*(k1 + 2*k2 + 2*k3 + k4)
module rk_fourth
    use range_kutta
    implicit none

    type, extends(range_kutta_base) :: rk_fourth_solve
        contains
        procedure :: solve => solve_fourth
    end type rk_fourth_solve
    
    contains

    function solve_fourth(self, prob) result(y_req)
        class(rk_fourth_solve), intent(inout) :: self
        type(rk_prob), intent(inout) :: prob
        real :: y_req
        integer :: steps, i

        steps = int((prob%x_p - prob%t) / prob%h + 0.5)

        do i = 1, steps
            prob%y_req = solve_fourth_once(self, prob)
            call prob%next()
        end do
        y_req = prob%y_req
    end function solve_fourth

    function solve_fourth_once(self, prob) result(y_req)
        class(rk_fourth_solve), intent(inout) :: self
        type(rk_prob), intent(inout) :: prob
        real :: y_req, k1, k2, k3, k4

        k1 = prob%h * prob%f(prob%t, prob%y)
        k2 = prob%h * prob%f((prob%t + (prob%h / 2.0)), (prob%y + (k1 / 2.0)))
        k3 = prob%h * prob%f((prob%t + (prob%h / 2.0)), (prob%y + (k2 / 2.0)))
        k4 = prob%h * prob%f(prob%t + prob%h, prob%y + k3)

        prob%is_first_iter = .false.
        y_req = prob%y + ((1.0 / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4))
    end function solve_fourth_once
end module rk_fourth