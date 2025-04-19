! Module for the Range Kutta 4th Order Method
! k1 = h * f(t, y)
! k2 = h * f(t + h/2, y + k1/2)
! y_req = y + k2

! OR=====================Heun's:

! k1 = h * f(t, y)
! k2 = h * f(t + h, y + k1)
! y_req = y + 0.5 * (k1 + k2)

module rk_two
    use range_kutta
    implicit none

    type, extends(range_kutta_base) :: rk_two_solve
        contains
        procedure :: solve => solve_two
    end type rk_two_solve

    contains

    function solve_two(self, prob) result(y_req)
        class(rk_two_solve), intent(inout) :: self
        type(rk_prob), intent(inout) :: prob
        real :: y_req
        integer :: steps, i

        steps = int((prob%x_p - prob%t) / prob%h + 0.5)

        do i = 1, steps
            prob%y_req = solve_two_once(self, prob)
            call prob%next()
        end do
        y_req = prob%y_req
    end function solve_two

    function solve_two_once(self, prob) result(y_req)
        class(rk_two_solve), intent(inout) :: self
        type(rk_prob), intent(inout) :: prob
        real :: y_req, k1, k2, k3, k4

        k1 = prob%h * prob%f(prob%t, prob%y)
        k2 = prob%h * prob%f((prob%t + (prob%h / 2.0)), (prob%y + (k1 / 2.0)))

        prob%is_first_iter = .false.
        y_req = prob%y + k2
    end function solve_two_once

    function solve_two_once_heun(self, prob) result(y_req)
        class(rk_two_solve), intent(inout) :: self
        type(rk_prob), intent(inout) :: prob
        real :: y_req, k1, k2

        k1 = prob%h * prob%f(prob%t, prob%y)
        k2 = prob%h * prob%f((prob%t + prob%h), (prob%y + k1 / 2.0))

        prob%is_first_iter = .false.
        y_req = prob%y + (0.5 * (k1 + k2))
    end function solve_two_once_heun
end module rk_two