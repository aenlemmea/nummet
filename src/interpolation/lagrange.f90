module lagrange
    use newton, only: read_table, read_table_from_array, x, y, free_arrays
    implicit none
    contains
    ! Lagrange's Interpolation
    ! @params n: size of x array
    ! @params x_p: Point to interpolate
    ! @result y_xp: The value of x_p post interpolation
    ! Unlike Newton's divided diff, this is not recursive
    ! Works on unequal intervals
    function lagrange_interpolate(n, x_p) result(y_xp)
        integer, intent(in) :: n
        real, intent(in) :: x_p
        real :: y_xp, side
        integer :: i, j

        y_xp = 0.0

        do i = 1, n
            side = y(i)
            do j = 1, n
                if (j /= i) then ! It is not equal to, not divide equal to
                    side = side * (x_p - x(j)) / (x(i) - x(j))
                end if
            end do
            y_xp = y_xp + side
        end do
        call free_arrays()
    end function lagrange_interpolate
end module lagrange