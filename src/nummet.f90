module nummet
   implicit none
   private

   public :: say_hello
   contains
   subroutine say_hello
      print *, "Hello, nummet!"
   end subroutine say_hello

   subroutine read_table()
      integer :: n, i
      real :: x_p

      is_init = .true.
      print *, "Enter number of data points: "
      read *, n

      allocate(x(n), y(n))

      print *, "Enter the x values: "
      do i = 1, n
          read *, x(i)
      end do

      print *, "Enter the y values: "
      do i = 1, n
          read *, y(i)
      end do

      allocate(diff_table(n, n))

      if (.not. allocated(x)) stop "Array X is not allocated"
      if (.not. allocated(y)) stop "Array Y is not allocated"
      if (.not. allocated(y)) stop "Array DIFF_TABLE is not allocated"

      print *, "Enter x to interpolate: "
      read *, x_p

  end subroutine read_table
end module nummet
