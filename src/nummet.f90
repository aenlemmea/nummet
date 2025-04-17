module nummet
   implicit none
   private

   public :: say_hello
   contains
   subroutine say_hello
      print *, "Hello, nummet!"
   end subroutine say_hello
end module nummet
