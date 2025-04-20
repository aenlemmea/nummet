module sysqn_test
    use lufac
    use testdrive, only : new_unittest, unittest_type, error_type, check
    implicit none
    private
  
    public :: collect_suite_sysqn
  contains
  
  subroutine collect_suite_sysqn(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
  
    testsuite = [ &
      new_unittest("valid", test_valid) &
    !   new_unittest("invalid", test_invalid, should_fail=.true.) &
      ]
  
  end subroutine collect_suite_sysqn
  
  subroutine test_valid(error)
    type(error_type), allocatable, intent(out) :: error
    call solve_linear_system()
  end subroutine test_valid
  
!   subroutine test_invalid(error)
!     type(error_type), allocatable, intent(out) :: error
!     ! ...
!   end subroutine test_invalid
  
  end module sysqn_test