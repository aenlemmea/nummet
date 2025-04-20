! Boilerplate from test_drive library
program check
    use, intrinsic :: iso_fortran_env, only : error_unit
    use testdrive, only : run_testsuite, new_testsuite, testsuite_type
    use interpolation_test, only: collect_suite_interpolation
    use diffyq_test, only: collect_suite_diffyq
    use algebqn_test, only: collect_suite_algebqn
    use sysqn_test, only: collect_suite_sysqn
    implicit none
    integer :: stat, is
    type(testsuite_type), allocatable :: testsuites(:)
    character(len = *), parameter :: fmt = '("#", *(1x, a))'

    stat = 0

    testsuites = [&
    new_testsuite("Interpolation Suite", collect_suite_interpolation), &
    new_testsuite("DiffyQ Suite", collect_suite_diffyq), &
    new_testsuite("Algebqn Suite", collect_suite_algebqn), &
    new_testsuite("Sysqn Suite", collect_suite_sysqn) &
    ! new_testsuite("suite2", collect_suite2) &
    ]

    do is = 1, size(testsuites)
        write(error_unit, fmt) "Testing:", testsuites(is)%name
        call run_testsuite(testsuites(is)%collect, error_unit, stat)
    end do

    if (stat > 0) then
        write(error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
        error stop
    end if

end program check