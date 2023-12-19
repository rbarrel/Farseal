module CoolingTest
  use iso_fortran_env, only: int32, real32
  use testdrive, only: new_unittest, unittest_type, error_type, &
    check, test_failed, skip_test
  !use Farseal, only: cooling
  implicit none

  private
  public :: collect_cooling_suite

  contains

    subroutine collect_cooling_suite(testsuite)
      type(unittest_type), allocatable, dimension(:), intent(out) :: testsuite

      testsuite = [ &
        new_unittest("Test Cooling Interface", test_cooling) &
      ]

    end subroutine collect_cooling_suite

    subroutine test_cooling(error)
      type(error_type), allocatable, intent(out) :: error
    end subroutine test_cooling

end module CoolingTest
