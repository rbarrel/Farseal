module CoolingTest
  use iso_fortran_env, only: int32, real32
  use testdrive, only: new_unittest, unittest_type, error_type, &
    check, test_failed, skip_test
  use Farseal, only: CoolingType
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
      type(CoolingType), allocatable :: Cooler
      real(kind=real32) :: starting_temp, ending_temp

      starting_temp = Cooler%temp
      call Cooler%cool()
      ending_temp = Cooler%temp

      call check(ending_temp < starting_temp)
      if (allocated(error)) return

    end subroutine test_cooling

end module CoolingTest
