module Cooling
  use iso_fortran_env, only: int32, real32
  use testdrive, only: new_unittest, unittest_type, error_type, &
    check, test_failed, skip_test
  use Farseal, only: CoolingType, CoolingMethods
  implicit none

  private
  public :: collect_cooling_suite

  contains

    subroutine collect_cooling_suite(testsuite)
      type(unittest_type), allocatable, dimension(:), intent(out) :: testsuite

      testsuite = [ &
        new_unittest("Test Cooling Interface", test_cooling), &
        new_unittest("Test Custom Cooling", test_custom_cooling) &
      ]

    end subroutine collect_cooling_suite

    subroutine test_cooling(error)
      type(error_type), allocatable, intent(out) :: error
      type(CoolingType), allocatable :: Cooler
      real(kind=real32) :: starting_temp, ending_temp

      Cooler = CoolingType()
      call Cooler%init()
      Cooler%method = CoolingMethods%ExpAdd
      starting_temp = Cooler%temp
      call Cooler%cool()
      ending_temp = Cooler%temp

      call check(error, ending_temp < starting_temp)
      if (allocated(error)) return

    end subroutine test_cooling

    subroutine test_custom_cooling(error)
      type(error_type), allocatable, intent(out) :: error
      type(CoolingType), allocatable :: Cooler
      real(kind=real32) :: starting_temp, ending_temp

      Cooler = CoolingType()
      call Cooler%init()
      Cooler%method = CoolingMethods%Custom
      Cooler%cool => FooBar
      starting_temp = Cooler%temp
      call Cooler%cool()
      ending_temp = Cool%temp

      call check(error, ending_temp < starting_temp)
      if (allocated(error)) return

      contains

        subroutine FooBar(self)
          class(CoolingType) :: self

          self%temp = self%temp - 1

        end subroutine FooBar

    end subroutine test_custom_cooling

end module Cooling
