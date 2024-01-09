module Objective
  use iso_fortran_env, only: int32, real32
  use testdrive, only: new_unittest, unittest_type, error_type, &
    check, test_failed, skip_test
  use Farseal, only: ObjectiveType, AnnealerType
  implicit none

  private
  public :: collect_objective_suite

  type, extends(ObjectiveType) :: TestObjectiveType
    real(kind=real32) :: expected_output = 10_real32
    contains
      procedure :: evaluate
  end type TestObjectiveType

  contains

    subroutine collect_objective_suite(testsuite)
      type(unittest_type), allocatable, dimension(:), intent(out) :: testsuite

      testsuite = [ &
        new_unittest("Test Objective Function Interface", test_objective) &
      ]

    end subroutine collect_objective_suite

    subroutine test_objective(error)
      type(error_type), allocatable, intent(out) :: error
      type(TestObjectiveType), allocatable :: objective
      integer(kind=real32), dimension(:), allocatable :: state
      real(kind=real32) :: expected_energy, actual_energy

      expected_energy = 10_real32
      allocate(state(1))
      state(:) = 0

      objective = TestObjectiveType()
      actual_energy = objective%evaluate(state)

      call check(error, expected_energy, actual_energy)
      if (allocated(error)) return

    end subroutine test_objective

    function evaluate(self, state) result(energy)
      class(TestObjectiveType), intent(inout) :: self
      integer(kind=real32), dimension(:), intent(in) :: state
      real(kind=real32) :: energy

      energy = self%expected_output

    end function evaluate

end module Objective
