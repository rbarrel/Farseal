module Objective
  use iso_fortran_env, only: int32, real32
  use testdrive, only: new_unittest, unittest_type, error_type, &
    check, test_failed, skip_test
  use Farseal, only: ObjectiveType, AnnealerType, objective_interface
  implicit none

  private
  public :: collect_objective_suite

  type, extends(AnnealerType) :: TestAnnealerType
  end type TestAnnealerType

  type, extends(ObjectiveType) :: TestObjectiveType
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
      type(TestObjectiveType), allocatable :: Objective
      type(TestAnnealerType), allocatable :: Annealer
      real(kind=real32) :: expected_energy, actual_energy

      expected_energy = 1.0_real32

      Annealer = TestAnnealerType()
      Objective = ObjectiveType(evaluate)
      actual_energy = Objective%evaluate(Annealer)

      call check(error, expected_energy, actual_energy)
      if (allocated(error)) return

    end subroutine test_objective

    function evaluate(self, Annealer) result(energy)
      class(TestObjectiveType), intent(inout) :: self
      class(TestAnnealerType), allocatable :: Annealer
      real(kind=real32) :: energy, expected_energy

      energy = expected_energy

    end function evaluate

end module Objective
