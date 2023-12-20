module Objective
  use iso_fortran_env, only: int32, real32
  use testdrive, only: new_unittest, unittest_type, error_type, &
    check, test_failed, skip_test
  use Farseal, only: ObjectiveType, AnnealerType
  implicit none

  private
  public :: collect_objective_suite

  contains

    subroutine collect_objective_suite(testsuite)
      type(unittest_type), allocatable, dimension(:), intent(out) :: testsuite

      testsuite = [ &
        new_unittest("Test Objective Function Interface", test_objective) &
      ]

    end subroutine collect_objective_suite

    subroutine test_objective(error)
      type(error_type), allocatable, intent(out) :: error
      type(ObjectiveType), allocatable :: Objective
      type(AnnealerType), allocatable :: Annealer
      real(kind=real32) :: expected_energy, actual_energy

      actual_energy = 1.0_real32

      Annealer = AnnealerType()
      Objective = ObjectiveType(FooBar)
      actual_energy = Objective%evaluate(Annealer)

      call check(error, expected_energy, actual_energy)
      if (allocated(error)) return

      contains

        function FooBar(Annealer) result(energy)
          type(AnnealerType), intent(in) :: Annealer
          real(kind=real32) :: energy

          energy = actual_energy

        end function FooBar

    end subroutine test_objective

end module Objective
