module Farseal
  use iso_fortran_env, only: real32, int32
  implicit none

  private
  public CoolingMethods, CoolingType, ObjectiveType, AnnealerType

  real(kind=real32), parameter :: pi = 4.0_real32 * atan(1.0_real32)

  type :: CoolingMethod_Values
    integer :: LinMult = 1
    integer :: ExpMult = 2
    integer :: LogMult = 3
    integer :: QuadMult = 4
    integer :: LinAdd = 5
    integer :: QuadAdd = 6
    integer :: ExpAdd = 7
    integer :: TrigAdd = 8
    integer :: Custom = 9
  end type CoolingMethod_Values

  type(CoolingMethod_Values), parameter :: CoolingMethods = CoolingMethod_Values()

  type :: CoolingType
    integer :: method = CoolingMethods%LinMult
    real(kind=real32) :: tmax = 100
    real(kind=real32) :: tmin = 0
    real(kind=real32) :: alpha = 0.01
    integer(kind=int32) :: k = 1
    integer(kind=int32) :: n = 100
    real(kind=real32) :: temp = -1
    procedure(cooling_subroutine), pointer :: cool => CoolingMethod_cool
    logical, private :: initialized = .false.

    contains

      procedure :: init => CoolingMethod_init

  end type CoolingType

  interface
    subroutine cooling_subroutine(self)
      import CoolingType
      class(CoolingType), intent(inout) :: self
    end subroutine cooling_subroutine
  end interface

  type :: AnnealerType
  end type AnnealerType

  type :: ObjectiveType
    procedure(objective_function), pointer :: evaluate
  end type ObjectiveType

  abstract interface
    real(kind=real32) function objective_function(self, Annealer)
      import real32, ObjectiveType, AnnealerType
      class(ObjectiveType) :: self
      type(AnnealerType), allocatable :: Annealer
    end function objective_function
  end interface

  contains

    subroutine CoolingMethod_init(self)
      class(CoolingType), intent(out) :: self

      self%temp = self%tmax
      self%initialized = .true.

    end subroutine CoolingMethod_init

    subroutine CoolingMethod_cool(self)
      class(CoolingType), intent(inout) :: self

      if (.not. self%initialized) then
        call self%init()
      end if

      self%k = self%k + 1

      select case(self%method)
      case (1) ! LinMult
        self%temp = self%tmax / &
          (1.0_real32 + self%alpha * self%k)
      case (2) ! ExpMult
        self%temp = self%tmax * self%alpha ** self%k
      case (3) ! LogMult
        self%temp = self%tmax / &
          (1.0_real32 + self%alpha * log10(self%k + 1.0_real32))
      case (4) ! QuadMult
        self%temp = self%tmax / &
          (1.0_real32 + self%alpha * (self%k * 1.0_real32) ** 2)
      case (5) ! LinAdd
        self%temp = self%tmin + (self%tmax - self%tmin) * &
          (self%n * 1.0_real32 - self%k) / (self%n * 1.0_real32)
      case (6) ! QuadAdd
        self%temp = self%tmin + (self%tmax - self%tmin) * &
          ((self%n * 1.0_real32 - self%k) / (self%n * 1.0_real32)) ** 2
      case (7) ! ExpAdd
        self%temp = self%tmin + (self%tmax - self%tmin) / &
          (1.0_real32 + exp(2.0_real32 * log(self%tmax - self%tmin) * &
          (self%k - 0.5_real32 * self%n) / (self%n * 1.0_real32)))
      case (8) ! TrigAdd
        self%temp = self%tmin + 0.5_real32 * (self%tmax - self%tmin) &
          * (1.0_real32 + cos(self%k * pi / (self%n * 1.0_real32)))
      case (9) ! Custom
        stop 'The user MUST select a valid cooling option or give a custom cooling function'
      end select

    end subroutine CoolingMethod_cool

end module Farseal
