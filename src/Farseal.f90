module Farseal
  use iso_fortran_env, only: real32, int32
  implicit none

  private

  real(kind=real32), parameter :: pi = 4.0_real32 * atan(1.0_real32)

  type, private :: CoolingMethod_Values
    integer :: LinMult
    integer :: ExpMult
    integer :: LogMult
    integer :: QuadMult
    integer :: LinAdd
    integer :: QuadAdd
    integer :: ExpAdd
    integer :: TrigAdd
  end type

  type(CoolingMethod_Values), public, parameter :: CoolingMethods = CoolingMethod_Values(1, 2, 3, 4, 5, 6, 7, 8)

  type, public :: CoolingType
    integer :: method = 1
    real(kind=real32) :: temp = 100
    real(kind=real32) :: tmin = 0
    real(kind=real32) :: tmax = 100
    real(kind=real32) :: alpha = 0.01
    integer(kind=int32) :: k = 1
    integer(kind=int32) :: n = 100
    contains
      procedure :: cool

  end type

  contains

    subroutine cool(self)
      class(CoolingType), intent(inout) :: self

      self%k = self%k + 1

      select case(self%method)
      case (1) ! LinMult
        self%temp = self%tmax / (1.0_real32 + self%alpha * self%k)
      case (2) ! ExpMult
        self%temp = self%tmax * self%alpha ** self%k
      case (3) ! LogMult
        self%temp = self%tmax / (1.0_real32 + self%alpha * log10(self%k + 1.0_real32))
      case (4) ! QuadMult
        self%temp = self%tmax / (1.0_real32 + self%alpha * (self%k * 1.0_real32) ** 2)
      case (5) ! LinAdd
        self%temp = self%tmin + (self%tmax - self%tmin) * (self%n * 1.0_real32 - self%k) / (self%n * 1.0_real32)
      case (6) ! QuadAdd
        self%temp = self%tmin + (self%tmax - self%tmin) * ((self%n * 1.0_real32 - self%k) / (self%n * 1.0_real32)) ** 2
      case (7) ! ExpAdd
        self%temp = self%tmin + (self%tmax - self%tmin) / &
          (1.0_real32 + exp(2.0_real32 * log(self%tmax - self%tmin) * (self%k - 0.5_real32 * self%n) / (self%n * 1.0_real32)))
      case (8) ! TrigAdd
        self%temp = self%tmin + 0.5_real32 * (self%tmax - self%tmin) * (1.0_real32 + cos(self%k * pi / (self%n * 1.0_real32)))
      end select

    end subroutine cool

end module Farseal
