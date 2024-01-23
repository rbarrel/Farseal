module Farseal
  use iso_fortran_env, only: real32, int32
  implicit none

  private
  public CoolingType, CoolingMethods, DiscreteAnnealType, ObjectiveType, AnnealerType

  real(kind=real32), parameter :: pi = 4.0_real32 * atan(1.0_real32)
  real(kind=real32), parameter :: prog_len = 91.0_real32

  !!!!!!!!!!!!!!!!!!!!!
  !!! Cooling Types !!!
  !!!!!!!!!!!!!!!!!!!!!

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
    integer :: method = CoolingMethods%LinAdd
    real(kind=real32) :: t_max = 100
    real(kind=real32) :: t_min = 0
    real(kind=real32) :: alpha = 0.01
    integer(kind=int32) :: k = 1
    integer(kind=int32) :: n = 100
    real(kind=real32) :: temp = -1
    logical :: mon_cool = .true.
    logical, private :: initialized = .false.
    procedure(cooling_subroutine), pointer :: cool => null()
    contains
      procedure, pass(self) :: init => CoolingMethod_init
  end type CoolingType

  interface
    subroutine cooling_subroutine(self)
      import CoolingType
      class(CoolingType), intent(inout) :: self
    end subroutine cooling_subroutine
  end interface


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! Objective Function Type !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  type :: ObjectiveType
    integer :: A = -1
    integer :: B = -1
    integer :: C = -1
    integer :: D = -1
    integer :: E = -1
    integer :: F = -1
    integer :: G = -1
    integer :: H = -1
    integer :: I = -1
    integer :: J = -1
    integer :: K = -1
    integer :: L = -1
    integer :: M = -1
    integer :: N = -1
    integer :: O = -1
    integer :: P = -1
    integer :: Q = -1
    integer :: R = -1
    integer :: S = -1
    integer :: T = -1
    integer :: U = -1
    integer :: V = -1
    integer :: W = -1
    integer :: X = -1
    integer :: Y = -1
    integer :: Z = -1
    procedure(objective_interface), pointer :: energy => null()
  end type ObjectiveType

  interface
    function objective_interface(self, state) result(energy)
      use iso_fortran_env, only: real32
      import ObjectiveType
      class(ObjectiveType), intent(inout) :: self
      integer, dimension(:), intent(in) :: state
      real(kind=real32) :: energy
    end function objective_interface
  end interface


  !!!!!!!!!!!!!!!!!!!!!!!
  !!! Annealing Types !!!
  !!!!!!!!!!!!!!!!!!!!!!!

  type :: AnnealerType
    integer :: size_states = 0
    real(kind=real32) :: max_step = 100.0_real32
    integer :: total_steps = 0
    type(CoolingType) :: cooler = CoolingType()
    real(kind=real32) :: e_best = 1.0e30_real32
    logical :: prog_bar = .false.
    real(kind=real32) :: resvar = 0.0_real32
    type(ObjectiveType) :: objective = ObjectiveType()
    ! @@@ TODO: Create a flag for verbose logging in here
    contains
      procedure, pass(self) :: optimize
  end type AnnealerType

  type, extends(AnnealerType) :: DiscreteAnnealType
    integer, dimension(:), pointer :: state_curr
    integer, dimension(:), allocatable :: state_neigh, state_best, var_values
    integer :: num_perturb = 0
    contains
      procedure, pass(self) :: get_neigh => get_neigh_disc
  end type DiscreteAnnealType

  ! @@@ TODO: ContinuousAnnealType
  ! @@@ TODO: CombinatorialAnnealType

  !!!!!!!!!!!!!!!!!!
  !!! Procedures !!!
  !!!!!!!!!!!!!!!!!!

  contains

    subroutine CoolingMethod_init(self)
      class(CoolingType), intent(out) :: self

      if (.not. self%initialized) then
        self%initialized = .true.
        self%temp = self%t_max
        self%cool => CoolingMethod_cool
      end if

    end subroutine CoolingMethod_init

    subroutine CoolingMethod_cool(self)
      class(CoolingType), intent(inout) :: self

      call self%init()

      self%k = self%k + 1

      select case(self%method)
      case (1) ! LinMult
        self%temp = self%t_max / &
          (1.0_real32 + self%alpha * self%k)
      case (2) ! ExpMult
        self%temp = self%t_max * self%alpha ** self%k
      case (3) ! LogMult
        self%temp = self%t_max / &
          (1.0_real32 + self%alpha * log10(self%k + 1.0_real32))
      case (4) ! QuadMult
        self%temp = self%t_max / &
          (1.0_real32 + self%alpha * (self%k * 1.0_real32) ** 2)
      case (5) ! LinAdd
        self%temp = self%t_min + (self%t_max - self%t_min) * &
          (self%n * 1.0_real32 - self%k) / (self%n * 1.0_real32)
      case (6) ! QuadAdd
        self%temp = self%t_min + (self%t_max - self%t_min) * &
          ((self%n * 1.0_real32 - self%k) / (self%n * 1.0_real32)) ** 2
      case (7) ! ExpAdd
        self%temp = self%t_min + (self%t_max - self%t_min) / &
          (1.0_real32 + exp(2.0_real32 * log(self%t_max - self%t_min) * &
          (self%k - 0.5_real32 * self%n) / (self%n * 1.0_real32)))
      case (8) ! TrigAdd
        self%temp = self%t_min + 0.5_real32 * (self%t_max - self%t_min) &
          * (1.0_real32 + cos(self%k * pi / (self%n * 1.0_real32)))
      case (9) ! Custom
        stop 'The user MUST select a valid cooling option or give a custom cooling function'
      end select

    end subroutine CoolingMethod_cool

    subroutine optimize(self)
      class(AnnealerType),intent(inout) :: self
      integer :: step
      real(kind=real32) :: e_neigh, temp_r, e_curr, t_curr

      e_neigh = 1000_real32
      e_curr = 1000_real32

      call self%cooler%init()

      select type(self)
        type is(DiscreteAnnealType)

          self%size_states=size(self%state_curr)

          if (.not. allocated(self%state_neigh)) then
            allocate(self%state_neigh(self%size_states))
          end if

          if (.not. allocated(self%state_best)) then
            allocate(self%state_best(self%size_states))
          end if

          if (.not. allocated(self%var_values)) then
            call init_var_values(self)
          end if

          self%state_neigh = self%state_curr
          self%state_best = self%state_curr

          e_curr = self%objective%energy(self%state_curr)
          self%e_best = e_curr

      end select

      t_curr = self%cooler%t_max
      step = 0
      self%total_steps = 0

      if (self%prog_bar) write(*, '(A)', advance='NO') 'PROGRESS:'

      ! Start Annealing
      do while(step < self%max_step .AND. t_curr > self%cooler%t_min)
        self%total_steps = self%total_steps + 1

        select type(self)
          type is(DiscreteAnnealType)
            self%state_neigh = self%get_neigh(self%state_curr)
            e_neigh = self%objective%energy(self%state_neigh)
        end select

        ! check and see if we accept the new temperature (lower temps always accepted)
        call random_number(temp_r)
        if (temp_r <= accept_prob(e_curr, e_neigh, t_curr)) then
          step = step + 1

          if (self%prog_bar .and. mod(step, nint(self%max_step/prog_len)) == 0) &
            write(*, '(A)', advance='NO') '*'

          select type(self)
            type is(DiscreteAnnealType)
              self%state_curr = self%state_neigh
          end select

          e_curr = e_neigh

        else
          ! 1% chance to accept a higher temp (critical for combinatorials)
          call random_number(temp_r)
          if (temp_r <= 0.01_real32) then
            step = step + 1

            if (self%prog_bar .and. mod(step, nint(self%max_step/prog_len)) == 0) &
              write(*, '(A)', advance='NO') '*'
          end if

        end if

        call self%cooler%cool()
        t_curr = self%cooler%temp

        if (e_curr < self%e_best) then
          self%e_best=e_curr

          select type(self)
            type is(DiscreteAnnealType)
              self%state_best=self%state_neigh
          end select
        end if

        if (.not. self%cooler%mon_cool) then
          t_curr=t_curr*(1.0_real32+(e_curr-self%e_best)/e_curr)
        end if

        if (abs(t_curr) .LE. self%resvar) then
          self%resvar=self%resvar/2.0_real32

          select type(self)
            type is(DiscreteAnnealType)
              e_curr=self%e_best
              self%state_curr=self%state_best
          end select
        end if

      end do

      if (self%prog_bar) write(*, '(A)') '*'

      select type(self)
        type is(DiscreteAnnealType)
          self%state_curr=self%state_best
      end select

    end subroutine optimize

    function get_neigh_disc(self, s_curr)
      class(DiscreteAnnealType), intent(inout) :: self
      integer, dimension(:), intent(in) :: s_curr
      integer, dimension(size(s_curr)) :: get_neigh_disc

!      real(kind=real32) :: temp_r
!      integer :: i, j, val, num_perturb, temp_i
!      real(kind=real32), dimension(size(self%var_values)) :: rand
!      integer, dimension(:), allocatable :: perturb_locs
!
!      ! this line is literally just to ensure that it doesn't complain about not using the variables
!      if (.false.) temp_r = self%e_best
!
!      ! perturb at least one parameter
!      num_perturb = max(1, self%num_perturb)
!
!      ! find which locations we are perturbing
!      allocate(perturb_locs(num_perturb))
!      perturb_locs = 0
!      i = 1
!      do while(i .le. num_perturb)
!        ! get a random index for the parameters
!        call random_number(temp_r)
!        temp_i = 1 + floor(temp_r * size(s_curr))
!        j = i
!        do j = 1, i - 1
!          if (perturb_locs(j) .eq. temp_i) exit
!        end do
!        ! if we exited early then we already found the index so we don't use that index and repeat
!        if (j .EQ. i) then
!          perturb_locs(j)=temp_i
!          i=i+1
!        end if
!      end do
!
!      ! now actually perturb the system at each location
      get_neigh_disc = s_curr
!      do i=1,num_perturb
!        j=perturb_locs(i)
!        ! pick new discrete value
!        call random_number(rand)
!        val=self%var_values(maxloc(rand,dim=1))
!        ! make sure that a new value is chosen
!        if (val .EQ. s_curr(j)) then
!          val=self%var_values(minloc(rand,dim=1))
!        end if
!        ! perturb the state
!        get_neigh_disc(j)=val
!      end do
!
    end function get_neigh_disc

  function accept_prob(e_current,e_neigh,t_current)
    real(kind=real32),intent(in) :: e_current,e_neigh,t_current
    real(kind=real32) :: accept_prob
    real(kind=real32) :: delta_e

    delta_e=e_neigh-e_current
    if(-delta_e/t_current .le. -700.0_real32)then
      accept_prob=0.0_real32
    elseif(-delta_e/t_current .GE. 700.0_real32)then
      accept_prob=10.0_real32
    else
      accept_prob=exp(-delta_e/t_current)
    endif
    if(delta_e .LE. 0.0_real32)accept_prob=10.0_real32
    if(ISNAN(accept_prob))accept_prob=0.0_real32
  end function accept_prob

  subroutine init_var_values(self)
    class(DiscreteAnnealType),intent(inout) :: self
    integer :: var_vals(self%size_states)
    integer :: i,num_unique

    if(any(var_vals==-999999999))stop 'if you are going to use -999999999 as a value, specify the &
        &var_values!'
    !traverse the initial state to see what the unique values are. uses -999999999
    !to indicate values not yet found
    var_vals=-999999999
    num_unique=0
    do i=1,self%size_states
      !check to see if we've already recorded this value, if not then add it to the list!
      if(.not. any(var_vals==self%state_curr(i)))then
        num_unique=num_unique+1
        var_vals(num_unique)=self%state_curr(i)
      endif
    enddo
    !num_unique is the number of unique values
    allocate(self%var_values(num_unique))
    do i=1,num_unique
      self%var_values(i)=var_vals(i)
    enddo
  endsubroutine

end module Farseal
