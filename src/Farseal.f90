module Farseal
  use iso_fortran_env, only: real32, int32
  implicit none

  private
  public CoolingMethods, DiscreteAnnealType, ObjectiveType

  real(kind=real32), parameter :: pi = 4.0_real32 * atan(1.0_real32)


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
    real(kind=real32) :: tmax = 100
    real(kind=real32) :: tmin = 0
    real(kind=real32) :: alpha = 0.01
    integer(kind=int32) :: k = 1
    integer(kind=int32) :: n = 100
    real(kind=real32) :: temp = -1
    logical :: mon_cool = .true.
    logical, private :: initialized = .false.
    contains
      procedure, pass(self) :: cool => CoolingMethod_cool
      procedure, pass(self) :: init => CoolingMethod_init
  end type CoolingType

  interface
    subroutine cooling_subroutine(self)
      import CoolingType
      class(CoolingType), intent(inout) :: self
    end subroutine cooling_subroutine
  end interface


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! Objective Function Types !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  type, abstract :: ObjectiveType
    contains
      procedure(objective_interface), deferred, pass(self) :: energy
  end type ObjectiveType

  interface
    function objective_interface(self, state) result(energy)
      use iso_fortran_env, only: real32
      import ObjectiveType
      class(ObjectiveType), intent(inout) :: self
      real(kind=real32), dimension(:), intent(in) :: state
      real(kind=real32) :: energy
    end function objective_interface
  end interface


  !!!!!!!!!!!!!!!!!!!!!!!
  !!! Annealing Types !!!
  !!!!!!!!!!!!!!!!!!!!!!!

  type, abstract :: AnnealerType
    integer :: max_step = 100
    integer :: total_steps = 0
    type(CoolingType) :: cooler = CoolingType()
    real(kind=real32) :: e_best = 1.0e30_real32
    logical :: prog_bar = .false.
    real(kind=real32) :: resvar = 0.0_real32
    contains
      procedure, pass(self) :: optimize
  end type AnnealerType

  type, abstract, extends(AnnealerType) :: DiscreteAnnealType
    real(kind=real32), pointer, dimension(:) :: state_curr
    integer, dimension(:), allocatable :: state_neigh, state_best, var_values
    integer :: num_perturb = 0
    contains
      procedure, pass(self) :: get_neigh => get_neigh_disc
  end type DiscreteAnnealType


  !!!!!!!!!!!!!!!!!!
  !!! Procedures !!!
  !!!!!!!!!!!!!!!!!!

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

    subroutine optimize(self)
      class(AnnealerType),intent(inout) :: self
!      integer :: step
!      real(kind=real32) :: e_neigh, temp_r, e_curr, t_curr
!
!      e_neigh = 1000_real32
!      e_curr = 1000_real32
!
!      ! @@@ TODO: Set Cooling
!
!      ! allocate the neighbor state variables and set the cooling, also set initial energy
!      self%size_states=size(self%state_curr)
!      if (.not. allocated(self%state_neigh)) then
!        allocate(self%state_neigh(self%size_states))
!      end if
!      if (.not. allocated(self%state_best)) then
!        allocate(self%state_best(self%size_states))
!      end if
!      ! if the var_values have not been given, then traverse the initial state and pick out each
!      ! unique value and assume those are the possible values
!      if (.not. allocated(self%var_values)) then
!        CALL init_var_values(self)
!      end if
!      self%state_neigh=self%state_curr
!      self%state_best=self%state_curr
!      !set energy to current energy
!      e_curr=self%energy(self%state_curr)
!      self%e_best=e_curr
!
!      t_curr=self%t_max
!      step=0
!      self%total_steps=0
!      if (self%prog_bar) write(*, '(A)', advance='NO') 'PROGRESS:'
!      ! actual simulated annealing happens here
!      do while(step .LT. self%max_step .AND. t_curr .gt. self%t_min)
!        self%total_steps=self%total_steps+1
!        ! get a new neighbor and compute energy
!        self%state_neigh=self%get_neigh(self%state_curr)
!        e_neigh=self%energy(self%state_neigh)
!        ! check and see if we accept the new temperature (lower temps always accepted)
!        CALL random_number(temp_r)
!        if (temp_r .LE. accept_prob(e_curr,e_neigh,t_curr)) then
!          ! if we accept then it always counts as a new step
!          step=step+1
!          if (self%prog_bar .AND. mod(step,nint(self%max_step/91_real32)) .eq. 0) &
!            write(*, '(A)', advance='NO') '*'
!          self%state_curr=self%state_neigh
!          e_curr=e_neigh
!        ELSE
!          ! otherwise, it has a 1% chance to count as a new step to finish the problem
!          ! especially important for combinatorials
!          CALL random_number(temp_r)
!          if (temp_r .LE. 0.01D0) then
!            step=step+1
!            if (self%prog_bar .AND. mod(step,nint(self%max_step/91_real32)) .eq. 0) &
!              write(*, '(A)', advance='NO') '*'
!          end if
!        end if
!        ! cool the temperature
!        t_curr=self%cool(self%t_min,self%t_max,self%alpha,step,self%max_step)
!        ! if it is the best energy, it's our new best value
!        if (e_curr .LT. self%e_best) then
!          self%e_best=e_curr
!          self%state_best=self%state_neigh
!        end if
!        if (.not. self%mon_cool)t_curr=t_curr*(1.0_real32+(e_curr-self%e_best)/e_curr)
!        ! rewind to best value
!        if (abs(t_curr) .LE. self%resvar) then
!          self%resvar=self%resvar/2.0_real32
!          e_curr=self%e_best
!          self%state_curr=self%state_best
!        end if
!      end do
!      if (self%prog_bar) write(*, '(A)') '*'
!
!      ! set to the best state we ended up finding.
!      self%state_curr=self%state_best
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

end module Farseal
