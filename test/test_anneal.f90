module Anneal
  use rsb
  use blas_sparse
  use iso_fortran_env, only: int32, real32, real64
  use testdrive, only: new_unittest, unittest_type, error_type, &
    check, test_failed, skip_test
  use Farseal, only: DiscreteAnnealType, ObjectiveType, CoolingMethods
  implicit none
private
  public :: collect_anneal_suite

  type, extends(ObjectiveType) :: IsingHamiltonianType
    integer :: J, H
    contains
      procedure :: evaluate => ising_hamiltonian
  end type IsingHamiltonianType

  contains

    subroutine collect_anneal_suite(testsuite)
      type(unittest_type), allocatable, dimension(:), intent(out) :: testsuite

      testsuite = [ &
        new_unittest("Test Discrete Anneal", test_discrete_anneal) &
      ]

    end subroutine collect_anneal_suite

    ! Test with Ising Model
    subroutine test_discrete_anneal(error)
      type(error_type), allocatable, intent(out) :: error
      type(DiscreteAnnealType), allocatable :: annealer
      type(IsingHamiltonianType), allocatable :: IsingHamiltonian
      integer :: n_spins, J, H, nnzJ, nnzH
      ! @@@ TODO: target was necessary for compilation to take place
      real(kind=real32), dimension(:), allocatable, target :: state
      type(c_ptr), parameter :: eo = c_null_ptr
      integer :: istat = 0
      integer, dimension(:), allocatable :: IJ, JJ, IH, JH
      real(kind=real32), dimension(:), allocatable :: VJ, VH
      real(kind=real32) :: energy_initial

      istat = rsb_lib_init(eo)
      if (istat .ne. 0) stop

      n_spins = 10
      allocate(state(n_spins))
      state = [-1, 1, -1, 1, -1, 1, -1, 1, -1, 1]

      annealer = DiscreteAnnealType(state_curr=state)
      annealer%max_step = 100
      annealer%alpha = 0.01
      annealer%t_max = 100.0
      annealer%t_min = 0.0
      annealer%cool_opt = CoolingMethods%QuadAdd
      annealer%mon_cool = .true.
      annealer%prog_bar = .true.
      annealer%resvar = 0.0
      !annealer%objective => IsingHamiltonian
      annealer%var_values = [1, -1]
      annealer%num_perturb = 1
      allocate(annealer%state_curr(n_spins))
      annealer%state_curr = state

      ! J matrix
      nnzJ = 3
      allocate(VJ(nnzJ))
      VJ = [10, 20, 30]
      allocate(IJ(nnzJ))
      IJ = [3, 2, 1]
      allocate(JJ(nnzJ))
      JJ = [2, 3, 1]

      call suscr_begin(n_spins, n_spins, J, istat)
      call ussp(J, blas_symmetric, istat)
      call suscr_insert_entries(J, nnzJ, VJ, IJ, JJ, istat)
      call uscr_end(J, istat)

      ! H vector
      nnzH = 5
      allocate(VH(nnzH))
      VH = [5, 4, 3, 2, 1]
      allocate(IH(nnzH))
      IH = [10, 8, 6, 2, 4]
      allocate(JH(nnzH))
      JH = [1, 1, 1, 1, 1]

      call suscr_begin(n_spins, 1, H, istat)
      call ussp(H, blas_general, istat)
      call suscr_insert_entries(H, nnzH, VH, IH, JH, istat)
      call uscr_end(J, istat)

      energy_initial = IsingHamiltonian%evaluate(annealer%state_curr)

      call annealer%optimize()

      ! TODO: Change Condition to Match the Best Hypothetical State Given IC
      if (.not. annealer%e_best < energy_initial) then
        call test_failed(error, "Blarg")
      end if

      call usds(J, istat)
      call usds(H, istat)
      istat = rsb_lib_exit(eo)

      if (allocated(error)) return

    end subroutine test_discrete_anneal

    !> Ising Model Hamiltonian (with Magnetic Moment, mu = 1)
    function ising_hamiltonian(self, state)
      class(IsingHamiltonianType), intent(inout) :: self
      real(kind=real32), dimension(:), intent(in) :: state
      real(kind=real32), dimension(:,:), allocatable :: state_mat, sigma_sigma
      real(kind=real32) :: ising_hamiltonian

      integer :: istat = 0
      integer :: trans = blas_no_trans
      integer :: incState = 1, inc_H_sigma = 1
      real(kind=real32) :: alpha = 1.0_real32
      real(kind=real32), dimension(:), allocatable :: H_sigma
      real(kind=real32), dimension(:,:), allocatable :: J_sigma_sigma

      state_mat = reshape(spread(state, 2, size(state)), [size(state), size(state)])

      sigma_sigma = state_mat * transpose(state_mat)
      allocate(J_sigma_sigma(size(state), size(state)))
      J_sigma_sigma(:,:) = 0
      ! @@@ TODO: Generic is usmm
      ! order(int), transA(int), nrhs(int), alpha(real), A(int), b(real(:,:)), ldb(int), c(real(:,:)), ldc(int), istat(int)
      call susmm( &
        order=blas_colmajor, &
        transA=trans, &
        nrhs=size(state), &
        alpha=alpha, &
        A=self%J, &
        b=sigma_sigma, &
        ldb=size(state), &
        c=J_sigma_sigma, &
        ldc=size(state), &
        istat=istat &
      )
      !call susmm(order=order, transA=trans, alpha=alpha, A=self%J, x=sigma_sigma, J_sigma_sigma, istat)
      J_sigma_sigma = sum(J_sigma_sigma * (-1))

      allocate(H_sigma(size(state)))
      H_sigma(:) = 0.0_real32
      call usmv( &
        transA=trans, &
        alpha=alpha, &
        A=self%H, &
        x=state, &
        incX=incState, &
        y=H_sigma, &
        incY=inc_H_sigma, &
        istat=istat &
      )

      ising_hamiltonian = sum(H_sigma)
      !ising_hamiltonian = J_sigma_sigma - H_sigma

    end function ising_hamiltonian

end module Anneal
