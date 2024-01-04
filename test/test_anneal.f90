module Anneal
  use rsb
  use blas_sparse
  use iso_fortran_env, only: int32, real32
  use testdrive, only: new_unittest, unittest_type, error_type, &
    check, test_failed, skip_test
  use Farseal, only: DiscreteAnnealType
  implicit none

  private
  public :: collect_anneal_suite

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
      integer :: n_spins, J, H, nnzJ, nnzH
      integer, dimension(:), allocatable :: state
      type(c_ptr), parameter :: eo = c_null_ptr
      integer :: istat = 0
      real(kind=real32), dimension(:), allocatable :: VJ, IJ, JH, VH, IH
      real(kind=real32) :: energy_initial

      istat = rsb_lib_init(eo)
      if (istat .ne. 0) stop

      n_spins = 10
      allocate(state(n_spins))
      state = [-1, 1, -1, 1, -1, 1, -1, 1, -1, 1]

      annealer%max_step = 100
      annealer%alpha = 0.01
      annealer%t_max = 100.0
      annealer%t_min = 0.0
      annealer%cool_opt = "QuadAdd"
      annealer%mon_cool = .true.
      annealer%prog_bar = .true.
      annealer%resvar = 0.0
      annealer%energy => ising_hamiltonian
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
      call uscr_insert_entries(J, nnzJ, VJ, IJ, JJ, istat)
      call uscr_end(J, istat)

      ! H vector
      nnzH = 5
      allocate(VH(nnzH))
      VH = [5, 4, 3, 2, 1]
      allocate(IH(nnzH))
      IH = [10, 8, 6, 2, 4]

      call suscr_begin(n_spins, 1, H, istat)
      call uscr_insert_entries(H, nnzH, VH, IH, JH, istat)
      call uscr_end(J, istat)

      energy_initial = ising_hamiltonian(annealer, annealer%state_curr)

      call annealer%optimize()

      ! TODO: Change Condition to Match the Best Hypothetical State Given IC
      if (.not. annealer%e_best < energy_initial) then
        call test_failed(error, "")
      end if

      call usds(J, istat)
      call usds(H, istat)
      istat = rsb_lib_exit(eo)

      if (allocated(error)) return

    end subroutine test_discrete_anneal

    !> Ising Model Hamiltonian (with Magnetic Moment, mu = 1)
    function ising_hamiltonian(sa, state)
      class(DiscreteAnnealType), intent(inout) :: sa
      integer, dimension(:), intent(in) :: state
      integer, dimension(:,:), allocatable :: state_mat, sigma_sigma
      real(8) :: ising_hamiltonian

      integer :: istat = 0
      integer :: trans = blas_no_trans
      integer :: incB = 1, inc_H_sigma = 1
      real(kind=real32) :: alpha = 1.0_real32
      real(kind=real32), dimension(:), allocatable :: H_sigma
      real(kind=real32), dimension(:,:), allocatable :: J_sigma_sigma

      state_mat = reshape(spread(state, 2, size(state)), [size(state), size(state)])

      ! TODO J_sigma_sigma looks too memory intensive (maybe?)
      ! What is the output of J * J_sigma_sigma?
      sigma_sigma = state_mat * transpose(state_mat)
      allocate(J_sigma_sigma(size(state), size(state)))
      J_sigma_sigma(:,:) = 0
      call usmm(J, sigma_sigma, J_sigma_sigma, istat)
      J_sigma_sigma = sum(J_sigma_sigma * -1)

      allocate(H_sigma(size(state)))
      H_sigma(:) = 0
      call usmv(trans, alpha, H, state, incB, H_sigma, inc_H_sigma, istat)

      ising_hamiltonian = J_sigma_sigma - H_sigma

    end function ising_hamiltonian

end module Anneal
