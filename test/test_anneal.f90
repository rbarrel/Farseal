module Anneal
  use rsb
  use blas_sparse
  use iso_c_binding
  use iso_fortran_env, only: int32, real32, real64
  use testdrive, only: new_unittest, unittest_type, error_type, &
    check, test_failed, skip_test
  use Farseal, only: DiscreteAnnealType, ObjectiveType, CoolingMethods
  implicit none
private
  public :: collect_anneal_suite

  type, extends(ObjectiveType) :: IsingHamiltonianType
    integer :: J = -1, H = -1
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
      !integer(kind=rsb_idx_kind) :: res
      !integer(kind=c_int), target :: ione = 1
      integer, dimension(:), allocatable :: IJ, JJ, IH, JH
      real(kind=real32), dimension(:), allocatable :: VJ, VH
      real(kind=real32) :: energy_initial

      istat = rsb_lib_init(eo)
      if (istat .ne. 0) stop
      !res = rsb_lib_set_opt(rsb_io_want_verbose_tuning,c_loc(ione))

      n_spins = 100
      allocate(state(n_spins))
      state = [ &
        1, 1, 1, 1, 1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, &
        -1, -1, 1, 1, 1, -1, 1, -1, -1, 1, 1, 1, 1, 1, 1, -1, -1, -1, &
        -1, -1, 1, 1, -1, 1, -1, 1, -1, -1, -1, 1, -1, 1, 1, -1, -1, 1, &
        1, -1, -1, -1, 1, 1, -1, -1, -1, 1, 1, 1, 1, -1, 1, 1, 1, -1, &
        1, -1, -1, -1, 1, 1, -1, 1, -1, 1, 1, -1, 1, 1, -1, 1, 1, 1, 1, &
        -1, -1, -1, 1, 1, 1, -1, -1, 1, 1, -1, -1 &
      ]
      !state(:) = 1

      IsingHamiltonian = IsingHamiltonianType()

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
      nnzJ = 10
      allocate(VJ(nnzJ))
      VJ = [ &
        -1.561563606294591, &
        -9.404055611238594, &
        -5.627240503927933, &
        0.10710576206724731, &
        -9.469280606322727, &
        -6.02324698626703, &
        2.997688755590463, &
        0.8988296120643327, &
        -5.591187559186066, &
        1.7853136775181753 &
      ]
      allocate(IJ(nnzJ))
      !IJ = [18, 4, 45, 40, 36, 22, 16, 89, 14, 96]
      IJ = [19, 5, 46, 41, 37, 23, 17, 90, 15, 97]
      allocate(JJ(nnzJ))
      !JJ = [24, 9, 6, 12, 57, 86, 79, 35, 24, 74]
      JJ = [25, 10, 7, 13, 58, 87, 80, 36, 25, 75]

      call suscr_begin(n_spins, n_spins, J, istat)
      call ussp(J, blas_symmetric, istat)
      call ussp(J, blas_one_base, istat)
      call ussp(J, blas_rowmajor, istat)
      call suscr_insert_entries(J, nnzJ, VJ, IJ, JJ, istat)
      call uscr_end(J, istat)

      IsingHamiltonian%J = J

      ! H vector
      nnzH = 10
      allocate(VH(nnzH))
      VH = [ &
        -8.145083132397042, &
        -8.06567246333072, &
        6.949887326949195, &
        2.0745206273378223, &
        6.142565465487603, &
        4.594635733876357, &
        0.7245618290940143, &
        9.462315279587411, &
        -2.42931245583293, &
        1.0408126254645396 &
      ]
      allocate(IH(nnzH))
      !IH = [0, 97, 20, 89, 54, 43, 35, 19, 27, 13]
      IH = [1, 98, 21, 90, 55, 44, 36, 20, 28, 14]
      allocate(JH(nnzH))
      JH(:) = 1

      call suscr_begin(n_spins, 1, H, istat)
      call ussp(H, blas_general, istat)
      call ussp(H, blas_one_base, istat)
      call ussp(H, blas_rowmajor, istat)
      call suscr_insert_entries(H, nnzH, VH, IH, JH, istat)
      call uscr_end(H, istat)

      IsingHamiltonian%H = H

      energy_initial = IsingHamiltonian%evaluate(annealer%state_curr)
      !print *, energy_initial

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
      real(kind=real32) :: ising_hamiltonian

      integer :: istat = 0
      integer :: incState = 1, inc_H_sigma = 1, inc_J_sigma = 1
      real(kind=real32) :: alpha = 1.0_real32
      real(kind=real32) :: J_sigma_sigma
      real(kind=real32), dimension(:), allocatable :: H_sigma, J_sigma

      allocate(J_sigma(size(state)))
      J_sigma(:) = 0.0_real32
      call usmv( &
        transA=blas_no_trans, &
        alpha=alpha, &
        A=self%J, &
        x=state, &
        incX=incState, &
        y=J_sigma, &
        incY=inc_J_sigma, &
        istat=istat &
      )

      J_sigma_sigma = dot_product(J_sigma, state)

      allocate(H_sigma(size(state)))
      H_sigma(:) = 0.0_real32
      call usmv( &
        transA=blas_trans, &
        alpha=alpha, &
        A=self%H, &
        x=state, &
        incX=incState, &
        y=H_sigma, &
        incY=inc_H_sigma, &
        istat=istat &
      )

      ising_hamiltonian = -1 * (J_sigma_sigma + sum(H_sigma))
      print *, ising_hamiltonian

    end function ising_hamiltonian

end module Anneal
