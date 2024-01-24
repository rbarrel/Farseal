module Anneal
  use rsb
  use blas_sparse
  use iso_c_binding
  use iso_fortran_env, only: int32, real32, real64
  use testdrive, only: new_unittest, unittest_type, error_type, &
    check, test_failed, skip_test
  use Farseal, only: CoolingMethods, ObjectiveType, DiscreteAnnealType
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
    ! Hardcoded Value Determined with ising_test_setup.sage
    subroutine test_discrete_anneal(error)
      type(error_type), allocatable, intent(out) :: error
      type(DiscreteAnnealType), allocatable :: annealer
      integer :: n_spins, J, H, nnzJ, nnzH
      integer, dimension(:), allocatable, target :: state
      type(c_ptr), parameter :: eo = c_null_ptr
      integer :: istat = 0
      integer, dimension(:), allocatable :: IJ, JJ, IH, JH
      real(kind=real32), dimension(:), allocatable :: VJ, VH
      real(kind=real32) :: initial_energy, expected_initial_energy

      istat = rsb_lib_init(eo)
      if (istat .ne. 0) stop

      n_spins = 16
      allocate(state(n_spins))
      state(:) = 1

      annealer = DiscreteAnnealType(state_curr=state)
      annealer%max_step = 100
      annealer%prog_bar = .false.
      annealer%resvar = 0.0
      annealer%var_values = [1, -1]
      annealer%num_perturb = 1
      annealer%cooler%mon_cool = .false.
      annealer%cooler%method = CoolingMethods%ExpMult
      allocate(annealer%state_curr(n_spins))
      annealer%state_curr = state
      annealer%objective%energy => ising_hamiltonian

      ! J matrix
      nnzJ = 100
      allocate(VJ(nnzJ))
      allocate(IJ(nnzJ))
      allocate(JJ(nnzJ))

      VJ = [ &
        9.790467012731906, 2.7999951970818575, 1.1389948754929247, &
        3.6922850197974917, 6.857038403796192, 5.519998230924896, &
        -5.419038560717913, -9.357995121919245, -3.690939038818362, &
        -4.645182480485945, -5.78034312827347, 8.858194286701089, &
        7.527352529453378, -3.706442384030442, 3.1087733058976, &
        -2.0873619787867153, 8.29095179481087, -0.8229629482520231, &
        -4.702396670038951, -5.067449846120331, 1.2273626832630153, &
        -4.745167829541295, 1.6917198044708108, 7.956457672049538, &
        -2.0119898971920547, -5.613584816854333, 9.950752129902206, &
        0.19052587352929073, -8.181811756524123, -9.057672491505308, &
        -7.8070173929868165, 2.5489208340618, 5.841587287259282, &
        -1.5568006640063192, -8.729445876960858, -2.3676142698692644, &
        9.922427604801936, 0.5822869019827408, 9.421567552272364, &
        7.215594044689961, -9.770379561143606, 4.414436387203892, &
        3.634207380531496, 0.7394066081759032, -4.663496200949144, &
        2.819235971596161, -7.768956528082471, -1.3046949866178998, &
        -0.9255258734158716, 9.076318550421604, 7.517058807563881, &
        -4.732218984978185, 0.011722261005965962, -6.426962389397373, &
        8.25255678689641, 7.410371396735339, -4.0311041710273425, &
        2.7789898973201037, 2.1794042287634454, -6.943214629007304, &
        5.250216001503025, 0.7875806023925147, 5.572529572611165, &
        0.6070734439035501, -9.98856207744113, -3.516878859906538, &
        -9.610465152283354, 8.581972325292341, 7.574437556463685, &
        6.6333105872235905, -3.8497174919467714, -8.84149667011625, &
        7.560191984080809, 8.938988905959881, -8.286930958642424, &
        -0.28019073366772496, -8.615749630632328, 5.212043305144631, &
        5.316688586139755, -7.432170710004744, -0.49435243802537343, &
        0.9960718698988771, -4.698867421198818, 7.448660821705147, &
        -1.5372411959822614, -5.764035891158359, 0.7859217755891663, &
        4.598621381799523, -5.9769787322060814, -3.7656741739821005, &
        9.902987133217891, 2.9975611527890695, -1.2379983217099184, &
        0.35151682071181156, -7.579916082634686, -5.506053259368853, &
        -3.238288757050894, 1.7661743691446663, -5.39770534806846, &
        -5.595652310968811 &
      ]

      IJ = [ &
        4, 1, 12, 5, 4, 16, 3, 12, 2, 11, 16, &
        15, 9, 2, 10, 7, 1, 1, 2, 4, 4, 9, &
        10, 16, 9, 4, 12, 11, 12, 16, 7, 4, 8, &
        10, 5, 13, 1, 13, 13, 3, 12, 16, 6, 14, &
        3, 15, 13, 6, 16, 15, 7, 2, 6, 6, 15, &
        5, 1, 12, 8, 9, 2, 7, 2, 9, 5, 11, &
        10, 6, 10, 4, 12, 2, 13, 11, 4, 5, 13, &
        15, 2, 13, 14, 8, 11, 6, 3, 6, 6, 4, &
        5, 11, 2, 10, 16, 3, 9, 16, 11, 8, 12, 5 &
      ]

      JJ = [ &
        10, 7, 14, 7, 15, 16, 4, 13, 11, 14, 14, &
        5, 12, 7, 8, 13, 9, 8, 8, 8, 12, 2, &
        11, 15, 16, 3, 8, 7, 4, 4, 12, 9, 3, &
        7, 8, 16, 2, 3, 15, 9, 3, 1, 8, 14, &
        8, 13, 4, 7, 8, 14, 2, 9, 12, 9, 10, &
        4, 12, 11, 6, 10, 16, 1, 5, 14, 12, 1, &
        15, 13, 4, 2, 5, 2, 8, 10, 11, 11, 2, &
        12, 10, 14, 5, 5, 3, 14, 10, 15, 11, 6, &
        5, 6, 3, 12, 7, 12, 9, 12, 12, 7, 1, 6 &
      ]

      call suscr_begin(n_spins, n_spins, J, istat)
      call ussp(J, blas_symmetric, istat)
      call suscr_insert_entries(J, nnzJ, VJ, IJ, JJ, istat)
      call uscr_end(J, istat)

      annealer%objective%J = J

      ! H vector
      nnzH = 12
      allocate(VH(nnzH))
      allocate(IH(nnzH))
      allocate(JH(nnzH))

      VH = [ &
        -9.404055611238594, -5.627240503927933, 0.10710576206724731, &
        -9.469280606322727, -6.02324698626703, 2.997688755590463, &
        0.8988296120643327, -5.591187559186066, 1.7853136775181753, &
        6.1886091335565325, -9.87002480643878, 6.116385036656158 &
      ]
      IH = [4, 1, 12, 5, 16, 14, 3, 2, 9, 13, 11, 15]
      JH(:) = 1

      call suscr_begin(n_spins, 1, H, istat)
      call ussp(H, blas_general, istat)
      call suscr_insert_entries(H, nnzH, VH, IH, JH, istat)
      call uscr_end(H, istat)

      annealer%objective%H = H

      expected_initial_energy = 27.17429127367062
      initial_energy = annealer%objective%energy(annealer%state_curr)

      if (abs(expected_initial_energy - initial_energy) > 1e-3) then
        call test_failed(error, "Ising Hamiltonian function error.")
      end if

      if (allocated(error)) return

      call annealer%optimize()

      if (annealer%e_best >= initial_energy) then
        call test_failed(error, "Discrete Annealer could not find a better state.")
      end if

      call usds(J, istat)
      call usds(H, istat)
      istat = rsb_lib_exit(eo)

      if (allocated(error)) return

    end subroutine test_discrete_anneal

    !> Ising Model Hamiltonian (with Magnetic Moment, mu = 1)
    function ising_hamiltonian(self, state)
      class(ObjectiveType), intent(inout) :: self
      integer, dimension(:), intent(in) :: state
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
        x=real(state, kind=real32), &
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
        x=real(state, kind=real32), &
        incX=incState, &
        y=H_sigma, &
        incY=inc_H_sigma, &
        istat=istat &
      )

      ising_hamiltonian = -1 * (J_sigma_sigma + sum(H_sigma))

    end function ising_hamiltonian

end module Anneal
