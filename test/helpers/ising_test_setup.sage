import random
from itertools import product

def J(size, nnzJ, seed=None):
    if seed is not None:
        random.seed(int(seed))

    IJ, JJ = zip(*random.sample(
        list(product(range(size), repeat=2)), nnzJ
    ))
    VJ = [random.uniform(-10, 10) for i in range(nnzJ)]
    return IJ, JJ, VJ

def H(size, nnzH, seed=None):
    if seed is not None:
        random.seed(int(seed))

    IH = random.sample(range(size), nnzH)
    VH = [random.uniform(-10, 10) for i in range(nnzH)]
    return IH, VH

def ising_hamiltonian(J, H, state):
    return -1 * J * state * state - H * state

def get_next_state(n, size):
    return vector(list(map(int, list(bin(n)[2:].zfill(size).replace("0", "3"))))) - vector([2] * size)

def solve(size, IJ, JJ, VJ, IH, VH):
    j_mat = matrix(
        {(ij, jj): vj for ij, jj, vj in zip(IJ, JJ, VJ)},
        nrows=size,
        ncols=size
    )

    h_mat = vector(RR, size, {ih: vh for ih, vh in zip(IH, VH)})

    best_energy, best_state = 1e8, []
    for i in range(2 ** size):
        state = get_next_state(i, size)
        energy = ising_hamiltonian(j_mat, h_mat, state)

        if energy < best_energy:
            best_energy = energy
            best_state = state.list()

    return best_energy, best_state

def main(seed=None):
    size, nnzJ, nnzH = 16, 100, 12
    IJ, JJ, VJ = J(size, nnzJ, seed)
    IH, VH = H(size, nnzH, seed)
    best_energy, best_state = solve(size, IJ, JJ, VJ, IH, VH)

    # +1 Because Fortran is 1-indexed
    print(f"Number of Spins: {size}")
    print(f"Number of Nonzero J Elements: {nnzJ}")
    print(f"Number of Nonzero H Elements: {nnzH}")
    print(f"Ith Index of J Matrix Nonzero Elements: {[v+1 for v in IJ]}")
    print(f"Jth Index of J Matrix Nonzero Elements: {[v+1 for v in JJ]}")
    print(f"Value of J Matrix Nonzero Elements: {VJ}")
    print(f"Ith Index of H Vector Nonzero Elements: {[v+1 for v in IH]}")
    print(f"Value of H Vector Nonzero Elements: {VH}")
    print(f"Optimal Ising Model Energy: {best_energy}")
    print(f"Optimal Ising Model State: {best_state}")

main(42)
