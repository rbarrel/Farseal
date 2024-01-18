import random
from pathlib import Path
from itertools import product

def state(size):
    return random.choices([-1, 1], k=size)

def J(size, nnzJ):
    IJ, JJ = zip(*random.sample(
        list(product(range(size), range(size))), nnzJ
    ))
    VJ = [random.uniform(-10, 10) for i in range(nnzJ)]
    return IJ, JJ, VJ

def H(size, nnzH):
    IH = random.sample(range(size), nnzH)
    VH = [random.uniform(-10, 10) for i in range(nnzH)]
    return IH, VH

def main(seed=42):
    random.seed(int(seed))
    size, nnzJ, nnzH = 100, 10, 10

    IJ, JJ, VJ = J(size, nnzJ)
    j_mat = matrix(
        {(ij, jj): vj for ij, jj, vj in zip(IJ, JJ, VJ)},
        nrows=size,
        ncols=size
    )

    IH, VH = H(size, nnzH)
    h_mat = vector(RR, size, {ih: vh for ih, vh in zip(IH, VH)})

    state_vec = vector(state(size))

    energy_init = -1 * j_mat * state_vec * state_vec - h_mat * state_vec
    print(f"Initial Energy is: {energy_init}")

main(42)
