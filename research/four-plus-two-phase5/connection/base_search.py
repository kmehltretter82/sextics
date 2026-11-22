"""Check rational independence of the inhomogeneous column by modular minors."""

from itertools import product
from random import Random
import numpy as np

PRIME = 1009
DEGREE = 12


def inverse(matrix):
    """Inverse over F_p, or None for a singular evaluation point."""
    size = len(matrix)
    a = np.concatenate((matrix, np.eye(size, dtype=np.int64)), axis=1) % PRIME
    for column in range(size):
        candidates = np.flatnonzero(a[column:, column])
        if not len(candidates):
            return None
        pivot = column + candidates[0]
        a[[column, pivot]] = a[[pivot, column]]
        a[column] = a[column] * pow(int(a[column, column]), -1, PRIME) % PRIME
        factors = a[:, column].copy()
        factors[column] = 0
        a = (a - factors[:, None] * a[column]) % PRIME
    return a[:, size:]


def evaluate(r, s):
    p = PRIME
    k = 5 * pow(7, -1, p) % p
    third = pow(3, -1, p)
    f = (1 + r**6 + k * r + 2 * k * r**5) % p
    g = (1 + s**6 + k * s + 3 * k * s**5) % p
    fr, gs = (6 * r**5 + k + 10 * k * r**4) % p, (6 * s**5 + k + 15 * k * s**4) % p
    b, br, bs = k * r * s % p, k * s % p, k * r % p
    coeff = np.array([f, b, 2 * b, g])
    mat = np.zeros((6, 6), dtype=np.int64)
    for j in range(3):
        mat[j : j + 4, j] = coeff
        mat[j : j + 4, j + 3] = coeff * (2 * j + 1 - 4 * third * np.arange(4)) % p
    inv = inverse(mat)
    if inv is None:
        return None
    red = inv[:3]

    def connection(derivative):
        columns = np.zeros((6, 3), dtype=np.int64)
        for j in range(3):
            columns[j : j + 4, j] = derivative
        return -2 * third * (red @ columns) % p

    ar = connection([fr, br, 2 * br, 0])
    ass = connection([0, bs, 2 * bs, gs])
    disc = (
        -27 * f * f * g * g
        - 32 * f * b**3
        + 36 * f * b * b * g
        + 4 * b**4
        - 4 * b**3 * g
    ) % p
    df = (-54 * f * g * g - 32 * b**3 + 36 * b * b * g) % p
    db = (-96 * f * b * b + 72 * f * b * g + 16 * b**3 - 12 * b * b * g) % p
    dg = (-54 * f * f * g + 36 * f * b * b - 4 * b**3) % p
    divisors = {
        "f": (f, fr, 0),
        "g": (g, 0, gs),
        "disc": (disc, (df * fr + db * br) % p, (dg * gs + db * bs) % p),
        "reciprocal": ((g - 8 * f) % p, -8 * fr % p, gs),
    }
    return ar, ass, red[:, 2], divisors


def rank(matrix):
    a = matrix.copy() % PRIME
    row = 0
    for column in range(a.shape[1]):
        candidates = np.flatnonzero(a[row:, column])
        if not len(candidates):
            continue
        pivot = row + candidates[0]
        a[[row, pivot]] = a[[pivot, row]]
        normalized = a[row, column:] * pow(int(a[row, column]), -1, PRIME) % PRIME
        a[row + 1 :, column:] = (
            a[row + 1 :, column:] - a[row + 1 :, column, None] * normalized
        ) % PRIME
        row += 1
        if row == len(a):
            break
    return row


def check(support):
    p = PRIME
    exps = np.array(list(product(range(DEGREE + 1), repeat=2)))
    count = 6 * len(exps)
    rng = Random(2000)
    rows = []
    while 3 * len(rows) < count + 18:
        r, s = [rng.randrange(1, p) for _ in range(2)]
        data = evaluate(r, s)
        if data is None:
            continue
        ar, ass, target, divisors = data
        D, dr, ds = 1, 0, 0
        for component in support:
            v, vr, vs = divisors[component]
            D, dr, ds = D * v % p, (dr * v + D * vr) % p, (ds * v + D * vs) % p
        if not D:
            continue
        inv_d = pow(D, -1, p)
        rp = np.array([pow(r, i, p) for i in range(DEGREE + 1)])
        sp = np.array([pow(s, j, p) for j in range(DEGREE + 1)])
        val = rp[exps[:, 0]] * sp[exps[:, 1]] % p * inv_d % p
        val_r = val * (exps[:, 0] * pow(r, -1, p) - dr * inv_d) % p
        val_s = val * (exps[:, 1] * pow(s, -1, p) - ds * inv_d) % p
        identity = np.eye(3, dtype=np.int64)[:, :, None]
        B = ar[:, :, None] * val + identity * val_r
        C = ass[:, :, None] * val + identity * val_s
        rows.append(
            np.hstack([B.reshape(3, -1), C.reshape(3, -1), target.reshape(3, 1)]) % p
        )
    matrix = np.concatenate(rows)
    augmented_rank = rank(matrix)
    assert augmented_rank == count + 1
    print(
        f"{'*'.join(support)}: {len(matrix)} rows, rank {count}, "
        f"augmented rank {augmented_rank} over F_{p}.",
        flush=True,
    )


if __name__ == "__main__":
    for divisor in [
        ("reciprocal",),
        ("disc",),
        ("disc", "reciprocal"),
        ("f", "g", "disc"),
    ]:
        check(divisor)
