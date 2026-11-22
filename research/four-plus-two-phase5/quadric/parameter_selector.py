"""Midpoint Bezout selector for the actual 70-equation smooth contact parameter slice."""

import sys
from pathlib import Path
from argparse import ArgumentParser
import numpy as np
from flint import nmod_mat
from fast_sweep import extend


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("upstream")
    parser.add_argument("--helpers", default=str(Path(__file__).resolve().parent))
    opts = parser.parse_args()
    sys.path.insert(0, opts.helpers)
    sys.argv[1] = opts.upstream
    import smooth_span as seed

    st = seed.upstream
    p = 11
    y0, y1 = seed.seed()
    J = st.differential(y0)
    _, cols = st.base_rref(J)
    J0 = J[:, cols]
    zero = np.zeros(70, dtype=np.int64)

    def det(A):
        return int(nmod_mat(A.shape[0], A.shape[1], A.flatten().tolist(), p).det()) % p

    norm = det(J0)
    assert norm == 8
    # Work in x=(z,u) with y=y0+z_piv+u*y1. Choose H=u+Lz through p,
    # so the tangent vector has H'=1. Replace u by H-Lz in the pencil of linear forms.
    L = np.arange(1, 71, dtype=np.int64) ** 2 % p
    eye = np.eye(90, dtype=np.int64)[:, cols]
    change = (eye - y1[:, None] * L) % p

    def coeff(y):
        midpoint = (6 * (y + y0)) % p
        full = st.differential(midpoint)
        return full @ change % p, full @ y1 % p

    E0, k0 = coeff(y0)
    assert det(E0) == norm and not np.any(k0)
    rng = np.random.default_rng(121)
    for _ in range(12):
        z = rng.integers(0, p, 70)
        H = int(rng.integers(p))
        y = (y0 + change @ z + H * y1) % p
        M, k = coeff(y)
        res = (seed.quadratic(y, y) - seed.quadratic(y0, y0)) % p
        assert not np.any((M @ z + k * H - res) % p)
    arr = extend(seed, 32, 8)
    Ms = []
    ks = []
    for n, y in enumerate(arr):
        D = st.differential(y if n == 0 else (y * 6) % p)
        Ms.append(D @ change % p)
        ks.append(D @ y1 % p)
    zs = np.array([zero, zero] + [y[cols] for y in arr[2:]], dtype=np.int64)
    Hs = zs @ L % p
    Hs[1] = 1
    assert all(
        not np.any(sum(Ms[i] @ zs[n - i] + ks[i] * Hs[n - i] for i in range(n + 1)) % p)
        for n in range(32)
    )
    print(
        "determinant at seed:",
        norm,
        "mod 11; twelve full 70-coordinate midpoint identity tests; all 70 series identities through order 31",
    )


if __name__ == "__main__":
    main()
