"""Compute scaled Cartier on the pole-two anti-sign Jacobian quotient."""

from itertools import product
from math import comb
import sys
import numpy as np


def powers(n):
    return [
        (a, b, c, n - a - b - c)
        for a in range(n + 1)
        for b in range(n - a + 1)
        for c in range(n - a - b + 1)
    ]


def rref(A, p):
    a = A.copy() % p
    row = 0
    piv = []
    for j in range(a.shape[1]):
        avail = np.flatnonzero(a[row:, j])
        if not len(avail):
            continue
        pick = row + avail[0]
        a[[pick, row]] = a[[row, pick]]
        a[row] = a[row] * pow(int(a[row, j]), -1, p) % p
        fac = a[:, j].copy()
        fac[row] = 0
        a = (a - fac[:, None] * a[row]) % p
        piv.append(j)
        row += 1
        if row == len(a):
            break
    return a[:row], piv


def run(p):
    mod = p * p
    k = 5 * pow(7, -1, mod) % mod
    H = {
        (6, 0, 0, 0): 1,
        (0, 6, 0, 0): 1,
        (0, 0, 6, 0): 1,
        (0, 0, 0, 6): 1,
        (5, 1, 0, 0): k,
        (1, 5, 0, 0): 2 * k,
        (0, 0, 5, 1): k,
        (0, 0, 1, 5): 3 * k,
        (3, 1, 1, 1): k,
        (1, 1, 3, 1): 2 * k,
    }
    # Four exponents here specify the polynomial numerator; e_polar includes the toric +1.
    mons = [m + (j,) for j in (2, 4) for m in powers(8 - j) if (m[2] + m[3]) % 2 == 1]
    terms = {m: i for i, m in enumerate(mons)}
    size = len(mons)
    jac = []
    for i, j in product(range(4), repeat=2):
        if (i // 2) == (j // 2):
            continue
        col = np.zeros(size, dtype=np.int64)
        for m, coef in H.items():
            if m[j]:
                v = tuple(m[l] + int(i == l) - int(j == l) for l in range(4)) + (2,)
                col[terms[v]] = coef * m[j] % p
        jac.append(col)
    J, piv = rref(np.array(jac), p)
    rest = [i for i in range(size) if i not in piv]
    assert len(piv) == 8 and len(rest) == 48
    # Reduce raw polynomial coefficients modulo Jac(H). Rows represent free quotient coordinates.
    Q = np.zeros((len(rest), size), dtype=np.int64)
    Q[:, rest] = np.eye(len(rest), dtype=np.int64)
    Q[:, piv] = -J[:, rest].T % p
    blocks = {}
    for j, jp in product((2, 4), repeat=2):
        if (p * jp - j) % 6:
            continue
        e = (p * jp - j) // 6
        power = 2 * p - 2 - e
        blocks[power] = (j, jp, comb(2 * p - 2, e) % mod)
    maxn = max(blocks)
    coef = np.ones((1, 1, 1), dtype=np.int64)
    matrix = np.zeros((size, size), dtype=np.int64)
    for n in range(1, maxn + 1):
        old = coef
        side = 6 * n + 1
        coef = np.zeros((side, side, side), dtype=np.int64)
        for (a, b, c, _), weight in H.items():
            s = len(old)
            coef[a : a + s, b : b + s, c : c + s] += old * weight
        coef %= mod
        if n not in blocks:
            continue
        j, jp, binomial = blocks[n]
        for v, w in product(mons, repeat=2):
            if v[4] != j or w[4] != jp:
                continue
            vals = tuple(p * (w[q] + 1) - (v[q] + 1) for q in range(4))
            if min(vals) < 0 or sum(vals) != 6 * n:
                continue
            amount = int(coef[vals[:3]]) * binomial % mod
            assert amount % p == 0
            matrix[terms[w], terms[v]] = -(amount // p) % p
    assert np.all(Q @ matrix @ J.T % p == 0)
    V = Q @ matrix[:, rest] % p
    target = Q[:, terms[(3, 0, 3, 0, 2)]]
    a = target.copy()
    for _ in range(48):
        a = V @ a % p
    stable = np.eye(len(rest), dtype=np.int64)
    for _ in range(48):
        stable = V @ stable % p
    _, cols = rref(stable, p)
    tr = int(np.trace(V @ V) % p)
    print(
        f"p={p}: raw dimension {size}, quotient {len(rest)}, trace(V²) {tr}, "
        f"stable rank {len(cols)}, first coordinate of V^48(A) {a[0]}",
        flush=True,
    )
    expected = {11: -36, 17: 16, 23: -48}
    if p in expected:
        assert tr == expected[p] % p
    assert np.any(a)
    receipt = {11: (30, 8), 17: (32, 7), 23: (32, 19)}
    if p in receipt:
        assert (len(cols), int(a[0])) == receipt[p]

    return mons, rest, V


if __name__ == "__main__":
    for a in sys.argv[1:]:
        run(int(a))
