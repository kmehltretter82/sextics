"""Independent binary-block multinomial calculation of the divided pole-two matrix."""

from functools import cache
from math import comb
from flint import nmod_mat
import sys


def grades(n):
    return [
        (a, b, c, n - a - b - c)
        for a in range(n + 1)
        for b in range(n - a + 1)
        for c in range(n - a - b + 1)
    ]


def check(p):
    mod = p * p
    k = 5 * pow(7, -1, mod) % mod
    mons = [(*v, j) for j in (2, 4) for v in grades(8 - j) if (v[2] + v[3]) % 2]
    idx = {v: i for i, v in enumerate(mons)}

    @cache
    def binary(n, g):
        if not n:
            return [1]
        old = binary(n - 1, g)
        new = [0] * (6 * n + 1)
        for i, c in enumerate(old):
            for j, w in ((0, 1), (1, k), (5, g * k), (6, 1)):
                new[i + j] = (new[i + j] + c * w) % mod
        return new

    @cache
    def Hpow(n, A, B, C, D):
        if min(A, B, C, D) < 0:
            return 0
        total = 0
        for i in range(min(n, A // 3, B, C, D) + 1):
            for j in range(min(n - i, A - 3 * i, B - i, (C - i) // 3, D - i) + 1):
                t = A + B - 4 * i - 2 * j
                if t < 0 or t % 6:
                    continue
                h = t // 6
                g = n - i - j - h
                if g < 0 or not 0 <= B - i - j <= 6 * h or not 0 <= D - i - j <= 6 * g:
                    continue
                coeff = (
                    comb(n, i)
                    * comb(n - i, j)
                    * comb(n - i - j, h)
                    * pow(k, i, mod)
                    * pow(2 * k, j, mod)
                )
                total += coeff * binary(h, 2)[B - i - j] * binary(g, 3)[D - i - j]
        return total % mod

    N = len(mons)
    L = []
    for w in mons:
        row = []
        for v in mons:
            if (p * w[4] - v[4]) % 6:
                row.append(0)
                continue
            e = (p * w[4] - v[4]) // 6
            n = 2 * p - 2 - e
            delta = [p * (x + 1) - (y + 1) for x, y in zip(w[:4], v[:4])]
            c = comb(2 * p - 2, e) * Hpow(n, *delta) % mod
            assert c % p == 0
            row.append((-c // p) % p)
        L.append(row)
    V = nmod_mat(L, p)
    f = {
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
    jac = []
    for a in range(4):
        for b in range(4):
            if a // 2 == b // 2:
                continue
            row = [0] * N
            for m, t in f.items():
                if m[b]:
                    v = tuple(y + (i == a) - (i == b) for i, y in enumerate(m)) + (2,)
                    row[idx[v]] = m[b] * t % p
            jac.append(row)
    J = nmod_mat(jac, p)
    assert J.rank() == 8

    def aug(X, Y):
        return nmod_mat(X.tolist() + Y.tolist(), p)

    assert aug(J, J * V.transpose()).rank() == 8
    unit = nmod_mat([[int(m == (3, 0, 3, 0, 2))] for m in mons], p)
    VV = V**48
    assert aug(J, (VV * unit).transpose()).rank() == 9
    stable = aug(J, VV.transpose()).rank() - 8
    # In a basis adapted to the invariant Jacobian rows, subtract the trace of its restriction.
    pivot = [
        next(j for j, x in enumerate(row) if x) for row in J.rref()[0].tolist()[:8]
    ]
    Jred = J.rref()[0]
    square = V**2
    JV = Jred * square.transpose()
    tr = (
        sum(int(square[i, i]) for i in range(N))
        - sum(int(JV[i, pivot[i]]) for i in range(8))
    ) % p
    print(
        p,
        "jacobian-rank",
        J.rank(),
        "stable quotient rank",
        stable,
        "matrix-square trace",
        tr,
        "rank with V^48A",
        9,
        flush=True,
    )


if __name__ == "__main__":
    for p in map(int, sys.argv[1:]):
        check(p)
