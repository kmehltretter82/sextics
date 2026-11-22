"""Plane-section interpolation on the component through the smooth 11-adic contact.

Usage: python smooth_sweep.py UPSTREAM_CONTACT_SPAN_DIR [degree [precision]]
Requirements: NumPy, SymPy, python-flint.
"""

import sys
from itertools import product
import numpy as np
from flint import nmod_poly, nmod_mat
import smooth_span as seed

p = 11


def S(coeff=0):
    return nmod_poly(coeff if isinstance(coeff, list) else [int(coeff)], p)


def extend(n):
    st = seed.upstream
    y0, y1 = seed.seed()
    D = st.differential(y0)
    _, cols = st.base_rref(D)
    inv = st.base_inverse(D[:, cols])
    z = np.zeros((n, 90), dtype=np.int64)
    z[0] = y0
    z[1] = y1
    zd = z.astype(np.float64)
    assert n * 100 * len(st.term_row) < 2**52
    for i in range(2, n):
        gram = zd[1:i].T @ zd[i - 1 : 0 : -1]
        rhs = (
            np.bincount(
                st.term_row,
                weights=st.term_sign * gram[st.term_i, st.term_j],
                minlength=70,
            ).astype(np.int64)
            % p
        )
        z[i, cols] = -inv @ rhs % p
        zd[i] = z[i]
        assert not np.any((D @ z[i] + rhs) % p)
    params = [S(a.tolist()) for a in z.T]
    terms = [S() for _ in st.plus_basis]
    for row, i, j, sign in zip(st.term_row, st.term_i, st.term_j, st.term_sign):
        terms[row] = terms[row] + int(sign) * params[i].mul_low(params[j], n)
    targets = seed.quadratic(y0, y0)
    assert all(t == S(e) for t, e in zip(terms, targets))
    return params


def restricted(params, s, N):
    st = seed.upstream
    offset = 0
    blocks = []
    for d, par in st.blocks:
        f = {}
        for i, j, _, b, u in st.monomials(d, par):
            key = i, j, u
            f[key] = f.get(key, S()) + pow(s, b, p) * params[offset]
            offset += 1
        blocks.append(f)
    Q, Q2, _, _, A, B = blocks
    h = {m: A.get(m, S()) + B.get(m, S()) for m in set(A) | set(B)}

    def elim(f):
        out = {}
        for (a, b, k), v in f.items():
            if k == 0:
                out[a, b] = out.get((a, b), S()) + v.mul_low(Q[0, 0, 1], N)
            else:
                assert k == 1
                for (i, j, t), c in Q.items():
                    if t == 0:
                        out[i + a, j + b] = out.get((i + a, j + b), S()) - v.mul_low(
                            c, N
                        )
        return out

    return elim(Q2), elim(h)


def eval_grad(f, a, b, n):
    out = [S(), S(), S()]
    powers = lambda x: [S(1), x, x.mul_low(x, n), x.mul_low(x, n).mul_low(x, n)]
    A, B = powers(a), powers(b)
    for (i, j), v in f.items():
        for k, c, a0, b0 in [(0, 1, i, j), (1, i, i - 1, j), (2, j, i, j - 1)]:
            if c:
                out[k] += c * A[a0].mul_low(B[b0], n).mul_low(v, n)
    return out


def points(q, h, N):
    seeds = []
    for a, b in product(range(p), repeat=2):
        if not any(eval_grad(f, S(a), S(b), 1)[0] for f in [q, h]):
            seeds.append((a, b))
    assert seeds == [(2, 8), (4, 7), (8, 7)]
    for a, b in seeds:
        a, b = S(a), S(b)
        n = 1
        while n < N:
            n = min(N, 2 * n)
            Q, Qa, Qb = eval_grad(q, a, b, n)
            H, Ha, Hb = eval_grad(h, a, b, n)
            v = (Qa.mul_low(Hb, n) - Qb.mul_low(Ha, n)).inverse_series_trunc(n)
            a, b = (
                a - v.mul_low(Hb.mul_low(Q, n) - Qb.mul_low(H, n), n),
                b - v.mul_low(Qa.mul_low(H, n) - Ha.mul_low(Q, n), n),
            )
        assert not eval_grad(q, a, b, n)[0] and not eval_grad(h, a, b, n)[0]
        yield a, b


def main():
    D = int(sys.argv[2]) if len(sys.argv) > 2 else 64
    n = int(sys.argv[3]) if len(sys.argv) > 3 else 760
    params = extend(n)
    q, h = restricted(params, 2, n)
    pts = list(points(q, h, n))
    columns = (D + 1) * (D + 2) // 2
    M = np.zeros((len(pts), n, columns), dtype=np.uint8)
    for sheet, (a, b) in enumerate(pts):
        bp = [S(1)]
        for j in range(D):
            bp.append(bp[-1].mul_low(b, n))
        for i in range(D + 1):
            for j in range(D - i + 1):
                col = (i + j) * (i + j + 1) // 2 + i
                vals = list(bp[j])
                M[sheet, : len(vals), col] = vals
                bp[j] = bp[j].mul_low(a, n)
    base = M.reshape(-1, columns)
    rank = nmod_mat(*base.shape, base.flatten().tolist(), p).rank()
    print(
        "smooth F11 slice d=2c; all three rational sheets; order",
        n - 1,
        "degree",
        D,
        "rank",
        rank,
        "of",
        columns,
        "nullity",
        columns - rank,
        flush=True,
    )


if __name__ == "__main__":
    main()
