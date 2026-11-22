"""Verify the smooth branch sextic and its disjoint contact divisor on the F11 seed."""

from itertools import permutations, product
from flint import nmod_mat, nmod_poly

p = 11


def add(*terms):
    h = {}
    for num, f in terms:
        for e, v in f.items():
            h[e] = (h.get(e, 0) + num * v) % p
    return {e: v for e, v in h.items() if v}


def mul(*factors):
    out = {(0, 0): 1}
    for f in factors:
        z = {}
        for (i, j), val in out.items():
            for (k, l), c in f.items():
                e = (i + k, j + l)
                z[e] = (z.get(e, 0) + val * c) % p
        out = {e: v for e, v in z.items() if v}
    return out


def data(lift=False):
    # m,n,u,v,e=(2,6,1,7,7). In the physical Segre coordinates K=5.
    m, n, u, v, e = (57, 83, 67, 84, 18) if lift else (2, 6, 1, 7, 7)
    K = e * (v - u) ** 2 * (m + n) * pow((n - m) ** 2, -1, p) % p
    a = {(0, 0): 1, (1, 1): K}
    b = {(0, 0): n, (1, 1): m * K % p}
    c = {(1, 0): 1, (0, 1): 1}
    d = {(1, 0): v, (0, 1): u}
    kk = 5 * pow(7, -1, p) % p
    core = [
        (1, [a] * 6),
        (1, [b] * 6),
        (1, [c] * 6),
        (1, [d] * 6),
        (kk, [a] * 5 + [b]),
        (2 * kk, [a] + [b] * 5),
        (kk, [c] * 5 + [d]),
        (3 * kk, [c] + [d] * 5),
        (kk, [a] * 3 + [b, c, d]),
        (2 * kk, [a, b] + [c] * 3 + [d]),
    ]
    H = add(*[(k, mul(*forms)) for k, forms in core])
    P, Q, R, T, S = (70, 19, 114, 36, 44) if lift else (4, 8, 4, 3, 0)
    q = add(
        (P, mul(a, a)), (Q, mul(b, b)), (R, mul(c, c)), (T, mul(c, d)), (S, mul(d, d))
    )
    hmat = [[4, 4, 9, 4], [6, 0, 6, 9], [9, 0, 10, 0], [7, 3, 9, 6]]
    if lift:
        hmat = [
            [37, 4, 119, 26],
            [17, 11, 116, 64],
            [75, 88, 54, 99],
            [62, 25, 53, 105],
        ]
    h = {(i, j): v for i, row in enumerate(hmat) for j, v in enumerate(row) if v}
    hb = {e: (-1) ** sum(e) * v % p for e, v in h.items()}
    assert H == add((1, mul(q, q, q)), (1, mul(h, hb)))
    return H, q, h, hb, K


def homogeneous_derivatives(f, deg):
    # Coordinates on the first P^1 are W0,W1 and the chart variable is W1/W0.
    for slot in (0, 1):
        for component in (0, 1):
            out = {}
            for term, v in f.items():
                k = term[slot] if component else deg[slot] - term[slot]
                key = list(term)
                key[slot] -= component
                if k % p:
                    out[tuple(key)] = v * k % p
            bideg = list(deg)
            bideg[slot] -= 1
            yield out, tuple(bideg)


def rows_for(gens, bidegree):
    n, m = bidegree
    rows = []
    for f, (a, b) in gens:
        for i, j in product(range(n - a + 1), range(m - b + 1)):
            row = [0] * ((n + 1) * (m + 1))
            for (k, l), v in f.items():
                row[(i + k) * (m + 1) + j + l] = v
            rows.append(row)
    return rows


def macaulay(gens, bidegree, minor=True):
    rows = rows_for(gens, bidegree)
    M = nmod_mat(rows, p)
    rk = M.rank()
    if minor and rk == M.ncols():
        R, rr = M.transpose().rref()
        selected = [next(i for i, c in enumerate(row) if c) for row in R.tolist()[:rr]]
        det = int(nmod_mat([rows[k] for k in selected], p).det())
        assert det
    else:
        det = None
    return len(rows), M.ncols(), rk, det


def lift_rank(rows):
    import numpy as np

    M = np.array(rows, dtype=np.int64) % 121
    rank = 0
    pivot = []
    for j in range(M.shape[1]):
        candidates = np.flatnonzero(M[rank:, j] % 11)
        if not len(candidates):
            continue
        i = rank + int(candidates[0])
        M[[rank, i]] = M[[i, rank]]
        M[rank] = M[rank] * pow(int(M[rank, j]), -1, 121) % 121
        coeff = M[rank + 1 :, j].copy()
        M[rank + 1 :] = (M[rank + 1 :] - coeff[:, None] * M[rank]) % 121
        pivot.append(j)
        rank += 1
    extra = [j for j in range(M.shape[1]) if j not in pivot]
    rest = M[rank:, extra]
    assert not np.any(rest % 11)
    N = nmod_mat((rest // 11).tolist(), 11)
    return rank, N.rank()


def ruling_resultant(q, h):
    # Bihomogenize in z, then take the degree-(2,3) Sylvester determinant.
    def coeff(f, deg):
        return [
            nmod_poly([f.get((i, j), 0) for i in range(deg + 1)], 11)
            for j in reversed(range(deg + 1))
        ]

    qz, hz = coeff(q, 2), coeff(h, 3)
    zero = nmod_poly([], 11)
    mat = []
    for row, times in [(qz, 3), (hz, 2)]:
        for j in range(times):
            mat.append([zero] * j + row + [zero] * (times - j - 1))
    det = zero
    for perm in permutations(range(5)):
        v = nmod_poly(
            [(-1) ** sum(perm[k] > perm[l] for k in range(5) for l in range(k + 1, 5))],
            11,
        )
        for k, l in enumerate(perm):
            v *= mat[k][l]
        det += v
    assert det.degree() == 12 and det.gcd(det.derivative()).degree() == 0
    return det, sorted((f.degree(), exp) for f, exp in det.factor()[1])


if __name__ == "__main__":
    H, q, h, hb, K = data()
    fixed = list(product((0, 6), repeat=2))
    corners = [H.get(e, 0) for e in fixed]
    print("Segre K", K, "fixed corner values", corners)
    assert corners == [1, 8, 0, 6]
    assert [H.get(e, 0) for e in ((4, 0), (5, 1), (6, 2))] == [5, 0, 5]
    common = macaulay([(q, (2, 2)), (h, (3, 3)), (hb, (3, 3))], (7, 7))
    print("q,h,conjugate h", common)
    assert common == (86, 64, 64, 7)
    resultant, factors = ruling_resultant(q, h)
    print("ruling resultant", resultant, "factor degrees and multiplicities", factors)
    assert factors == [(1, 1), (5, 1), (6, 1)]
    # A divisor from B0=W0^3 Z0^3, i.e. affine chart coefficient 1.
    B0 = {(0, 0): 1}
    clean = macaulay([(q, (2, 2)), (H, (6, 6)), (B0, (3, 3))], (10, 10))
    print("auxiliary divisor avoids q=0", clean)
    assert clean[1] == clean[2]
    # K_C tau=3L+D = 5L-sigma D via the denominator q.
    # V is the (5,5) numerator space, while W vanishes on C intersect B0.
    V = rows_for([(q, (2, 2)), (hb, (3, 3))], (5, 5))
    W = rows_for([(B0, (3, 3))], (5, 5))
    ranks = [nmod_mat(M, 11).rank() for M in (V, W, V + W)]
    print(
        "degree-five restriction spaces V,W,V+W",
        ranks,
        "evaluation rank",
        ranks[2] - ranks[1],
    )
    assert ranks == [24, 9, 32] and ranks[2] - ranks[1] == 23
    p = 121
    H, q, h, hb, K = data(lift=True)
    print("lifted fixed corner values", [H.get(e, 0) for e in fixed])
    print("lifted moving sheet norm holds in all 49 coefficients")
    der = list(homogeneous_derivatives(H, (6, 6)))
    rows = rows_for(der, (19, 19))
    ranks = lift_rank(rows)
    print(
        "gradient matrix",
        len(rows),
        len(rows[0]),
        "unit and valuation-one pivots",
        ranks,
        flush=True,
    )
    assert sum(ranks) == len(rows[0])
