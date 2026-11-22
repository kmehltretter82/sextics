"""Continue the positive characteristic-11 contact and compute its obstruction span.

Argument: the directory containing the upstream contact_series_span.py.
"""

import sys
from pathlib import Path
import numpy as np
from fiber_current import QBASE, AE, BO, AB, mon, table, times, add
from fiber_point import rank, det, xbase, Dual, chart, cubic
from itertools import product
from fiber_current import ev, vec, kernel

sys.path.insert(0, str(Path(sys.argv[1]).resolve()))
import contact_series_span as upstream

p = upstream.P = 11
order = 60


def app(f):
    return {m + (0,): int(c) % p for m, c in f.items()}


def seed():
    Q = table(QBASE, [4, 8, 4, 3, 0])
    q = table(QBASE + [AB], [4, 4, 5, 10, 7, 1])
    A = table(AE, [5, 3, 1, 8, 9, 5, 6, 4])
    B = table(BO, [7, 9, 5, 2, 4, 5, 0, 6])
    W = table(mon(4, 0), [10, 6, 7, 6, 6, 8, 1, 5, 2, 4, 5, 10, 3, 6, 3, 3, 6, 10, 0])
    H = add(
        (times(W, q), 1),
        (times(times(Q, Q), Q), 1),
        (times(A, A), 1),
        (times(B, B), -1),
    )
    assert H == table(
        [
            (6, 0, 0, 0),
            (0, 6, 0, 0),
            (0, 0, 6, 0),
            (0, 0, 0, 6),
            (5, 1, 0, 0),
            (1, 5, 0, 0),
            (0, 0, 5, 1),
            (0, 0, 1, 5),
            (3, 1, 1, 1),
            (1, 1, 3, 1),
        ],
        [1, 1, 1, 1, 7, 3, 7, 10, 7, 3],
    )
    for w, z in product(range(7), repeat=2):
        pt = tuple(getattr(e, "val", e) for e in chart([Dual(y) for y in xbase], w, z))
        assert ev(add((A, 1), (B, 1)), pt) % p == cubic(w, z) % p
    tau = [
        0,
        6,
        3,
        3,
        1,
        4,
        6,
        7,
        3,
        9,
        0,
        3,
        5,
        0,
        8,
        6,
        2,
        0,
        2,
        5,
        5,
        5,
        3,
        0,
        8,
        0,
        6,
        1,
        5,
        0,
        5,
        2,
        10,
        9,
        7,
        10,
        1,
        10,
        2,
        1,
        4,
        5,
        7,
        1,
        0,
    ]
    Qd = table(QBASE, tau[:5])
    qd = table(QBASE, tau[5:10])
    Ad = table(AE, tau[10:18])
    Bd = table(BO, tau[18:26])
    Wd = table(mon(4, 0), tau[26:])
    Z = {(0, 0, 0, 0, 1): 1}

    def plus(a, b, s=1):
        return upstream.poly_add(a, b, s)

    mul = upstream.poly_mul
    Q = app(Q)
    Qd = app(Qd)
    g1 = plus(plus(mul(Z, Z), mul(Z, Q), -1), mul(Q, Q))
    dg1 = plus(mul(Q, Qd), mul(Q, Qd))
    dg1 = plus(dg1, mul(Z, Qd), -1)
    return upstream.encode(
        [plus(Z, Q), app(q), g1, app(W), app(A), app(B)]
    ), upstream.encode([Qd, app(qd), dg1, app(Wd), app(Ad), app(Bd)])


def quadratic(a, b):
    return (
        np.bincount(
            upstream.term_row,
            weights=upstream.term_sign * a[upstream.term_i] * b[upstream.term_j],
            minlength=70,
        ).astype(np.int64)
        % p
    )


def matrix(coeff):
    M = np.zeros((60, 72), dtype=np.int64)
    M[upstream.ideal_rows, upstream.ideal_cols] = coeff[upstream.ideal_parameters]
    return M


def main():
    initial, tangent = seed()
    diff = upstream.differential(initial)
    _, cols = upstream.base_rref(diff)
    assert len(cols) == 70 and not np.any(diff @ tangent % p)
    inverse = upstream.base_inverse(diff[:, cols])
    jets = np.zeros((order + 1, 90), dtype=np.int64)
    jets[0] = initial
    jets[1] = tangent
    for t in range(2, order + 1):
        rhs = sum(quadratic(jets[i], jets[t - i]) for i in range(1, t))
        jets[t, cols] = -inverse @ rhs % p
        assert not np.any(
            sum(quadratic(jets[i], jets[t - i]) for i in range(t + 1)) % p
        )
    M = np.stack([matrix(z) for z in jets])
    _, cols = upstream.base_rref(M[0])
    assert len(cols) == 59
    _, rows = upstream.base_rref(M[0][:, cols].T)
    (free,) = set(range(60)) - set(rows)
    inv = upstream.base_inverse(M[0][np.ix_(rows, cols)].T)
    ell = np.zeros((order + 1, 60), dtype=np.int64)
    ell[0, free] = 1
    ell[0, rows] = -inv @ M[0, free, cols] % p
    for t in range(1, order + 1):
        rhs = sum(M[i][:, cols].T @ ell[t - i] for i in range(1, t + 1))
        ell[t, rows] = -inv @ rhs % p
    for t in range(order + 1):
        assert not np.any(sum(ell[t - i] @ M[i] for i in range(t + 1)) % p)
    J = upstream.fixed_jacobian_directions(initial)
    assert not np.any(ell @ J % p)
    for u in [None, 0, 1, 2]:
        selected = [
            i for i, m in enumerate(upstream.minus_basis) if u is None or m[-1] == u
        ]
        mat = ell[:, selected]
        rk, rows, cols = rank(mat)
        expected = {None: (48, 9), 0: (32, 1), 1: (16, 7), 2: (0, 1)}
        assert (rk, det(mat[np.ix_(rows, cols)])) == expected[u]
        print(
            "sector",
            u,
            "rank",
            rk,
            "minor",
            det(mat[np.ix_(rows, cols)]),
            "last jet row",
            max(rows, default=None),
        )
    # Check the scalar between the obstruction row and the canonical
    # residue on all 60 anti-invariant weighted sextics at the new seed.
    f1, f2, g1, g2, A, B = upstream.decode(initial)
    ds = upstream.decode(tangent)
    physical = lambda f: {m[:4]: v for m, v in f.items() if not m[-1]}
    Q = physical(f1)
    h = physical(upstream.poly_add(A, B))
    dh = physical(upstream.poly_add(ds[4], ds[5]))
    forms = [physical(f2), h, add((times(Q, Q), 3)), physical(g2)]
    b9 = mon(9)
    mac = np.array(
        [
            vec(times(f, {m: 1}), b9)
            for f, d in zip(forms, [2, 3, 4, 4])
            for m in mon(9 - d)
        ]
    )
    b10 = mon(10)
    mac10 = np.array(
        [
            vec(times(f, {m: 1}), b10)
            for f, d in zip(forms, [2, 3, 4, 4])
            for m in mon(10 - d)
        ]
    )
    assert len(upstream.base_rref(mac10)[1]) == len(b10) == 286
    soc = kernel(mac)
    minusQ = add((Q, -1))
    top = []
    for m in upstream.minus_basis:
        poly = {m[:4]: 1}
        for _ in range(m[4]):
            poly = times(poly, minusQ)
        top.append(int(soc @ vec(times(poly, dh), b9)) % p)
    j = upstream.minus_basis.index((3, 0, 3, 0, 0))
    assert top[j] == 1 and ell[0, j]
    ratio = top[j] * pow(int(ell[0, j]), -1, p) % p
    assert np.array_equal(top, ratio * ell[0] % p)
    print("seed raw-socle/obstruction scalar", ratio, "target obstruction", ell[0, j])
    return ell, jets


if __name__ == "__main__":
    main()
