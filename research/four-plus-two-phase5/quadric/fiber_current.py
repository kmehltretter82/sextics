"""Canonical single-sheet current in the 11-adic lift of a coupled non-nodal contact."""

from functools import reduce
from itertools import permutations, product
from operator import mul
import numpy as np
import fiber_point as fp

p = 11
inv = lambda x: pow(int(x) % p, -1, p)
QBASE = [(2, 0, 0, 0), (0, 2, 0, 0), (0, 0, 2, 0), (0, 0, 1, 1), (0, 0, 0, 2)]
AB = (1, 1, 0, 0)
AE = [
    (3, 0, 0, 0),
    (0, 3, 0, 0),
    (1, 0, 2, 0),
    (1, 0, 1, 1),
    (1, 0, 0, 2),
    (0, 1, 2, 0),
    (0, 1, 1, 1),
    (0, 1, 0, 2),
]
BO = [
    (2, 0, 1, 0),
    (2, 0, 0, 1),
    (0, 2, 1, 0),
    (0, 2, 0, 1),
    (0, 0, 3, 0),
    (0, 0, 2, 1),
    (0, 0, 1, 2),
    (0, 0, 0, 3),
]


def mon(n, parity=None):
    return [
        t
        for t in product(range(n + 1), repeat=4)
        if sum(t) == n and (parity is None or (t[2] + t[3]) % 2 == parity)
    ]


def table(basis, values):
    return dict(zip(basis, map(int, values)))


def times(f, g):
    z = {}
    for m, v in f.items():
        for n, w in g.items():
            u = tuple(a + b for a, b in zip(m, n))
            z[u] = (z.get(u, 0) + v * w) % p
    return {e: v for e, v in z.items() if v}


def add(*terms):
    z = {}
    for f, s in terms:
        for e, v in f.items():
            z[e] = (z.get(e, 0) + s * v) % p
    return {e: v for e, v in z.items() if v}


def ev(f, x):
    return (
        sum(
            c * reduce(mul, (pow(int(a), i, p) for a, i in zip(x, m)), 1)
            for m, c in f.items()
        )
        % p
    )


def deriv(f, j):
    return {
        tuple(v - int(i == j) for i, v in enumerate(m)): m[j] * c % p
        for m, c in f.items()
        if m[j]
    }


def vec(f, base):
    return [f.get(m, 0) for m in base]


def rref(a):
    a = np.array(a, dtype=np.int64) % p
    row = 0
    cols = []
    for j in range(a.shape[1]):
        i = next((i for i in range(row, len(a)) if a[i, j] % 11), None)
        if i is None:
            continue
        a[[row, i]] = a[[i, row]]
        a[row] = a[row] * inv(a[row, j]) % p
        tail = a[:, j].copy()
        tail[row] = 0
        a = (a - tail[:, None] * a[row]) % p
        row += 1
        cols.append(j)
    return a, cols


def solve(a, b):
    n = len(a[0])
    aug, cols = rref(np.column_stack((a, b)))
    assert cols == list(range(n))
    v = aug[:n, -1]
    assert not np.any((np.array(a) @ v - b) % p)
    return v


def kernel(a):
    rr, cols = rref(a)
    (free,) = set(range(len(a[0]))) - set(cols)
    v = np.zeros(len(a[0]), dtype=np.int64)
    v[free] = 1
    v[cols] = -rr[: len(cols), free] % p
    assert not np.any((np.array(a) @ v) % p)
    return v


def main(sign=1, lifted=None, modulus=11):
    global p
    p = modulus
    kappa = 5 * inv(7) % p
    if lifted is None:
        lifted = fp.xbase + fp.hbase.flatten().tolist() + [2, 10]
    vals = list(map(int, lifted[:10]))
    coeff = list(map(int, lifted[10:26]))
    wp, zp = map(int, lifted[-2:])
    m, n, u, v, e = vals[:5]
    K = e * (v - u) ** 2 * (m + n) * inv((n - m) ** 2) % p

    def physical(w, z):
        return tuple(
            a % p for a in (1 + K * w * z, n + m * K * w * z, w + z, v * w + u * z)
        )

    def cubic(w, z):
        return (
            sum(
                c * w**i * z**j for c, (i, j) in zip(coeff, product(range(4), repeat=2))
            )
            % p
        )

    H = table(
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
        [1, 1, 1, 1, kappa, 2 * kappa, kappa, 3 * kappa, kappa, 2 * kappa],
    )
    Q1 = table(QBASE, vals[5:])
    D = -inv(m + n)
    Q2 = table(QBASE + [AB], [D * m * n, D, e * u * v, -e * (u + v), e, 1])
    base = AE + BO
    samples = list(product(range(4), repeat=2))
    fpoints = [physical(w, z) for w, z in samples]
    assert all(ev(Q2, pt) % p == 0 for pt in fpoints)
    hv = solve(
        [[int(ev({v: 1}, pt)) for v in base] for pt in fpoints],
        [cubic(sign * w, sign * z) % p for w, z in samples],
    )
    A, B = table(AE, hv[:8]), table(BO, hv[8:])
    h = add((A, 1), (B, 1))
    assert all(
        ev(h, physical(w, z)) % p == cubic(sign * w, sign * z) % p
        for w, z in product(range(7), repeat=2)
    )
    C2 = times(Q1, Q1)
    G1 = add((C2, 3))
    R = add((H, 1), (times(Q1, C2), -1), (times(A, A), -1), (times(B, B), 1))
    mb4 = mon(4, 0)
    mb6 = mon(6, 0)
    mult = np.array([vec(times(Q2, {m: 1}), mb6) for m in mb4]).T
    G = table(mb4, solve(mult, vec(R, mb6)))
    J = []
    for poly, bs, k in [
        (C2, QBASE, -3),
        (G, QBASE, -1),
        (A, AE, -2),
        (B, BO, 2),
        (Q2, mb4, -1),
    ]:
        J.extend([vec(times(poly, {m: k}), mb6) for m in bs])
    tangent = kernel(np.array(J).T)
    dq = table(QBASE, tangent[5:10])
    dh = table(base, tangent[10:26])
    vr = times(dh, {(3, 0, 3, 0): 1})
    forms = [Q2, h, G1, G]
    degree = [2, 3, 4, 4]
    b9 = mon(9)
    mac = np.array(
        [vec(times(f, {m: 1}), b9) for f, d in zip(forms, degree) for m in mon(9 - d)]
    )
    soc = kernel(mac)
    res = lambda f: int(soc @ vec(f, b9)) % p
    jac = {}
    for perm in permutations(range(4)):
        term = reduce(times, [deriv(f, i) for i, f in zip(perm, forms)])
        s = (-1) ** sum(i > j for v, i in enumerate(perm) for j in perm[v + 1 :])
        jac = add((jac, 1), (term, s))
    norm = res(jac)
    ell = 96 * inv(norm) * res(vr) % p
    a, _, c, _ = physical(wp, zp)
    r, s, alpha = 4, 2, sign * c * inv(a) % p
    pt = [1, r, alpha, s * alpha]
    assert ev(Q2, pt) % p == ev(h, pt) % p == 0

    def partial(f):
        db, dc, dd = (ev(deriv(f, i), pt) % p for i in [1, 2, 3])
        return np.array([db, alpha * dd, dc + s * dd], dtype=np.int64) % p

    qr, qs, qa = map(int, partial(Q2))
    hr, hs, ha = map(int, partial(h))
    br, bs, ba = map(int, -partial(Q1) % p)
    j = (ev(dq, pt) * ha - ev(dh, pt) * qa) % p
    er = (qa * hr - ha * qr) % p
    es = (qa * hs - ha * qs) % p
    ez = (qr * hs - qs * hr) % p
    ey = (er * bs - es * br + ba * ez) % p
    beta = -ev(Q1, pt) % p
    fy = 1 + s**6 + kappa * (s + 3 * s**5)
    fys = 6 * s**5 + kappa * (1 + 15 * s**4)
    ps = [
        (kappa * r * s, kappa * s, kappa * r),
        (2 * kappa * r * s, 2 * kappa * s, 2 * kappa * r),
        (fy, 0, fys),
    ]
    entries = [
        (
            alpha * (3 * v * beta * ey - beta**2 * (er * ds - es * dr))
            - 2 * i * v * beta**2 * ez
        )
        % p
        for i, (v, dr, ds) in enumerate(ps, 1)
    ]
    c0, c2, c4 = entries
    c2 += alpha**2 * c4
    c0 += alpha**2 * c2
    triple = [4 * ell * inv(3 * j) * x % p for x in (c0, c2, c4)]
    # Check against the current before the triangular division in t.
    U, V, Z = [ell * c * inv(j) % p for c in (er, es, ez)]
    pall = [
        (1 + r**6 + kappa * (r + 2 * r**5), 6 * r**5 + kappa * (1 + 10 * r**4), 0)
    ] + ps

    def pjets(t):
        return tuple(
            sum(row[i] * t ** (2 * q) for q, row in enumerate(pall)) % p
            for i in range(3)
        ) + (
            sum(2 * q * row[0] * t ** (2 * q - 1) for q, row in enumerate(pall) if q)
            % p,
        )

    P0, pr, ps0, pa = pjets(alpha)
    assert (P0 + beta**3) % p == 0 and P0
    T = (pa * Z + U * ps0 - V * pr) * inv(P0) % p
    for t in range(11):
        den = (t * t - alpha * alpha) % p
        if not den % 11:
            continue
        scalar = 2 * beta**2 * inv(den) % p
        qt, qr, qs = (-scalar * t * Z, scalar * alpha * V, -scalar * alpha * U)
        divq = -4 * beta**2 * alpha * T * inv(3 * den) % p
        P, Pr, Ps, Pt = pjets(t)
        density = (-P * divq + 2 * inv(3) * (qt * Pt + qr * Pr + qs * Ps)) % p
        assert density == (triple[0] + t * t * triple[1] + t**4 * triple[2]) % p
    print(
        "sign, alpha, jac-socle, target-socle, ell, j, E_r E_s E_Z E_Y, S, triple",
        sign,
        alpha,
        norm,
        res(vr),
        ell,
        j,
        (er, es, ez, ey),
        entries,
        triple,
    )
    return triple


if __name__ == "__main__":
    assert main(1) == main(-1) == [9, 9, 6]
    from fiber_lift import main as lift

    point = lift()
    modulus = 11**7
    a = main(1, point, modulus)
    b = main(-1, point, modulus)
    assert a == b == [13036912, 13328335, 93231] and [q % 11 for q in a] == [9, 9, 6]
    print("11-adic single-sheet current modulo", modulus, a)
