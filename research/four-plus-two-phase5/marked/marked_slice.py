"""Eliminate two axis coordinates and lambda on an entire coupled one-parameter slice."""

from itertools import product
from random import Random
from collections import Counter
import numpy as np
from marked_det import multiplication, inv_det
from marked_axis import geometry, seed, QMON
from marked_axis_fiber import axis


p = 11


# F121=F11[u]/(u^2+1), element a+b*u stored as a+11*b, 0 <= a,b < 11.
def elt(a, b=0):
    return a % p + p * (b % p)


def add(x, y):
    return elt(x % p + y % p, x // p + y // p)


def neg(x):
    return elt(-(x % p), -(x // p))


def mul(x, y):
    return elt(
        (x % p) * (y % p) - (x // p) * (y // p), (x % p) * (y // p) + (x // p) * (y % p)
    )


def power(x, n):
    z = 1
    while n:
        if n % 2:
            z = mul(z, x)
        x = mul(x, x)
        n //= 2
    return z


def value(coeffs, t):
    ans = 0
    for c in reversed(coeffs):
        ans = add(mul(ans, t), c)
    return ans


def strip(f):
    while f and not f[-1]:
        f.pop()
    return f


def monic(f):
    if not f:
        return []
    inv = power(f[-1], 119)
    return [mul(c, inv) for c in f]


def gcd(f, g):
    while g:
        g = monic(g)
        f = f.copy()
        while len(f) >= len(g):
            lead = f[-1]
            offset = len(f) - len(g)
            for i, c in enumerate(g):
                f[offset + i] = add(f[offset + i], neg(mul(lead, c)))
            strip(f)
        f, g = g, f
    return monic(f)


def slice_data(params):
    vals = np.array(
        [
            multiplication(
                p,
                p,
                4,
                2,
                [v + (i == 6) * t - (i == 8) * t for i, v in enumerate(params)],
            )
            for t in range(7)
        ]
    )[:, :, :, 0]
    V = np.array([[pow(t, n, p) for n in range(7)] for t in range(7)])
    iv, _ = inv_det(V, p, p)
    coeffs = np.tensordot(iv, vals, axes=1) % p
    for t in (7, 8, 9, 10):
        held = np.array(
            multiplication(
                p,
                p,
                4,
                2,
                [v + (i == 6) * t - (i == 8) * t for i, v in enumerate(params)],
            )
        )[:, :, 0]
        row = np.array([pow(t, n, p) for n in range(7)])
        assert np.array_equal(np.tensordot(row, coeffs, axes=1) % p, held)
    assert not np.any(coeffs[4:, 8])
    return coeffs


def main():
    b, Q, h = seed()
    params = b + [Q[k] for k in QMON]
    F, P = geometry(params, 4, 2, p)
    g = {k: v * pow(F[0, 0], -1, p) % p for k, v in F.items()}
    roots = []
    for w in (True, False):
        f, _ = axis(p, g, w)
        cf = [int(f.nth(n)) % p for n in range(9)]
        rr = [t for t in range(p * p) if value(cf, t) == 0]
        assert len(rr) == 8
        roots.append(rr)
        print(("w" if w else "z") + " axis roots as a+11b:", rr)
    cf = slice_data(params)
    sol = []
    deg = Counter()
    u_deg = Counter()
    zero_counts = Counter()
    polys = []
    for w, z in product(*roots):
        basis = np.array(
            [mul(power(w, i), power(z, j)) for i, j in product(range(8), repeat=2)]
        )
        vals = (cf @ (basis % p)) % p + p * ((cf @ (basis // p)) % p)
        row = [strip(list(vals[:, i])) for i in range(9)]
        polys.append(row)
        d = []
        for f in row:
            d = gcd(d, f)
        assert d
        deg[len(d) - 1] += 1
        if len(d) > 1:
            sol.append((w, z, d))
        U = []
        for k in range(7):
            a = 0
            for f in row:
                a = add(a, f[k] if k < len(f) else 0)
            U.append(a)
        U = strip(U)
        u_deg[len(U) - 1] += 1
        zero_counts[str(value(U, 0) == 0)] += 1
    print("common-gcd degree distribution", dict(deg), "sheets with positive gcd", sol)
    print(
        "U degree distribution", dict(u_deg), "U(0)=0 distribution", dict(zero_counts)
    )
    assert deg == {0: 63, 1: 1} and sol == [(6, 9, [0, 1])]
    assert zero_counts["True"] == 1

    # Univariate norm in a monic algebra over F121. The Euclidean recursion
    # uses Res(f,g)=(-1)^(mn) lc(g)^(m-r) Res(g,f mod g).
    def resultant(f, g):
        if not g:
            return 0
        m, n = len(f) - 1, len(g) - 1
        if n == 0:
            return power(g[0], m)
        if m < n:
            return mul(elt((-1) ** (m * n)), resultant(g, f))
        h = f.copy()
        base = monic(g)
        while len(h) >= len(g):
            offset = len(h) - len(g)
            lead = h[-1]
            for i, c in enumerate(base):
                h[i + offset] = add(h[i + offset], neg(mul(lead, c)))
            strip(h)
        if not h:
            return 0
        scalar = mul(elt((-1) ** (m * n)), power(g[-1], m - (len(h) - 1)))
        return mul(scalar, resultant(g, h))

    def lin(row, weights):
        out = []
        for k in range(7):
            y = 0
            for f, w in zip(row, weights):
                y = add(y, mul(w, f[k] if k < len(f) else 0))
            out.append(y)
        return strip(out)

    chosen = list(product(*roots)).index((6, 9))
    rng = Random(42)
    for trial in range(100):
        weights = [1] * 8 if trial == 0 else [rng.randrange(11) for _ in range(7)] + [1]
        pr = 1
        for i, row in enumerate(polys):
            U = monic(lin(row, [1] * 9))
            W = lin(row, weights)
            if i == chosen:
                assert U[0] == W[0] == 0 and U[1]
                U = U[1:]
            pr = mul(pr, resultant(U, W))
            if not pr:
                break
        if pr:
            break
    else:
        raise ValueError("no separating W found")
    assert 1 <= pr <= 10
    print(
        "second determinant: weights",
        weights,
        "rank 383, tr(adj(W))",
        int(pr),
        "mod 11",
    )


if __name__ == "__main__":
    main()
