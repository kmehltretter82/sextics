"""Eleven-variable marked norm-contact equations. Uses only Python's standard library."""

from itertools import product
from math import prod

QMON = [(0, 2), (1, 1), (2, 0), (2, 2)]
HMON = [m for m in product(range(4), repeat=2) if m != (0, 0)]
RES = [
    (i, j)
    for i in range(4, 7)
    for j in range(1, 7)
    if (i + j) % 2 == 0 and (i, j) != (6, 6)
]


def plus(mod, *pieces):
    ans = {}
    for coeff, poly in pieces:
        for m, value in poly.items():
            ans[m] = (ans.get(m, 0) + coeff * value) % mod
    return {m: v for m, v in ans.items() if v}


def mul(mod, *factors):
    ans = {(0, 0): 1}
    for factor in factors:
        new = {}
        for (i, j), a in ans.items():
            for (k, col), b in factor.items():
                pos = (i + k, j + col)
                new[pos] = (new.get(pos, 0) + a * b) % mod
        ans = {m: v for m, v in new.items() if v}
    return ans


def geometry(values, r, s, mod):
    x, y, v, j, Z, *q = values
    a = {(0, 0): x, (1, 1): 1 - x}
    b = {(0, 0): y, (1, 1): r - y}
    c = {(1, 0): v, (0, 1): 1 - v}
    d = {(1, 0): j, (0, 1): s - j}
    k = 5 * pow(7, -1, mod) % mod

    def binary(a, b, factor):
        return plus(
            mod,
            (1, mul(mod, *[a] * 6)),
            (1, mul(mod, *[b] * 6)),
            (k, mul(mod, *[a] * 5, b)),
            (factor * k, mul(mod, *[b] * 5, a)),
        )

    H = plus(
        mod,
        (1, binary(a, b, 2)),
        (Z**3, binary(c, d, 3)),
        (k * Z, mul(mod, a, a, a, b, c, d)),
        (2 * k * Z**2, mul(mod, a, b, c, c, c, d)),
    )
    P = sum(H.values()) % mod
    Q = {(0, 0): 1 - sum(q), **dict(zip(QMON, q))}
    return plus(mod, (1, H), (-P, mul(mod, Q, Q, Q))), P


def norm(h, i, j, mod):
    return (
        sum(
            a * (-1) ** (i - k + j - col) * h.get((i - k, j - col), 0)
            for (k, col), a in h.items()
            if k <= i and col <= j
        )
        % mod
    )


def determinant(matrix, mod):
    a = [row.copy() for row in matrix]
    result = 1
    for j in range(len(a)):
        k = next((i for i in range(j, len(a)) if a[i][j] % mod), None)
        if k is None:
            return 0
        if k != j:
            a[j], a[k] = a[k], a[j]
            result = -result
        result = result * a[j][j] % mod
        inv = pow(a[j][j], -1, mod)
        for i in range(j + 1, len(a)):
            q = a[i][j] * inv % mod
            a[i] = [(x - q * y) % mod for x, y in zip(a[i], a[j])]
    return result % mod


def det3(m):
    (a, b, c), (d, e, f), (g, h, i) = m
    return a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)


def solve3(matrix, rhs, mod):
    inverse = pow(det3(matrix), -1, mod)
    ans = []
    for col in range(3):
        m = [row.copy() for row in matrix]
        for i in range(3):
            m[i][col] = rhs[i]
        ans.append(det3(m) * inverse % mod)
    return ans


def reconstruct(values, r, s, mod):
    F, P = geometry(values[:9], r, s, mod)
    F0 = F.get((0, 0), 0)
    g = plus(mod, (pow(F0, -1, mod), F))
    half = pow(2, -1, mod)
    x, y, v, j, Z = values[:5]
    units = [Z, P, y - r * x, j - s * v, F0]
    h = {(0, 0): 1}
    res = []
    for direction, T in enumerate(values[9:]):

        def key(t):
            return (t, 0) if direction == 0 else (0, t)

        # Primitive coordinate T=V+E stays unramified at all four seed orderings.
        L = -g.get(key(6), 0) % mod
        A = (T * T + L + g.get(key(2), 0)) * half % mod
        D = 2 * T * (A + 1) % mod
        E = (A * A + T * T * L + 2 * L - g.get(key(4), 0)) * pow(D, -1, mod) % mod
        V = (T - E) % mod
        Y = (A - T * E) % mod
        res.append((E * E - L) % mod)
        h[key(1)], h[key(2)], h[key(3)] = V, Y, E
        units.append(D)
        if direction == 1:
            units.extend([E, E - V * Y])
    for i in [1, 2, 3]:
        degrees = [1, 3, 5] if i % 2 else [2, 4, 6]
        base = [norm(h, i, j, mod) for j in degrees]
        columns = []
        for col in [1, 2, 3]:
            h[i, col] = 1
            columns.append(
                [(norm(h, i, j, mod) - t) % mod for j, t in zip(degrees, base)]
            )
            del h[i, col]
        rhs = [(g.get((i, j), 0) - t) % mod for j, t in zip(degrees, base)]
        for col, t in enumerate(solve3(list(map(list, zip(*columns))), rhs, mod), 1):
            h[i, col] = t
    res.extend((norm(h, i, j, mod) - g.get((i, j), 0)) % mod for i, j in RES)
    res.append(sum(h.values()) % mod)
    return res, h, units, g


def seed():
    # The archived physical contact has kappa=7, alpha=6, and root directions
    # (2,6) and (1,7). Its old Segre mark is (2,10). Normalize a(mark)=1.
    orig = [4, 4, 9, 4, 6, 0, 6, 9, 9, 0, 10, 0, 7, 3, 9, 6]
    wp, zp, scale = 2, 10, 2
    h = {
        (i, j): orig[4 * i + j]
        * pow(wp, i, 11)
        * pow(zp, j, 11)
        * pow(scale**3, -1, 11)
        % 11
        for i, j in product(range(4), repeat=2)
    }
    q = [4, 7, 10, 8]
    Q = {(0, 0): (1 - sum(q)) % 11, **dict(zip(QMON, q))}
    return [6, 3, 2, 3, 3], Q, h


def swap(params, Q, h, xray, yray, p=11, r=4, s=2):
    x, y, v, j, Z = params
    if xray:
        x, y = 1 - x, r - y
    if yray:
        v, j = 1 - v, s - j

    def transform(poly, deg):
        out = {}
        for (i, j), t in poly.items():
            if xray:
                i, j = deg - j, deg - i
            if yray:
                i, j = j, i
            out[i, j] = t % p
        return out

    return [x % p, y % p, v % p, j % p, Z], transform(Q, 2), transform(h, 3)


def jet_det(values, r, s, p):
    mod = p * p
    baseline = reconstruct(values, r, s, mod)[0]
    columns = []
    for j in range(len(values)):
        moved = values.copy()
        moved[j] += p
        residual = reconstruct(moved, r, s, mod)[0]
        cols = [(a - b) % mod for a, b in zip(residual, baseline)]
        assert all(t % p == 0 for t in cols)
        columns.append([t // p for t in cols])
    return determinant(list(map(list, zip(*columns))), p)


def verify():
    p, r, s = 11, 4, 2
    args = seed()
    seen = set()
    for xray, yray in product([False, True], repeat=2):
        params, Q, h = swap(*args, xray, yray)
        u = h[0, 0]
        inv = pow(u, -1, p)
        values = (
            params
            + [Q[m] for m in QMON]
            + [(h[1, 0] + h[3, 0]) * inv % p, (h[0, 1] + h[0, 3]) * inv % p]
        )
        assert tuple(values) not in seen
        seen.add(tuple(values))
        residual, built, units, g = reconstruct(values, r, s, p)
        assert prod(units) % p and residual == [0] * 11
        assert all((built[m] * u - t) % p == 0 for m, t in h.items())
        assert all(
            (g.get((i, j), 0) - norm(built, i, j, p)) % p == 0
            for i, j in product(range(7), repeat=2)
        )
        det = jet_det(values, r, s, p)
        assert det
        print(f"order {int(xray)}{int(yray)}: variables {values}; h00 {u}; det {det}")
    print(
        "Four distinct simple roots; every reconstructed coefficient and all 49 norm slots agree."
    )


if __name__ == "__main__":
    verify()
