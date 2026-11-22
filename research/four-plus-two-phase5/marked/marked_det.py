"""Nine-equation determinant oracle on the full, unreduced 8-by-8 axis chart."""

from argparse import ArgumentParser
import numpy as np
import sympy as sp
from marked_axis import geometry, seed, QMON, RES, determinant


def edge(p, mod, g, first):
    t = sp.Symbol("t")

    def poly(a):
        return sp.Poly(a, t, domain=sp.ZZ)

    T = poly(t)

    def key(n):
        return (n, 0) if first else (0, n)

    L = -g.get(key(6), 0) % mod
    A = (T * T + L + g.get(key(2), 0)) * pow(2, -1, mod)
    C = A * A + T * T * L + 2 * L - g.get(key(4), 0)
    D = 2 * T * (A + 1)

    def reduce_coeff(f):
        return sp.Poly.from_dict(
            {e: int(c) % mod for e, c in f.terms()}, (t,), domain=sp.ZZ
        )

    f = reduce_coeff(16 * (C * C - L * D * D))
    assert f.degree() == 8 and f.LC() == 1

    def rem(a):
        return reduce_coeff(a.rem(f))

    ff = sp.Poly(f, t, modulus=p)
    if ff.gcd(ff.diff()).degree() > 0:
        raise ValueError("axis is ramified")

    def inv(a):
        b = poly(sp.Poly(a, t, modulus=p).invert(ff).as_expr())
        b = rem(b * (2 - a * b))
        assert rem(b * a) == 1
        return b

    E = rem(C * inv(D))
    V = rem(T - E)
    Y = rem(A - T * E)
    vals = {"V": V, "Y": Y, "E": E}
    if not first:
        vals.update(ie=inv(E), idiff=inv(E - V * Y))

    def arr(h):
        return np.array([int(h.nth(i)) % mod for i in range(8)], dtype=np.int64)

    return arr(f), {k: arr(v) for k, v in vals.items()}


def multiplication(p, mod, r, s, params):
    F, P = geometry(params, r, s, mod)
    f0 = F.get((0, 0), 0)
    x, y, v, j, Z = params[:5]
    if not (Z * P * (y - r * x) * (j - s * v) * f0) % p:
        raise ValueError("outside moving mark/open")
    g = {m: u * pow(f0, -1, mod) % mod for m, u in F.items()}
    a, wx = edge(p, mod, g, True)
    b, zx = edge(p, mod, g, False)
    zero = np.zeros((8, 8), dtype=np.int64)
    one = zero.copy()
    one[0, 0] = 1

    def W(f):
        return (
            np.vstack((np.zeros((1, 8), dtype=np.int64), f[:-1])) - a[:, None] * f[-1]
        ) % mod

    def T(f):
        return (
            np.column_stack((np.zeros(8, dtype=np.int64), f[:, :-1]))
            - f[:, -1, None] * b
        ) % mod

    def mul(f, h):
        out = zero.copy()
        fh = f
        for row in h:
            ft = fh
            for coeff in row:
                if coeff:
                    out += coeff * ft
                ft = T(ft)
            fh = W(fh)
        return out % mod

    def regular(f):
        cols = []
        fw = f
        for i in range(8):
            ft = fw
            for j in range(8):
                cols.append(ft.flatten())
                ft = T(ft)
            fw = W(fw)
        return np.array(cols).T

    def embed(v, dir):
        x = zero.copy()
        if dir == 0:
            x[:, 0] = v
        else:
            x[0, :] = v
        return x

    h = {(0, 0): one}
    for n, k in enumerate(("V", "Y", "E"), 1):
        h[n, 0] = embed(wx[k], 0)
        h[0, n] = embed(zx[k], 1)
    V, Y, E, idiff, ie = [embed(zx[k], 1) for k in ("V", "Y", "E", "idiff", "ie")]

    def norm(i, j):
        return (
            sum(
                (
                    (-1) ** (i - k + j - col) * mul(a, h[i - k, j - col])
                    for (k, col), a in h.items()
                    if (i - k, j - col) in h
                ),
                zero,
            )
            % mod
        )

    half = pow(2, -1, mod)
    for i in (1, 2, 3):
        js = [1, 3, 5] if i % 2 else [2, 4, 6]
        d, e, f = [(half * (g.get((i, j), 0) * one - norm(i, j))) % mod for j in js]
        if i % 2:
            h1 = d
            m = (e - mul(Y, h1)) % mod
            h2 = mul(mul(Y, m) - f, idiff)
            h3 = (m + mul(V, h2)) % mod
        else:
            h3 = -mul(f, ie) % mod
            h1 = mul(mul(Y, d) - mul(V, h3) - e, idiff)
            h2 = (d + mul(V, h1)) % mod
        h[i, 1], h[i, 2], h[i, 3] = h1, h2, h3
    return [regular((norm(i, j) - g.get((i, j), 0) * one) % mod) for i, j in RES] + [
        regular(sum(h.values(), zero) % mod)
    ]


def stable_minor(A, p):
    U = A.copy() % p
    rows = list(range(len(U)))
    cols = []
    n = 0
    for j in range(len(U)):
        k = next((k for k in range(n, len(U)) if U[k, j] % p), None)
        if k is None:
            continue
        U[[n, k]] = U[[k, n]]
        rows[n], rows[k] = rows[k], rows[n]
        factor = (U[n + 1 :, j] * pow(int(U[n, j]), -1, p)) % p
        U[n + 1 :] = (U[n + 1 :] - factor[:, None] * U[n]) % p
        cols.append(j)
        n += 1
    return rows[:n], cols


def inv_det(A, p, mod):
    n = len(A)
    a = np.concatenate((A.copy(), np.eye(n, dtype=np.int64)), axis=1)
    det = 1
    for j in range(n):
        k = next(i for i in range(j, n) if a[i, j] % p)
        pivot = int(a[k, j])
        if k != j:
            a[[k, j]] = a[[j, k]]
            det = -det
        det = det * pivot % mod
        a[j] = a[j] * pow(pivot, -1, mod) % mod
        coef = a[:, j].copy()
        coef[j] = 0
        a = (a - coef[:, None] * a[j]) % mod
    return a[:, n:], det


def oracle(mats, weights, p, mod, minor):
    n = len(mats[0])
    U = (
        sum((c * M for c, M in zip(weights, mats)), np.zeros((n, n), dtype=np.int64))
        % mod
    )
    rows, cols = minor
    if len(rows) != n - 1:
        raise ValueError("rank is not n-1 on the special fiber")
    ri = next(i for i in range(n) if i not in rows)
    ci = next(i for i in range(n) if i not in cols)
    # Stable row pivots are not necessarily sorted. Sorting fixes the inverse permutation sign.
    rows = sorted(rows) + [ri]
    cols = sorted(cols) + [ci]
    sgn = (-1) ** (ri + ci)

    def reorder(M):
        return M[np.ix_(rows, cols)]

    U = reorder(U)
    inv, det = inv_det(U[:-1, :-1], p, mod)
    left = U[-1, :-1] @ inv % mod
    right = inv @ U[:-1, -1] % mod
    gamma = (U[-1, -1] - left @ U[:-1, -1]) % mod

    def coefficient(M):
        M = reorder(M)
        C = M[:-1, :-1]
        derivative = (
            M[-1, -1] - left @ M[:-1, -1] - M[-1, :-1] @ right + left @ C @ right
        ) % mod
        return (sgn * det * (gamma * np.sum(inv * C.T) + derivative)) % mod

    ell = coefficient(np.eye(n, dtype=np.int64))
    degree0 = sgn * det * gamma % mod
    # U has nonzero last weight, hence U,R0,...,R7 form an equivalent set of generators.
    return [int(degree0)] + [int(coefficient(M)) for M in mats[:-1]], int(ell)


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--weights", type=int, nargs=9)
    args = parser.parse_args()
    p, r, s = 11, 4, 2
    b, Q, h = seed()
    params = b + [Q[k] for k in QMON]
    M = multiplication(p, p, r, s, params)
    trial = (
        [args.weights]
        if args.weights
        else [[1] * 9, list(range(1, 10)), [1, 0, 0, 0, 0, 0, 0, 0, 1]]
    )
    for weights in trial:
        if not weights[-1] % p:
            continue
        U = sum(w * f for w, f in zip(weights, M)) % p
        minor = stable_minor(U, p)
        if len(minor[0]) == 63:
            vals, t = oracle(M, weights, p, p, minor)
            if t % p:
                break
    else:
        raise ValueError("no proposed weight gives the required simple zero")
    assert vals == [0] * 9
    print("weights", weights, "rank 63 trace(adj(U))", t, "equations", vals, flush=True)
    mod = p * p
    values, _ = oracle(multiplication(p, mod, r, s, params), weights, p, mod, minor)
    cols = []
    for i in range(9):
        moved = params.copy()
        moved[i] += p
        val, _ = oracle(multiplication(p, mod, r, s, moved), weights, p, mod, minor)
        diff = [(a - b) % mod for a, b in zip(val, values)]
        assert all(y % p == 0 for y in diff)
        cols.append([y // p for y in diff])
    jac = list(map(list, zip(*cols)))
    det = determinant(jac, p)
    assert det
    print("geometry point", params, "determinant Jacobian", det, "mod", p)
    for row in jac:
        print(*row)


if __name__ == "__main__":
    main()
