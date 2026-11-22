"""Eliminate the two axis parameters exactly over a fixed marked geometry, over F_p."""

from argparse import ArgumentParser
from itertools import product
import sympy as sp
from sympy.polys.rings import ring
from marked_axis import geometry, QMON, RES, seed, swap


def remove_nonunits(f, bad):
    while f.degree() > 0:
        g = f.gcd(bad)
        if g.degree() == 0:
            break
        f = f.exquo(g)
    return f.monic()


def axis(p, g, first):
    t = sp.Symbol("t")

    def poly(a):
        return sp.Poly(a, t, modulus=p)

    def key(n):
        return (n, 0) if first else (0, n)

    T = poly(t)
    L = -int(g.get(key(6), 0)) % p
    A = (T * T + L + g.get(key(2), 0)) * pow(2, -1, p)
    C = A * A + T * T * L + 2 * L - g.get(key(4), 0)
    D = 2 * T * (A + 1)
    f = remove_nonunits((C * C - L * D * D).monic(), D)
    if f.degree() == 0:
        return f, {}
    E = (C * D.invert(f)) % f
    V = (T - E) % f
    Y = (A - T * E) % f
    vals = {"E": E, "V": V, "Y": Y}
    if not first:
        f = remove_nonunits(f, E * (E - V * Y))
        if f.degree() > 0:
            vals.update(recip_E=E.invert(f), recip_diff=(E - V * Y).invert(f))
    return f, vals


def solve(p, r, s, params):
    F, P = geometry(params, r, s, p)
    x, y, v, j, Z = params[:5]
    F0 = F.get((0, 0), 0)
    if not (Z * P * (y - r * x) * (j - s * v) * F0) % p:
        raise ValueError("geometry outside marked axes chart")
    g = {pos: (f * pow(F0, -1, p)) % p for pos, f in F.items()}
    aw, wx = axis(p, g, True)
    az, ax = axis(p, g, False)
    dw, dz = int(aw.degree()), int(az.degree())
    if not dw or not dz:
        return (dw, dz), 0, None
    R, W, T = ring("W,T", sp.GF(p))
    coord = [W, T]

    def as_axis(poly, dim):
        return sum((int(c) * coord[dim] ** m[0] for m, c in poly.terms()), R.zero)

    rels = [as_axis(aw, 0), as_axis(az, 1)]

    def red(val):
        return R(val).rem(rels)

    ax = {k: as_axis(v, 1) for k, v in ax.items()}
    wx = {k: as_axis(v, 0) for k, v in wx.items()}
    h = {(0, 0): R.one}
    for n, name in enumerate(("V", "Y", "E"), 1):
        h[n, 0] = wx[name]
        h[0, n] = ax[name]
    V, Y, E = [ax[c] for c in ("V", "Y", "E")]
    half = pow(2, -1, p)

    def norm(i, j):
        return red(
            sum(
                (
                    a * (-1) ** (i - k + j - col) * h.get((i - k, j - col), 0)
                    for (k, col), a in h.items()
                    if k <= i and col <= j
                ),
                R.zero,
            )
        )

    for i in (1, 2, 3):
        cols = [1, 3, 5] if i % 2 else [2, 4, 6]
        b, c, d = [red(half * (g.get((i, j), 0) - norm(i, j))) for j in cols]
        if i % 2:
            h1 = b
            t = red(c - Y * h1)
            h2 = red((Y * t - d) * ax["recip_diff"])
            h3 = red(t + V * h2)
        else:
            h3 = red(-d * ax["recip_E"])
            h1 = red((Y * b - V * h3 - c) * ax["recip_diff"])
            h2 = red(b + V * h1)
        h[i, 1], h[i, 2], h[i, 3] = h1, h2, h3
    residual = [red(norm(i, j) - g.get((i, j), 0)) for i, j in RES] + [
        red(sum(h.values(), R.zero))
    ]
    basis = [(i, j) for i in range(dw) for j in range(dz)]
    piv = {}
    for f in residual:
        for i, j in basis:
            v = red(f * W**i * T**j)
            row = [int(v.get(m, 0)) % p for m in basis]
            for k in sorted(piv):
                if row[k]:
                    row = [(a - row[k] * b) % p for a, b in zip(row, piv[k])]
            k = next((k for k, c in enumerate(row) if c), None)
            if k is not None:
                piv[k] = [(a * pow(row[k], -1, p)) % p for a in row]
        if len(piv) == len(basis):
            break
    dim = len(basis) - len(piv)
    point = None
    if dim == 1:

        def normal(f):
            row = [int(f.get(m, 0)) % p for m in basis]
            for k in sorted(piv):
                if row[k]:
                    row = [(a - row[k] * b) % p for a, b in zip(row, piv[k])]
            return next(row[k] for k in range(len(basis)) if k not in piv)

        unit = pow(normal(R.one), -1, p)
        point = (normal(red(W)) * unit % p, normal(red(T)) * unit % p)
    return (dw, dz), dim, point


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--p", type=int)
    parser.add_argument("--r", type=int)
    parser.add_argument("--s", type=int)
    parser.add_argument("--geometry", type=int, nargs=9)
    args = parser.parse_args()
    if args.p:
        assert args.r is not None and args.s is not None and args.geometry
        print(solve(args.p, args.r, args.s, args.geometry))
        return
    for e, f in product([False, True], repeat=2):
        b, Q, h = swap(*seed(), e, f)
        params = b + [Q[k] for k in QMON]
        deg, dim, point = solve(11, 4, 2, params)
        truth = tuple(
            (
                h[(i, 0)]
                + h[(3 * i, 0)]
                + h[(0, 1 - i)]
                + h[(0, 3 * (1 - i))]
                - 2 * h[(0, 0)]
            )
            * pow(h[0, 0], -1, 11)
            % 11
            for i in (1, 0)
        )
        assert (dim, point) == (1, truth)
        print(
            f"ordering {int(e)}{int(f)}: axis degrees {deg}, quotient dimension {dim}, reduced point {point}"
        )
    b, Q, h = seed()
    params = b + [Q[k] for k in QMON]
    params[7] = (params[7] + 1) % 11
    deg, dim, point = solve(11, 4, 2, params)
    assert dim == 0
    print("altered Q20 control:", deg, "quotient dimension", dim)


if __name__ == "__main__":
    main()
