"""Certify isolated complex solutions of the lifted (44+2)-equation system."""

import argparse
from itertools import product
import json
from pathlib import Path
from flint import acb, acb_mat, arb, ctx


def exps(d):
    return [
        e
        for e in product(range(d + 1), repeat=4)
        if sum(e) == d and sum(e[2:]) % 2 == 0
    ]


q = [(2, 0, 0, 0), (0, 2, 0, 0), (0, 0, 2, 0), (0, 0, 1, 1), (0, 0, 0, 2)]
ab = (1, 1, 0, 0)
aa = [
    (3, 0, 0, 0),
    (0, 3, 0, 0),
    (1, 0, 2, 0),
    (1, 0, 1, 1),
    (1, 0, 0, 2),
    (0, 1, 2, 0),
    (0, 1, 1, 1),
    (0, 1, 0, 2),
]
bb = [
    (2, 0, 1, 0),
    (2, 0, 0, 1),
    (0, 2, 1, 0),
    (0, 2, 0, 1),
    (0, 0, 3, 0),
    (0, 0, 2, 1),
    (0, 0, 1, 2),
    (0, 0, 0, 3),
]
quart = exps(4)
six = exps(6)
index = {e: i for i, e in enumerate(six)}
piv = [i for i, e in enumerate(six) if e[0] and e[1]]
assert len(piv) == 19


def times(f, g):
    out = {}
    for m, x in f.items():
        for n, y in g.items():
            e = tuple(i + j for i, j in zip(m, n))
            out[e] = out.get(e, 0) + x * y
    return out


def system(z, r, s):
    C = dict(zip(q, z[:5]))
    Q = {ab: 1, **dict(zip(q, z[5:10]))}
    A = dict(zip(aa, z[10:18]))
    B = dict(zip(bb, z[18:26]))
    G = dict(zip(quart, z[26:45]))
    alpha = z[-1]
    k = acb(5) / 7
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
    square = times(C, C)
    val = [H.get(e, 0) for e in six] + [acb(0), acb(0)]
    for h, sign in [
        (times(C, square), -1),
        (times(A, A), -1),
        (times(B, B), 1),
        (times(Q, G), -1),
    ]:
        for e, x in h.items():
            val[index[e]] += sign * x
    J = acb_mat(46, 46)
    col = 0
    for mons, h, sign in [
        (q, square, -3),
        (q, G, -1),
        (aa, A, -2),
        (bb, B, 2),
        (quart, Q, -1),
    ]:
        for m in mons:
            for n, x in h.items():
                e = tuple(i + j for i, j in zip(m, n))
                J[index[e], col] = sign * x
            col += 1
    for i, form, mons, offset in [(44, Q, q, 5), (45, {**A, **B}, aa + bb, 10)]:
        for m, coeff in form.items():
            _, b, c, d = m
            base = r**b * s**d
            value = base * alpha ** (c + d)
            if m in mons:
                J[i, offset + mons.index(m)] = value
            val[i] += coeff * value
            if c + d:
                J[i, 45] += coeff * (c + d) * base * alpha ** (c + d - 1)
    return acb_mat([[v] for v in val]), J


def center(raw, r, s):
    y = [acb(arb(a), arb(b)).mid() for a, b in raw]
    z = y[:-1] + [acb(0)] * 19 + y[-1:]
    f, j = system(z, r, s)
    mat = acb_mat([[j[i, k] for k in range(26, 45)] for i in piv])
    g = mat.solve(acb_mat([[-f[i, 0]] for i in piv]))
    return z[:26] + [g[i, 0].mid() for i in range(19)] + z[-1:]


def rowsum(m):
    return max(
        sum(m[i, j].abs_upper() for j in range(m.ncols())).upper()
        for i in range(m.nrows())
    )


def root_box(record, r, s):
    """Enclose the unique fixed point in the complex sup-norm ball."""
    z = center(record["coordinates"], r, s)
    return contraction(z, r, s, arb("1e-44"), r)


def contraction(z, r, s, radius, r_span):
    margin = acb(arb(0, radius), arb(0, radius))
    box = [v + margin for v in z]
    I = acb_mat([[int(i == j) for j in range(46)] for i in range(46)])
    f, j = system(z, r, s)
    R = j.inv().mid()
    f, _ = system(z, r_span, s)
    _, boxj = system(box, r_span, s)
    L = rowsum(I - R * boxj)
    eps = rowsum(R * f)
    # The enclosing product of real intervals includes the complex balls.
    # ||I-RJ||<1 proves invertibility of both the midpoint R and each J.
    assert L < 1 and eps + L * radius < radius
    return box, eps, L


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("endpoints", type=Path)
    data = json.loads(p.parse_args().endpoints.read_text())
    ctx.prec = 320
    r, s = acb(1) / 10, acb(1)
    alphas = []
    for row in data:
        box, eps, L = root_box(row, r, s)
        assert all(not box[-1].overlaps(a) for a in alphas)
        alphas.append(box[-1])
        print(
            row["alpha"][0][:13],
            row["alpha"][1][:13],
            "Newton displacement",
            float(eps),
            "Lipschitz bound",
            float(L),
            flush=True,
        )
    print(len(alphas), "pairwise disjoint certified roots, complex sup radius 1e-44")


if __name__ == "__main__":
    main()
