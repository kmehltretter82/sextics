"""Enclose the source of the canonical current separately at each certified sheet."""

import argparse
from functools import reduce
from itertools import permutations
import json
from pathlib import Path
import numpy as np
from flint import acb, acb_mat, arb, ctx
import certify_roots as cr


def mons(n):
    return [
        (i, j, k, n - i - j - k)
        for i in range(n + 1)
        for j in range(n + 1 - i)
        for k in range(n + 1 - i - j)
    ]


def deriv(p, i):
    return {
        tuple(m[k] - int(i == k) for k in range(4)): m[i] * v
        for m, v in p.items()
        if m[i]
    }


def value(p, x):
    return sum(
        c * reduce(lambda a, b: a * b, (y**j for y, j in zip(x, m)), 1)
        for m, c in p.items()
    )


def det_poly(forms):
    out = {}
    for perm in permutations(range(4)):
        sign = (-1) ** sum(x > y for i, x in enumerate(perm) for y in perm[i + 1 :])
        for m, c in reduce(
            cr.times, [deriv(f, i) for f, i in zip(forms, perm)]
        ).items():
            out[m] = out.get(m, 0) + sign * c
    return out


def pivot_rows(rows):
    # Numerical row choices are only a heuristic. The interval solve below
    # verifies their invertibility, independently of all floating point choices.
    a = np.array([[complex(z.mid()) for z in row] for row in rows], dtype=complex)
    ids = np.arange(len(a))
    ip = 0
    for col in range(a.shape[1]):
        if ip == a.shape[1] - 1:
            break
        piv = ip + int(np.argmax(abs(a[ip:, col])))
        if abs(a[piv, col]) < 1e-10:
            continue
        a[[ip, piv]] = a[[piv, ip]]
        ids[[ip, piv]] = ids[[piv, ip]]
        a[ip + 1 :] -= (a[ip + 1 :, col] / a[ip, col])[:, None] * a[ip]
        ip += 1
    assert ip == a.shape[1] - 1
    return ids[:ip]


def source(box, r, s):
    _, J = cr.system(box, r, s)
    C = dict(zip(cr.q, box[:5]))
    Q = {cr.ab: acb(1), **dict(zip(cr.q, box[5:10]))}
    H = dict(zip(cr.aa + cr.bb, box[10:26]))
    G = dict(zip(cr.quart, box[26:45]))
    alpha = box[-1]
    x = (acb(1), r, alpha, s * alpha)

    def partial(f):
        db, dc, dd = (value(deriv(f, i), x) for i in (1, 2, 3))
        return db, alpha * dd, dc + s * dd

    qr, qs, qa = partial(Q)
    hr, hs, ha = partial(H)
    br, bs, ba = (-i for i in partial(C))
    tangent = J.solve(acb_mat([[0]] * 44 + [[-qr], [-hr]]))
    v2 = dict(zip(cr.q, [tangent[k, 0] for k in range(5, 10)]))
    v3 = dict(zip(cr.aa + cr.bb, [tangent[k, 0] for k in range(10, 26)]))
    forms = [Q, H, {k: 3 * v for k, v in cr.times(C, C).items()}, G]
    basis = mons(9)
    index = {m: i for i, m in enumerate(basis)}
    mac = []
    for form, d in zip(forms, (2, 3, 4, 4)):
        for m in mons(9 - d):
            row = [acb(0)] * len(basis)
            for e, c in form.items():
                row[index[tuple(i + j for i, j in zip(m, e))]] = c
            mac.append(row)
    rows = [mac[i] for i in pivot_rows(mac)]
    det = det_poly(forms)
    rows.append([det.get(e, acb(0)) for e in basis])
    soc = acb_mat(rows).solve(acb_mat([[0]] * (len(basis) - 1) + [[96]]))
    beta = -value(C, x)
    j = value(v2, x) * ha - value(v3, x) * qa
    er = qa * hr - ha * qr
    es = qa * hs - ha * qs
    ez = qr * hs - qs * hr
    ey = er * bs - es * br + ba * ez
    k = acb(5) / 7
    g = 1 + s**6 + k * (s + 3 * s**5)
    dg = 6 * s**5 + k * (1 + 15 * s**4)
    coeffs = [
        (k * r * s, k * s, k * r),
        (2 * k * r * s, 2 * k * s, 2 * k * r),
        (g, 0, dg),
    ]
    terms = [
        alpha * (3 * p * beta * ey - beta**2 * (er * ds - es * dr))
        - 2 * i * p * beta**2 * ez
        for i, (p, dr, ds) in enumerate(coeffs, 1)
    ]
    Pstar = (
        1
        + r**6
        + k * (r + 2 * r**5)
        + sum(p * alpha ** (2 * i) for i, (p, _, _) in enumerate(coeffs, 1))
    )
    assert (Pstar + beta**3).contains(0)
    raw = [
        terms[0] + alpha**2 * terms[1] + alpha**4 * terms[2],
        terms[1] + alpha**2 * terms[2],
        terms[2],
    ]
    sources = {}
    for m in mons(6):
        if (m[2] + m[3]) % 2 == 0:
            continue
        ell = sum(
            soc[index[tuple(i + j for i, j in zip(m, e))], 0] * c for e, c in v3.items()
        )
        sources[m] = [4 * ell * x / (3 * j) for x in raw]
        assert all(y.is_finite() for y in sources[m])
    return sources


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("endpoints", type=Path)
    p.add_argument("--limit", type=int, default=None)
    p.add_argument("--transport-test", action="store_true")
    args = p.parse_args()
    rows = json.loads(args.endpoints.read_text())[: args.limit]
    ctx.prec = 320
    r, s = acb(1) / 10, acb(1)
    step = acb(1) / 10**14
    direction = acb(arb(0, step.real))
    final, shifted = {}, {}
    for i, row in enumerate(rows):
        z = cr.center(row["coordinates"], r, s)
        box, _, _ = cr.root_box(row, r, s)
        results = source(box, r, s)
        for m, col in results.items():
            final[m] = [a + b for a, b in zip(final.get(m, (0, 0, 0)), col)]
        print("sheet", i, [v.str(9) for v in results[3, 0, 3, 0]], flush=True)
        if args.transport_test:
            outer = arb("1e-8")
            cr.contraction(z, r, s, outer, r + direction)
            trial = list(z)
            for _ in range(3):
                f, J = cr.system(trial, r + step, s)
                correction = J.solve(f)
                trial = [(v - correction[i, 0]).mid() for i, v in enumerate(trial)]
            tight = arb("1e-44")
            assert all((a - b).abs_upper() + tight < outer for a, b in zip(z, trial))
            moved, _, _ = cr.contraction(trial, r + step, s, tight, r + step)
            for m, col in source(moved, r + step, s).items():
                shifted[m] = [a + b for a, b in zip(shifted.get(m, (0, 0, 0)), col)]
    print(
        "unnormalized finite subsum for A:",
        [v.str(9) for v in final[3, 0, 3, 0]],
        flush=True,
    )
    from itertools import combinations

    for keys in combinations(final, 3):
        matrix = acb_mat([final[key] for key in keys]).transpose()
        det = matrix.det()
        if det.is_finite() and not det.contains(0):
            print("independent input exponents and subsum minor:", keys, det.str(9))
            target = acb_mat([[0], [0], [1]])
            weights = matrix.solve(target)
            print(
                "unique weights at first specialization:",
                [weights[k, 0].str(9) for k in range(3)],
            )
            if args.transport_test:
                error = (
                    acb_mat([shifted[key] for key in keys]).transpose() * weights
                    - target
                )
                assert any(not error[i, 0].contains(0) for i in range(3))
                print(
                    "same weights, continued sheets, defect at r+10^-14:",
                    [error[i, 0].str(9) for i in range(3)],
                )
            break


if __name__ == "__main__":
    main()
