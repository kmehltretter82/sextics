"""Construct the marked split-quadric contact system and verify the coupled F11 seed."""

from argparse import ArgumentParser
from collections import Counter
from itertools import product
from pathlib import Path
from sympy.polys.rings import ring
from sympy import GF

QMON = [(0, 2), (1, 1), (2, 0), (2, 2)]
HMON = [m for m in product(range(4), repeat=2) if m != (0, 0)]
ROWS = [m for m in product(range(7), repeat=2) if sum(m) % 2 == 0 and m != (6, 6)]
NAMES = (
    ["x", "y", "v", "j", "Z"]
    + [f"Q{i}{k}" for i, k in QMON]
    + [f"h{i}{k}" for i, k in HMON]
)


def build(p, r, s):
    R, *gens = ring(",".join(NAMES + ["w", "z"]), GF(p))
    x, y, v, j, Z, *extra = gens
    w, z = extra[-2:]
    k = 5 * pow(7, -1, p) % p
    a = x + (1 - x) * w * z
    b = y + (r - y) * w * z
    c = v * w + (1 - v) * z
    d = j * w + (s - j) * z

    def f(a, b, t):
        return a**6 + b**6 + k * (a**5 * b + t * a * b**5)

    marked = (
        f(a, b, 2) + Z**3 * f(c, d, 3) + k * Z * a * b * c * d * (a * a + 2 * Z * c * c)
    )
    P = f(1, r, 2) + Z**3 * f(1, s, 3) + k * r * s * Z * (1 + 2 * Z)
    Q = 1 + sum(q * (w**i * z**k - 1) for q, (i, k) in zip(gens[5:9], QMON))
    h = sum(q * (w**i * z**k - 1) for q, (i, k) in zip(gens[9:24], HMON))
    conjugate = sum(
        q * ((-1) ** (i + k) * w**i * z**k - 1) for q, (i, k) in zip(gens[9:24], HMON)
    )
    identity = marked - P * Q**3 - h * conjugate
    by_pair = {m: {} for m in ROWS + [(6, 6)]}
    for powers, coeff in identity.terms():
        label = powers[-2:]
        assert label in by_pair
        by_pair[label][powers[:-2]] = int(coeff) % p
    accum = {}
    for row in by_pair.values():
        for mon, c in row.items():
            accum[mon] = (accum.get(mon, 0) + c) % p
    assert not any(accum.values())
    return [by_pair[m] for m in ROWS]


def determinant(matrix, p):
    a = [row.copy() for row in matrix]
    result = 1
    for j in range(len(a)):
        k = next((i for i in range(j, len(a)) if a[i][j] % p), None)
        if k is None:
            return 0
        if k != j:
            a[j], a[k] = a[k], a[j]
            result = -result
        result = result * a[j][j] % p
        inv = pow(a[j][j], -1, p)
        for i in range(j + 1, len(a)):
            q = a[i][j] * inv % p
            if q:
                a[i] = [(x - q * y) % p for x, y in zip(a[i], a[j])]
    return result % p


def seed():
    p = 11
    r = 4
    s = 2
    k = 7
    m, n, u, v, e = 2, 6, 1, 7, 7
    wp, zp = 2, 10
    K = e * (v - u) ** 2 * (m + n) * pow((n - m) ** 2, -1, p) % p
    scale = (1 + K * wp * zp) % p
    cs = (wp + zp) % p
    alpha = cs * pow(scale, -1, p) % p
    Z = alpha * alpha % p
    x = pow(scale, -1, p)
    y = n * x % p
    vnew = wp * pow(cs, -1, p) % p
    j = v * wp * pow(cs, -1, p) % p
    P = (
        1
        + r**6
        + k * r
        + 2 * k * r**5
        + k * r * s * Z * (1 + 2 * Z)
        + (1 + s**6 + k * s + 3 * k * s**5) * Z**3
    ) % p
    assert Z == 3 and P == 6 and (y - r * x) * (j - s * vnew) % p
    R, w, z = ring("w,z", GF(p))
    a = x + (1 - x) * w * z
    b = y + (r - y) * w * z
    c = alpha * (vnew * w + (1 - vnew) * z)
    d = alpha * (j * w + (s - j) * z)
    Q = (4 * a * a + 8 * b * b + 4 * c * c + 3 * c * d) * pow(8, -1, p)
    coeffs = dict(Q.terms())
    assert sum(map(int, coeffs.values())) % p == 1
    qco = [int(coeffs.get(pair, 0)) % p for pair in QMON]
    old = [4, 4, 9, 4, 6, 0, 6, 9, 9, 0, 10, 0, 7, 3, 9, 6]
    hco = [
        old[4 * i + j] * pow(wp, i, p) * pow(zp, j, p) * pow(scale**3, -1, p) % p
        for i, j in HMON
    ]
    assert (sum(hco) + old[0] * pow(scale**3, -1, p)) % p == 0
    return [x, y, vnew, j, Z] + qco + hco


def val_jac(rows, values, p):
    vals = []
    jac = []
    for row in rows:
        val = 0
        gradient = [0] * len(values)
        for powers, c in row.items():
            t = c
            for a, i in zip(values, powers):
                t = t * pow(a, i, p) % p
            val = (val + t) % p
            for k, e in enumerate(powers):
                if not e:
                    continue
                t = c * e
                for j, (a, i) in enumerate(zip(values, powers)):
                    t = t * pow(a, i - int(k == j), p) % p
                gradient[k] = (gradient[k] + t) % p
        vals.append(val)
        jac.append(gradient)
    return vals, jac


def m4row(row):
    output = []
    for exps, c in sorted(row.items()):
        terms = [str(c)] if c != 1 or not any(exps) else []
        terms.extend(
            f"X{i}" + (f"^{n}" if n > 1 else "") for i, n in enumerate(exps) if n
        )
        output.append("*".join(terms))
    return "+".join(output) or "0"


def emit(rows, p, r, s, where):
    # One inverse each for Z, P(Z), and the two basis-change determinants.
    k = 5 * pow(7, -1, p) % p
    coef = [
        (1 + r**6 + k * r + 2 * k * r**5) % p,
        k * r * s % p,
        2 * k * r * s % p,
        (1 + s**6 + k * s + 3 * k * s**5) % p,
    ]

    def unit(terms, idx):
        terms = [(([(i, d) for i, d in m] + [(idx, 1)]), c) for m, c in terms]
        row = {
            tuple(dict(m).get(i, 0) for i in range(28)): c % p
            for m, c in terms
            if c % p
        }
        row[(0,) * 28] = p - 1
        return m4row(row)

    inv_rows = [
        unit([([(4, 1)], 1)], 24),
        unit([([] if n == 0 else [(4, n)], c) for n, c in enumerate(coef)], 25),
        unit([([(1, 1)], 1), ([(0, 1)], -r)], 26),
        unit([([(3, 1)], 1), ([(2, 1)], -s)], 27),
    ]
    where.write_text(
        "\n".join(
            ["$fieldsize " + str(p), "$vars 28 X"]
            + [m4row(row) for row in rows]
            + inv_rows
        )
        + "\n"
    )


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--p", type=int, default=11)
    parser.add_argument("--r", type=int, default=4)
    parser.add_argument("--s", type=int, default=2)
    parser.add_argument("--m4", type=Path)
    args = parser.parse_args()
    rows = build(args.p, args.r, args.s)
    stats = Counter(max(map(sum, row)) for row in rows)
    print(
        "equations",
        len(rows),
        "terms",
        sum(map(len, rows)),
        "degree histogram",
        dict(sorted(stats.items())),
        flush=True,
    )
    if (args.p, args.r, args.s) == (11, 4, 2):
        values = seed()
        res, jac = val_jac(rows, values, args.p)
        det = determinant(jac, args.p)
        assert res == [0] * 24 and det != 0
        print("seed coordinates", values, "Jacobian determinant", det, flush=True)
    if args.m4:
        emit(rows, args.p, args.r, args.s, args.m4)


if __name__ == "__main__":
    main()
