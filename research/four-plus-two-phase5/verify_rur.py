"""Verify Singular's binary contact algebra by direct substitution over Q."""

import re
from pathlib import Path
import sys
from flint import arb_poly, ctx, fmpq, fmpq_poly, fmpz_poly
from sympy import Matrix


def univariate(expr):
    out = [0] * 41
    for term in re.findall(r"[+-]?[^+-]+", expr):
        match = re.fullmatch(r"([+-]?)([0-9]*)(\*?m(?:\^([0-9]+))?)?", term)
        if match is None:
            raise ValueError(term[:100])
        sign, coef, var, pow_ = match.groups()
        out[int(pow_ or 1) if var else 0] += (-1 if sign == "-" else 1) * int(coef or 1)
    return fmpz_poly(out)


def read(path):
    raw = Path(path).read_text().splitlines()
    if raw[:4] != ["dimension and length", "0", "40", "lex basis"]:
        raise ValueError("Expected a zero-dimensional quotient of length 40")
    P = univariate(raw[4].split("=", 1)[1])
    polys = []
    for row, name in zip(raw[5:], ("v", "w", "iv"), strict=True):
        expr = row.split("=", 1)[1]
        c, rhs = expr.split("*" + name, 1)
        polys.append(fmpq_poly(-univariate(rhs), int(c)))
    return P, *polys


def verify_place(P, v, w, digits=32):
    """Lift the published unnormalized C,B seed and compare its normalized shape."""
    prime = 23
    seed = [15, 18, 16, 4, 6, 17, 22]

    def column(poly, shift):
        return ([0] * shift + list(poly) + [0] * 7)[:7]

    C, B = fmpz_poly(seed[:3]), fmpz_poly(seed[3:])
    jac = Matrix(
        [column(3 * C * C, i) for i in range(3)] + [column(-2 * B, i) for i in range(4)]
    ).T
    assert int(jac.det()) % prime != 0
    inverse = jac.inv_mod(prime)
    modulus = prime
    g7 = fmpz_poly([7, 5, 0, 0, 0, 15, 7])
    for _ in range(1, digits):
        C, B = fmpz_poly(seed[:3]), fmpz_poly(seed[3:])
        residual = 7 * (C**3 - B**2) - g7
        assert all(int(residual[i]) % modulus == 0 for i in range(7))
        error = Matrix(
            [
                (-int(residual[i]) // modulus) * pow(7, -1, prime) % prime
                for i in range(7)
            ]
        )
        seed = [x + modulus * (int(d) % prime) for x, d in zip(seed, inverse * error)]
        modulus *= prime
    C, B = fmpz_poly(seed[:3]), fmpz_poly(seed[3:])
    assert all(int(e) % modulus == 0 for e in 7 * (C**3 - B**2) - g7)
    m = pow(seed[0], 3, modulus)
    assert int(P(m)) % modulus == 0
    for shape, raw in zip((v, w), seed[1:3]):
        lhs = int(shape.numer()(m)) * seed[0]
        rhs = int(shape.denom()) * raw
        assert (lhs - rhs) % modulus == 0
    print(
        "verified selected 23-adic place through 23^%d; m mod 23^4 =" % digits,
        m % prime**4,
    )


def verify(path):
    P, v, w, iv = read(path)
    factors = P.factor()[1]
    print("shape degrees", [(x.degree(), n) for x, n in factors], flush=True)
    assert [(x.degree(), n) for x, n in factors] == [(40, 1)]
    verify_place(P, v, w)
    F = fmpq_poly(P)
    m = fmpq_poly([0, 1])
    K = lambda x: x % F
    one = fmpq_poly([1])
    assert K(iv * m * (m - 1) * w) == one
    cubic = [one]

    def product(a, b):
        res = [fmpq_poly() for _ in range(len(a) + len(b) - 1)]
        for i, val in enumerate(a):
            for j, val2 in enumerate(b):
                res[i + j] += val * val2
        return [K(y) for y in res]

    g = [one, one * fmpq(5, 7), 0, 0, 0, one * fmpq(15, 7), one]
    q = [one, v, w]
    f = [K(m * e - y) for e, y in zip(product(product(q, q), q), g)]
    content, inv, _ = (m - 1).xgcd(F)
    assert content == one
    for i in range(1, 4):
        prev = sum((cubic[j] * cubic[i - j] for j in range(1, i)), fmpq_poly())
        cubic.append(K((K(f[i] * inv) - prev) / 2))
    assert all(not K(f[i] - (m - 1) * p) for i, p in enumerate(product(cubic, cubic)))
    print(
        "verified all seven coefficients of C^3-B^2=g and the open-chart inverse",
        flush=True,
    )
    ctx.prec = 1024
    roots = P.complex_roots()
    real = [a for a, n in roots if a.imag.is_zero()]
    print("real embeddings:", len(real), "of", len(roots), flush=True)
    print("real m:", [float(a.real) for a in real], flush=True)
    (selected,) = [a.real for a in real if 1 < a.real < 2]
    u = selected.root(3)
    b = (selected - 1).sqrt()
    for name, values in [
        ("C", [u, u * arb_poly(v)(selected), u * arb_poly(w)(selected)]),
        ("B", [b * arb_poly(coefficient)(selected) for coefficient in cubic]),
    ]:
        print(name, [value.str(55) for value in values])
    return P, v, w, cubic


if __name__ == "__main__":
    verify(sys.argv[1])
