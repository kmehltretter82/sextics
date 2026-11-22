"""Classify cube-contact lines meeting the line a=b=0 on the catalogue core."""

from math import factorial

import sympy as sp

s, r = sp.symbols("s r")
k = sp.Rational(5, 7)
g = 1 + s**6 + k * s + 3 * k * s**5
f = 1 + r**6 + k * r + 2 * k * r**5


def poly(expr, var=s):
    return sp.Poly(expr, var, domain=sp.QQ)


def numerator(expr, var=s):
    return poly(sp.cancel(expr).as_numer_denom()[0], var)


def noninvariant():
    q = [g.diff(s, j) / factorial(j) for j in range(7)]
    alpha = sp.cancel(q[1] / (3 * g))
    beta0 = sp.cancel(q[2] / (3 * g) - alpha**2)
    u = sp.cancel(1 - 2 * s * alpha)
    rhs = sp.cancel(g * alpha**3 + 6 * g * alpha * beta0 - q[3])
    assert numerator(u).gcd(numerator(rhs)).degree() == 0
    B = sp.cancel(rhs / u)
    beta = sp.cancel(beta0 + s * B / (3 * g))
    C = sp.cancel(3 * g * alpha * beta**2 - q[5])
    D = sp.cancel(g * beta**3 - 1)
    remaining = numerator(3 * g * (beta**2 + alpha**2 * beta) - q[4] - s * C)
    for exception in (numerator(g), numerator(u)):
        common = remaining.gcd(exception)
        while common.degree() > 0:
            remaining = remaining.exquo(common)
            common = remaining.gcd(exception)
    remaining = remaining.monic()
    assert remaining.degree() == 35
    assert all(remaining.gcd(numerator(p)).degree() == 0 for p in (B, C))

    def residue(expr):
        top, bottom = sp.cancel(expr).as_numer_denom()
        return (poly(top) * poly(bottom).invert(remaining)) % remaining

    bv, cv, dv = map(residue, (B, C, D))
    ratio = bv**2 * (4 * k * cv).invert(remaining) % remaining
    f_ratio = (poly(1) + ratio**6 + k * ratio + 2 * k * ratio**5) % remaining
    compatibility = (dv * bv**3 - 8 * cv**3 * f_ratio) % remaining
    assert remaining.gcd(compatibility).degree() == 0
    print(
        "Noninvariant lines: degree 35 after saturation, all compatibility gcds are one."
    )


def invariant():
    for binary, var in ((f, r), (g, s)):
        p = poly(binary, var)
        assert p.gcd(p.diff()).degree() == 0 and p.eval(0) == 1
    sigma = sp.cancel(6 * f / (k * r))
    fixed = numerator(g.subs(s, sigma) - 8 * f, r)
    assert fixed.degree() == 36 and fixed.gcd(fixed.diff()).degree() == 0
    assert fixed.gcd(poly(r * f, r)).degree() == 0
    odd = sp.Matrix(
        [
            [sp.diff(f, r), k * r, 1],
            [k * sigma, 2 * k * r, 4],
            [2 * k * sigma, sp.diff(g, s).subs(s, sigma), 4],
        ]
    )
    assert fixed.gcd(numerator(odd.det(), r)).degree() == 0
    print("Mixed invariant lines: all 36 points are reduced in the full Grassmannian.")


if __name__ == "__main__":
    noninvariant()
    invariant()
