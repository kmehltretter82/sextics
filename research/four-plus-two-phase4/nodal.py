"""Check the three local sheets' boundary functional over F_23[v]/(v²-15)."""

import sympy as sp

v, s = sp.symbols("v s")
p = 23
d = 15


def pair(expr):
    """Evaluate a quotient in F_p[v]/(v²-d)."""
    numerator, denominator = sp.fraction(sp.cancel(expr))

    def reduce(poly):
        a = sp.rem(poly, v * v - d, v).expand()
        assert sp.degree(a, v) <= 1
        return tuple(int(a.coeff(v, i).p) * pow(int(a.coeff(v, i).q), -1, p) % p for i in range(2))

    a, b = reduce(numerator)
    c, e = reduce(denominator)
    norm = (c * c - d * e * e) % p
    assert norm
    return ((a * c - d * b * e) * pow(norm, -1, p) % p, (b * c - a * e) * pow(norm, -1, p) % p)


B = 4 + 6 * s + 17 * s * s + 22 * s**3
C = 15 + 18 * s + 16 * s * s
g = 1 + s**6 + 4 * s + 12 * s**5
U = 8 + 21 * v
V = (18 + 5 * v) + (17 + 12 * v) * s + (20 + 6 * v) * s * s


def main():
    assert all(c % p == 0 for c in sp.Poly(C**3 - B**2 - g, s).coeffs())
    shape = -sp.Rational(4, 3) * C * (C * g.diff(s) - 3 * g * C.diff(s)) * V * (U + V / 2)
    denom = B * (B**2 * U**3 + V**3)
    results = []
    for S in [1, 2, 3, 4]:

        def val(f):
            return f.subs(s, S)

        # An independent direct companion-matrix evaluation of the geometric source S3/j:
        # alpha^3=-1/B, and, at r=0, q_alpha=0, q_r=1, q_s=0, j=v2*h_alpha.
        A = sp.Matrix([[0, 0, -1 / val(B)], [1, 0, 0], [0, 1, 0]])
        E = sp.eye(3)
        beta = -(A**2) * val(C)
        ha = 3 * A**2 * val(B)
        Er = -ha
        EZ = A**3 * val(B.diff(s))
        EY = Er * (-(A**2) * val(C.diff(s))) + (-2 * A * val(C)) * EZ
        S3 = (
            A * (3 * val(g) * beta * EY - beta**2 * Er * val(g.diff(s))) - 6 * val(g) * beta**2 * EZ
        )
        j = (U * E + val(V) * A**2) * ha
        dtrace = sp.Rational(4, 9) * ((E / 2 - A**2) * S3 * j.inv()).trace()
        actual = pair(dtrace)
        assert actual == pair(val(shape) / val(denom))
        results.append(actual)
    print("D=(k2/2-k1)/ell(seed) on the three local sheets for s=1,2,3,4:", results)
    assert results[0] != results[1]
    return results


if __name__ == "__main__":
    main()
