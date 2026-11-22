"""Symbolic checks on degree-four and degree-five elliptic quotient loci."""

import sympy as sp

h, t, v = sp.symbols("h t v")

# A ± iB may be interchanged. Their product specifies the sextic as soon as
# the extra cube and B's factors of multiplicity divisible by three are removed.
for degree in (4, 5):
    if degree == 4:
        first = (t - h * v) * (t - v) ** 3
        second = (t + h * v) * (t + v) ** 3
        denominator = (first + second) / 2
        surviving = t**2 - h * h * v * v
        crossing = 1 - 3 * h - 3 * h * h
    else:
        first = (t * t + 3 * v * t + h * v * v) * (t - v) ** 3
        second = (t * t - 3 * v * t + h * v * v) * (t + v) ** 3
        denominator = (first - second) / 2
        surviving = t**4 + (2 * h - 9) * v * v * t * t + h * h * v**4
        crossing = h + 1
    assert sp.expand(first * second - surviving * (t * t - v * v) ** 3) == 0
    model = sp.Poly((surviving * denominator).subs(v, 1), t).monic()
    U, V, W = [sp.factor(model.nth(j)) for j in (2, 4, 0)]
    disc = sp.factor(sp.discriminant(model.as_expr(), t))
    assert sp.rem(sp.fraction(U)[0], crossing, h) == 0
    for unit in [sp.diff(U, h), V, W, disc]:
        numerator, den = sp.fraction(sp.factor(unit))
        assert sp.degree(sp.gcd(numerator, crossing, h), h) == 0
        assert sp.degree(sp.gcd(den, crossing, h), h) == 0
    print(f"d={degree}: U={U}, V={V}, W={W}, ordinary crossing={crossing}")

# y=g, epsilon=h-h0, c1 and c2 analytic units, b/g=c1/epsilon,
# f/g=c2/epsilon^3. Substitute r=c1*y/(k*s*epsilon) in the binary f.
points = [(0, 6), (1, 5), (5, 1), (6, 0), (1, 3)]
for valuation, expected in [(sp.Rational(3), [0, 1]), (sp.Rational(3, 5), [1, 6])]:
    weights = [n * valuation + d for n, d in points]
    edge = sorted(
        n for (n, d), weight in zip(points, weights) if weight == min(weights)
    )
    assert edge == expected
    print(f"ord_epsilon(g)={valuation}: first binomial powers {edge}")
