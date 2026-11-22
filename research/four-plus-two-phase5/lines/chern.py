"""Degree of the top Chern class on P(Sym² U*) over Gr(2,4), for lines convention."""

import sympy as sp

x, y, t = sp.symbols("x y t")
roots = [(6 - i) * x + i * y for i in range(7)]
chern = sp.Poly(sp.prod(1 + t * z for z in roots), t)
quadratic = [2 * x, x + y, 2 * y]
segre = [1]
denom = sp.Poly(sp.prod(1 + t * z for z in quadratic), t)
for n in range(1, 5):
    segre.append(-sum(denom.nth(i) * segre[n - i] for i in range(1, min(n, 3) + 1)))
pushed = sp.expand(sum(3**j * chern.nth(6 - j) * segre[j - 2] for j in range(2, 7)))
degree = sp.Poly((x - y) * pushed, x, y).coeff_monomial(x**3 * y**2)
assert degree == 2628
print("Pushforward:", pushed)
print("Degree:", degree)
