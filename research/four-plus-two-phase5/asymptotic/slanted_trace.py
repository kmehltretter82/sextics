"""Verify a nonzero normalized source projection by eliminating a sextic direction."""

import sympy as sp


A, B, D, G, lam = sp.symbols("A B D G lam")
norm = D**2 * (D**2 + 3 * A * B * G) ** 2 + G * (B**3 + G**2 * A**3) ** 2
assert sp.resultant(lam**6 + G, D + A * lam**5 + B * lam, lam) == sp.expand(norm)
log_ab = norm * norm.diff(A, B) - norm.diff(A) * norm.diff(B)
log_dd = norm * norm.diff(D, D) - norm.diff(D) ** 2
assert sp.expand(log_ab + G * log_dd) == 0
p, s, t, k = sp.symbols("p s t k")
r = p + lam * t
P = 1 + r**6 + k * r + 2 * k * r**5 + k * r * s * (t**2 + 2 * t**4) + G * t**6
coeff = sp.Poly(P, t)
L = (6 * p + 2 * k) * lam**5 + 2 * k * lam * s
assert sp.expand(coeff.nth(5) - L) == 0 and coeff.nth(6) == G + lam**6
residual = sp.Poly(5 * L * t**4 - sp.diff(P, t), t)
assert residual.nth(4) == 0 and sp.expand(residual.nth(5) + 6 * (G + lam**6)) == 0
assert residual.nth(3) != 0
print(
    "Sextic norm, zero divergence, degree-five leading coefficient, and cubic residual verified"
)
