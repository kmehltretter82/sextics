"""Coefficient identities for the rank-three fiber connection and base supports.

Run with Python 3 and SymPy. Every calculation is in characteristic zero.
"""

from itertools import combinations
import sympy as sp

z, a, b, d = sp.symbols("z a b d")
r, s, c, lam = sp.symbols("r s c lam")
channels = [a, b, d]
kappa = sp.Rational(5, 7)
f = 1 + r**6 + kappa * r + 2 * kappa * r**5
g = 1 + s**6 + kappa * s + 3 * kappa * s**5
P = a + b * z + 2 * b * z**2 + d * z**3
disc = sp.discriminant(P, z)


def coefficients(poly):
    polynomial = sp.Poly(poly, z)
    return sp.Matrix([polynomial.nth(i) for i in range(6)])


def polynomial_columns(poly):
    return sp.Matrix.hstack(*[coefficients(z**j * poly) for j in range(3)])


def factor(matrix):
    return matrix.applyfunc(sp.factor)


def at(matrix, x, value):
    return matrix.applyfunc(lambda entry: sp.limit(entry, x, value))


M = polynomial_columns(P).row_join(
    sp.Matrix.hstack(
        *[
            coefficients(z**j * ((2 * j + 1) * P - sp.Rational(4, 3) * z * P.diff(z)))
            for j in range(3)
        ]
    )
)
assert sp.factor(M.det(method="domain-ge") / (a * d * disc)) == sp.Rational(64, 27)
adj, denominator = M.to_DM().inv_den()
red = adj.to_Matrix()[:3, :] / denominator.as_expr()
conn = {
    x: factor(-sp.Rational(2, 3) * red * polynomial_columns(P.diff(x)))
    for x in channels
}
for x in channels:
    log_derivative = sp.diff(a, x) / a + 7 * sp.diff(d, x) / d + sp.diff(disc, x) / disc
    assert sp.cancel(6 * conn[x].trace() + log_derivative) == 0
for x, y in combinations(channels, 2):
    curvature = (
        conn[y].diff(x) - conn[x].diff(y) + conn[x] * conn[y] - conn[y] * conn[x]
    )
    assert curvature.applyfunc(sp.cancel) == sp.zeros(3)
print("Universal determinant, logarithmic trace and flatness verified.")
# An ordinary member of the transverse class. Deg 36 counts all finite collisions.
line = {a: f, b: kappa * r * (r + 1), d: g.subs(s, r + 1)}
boundary = sp.Poly((a * d * disc).subs(line), r)
assert boundary.degree() == 36 and boundary.gcd(boundary.diff()).degree() == 0
print(
    "The line s=r+1 has six simple f zeros, six simple g zeros and 24 simple discriminant zeros."
)

A = factor(conn[b].subs({a: 1, b: lam, d: c}))
Delta = disc.subs({a: 1, b: lam, d: c})
assert sp.factor(sp.discriminant(Delta, lam)) == -5038848 * c**4 * (c - 8) ** 4
assert sp.cancel(A[0, 1] * Delta - (2 * lam**2 - 3 * c * lam + 9 * c)) == 0
joint = at(-lam * A, lam, sp.oo) - sp.Rational(2, 3) * sp.eye(3)
expected = sp.Matrix([[-2, 0, 0], [4, 0, -5], [0, 0, -14]]) / 12
assert joint == expected
print("Pencil g=c f: joint residue", joint.tolist())
reciprocal = factor(A.subs(c, 8))
assert sp.expand(Delta.subs(c, 8) - 4 * (lam + 2) * (lam - 6) ** 3) == 0
for root, eigenvalue in [(6, -sp.Rational(1, 2)), (-2, -sp.Rational(1, 6))]:
    residue = at((lam - root) * reciprocal, lam, root)
    assert residue.eigenvals() == {0: 2, eigenvalue: 1}
print("The reciprocal member has finite residues (-1/2,0,0) and (-1/6,0,0).")
assert sp.Poly(f, r).gcd(sp.Poly(f.diff(r), r)).degree() == 0
assert sp.Poly(g, s).gcd(sp.Poly(g.diff(s), s)).degree() == 0

# Rational ODEs on the generic r- and s-lines. A nonzero minor in the
# augmented system after specialization proves independence on the generic line.
from itertools import product

target = red * coefficients(z**2)
core_coeffs = {a: f, b: kappa * r * s, d: g}
for x, other, bounds, expected_minor in [
    (r, s, (3, 1, -1), 643),
    (s, r, (1, 3, 5), 563),
]:
    line = {key: value.subs(other, 1) for key, value in core_coeffs.items()}
    E = sp.expand((a * d * disc).subs(line))
    assert sp.degree(E, x) == 18 and sp.gcd(E, E.diff(x), x) == 1
    A_line = sum(
        (conn[key].subs(line) * value.diff(x) for key, value in line.items()),
        start=sp.zeros(3),
    )
    A_num = (E * A_line).applyfunc(sp.cancel)
    target_num = (E * target.subs(line)).applyfunc(sp.cancel)
    columns = []
    for j, bound in enumerate(bounds):
        for i in range(bound + 1):
            v = sp.zeros(3, 1)
            v[j] = x**i
            columns.append(A_num * v + E * v.diff(x))
    columns.append(target_num)
    polynomials = [[sp.Poly(entry, x) for entry in column] for column in columns]
    degree = max(entry.degree() for column in polynomials for entry in column)
    coefficient_matrix = sp.Matrix(
        [
            [column[j].nth(i) for column in polynomials]
            for i, j in product(range(degree + 1), range(3))
        ]
    )
    rows = coefficient_matrix.T.rref()[1]
    minor = coefficient_matrix.extract(rows, range(len(columns)))
    assert minor.is_square
    det = minor.det()
    determinant_mod = int(det.p) * pow(int(det.q), -1, 1009) % 1009
    assert determinant_mod == expected_minor
    print(
        f"{x} coefficient matrix {coefficient_matrix.shape}: rows {rows}, "
        f"augmented determinant {determinant_mod} mod 1009."
    )
