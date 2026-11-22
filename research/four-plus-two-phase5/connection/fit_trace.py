"""Recover the scalar and a bounded vertical correction from an exact trace triple.

Call fit_trace(source, denom, r_degree, s_degree) with three rational functions
of r,s. The coupled P is fixed. A result includes a symbolic identity check.
"""

from itertools import product
import sympy as sp
from sympy.polys.matrices import DomainMatrix


r, s, t = sp.symbols("r s t")
coordinates = (r, s, t)
k = sp.Rational(5, 7)
f = 1 + r**6 + k * r + 2 * k * r**5
g = 1 + s**6 + k * s + 3 * k * s**5
p = f + k * r * s * (t**2 + 2 * t**4) + g * t**6


def source_of(flux):
    action = sum(q * p.diff(x) for q, x in zip(flux, coordinates))
    divergence = sum(q.diff(x) for q, x in zip(flux, coordinates))
    return sp.cancel(2 * action / 3 - p * divergence)


def fit_trace(source, denom, r_degree, s_degree):
    """Return (scalar, flux) for raw_source = scalar*t^4 + D_P(flux), or None."""
    if len(source) != 3 or min(r_degree, s_degree) < 0:
        raise ValueError("Expected three source coefficients and nonnegative bounds")
    denom = sp.Poly(denom, r, s, domain=sp.QQ).as_expr()
    if denom == 0:
        raise ValueError("The correction requires a nonzero polynomial denominator")
    raw = sum(sp.cancel(value) * t ** (2 * i) for i, value in enumerate(source))
    raw_num, raw_den = sp.fraction(sp.cancel(raw))
    if raw_den.has(t):
        raise ValueError("A trace coefficient must lie in the ratio base")

    # Numerators for the three flux components range over the normal forms
    # of the relative t-reduction, with even base components and odd t component.
    basis = []
    columns = []
    for axis, coordinate in enumerate(coordinates):
        for i, j, power in product(range(r_degree + 1), range(s_degree + 1), range(3)):
            monomial = r**i * s**j * t ** (2 * power + int(coordinate == t))
            basis.append((axis, monomial))
            numerator = sp.Rational(2, 3) * monomial * p.diff(
                coordinate
            ) * denom - p * (
                monomial.diff(coordinate) * denom - monomial * denom.diff(coordinate)
            )
            columns.append(sp.Poly(raw_den * numerator, r, s, t, domain=sp.QQ))
    scale_column = len(columns)
    columns.append(sp.Poly(raw_den * denom**2 * t**4, r, s, t, domain=sp.QQ))
    columns.append(sp.Poly(raw_num * denom**2, r, s, t, domain=sp.QQ))
    rows = sorted(set().union(*(column.monoms() for column in columns)))
    row_index = {exp: index for index, exp in enumerate(rows)}
    entries = {
        (row_index[exp], j): value
        for j, column in enumerate(columns)
        for exp, value in column.terms()
    }
    augmented = sp.SparseMatrix(len(rows), len(columns), entries)
    reduced, pivots = DomainMatrix.from_Matrix(augmented).rref()
    rhs_index = len(columns) - 1
    if rhs_index in pivots:
        return None
    reduced = reduced.to_Matrix()
    solution = [sp.S.Zero] * rhs_index
    for i, pivot in enumerate(pivots):
        solution[pivot] = reduced[i, rhs_index]
    scalar = solution[scale_column]
    flux = [sp.S.Zero] * 3
    for coefficient, (axis, monomial) in zip(solution, basis):
        flux[axis] += coefficient * monomial
    flux = tuple(sp.cancel(component / denom) for component in flux)
    assert sp.cancel(raw - scalar * t**4 - source_of(flux)) == 0
    return scalar, flux


def main():
    # A manufactured base correction on the literal coupled P. The true full
    # trace can replace this control once its three coefficients are known.
    denom = 1 + r
    known = (1 / denom, 0 * t, -t / (3 * denom**2))
    raw = source_of(known) + sp.Rational(7, 5) * t**4
    source = sp.Poly(raw, t)
    assert source.degree() == 4
    triple = [source.nth(2 * i) for i in range(3)]
    assert fit_trace(triple, 1, 0, 0) is None
    scalar, correction = fit_trace(triple, denom**2, 1, 0)
    assert scalar == sp.Rational(7, 5)
    print("Recovered scalar", scalar, "and exact coupled-chart correction", correction)


if __name__ == "__main__":
    main()
