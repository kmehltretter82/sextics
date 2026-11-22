"""Obstructions for the twisted connection on the coupled energy plane, over Q."""

import argparse
from itertools import product

import sympy as sp
from sympy.polys.matrices import DomainMatrix

from energy_connection import R, critical_check, energy, f, matrices, r


def energy_only(basis, multiplication, divergence):
    """Check the finite residues, exponents at infinity and forced candidate at r=1."""
    target = basis.index((3, 0))
    dim = len(basis)
    tt = multiplication.subs(r.as_expr(), 1)
    bb = divergence.subs(r.as_expr(), 1)
    identity = sp.eye(dim)
    target_vector = identity[:, target]
    cp = critical_check(1, tt)
    assert cp.eval(0) and not bb.free_symbols and not tt.free_symbols

    # The spectral projectors of tt are spanned by its first dim powers.
    power = identity
    for _ in range(dim):
        assert ((bb - identity) * power).trace() == 0
        power = power * tt

    spectrum = {
        sp.Rational(1, 2): 2,
        sp.Rational(5, 6): 4,
        sp.Rational(7, 6): 4,
        sp.Rational(3, 2): 2,
    }
    assert bb.eigenvals() == spectrum
    # A polynomial solution of degree n>1 needs 8/3-n in this spectrum.
    assert all(not (sp.Rational(8, 3) - b).is_integer for b in spectrum)
    # At zero the only linear candidate has slope -3/2 times the target.
    slope = -sp.Rational(3, 2) * target_vector
    assert (1 - sp.Rational(5, 3)) * slope == target_vector
    assert cp.eval(2) != 0
    gamma_at_two = (bb - identity) * (2 * identity - tt).inv()
    residual = 2 * gamma_at_two * slope
    assert any(residual)
    print(
        "Pure energy direction: 12 nilpotent critical residues, twist residue -5/3, "
        "forced linear candidate fails at E=2.",
        flush=True,
    )


def search(
    basis, multiplication, divergence, coupling, denominator, r_degree, energy_degree
):
    """Construct and row-reduce the coefficient equations for this numerator box."""
    dim = len(basis)
    monomials = list(product(range(r_degree + 1), range(energy_degree + 1)))
    columns = list(product(("X", "Z"), range(dim), monomials))
    identity = sp.eye(dim)
    delta = identity * energy.as_expr() - multiplication

    def convert(mat):
        return [[R.from_expr(mat[i, j]) for j in range(dim)] for i in range(dim)]

    C, J, B = map(convert, (coupling, delta, divergence))
    d = denominator
    dr, de = d.diff(r), d.diff(energy)
    indices = {}
    data = {}

    def append(component, column, polynomial):
        for exp, coefficient in polynomial.items():
            row = indices.setdefault((component, *exp), len(indices))
            data.setdefault(row, {})[column] = coefficient

    for col, (field, i, (j, k)) in enumerate(columns):
        m = r**j * energy**k
        twisted = 3 * energy * (d * m.diff(energy) - de * m) - 5 * d * m
        if field == "X":
            horizontal = 3 * energy * (d * m.diff(r) - dr * m)
            for row in range(dim):
                append(
                    row, col, C[row][i] * twisted + (horizontal if row == i else R.zero)
                )
        else:
            for row in range(dim):
                append(row, col, J[row][i] * twisted + 3 * energy * d * B[row][i] * m)
    append(basis.index((3, 0)), len(columns), 3 * energy * d**2)
    matrix = DomainMatrix.from_dod(data, (len(indices), len(columns) + 1), sp.QQ)
    _, pivots = matrix.rref()
    coefficient_rank = len(pivots) - int(len(columns) in pivots)
    print(
        f"deg(r,E) <= ({r_degree},{energy_degree}), denominator {d}, "
        f"{len(indices)} rows, {len(columns)} unknowns, "
        f"ranks {coefficient_rank}/{len(pivots)}",
        flush=True,
    )
    return coefficient_rank, len(pivots)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--r-degree", type=int, default=2)
    parser.add_argument("--E-degree", type=int, default=1)
    denominators = {"unit": R.one, "cube": energy - r**3, "central": energy - f}
    parser.add_argument("--denom", choices=["all", *denominators], default="all")
    args = parser.parse_args()
    if min(args.r_degree, args.E_degree) < 0:
        parser.error("Numerator degree bounds must be nonnegative.")
    basis, T, B, C = matrices()
    energy_only(basis, T, B)
    selected = (
        denominators.values() if args.denom == "all" else [denominators[args.denom]]
    )
    for d in selected:
        ranks = search(basis, T, B, C, d, args.r_degree, args.E_degree)
        if (args.r_degree, args.E_degree) == (2, 1):
            assert ranks == (144, 145)


if __name__ == "__main__":
    main()
