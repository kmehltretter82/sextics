"""Exact H-only pole-three spanning test for the four-plus-two weighted sextic.

The columns are degree-(n-5) Jacobian multiples, then products of sigma-odd
degree-(n-6) monomials with b^6 and d^6. A full-rank square submatrix over
the integers proves the claim over Q. GF(23) is used only to choose columns.
"""

import argparse
import json
from pathlib import Path
from flint import fmpz_mat, nmod_mat

# Coefficients of 7H in variables (a,b,c,d). The nonzero scalar does not
# affect its Jacobian ideal over Q.
H7 = {
    (6, 0, 0, 0): 7,
    (0, 6, 0, 0): 7,
    (0, 0, 6, 0): 7,
    (0, 0, 0, 6): 7,
    (5, 1, 0, 0): 5,
    (1, 5, 0, 0): 10,
    (0, 0, 5, 1): 5,
    (0, 0, 1, 5): 15,
    (3, 1, 1, 1): 5,
    (1, 1, 3, 1): 10,
}
DIRECTIONS = {"b6": (0, 6, 0, 0), "d6": (0, 0, 0, 6)}


def monomials(n, parity):
    return [
        (a, b, c, n - a - b - c)
        for a in range(n + 1)
        for b in range(n - a + 1)
        for c in range(n - a - b + 1)
        if (n - a - b) % 2 == parity
    ]


def plus(x, y):
    return tuple(i + j for i, j in zip(x, y))


def jacobian_columns(n):
    target = monomials(n, 1)
    loc = {m: k for k, m in enumerate(target)}
    data = []
    labels = []
    for var in range(4):
        for m in monomials(n - 5, 1 - (var >= 2)):
            vector = [0] * len(target)
            for e, coeff in H7.items():
                if e[var]:
                    exponent = list(plus(m, e))
                    exponent[var] -= 1
                    vector[loc[tuple(exponent)]] += e[var] * coeff
            data.append(vector)
            labels.append(dict(derivative="abcd"[var], monomial=m))
    return target, data, labels


def degree_certificate(n, p):
    target, data, labels = jacobian_columns(n)
    loc = {m: i for i, m in enumerate(target)}
    jacobian_cols = len(data)
    for name, t in DIRECTIONS.items():
        for m in monomials(n - 6, 1):
            vector = [0] * len(target)
            vector[loc[plus(m, t)]] = 1
            data.append(vector)
            labels.append(dict(direction=name, monomial=m))
    # The transpose of the list of column vectors is the raw coefficient matrix.
    matrix = nmod_mat(data, p).transpose()
    red, rank = matrix.rref()
    assert rank == len(target)
    chosen = [next(i for i, x in enumerate(row) if x) for row in red.tolist()[:rank]]
    jac_selected = [i for i in chosen if i < jacobian_cols]
    chosen_tests = [labels[i] for i in chosen if i >= jacobian_cols]
    integral = fmpz_mat([data[i] for i in chosen]).transpose()
    det = int(integral.det())
    assert det % p != 0
    print(
        f"core degree {n}, {len(target)} raw coordinates, "
        f"{len(jac_selected)} Jacobian columns, {len(chosen_tests)} products, "
        f"integer determinant mod {p} = {det % p}",
        flush=True,
    )
    return {
        "core_target_degree": n,
        "u_exponent": (12 - n) // 2,
        "raw_dimension": len(target),
        "jacobian_columns_selected": len(jac_selected),
        "product_columns_selected": len(chosen_tests),
        "selected_column_indices": chosen,
        "selected_test_numerators": chosen_tests,
        "determinant": str(det),
        "determinant_mod_p": det % p,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prime", type=int, default=23)
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()
    certs = [degree_certificate(n, args.prime) for n in (12, 10)]
    assert [c["product_columns_selected"] for c in certs] == [16, 32]
    if args.json:
        args.json.write_text(
            json.dumps(
                dict(
                    column_selection_prime=args.prime,
                    variables=["a", "b", "c", "d"],
                    certificates=certs,
                ),
                indent=2,
            )
            + "\n"
        )


if __name__ == "__main__":
    main()
