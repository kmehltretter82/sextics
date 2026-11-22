"""Compute the full weighted-family infinitesimal intertwiner at the Fermat reduction."""

from itertools import product
from flint import nmod_mat


PRIME = 5
DIMENSION = 48
# Each element is a coefficient and an exponent in (a,b,c,d,u).
DIRECTIONS = [
    [
        (2, (0, 0, 0, 4, 1)),
        (1, (0, 0, 2, 2, 1)),
        (3, (1, 1, 3, 1, 0)),
        (3, (2, 0, 2, 2, 0)),
        (2, (3, 3, 0, 0, 0)),
    ],
    [
        (2, (0, 0, 3, 1, 1)),
        (2, (1, 3, 0, 0, 1)),
        (3, (3, 1, 1, 1, 0)),
        (3, (3, 1, 2, 0, 0)),
    ],
    [
        (1, (0, 0, 3, 1, 1)),
        (1, (0, 2, 4, 0, 0)),
        (2, (1, 1, 0, 2, 1)),
        (2, (1, 1, 1, 3, 0)),
    ],
]


def negative_basis(degree):
    return [
        exp
        for exp in product(range(5), range(5), range(5), range(5), range(2))
        if sum(exp) + exp[4] == degree and (exp[2] + exp[3]) % 2
    ]


def multiplication(direction):
    source = negative_basis(6)
    index = {exp: i for i, exp in enumerate(negative_basis(12))}
    assert len(source) == len(index) == DIMENSION
    result = nmod_mat(DIMENSION, DIMENSION, PRIME)
    for scalar, term in direction:
        assert sum(term) + term[4] == 6 and (term[2] + term[3]) % 2 == 0
        for j, monomial in enumerate(source):
            exp = tuple(a + b for a, b in zip(monomial, term))
            if exp in index:
                result[index[exp], j] += scalar
    return result


def main():
    first, second, third = map(multiplication, DIRECTIONS)
    assert int(first.det()) == 1
    x, y = first.inv() * second, first.inv() * third
    identity = nmod_mat(
        [[int(i == j) for j in range(DIMENSION)] for i in range(DIMENSION)], PRIME
    )
    kernel, nullity = (x - 2 * identity).nullspace()
    assert nullity == 1
    basis = nmod_mat([[kernel[i, 0]] for i in range(DIMENSION)], PRIME)
    assert x * basis == 2 * basis
    ranks = [1]
    while basis.ncols() < DIMENSION:
        candidates = nmod_mat(
            [
                a + b + c
                for a, b, c in zip(
                    basis.tolist(), (x * basis).tolist(), (y * basis).tolist()
                )
            ],
            PRIME,
        )
        reduced, rank = candidates.rref()
        assert rank > basis.ncols()
        ranks.append(rank)
        pivots = [
            next(i for i, value in enumerate(row) if value)
            for row in reduced.tolist()[:rank]
        ]
        basis = nmod_mat(
            [[row[i] for i in pivots] for row in candidates.tolist()], PRIME
        )
    assert int(basis.det()) == 4
    print("det(M0) = 1, dim ker(X - 2I) = 1, orbit ranks", ranks, "det(basis) = 4")


if __name__ == "__main__":
    main()
