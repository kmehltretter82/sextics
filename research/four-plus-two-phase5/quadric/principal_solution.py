"""Formal polar primitive matching the 96 compact test values at the smooth point.

Arguments: upstream contact-span directory and the two-direction span.json.
Uses the Grothendieck socle convention for the common cup-product scalar.
"""

from pathlib import Path
from itertools import permutations
from functools import reduce
import json
import sys
import numpy as np
import canonical_variation as cv
from fiber_current import mon, vec, deriv, times, add, table, kernel, solve

p = 11


def pairings(selection):
    H = table(
        [
            (6, 0, 0, 0),
            (0, 6, 0, 0),
            (0, 0, 6, 0),
            (0, 0, 0, 6),
            (5, 1, 0, 0),
            (1, 5, 0, 0),
            (0, 0, 5, 1),
            (0, 0, 1, 5),
            (3, 1, 1, 1),
            (1, 1, 3, 1),
        ],
        [1, 1, 1, 1, 7, 14, 7, 21, 7, 14],
    )
    bs = mon(16, 0)
    grad = [deriv(H, i) for i in range(4)]
    mac = [
        vec(times(f, {e: 1}), bs)
        for i, f in enumerate(grad)
        for e in mon(11, int(i >= 2))
    ]
    soc = kernel(mac)
    hess = {}
    for perm in permutations(range(4)):
        a = reduce(times, [deriv(f, j) for f, j in zip(grad, perm)])
        sign = (-1) ** sum(x > y for i, x in enumerate(perm) for y in perm[i + 1 :])
        hess = add((hess, 1), (a, sign))
    raw_hess = int(soc @ vec(hess, bs)) % p
    # Hess(u^3+H)=6u Hess(H), full weighted Jacobian CI length 2*5^4.
    norm = 1250 * pow(6 * raw_hess, -1, p) % p
    y = [0] * 48
    for record in selection["certificates"]:
        for test in record["selected_test_numerators"]:
            if record["u_exponent"] == 0:
                y.append(0)
                continue
            a, b, c, d = test["monomial"]
            m = (
                a + 3,
                b + 6 * int(test["direction"] == "b6"),
                c + 3,
                d + 6 * int(test["direction"] == "d6"),
            )
            y.append((-2 * norm * int(soc[bs.index(m)])) % p)
    return np.array(y, dtype=np.int64), raw_hess


def main():
    filename = sys.argv[2]
    X = cv.complete_basis(filename)
    y, hess = pairings(json.loads(Path(filename).read_text()))
    c = solve(X.T, y)
    assert np.array_equal(c @ X % p, y)
    assert c[-1] == 8 and c[10::11].tolist() == [6, 9, 1, 8, 2, 8, 3, 7] and c[0] == 4
    hprime = (-np.arange(1, 97, dtype=np.int64) * c) % p
    print(
        "raw core Hessian socle",
        hess,
        "weighted-normalized target pairings",
        y.tolist(),
    )
    print("PP primitive tau^(-1)..tau^(-96)", c.tolist())
    print("last coefficient, primitive/differential", int(c[-1]), int(hprime[-1]))
    print("coefficients at multiples of 11 in primitive", c[10::11].tolist())
    return c, y


if __name__ == "__main__":
    main()
