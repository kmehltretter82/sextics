"""Direct interpolation from the four Segre axes at the exact coupled sample."""

from itertools import product
import numpy as np
from fiber_point import Dual, hbase, xbase, potential
from fiber_current import solve

p = 11


def main():
    # Independent physical evaluations recover all 49 sextic coefficients.
    grid = list(product(range(7), repeat=2))
    vand = [[pow(w, i, p) * pow(z, j, p) % p for i, j in grid] for w, z in grid]
    f = solve(
        vand, [potential([Dual(k) for k in xbase], w, z).val for w, z in grid]
    ).reshape(7, 7)
    for sign in (1, -1):
        known = np.array(
            [[sign ** (i + j) * hbase[i, j] % p for j in range(4)] for i in range(4)]
        )
        h = known.copy()
        h[1:3, 1:3] = 0
        half = pow(2, -1, p)
        h[1, 1] = (half * f[1, 1] + h[0, 1] * h[1, 0]) * pow(int(h[0, 0]), -1, p) % p
        h[1, 2] = (h[0, 2] * h[1, 3] - half * f[1, 5]) * pow(int(h[0, 3]), -1, p) % p
        h[2, 1] = (h[2, 0] * h[3, 1] - half * f[5, 1]) * pow(int(h[3, 0]), -1, p) % p
        h[2, 2] = (half * f[5, 5] + h[2, 3] * h[3, 2]) * pow(int(h[3, 3]), -1, p) % p
        assert np.array_equal(h, known)
        for i, j in grid:
            c = (
                sum(
                    h[k, l] * h[i - k, j - l] * (-1) ** (i + j - k - l)
                    for k, l in product(range(4), repeat=2)
                    if 0 <= i - k <= 3 and 0 <= j - l <= 3
                )
                % p
            )
            assert c == f[i, j]
    print(
        "both cubic factors reconstructed from four boundary cubics and four scalar equations"
    )


if __name__ == "__main__":
    main()
