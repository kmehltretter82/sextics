"""Hensel-lift the positive coupled sheet with the exact rational kappa=5/7."""

from itertools import product
import numpy as np
from fiber_point import xbase, hbase, rank, det


def solve(matrix, rhs):
    """Solve a nonsingular square system over the residue field."""
    work = np.column_stack((matrix, rhs)) % 11
    for j in range(len(matrix)):
        i = next(i for i in range(j, len(matrix)) if work[i, j])
        work[[i, j]] = work[[j, i]]
        work[j] = work[j] * pow(int(work[j, j]), -1, 11) % 11
        for i in range(len(matrix)):
            if i != j:
                work[i] = (work[i] - work[i, j] * work[j]) % 11
    answer = work[:, -1]
    assert not np.any((matrix @ answer - rhs) % 11)
    return answer


p = 11


def equations(point, modulus):
    m, n, u, v, e, P, Q, R, T, S, *tail = map(int, point)
    signs = [(-1) ** (i + j) for i, j in product(range(4), repeat=2)]
    coeff = tail[:-2]
    wp, zp = tail[-2:]
    K = e * (v - u) ** 2 * (m + n) * pow((n - m) ** 2, -1, modulus) % modulus
    k = 5 * pow(7, -1, modulus) % modulus

    def phys(w, z):
        return tuple(
            y % modulus
            for y in (1 + K * w * z, n + m * K * w * z, w + z, v * w + u * z)
        )

    def cub(w, z, sign):
        return (
            sum(
                c * t * w**i * z**j
                for (i, j), c, t in zip(product(range(4), repeat=2), coeff, sign)
            )
            % modulus
        )

    rows = []
    for w, z in product(range(7), repeat=2):
        a, b, c, d = phys(w, z)
        H = (
            a**6
            + b**6
            + c**6
            + d**6
            + k
            * (
                a**5 * b
                + 2 * a * b**5
                + c**5 * d
                + 3 * c * d**5
                + a * b * c * d * (a * a + 2 * c * c)
            )
        )
        Q1 = P * a * a + Q * b * b + R * c * c + T * c * d + S * d * d
        rows.append((H - Q1**3 - cub(w, z, [1] * 16) * cub(w, z, signs)) % modulus)
    a, b, c, d = phys(wp, zp)
    rows.extend([(b - 4 * a) % modulus, (d - 2 * c) % modulus, cub(wp, zp, [1] * 16)])
    return np.array(rows, dtype=np.int64)


def main():
    point = np.array(xbase + hbase.flatten().tolist() + [2, 10], dtype=np.int64)
    initial = equations(point, p * p)
    eye = np.eye(len(point), dtype=np.int64)
    jac = np.stack(
        [
            (equations(point + p * direction, p * p) - initial) // p % p
            for direction in eye
        ],
        axis=1,
    )
    k, pivots, _ = rank(jac)
    assert k == 28 and det(jac[pivots]) == 8
    residue = point.copy()
    for j in range(1, 8):
        mod = p**j
        values = equations(point, mod * p)
        assert not np.any(values % mod)
        correction = solve(jac[pivots], -(values[pivots] // mod))
        point += mod * correction
        assert not np.any(equations(point, mod * p))
        if j == 1:
            print("lift modulo 121:", point.tolist())
    assert np.array_equal(point % p, residue)
    print("all 49 grid equations and all three moving equations vanish modulo", p**8)

    return point


if __name__ == "__main__":
    main()
