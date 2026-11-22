"""Certify ranks of the coupled-core ratio resultant over F_23[v]/(v²−15).

Usage:
    python fiber_degree.py path/to/contact_span_coefficients.npz
"""

from itertools import product
import sys

import numpy as np


PRIME = 23
V_SQUARED = 15


def monomials(degree, parity):
    return sorted(
        powers + (u,)
        for powers in product(range(degree + 1), repeat=4)
        for u in range(degree // 2 + 1)
        if sum(powers) + 2 * u == degree and sum(powers[2:]) % 2 == parity
    )


def multiply(x, y):
    """Multiply truncated series with coefficients stored as real/imaginary pairs."""

    def convolution(i, j):
        return np.convolve(x[:, i], y[:, j])[: len(x)]

    real = convolution(0, 0) + V_SQUARED * convolution(1, 1)
    imag = convolution(0, 1) + convolution(1, 0)
    return np.stack((real % PRIME, imag % PRIME), axis=1)


def contact_equations(jets, r, s):
    """Evaluate a=1,b=r,c=alpha,d=s*alpha and eliminate the coordinate u using f1."""
    body = []
    offset = 0
    for degree, parity in ((2, 0), (2, 0), (4, 0), (4, 0), (3, 0), (3, 1)):
        basis = monomials(degree, parity)
        slots = np.zeros((degree // 2 + 1, degree + 1, len(jets), 2), dtype=np.int64)
        for j, (_, b, c, d, u) in enumerate(basis):
            value = pow(r, b, PRIME) * pow(s, d, PRIME) % PRIME
            slots[u, c + d] += jets[:, offset + j] * value
        offset += len(basis)
        body.append(slots % PRIME)

    assert offset == jets.shape[1]
    f1, f2, _, _, even_cubic, odd_cubic = body
    h = even_cubic + odd_cubic
    # The coefficient f1_u is an invertible series at this seed.
    assert f1[1, 0, 0].tolist() == [1, 0]
    conic = [multiply(f1[1, 0], f2[0, j]) - multiply(f2[1, 0], f1[0, j]) for j in (0, 2)]
    cubic = [multiply(f1[1, 0], h[0, i]) for i in range(4)]
    for i in (0, 1):
        for j in (0, 2):
            cubic[i + j] -= multiply(h[1, i], f1[0, j])
    return conic, cubic


def resultant_jets(jets, r, s):
    (u, v), h = contact_equations(jets, r, s)
    even = multiply(v, h[0]) - multiply(u, h[2])
    odd = multiply(v, h[1]) - multiply(u, h[3])
    return (multiply(v, multiply(even, even)) + multiply(u, multiply(odd, odd))) % PRIME


def determinant(matrix):
    """Field Gaussian elimination with sign in the supplied ordering."""
    a = matrix.copy() % PRIME
    n = len(a)
    assert a.shape == (n, n, 2)
    det = np.array([1, 0])
    for j in range(n):
        nonzero = np.flatnonzero(np.any(a[j:, j] != 0, axis=1))
        if not len(nonzero):
            return [0, 0]
        pivot = j + int(nonzero[0])
        if pivot != j:
            a[[j, pivot]] = a[[pivot, j]]
            det = -det
        x, y = a[j, j]
        det = np.array([det[0] * x + V_SQUARED * det[1] * y, det[0] * y + det[1] * x]) % PRIME
        inv_norm = pow(int((x * x - V_SQUARED * y * y) % PRIME), -1, PRIME)
        b, c = x * inv_norm % PRIME, -y * inv_norm % PRIME
        row0, row1 = a[j, :, 0].copy(), a[j, :, 1].copy()
        row_real = (b * row0 + V_SQUARED * c * row1) % PRIME
        row_imag = (b * row1 + c * row0) % PRIME
        factors_real = a[j + 1 :, j, 0, None].copy()
        factors_imag = a[j + 1 :, j, 1, None].copy()
        a[j + 1 :, :, 0] -= factors_real * row_real + V_SQUARED * factors_imag * row_imag
        a[j + 1 :, :, 1] -= factors_real * row_imag + factors_imag * row_real
        a[j + 1 :] %= PRIME
    return det.tolist()


def infinity_check(jets):
    """Check that the conic and cubic leading coefficients are generically coprime in s."""
    # Reconstruct the conic's first tau coefficient by interpolation at s=0,1,-1.
    vals = [contact_equations(jets, 0, s)[0][1][1] % PRIME for s in (0, 1, PRIME - 1)]
    v0, vp, vm = vals
    v1 = (vp - vm) * pow(2, -1, PRIME) % PRIME
    v2 = ((vp + vm) * pow(2, -1, PRIME) - v0) % PRIME
    # Verify the cubic leading coefficient at the seed on four distinct points.
    points = range(4)
    base = [contact_equations(jets, 0, j)[1][3][0] % PRIME for j in points]
    cubic = [(4, 0), (6, 0), (17, 0), (22, 0)]
    for s, value in zip(points, base):
        expected = sum(np.array(coeff) * pow(s, i, PRIME) for i, coeff in enumerate(cubic)) % PRIME
        assert np.array_equal(value, expected)
    sylvester = np.zeros((5, 5, 2), dtype=np.int64)
    for i in range(3):
        sylvester[i, i : i + 3] = [v2, v1, v0]
    for i in range(2):
        sylvester[3 + i, i : i + 4] = cubic[::-1]
    return determinant(sylvester)


if __name__ == "__main__":
    data = np.load(sys.argv[1], allow_pickle=False)
    jets = data["branch"][:49]
    assert jets.shape == (49, 90, 2)
    grid = list(product(range(7), repeat=2))
    matrix = np.stack([resultant_jets(jets, r, s) for r, s in grid], axis=1)
    field_det = determinant(matrix)
    line_det = determinant(matrix[3:10, :7])
    infinity_res = infinity_check(jets)
    assert field_det == [15, 5]
    assert not np.any(matrix[:3, :7]) and line_det == [14, 1]
    assert infinity_res == [8, 6]
    print("49 by 49 determinant on the 7 by 7 grid:", field_det)
    print("7 by 7 r=0 determinant, jet orders 3 through 9:", line_det)
    print("conic/cubic Sylvester determinant at alpha infinity:", infinity_res)
