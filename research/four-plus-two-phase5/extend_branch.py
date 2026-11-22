"""Extend the fixed contact branch with exact quadratic-convolution coefficient solves."""

import argparse
from pathlib import Path
import sys
import numpy as np


def generate(span, degree):
    P = span.P
    z0, ds, receipt = span.source_data()
    D = span.differential(z0)
    rows = receipt["incidence_derivative"]["minor_rows"]
    cols = receipt["incidence_derivative"]["minor_columns"]
    (missing,) = set(range(len(D))) - set(rows)
    inv = span.base_inverse(D[np.ix_(rows, cols)])
    ell = np.zeros(len(D), dtype=np.int64)
    ell[missing] = 1
    ell[rows] = -D[missing, cols] @ inv % P
    d1, d2 = [np.column_stack((u, np.zeros_like(u))) for u in ds]
    z = np.zeros((degree + 1, 90, 2), dtype=np.int64)
    z[0, :, 0] = z0
    z[1] = (d1 + span.mul([19, 8], d2)) % P
    denom = ell @ (span.quadratic(z[1], d2) + span.quadratic(d2, z[1])) % P
    invden = span.inv(denom)
    # Each product here is a dot product of nonnegative residue representatives.
    # Keep the integer totals strictly below the exact binary64 integer threshold.
    assert degree * 16 * 23**2 * len(span.term_sign) < 2**52
    gram = np.zeros((90, 90, 2), dtype=np.int64)
    blocks = [(0, 2), (1, 3), (4, 4), (5, 5)]
    offsets = list(span.offsets) + [90]
    z_float = np.zeros(z.shape, dtype=np.float64)
    z_float[:2] = z[:2]

    def convolution(n, lo):
        for left, right in blocks:
            a, b = (
                slice(offsets[left], offsets[left + 1]),
                slice(offsets[right], offsets[right + 1]),
            )
            end = n - lo + 1
            x = z_float[lo:end, a, 0]
            u = z_float[lo:end, a, 1]
            y = z_float[lo:end, b, 0][::-1]
            v = z_float[lo:end, b, 1][::-1]
            gram[a, b, 0] = x.T @ y + 15 * (u.T @ v)
            gram[a, b, 1] = x.T @ v + u.T @ y
        out = np.zeros((70, 2), dtype=np.int64)
        np.add.at(
            out, span.term_row, span.term_sign[:, None] * gram[span.term_i, span.term_j]
        )
        return out % P

    forcing = span.quadratic(z[1], z[1])
    for n in range(2, degree + 1):
        assert not np.any(ell @ forcing % P)
        normal = np.zeros((90, 2), dtype=np.int64)
        normal[cols] = -inv @ forcing[rows] % P
        assert not np.any((D @ normal + forcing) % P)
        inner = convolution(n + 1, 2)
        cross = span.quadratic(z[1], normal) + span.quadratic(normal, z[1])
        correction = -span.mul(ell @ (inner + cross) % P, invden) % P
        z[n] = (normal + span.mul(correction, d2)) % P
        z_float[n] = z[n]
        forcing = (inner + span.quadratic(z[1], z[n]) + span.quadratic(z[n], z[1])) % P
    return z


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("packet", type=Path)
    parser.add_argument("degree", type=int)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    sys.path.insert(0, str(args.packet / "contact-span"))
    import contact_series_span as span

    z = generate(span, args.degree)
    archive = np.load(args.packet / "contact-span/contact_span_coefficients.npz")[
        "branch"
    ]
    n = min(len(z), len(archive))
    assert np.array_equal(z[:n], archive[:n])
    np.savez_compressed(args.output, branch=z)
    print(
        "Verified original",
        n,
        "coefficients and fixed equation through order",
        args.degree,
    )


if __name__ == "__main__":
    main()
