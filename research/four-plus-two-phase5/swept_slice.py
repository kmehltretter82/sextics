"""Interpolate equations of the swept surface on a plane d=s*c from all six moving points."""

import argparse
import sys
from pathlib import Path

import numpy as np
from flint import nmod_mat, nmod_poly


class Series:
    """Truncated series in the branch parameter, over F_23[v]/(v^2-15)."""

    precision = 1

    def __init__(self, a=0, b=0):
        self.a = (
            a
            if isinstance(a, nmod_poly)
            else nmod_poly(a if isinstance(a, list) else [int(a)], 23)
        )
        self.b = (
            b
            if isinstance(b, nmod_poly)
            else nmod_poly(b if isinstance(b, list) else [int(b)], 23)
        )

    def __add__(self, y):
        if not isinstance(y, Series):
            y = Series(y)
        return Series(self.a + y.a, self.b + y.b)

    __radd__ = __add__

    def __sub__(self, y):
        return self + (-y)

    def __neg__(self):
        return Series(-self.a, -self.b)

    def __mul__(self, y):
        if not isinstance(y, Series):
            return Series(self.a * int(y), self.b * int(y))
        n = self.precision
        ac = self.a.mul_low(y.a, n)
        bd = self.b.mul_low(y.b, n)
        cross = (self.a + self.b).mul_low(y.a + y.b, n) - ac - bd
        return Series(ac + 15 * bd, cross)

    __rmul__ = __mul__

    def inv(self):
        n = self.precision
        norm = self.a.mul_low(self.a, n) - 15 * self.b.mul_low(self.b, n)
        reciprocal = norm.inverse_series_trunc(n)
        return Series(self.a.mul_low(reciprocal, n), -self.b.mul_low(reciprocal, n))

    def coefficients(self):
        result = np.zeros((self.precision, 2), dtype=np.int16)
        for j, poly in enumerate((self.a, self.b)):
            vals = list(poly.truncate(self.precision))
            result[: len(vals), j] = vals
        return result


def verify_fixed_core(jets, span):
    """Re-multiply the six defining forms using independent FLINT polynomials."""
    Series.precision = len(jets)
    coeffs = [Series(j[:, 0].tolist(), j[:, 1].tolist()) for j in jets.swapaxes(0, 1)]
    result = [Series() for _ in span.plus_basis]
    for row, i, j, sign in zip(span.term_row, span.term_i, span.term_j, span.term_sign):
        result[row] = result[row] + int(sign) * coeffs[i] * coeffs[j]
    central, _, _ = span.source_data()
    target = span.quadratic(
        np.stack([central, np.zeros_like(central)], axis=1),
        np.stack([central, np.zeros_like(central)], axis=1),
    )
    for value, expected in zip(result, target):
        assert not (value - Series(*expected)).coefficients().any()


def equations(jets, monomials, s):
    """Eliminate the monic weighted variable first, keeping a,b and setting c=1,d=s."""
    params = [Series(x[:, 0].tolist(), x[:, 1].tolist()) for x in jets.swapaxes(0, 1)]
    blocks = []
    offset = 0
    for deg, parity in ((2, 0), (2, 0), (4, 0), (4, 0), (3, 0), (3, 1)):
        form = {}
        for a, b, _, d, u in monomials(deg, parity):
            key = a, b, u
            form[key] = form.get(key, 0) + int(pow(s, d, 23)) * params[offset]
            offset += 1
        blocks.append(form)
    assert offset == 90
    f1, f2, _, _, even, odd = blocks
    h = dict(even)
    for m, c in odd.items():
        h[m] = h.get(m, 0) + c

    def eliminate(form):
        result = {(a, b): c * f1[0, 0, 1] for (a, b, u), c in form.items() if u == 0}
        for (a, b, u), factor in form.items():
            if u == 1:
                for (x, y, k), c in f1.items():
                    if k == 0:
                        m = a + x, b + y
                        result[m] = result.get(m, Series()) - factor * c
        return result

    return eliminate(f2), eliminate(h)


def evaluate(form, a, b):
    zero, one = Series(), Series(1)
    ap = [one, a, a * a, a * a * a]
    bp = [one, b, b * b, b * b * b]
    value, da, db = zero, zero, zero
    for (i, j), weight in form.items():
        value += weight * ap[i] * bp[j]
        if i:
            da += i * weight * ap[i - 1] * bp[j]
        if j:
            db += j * weight * ap[i] * bp[j - 1]
    return value, da, db


def roots_mod_seed(s):
    B = (4 + 6 * s + 17 * s**2 + 22 * s**3) % 23
    if not B:
        raise ValueError("The chosen binary cubic vanishes at this slice.")
    for x in range(23):
        for y in range(23):
            if (x**3 + 45 * x * y * y + B) % 23 == 0 and (
                3 * x * x * y + 15 * y**3
            ) % 23 == 0:
                yield Series(x, y)


def plane_branches(jets, monomials, s):
    n = len(jets)
    Series.precision = n
    conic, cubic = equations(jets, monomials, s)
    roots = list(roots_mod_seed(s))
    assert len(roots) == 3
    for root in roots:
        for a, b in ((root, Series()), (Series(), root)):
            solved = 1
            while solved < n:
                Series.precision = min(2 * solved, n)
                q, qa, qb = evaluate(conic, a, b)
                h, ha, hb = evaluate(cubic, a, b)
                invdet = (qa * hb - qb * ha).inv()
                a, b = a - invdet * (hb * q - qb * h), b - invdet * (qa * h - ha * q)
                solved = Series.precision
            assert not evaluate(conic, a, b)[0].coefficients().any()
            assert not evaluate(cubic, a, b)[0].coefficients().any()
            yield a, b


def evaluate_monomials(branches, degree):
    # Columns are indexed by (i,j), i+j <= degree. Each plane point has c=1.
    exps = [(i, n - i) for n in range(degree + 1) for i in range(n + 1)]
    matrix = np.empty((6, Series.precision, len(exps), 2), dtype=np.int16)
    for sheet, (a, b) in enumerate(branches):
        bpowers = [Series(1)]
        for j in range(degree):
            bpowers.append(bpowers[-1] * b)
        for i in range(degree + 1):
            for j in range(degree - i + 1):
                col = (i + j) * (i + j + 1) // 2 + i
                matrix[sheet, :, col] = bpowers[j].coefficients()
                bpowers[j] *= a
    return exps, matrix.reshape(-1, len(exps), 2)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("packet", type=Path)
    parser.add_argument("--jets", type=Path)
    parser.add_argument("--degree", type=int, default=26)
    parser.add_argument("--s", type=int, default=1)
    args = parser.parse_args()
    sys.path.insert(0, str(args.packet / "contact-span"))
    import contact_series_span as span

    path = args.jets or args.packet / "contact-span/contact_span_coefficients.npz"
    jets = np.load(path)["branch"]
    verify_fixed_core(jets, span)
    branches = list(plane_branches(jets, span.monomials, args.s))
    exps, m = evaluate_monomials(branches, args.degree)
    # Restrict scalars. Full prime-field rank 2*N is equivalent to extension rank N.
    real, imag = m[:, :, 0], m[:, :, 1]
    base = np.block([[real, 15 * imag % 23], [imag, real]])
    mat = nmod_mat(*base.shape, base.flatten().tolist(), 23)
    rank = mat.rank()
    print(
        "slice s",
        args.s,
        "precision",
        len(jets),
        "degree",
        args.degree,
        "extension rank",
        rank // 2,
        "of",
        len(exps),
        "kernel dimension",
        len(exps) - rank // 2,
        flush=True,
    )
    assert rank % 2 == 0


if __name__ == "__main__":
    main()
