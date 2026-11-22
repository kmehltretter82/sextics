#!/usr/bin/env python3
"""Four terms of the normalized local trace for the canonical input a³c³.

Argument: extracted four-plus-two-current directory from Discussion 120.
Uses the original packet's polynomial arithmetic and its two archived jet tables.
"""

import argparse
from itertools import permutations
from pathlib import Path
import sys

import numpy as np


def socle_coefficients(packet, chart, span):
    """Compute Res(a³c³ ∂τh) using the original presentation before eliminating u."""
    n = chart.N
    archive = np.load(packet / "contact-span/contact_span_coefficients.npz")
    original = archive["branch"]
    generators = []
    offset = 0
    for degree, parity in span.blocks:
        basis = span.monomials(degree, parity)
        generators.append({m: original[:n, offset + i].copy() for i, m in enumerate(basis)})
        offset += len(basis)
    f1, f2, g1, g2, even, odd = generators
    lead = f1[(0, 0, 0, 0, 1)]
    quad = chart.pscale({m[:4]: c for m, c in f1.items() if not m[4]}, chart.sinv(lead))
    h = chart.padd((even, 1), (odd, 1))
    forms = [chart.subst_u(f, quad) for f in (f2, h, g1, g2)]
    # h is differentiated in the original five-variable presentation.
    derivative = {}
    offset = len(span.monomials(2, 0)) * 2 + len(span.monomials(4, 0)) * 2
    for parity in (0, 1):
        for m in span.monomials(3, parity):
            derivative[m] = original[1 : n + 1, offset].copy() * np.arange(1, n + 1)[:, None] % 23
            offset += 1
    dh = chart.subst_u(derivative, quad)

    # In degree nine, normalize a left null functional of the ideal to have
    # value one on a free monomial. The Jacobian fixes its remaining scale.
    basis = chart.mon(9)
    index = {m: i for i, m in enumerate(basis)}
    rows = []
    for form, degree in zip(forms, (2, 3, 4, 4)):
        for monomial in chart.mon(9 - degree):
            row = np.zeros((len(basis), n, 2), dtype=np.int64)
            for m, coefficient in form.items():
                row[index[chart.madd(m, monomial)]] = coefficient
            rows.append(row)
    matrix = np.stack(rows)
    assert not np.any(matrix[:, :, 0, 1])
    _, pivots = span.base_rref(matrix[:, :, 0, 0].T)
    assert len(pivots) == len(basis) - 1
    _, pivot_monomials = span.base_rref(matrix[pivots, :, 0, 0])
    (free,) = set(range(len(basis))) - set(pivot_monomials)
    extra = np.zeros_like(matrix[0])
    extra[free] = chart.ONE
    square = np.concatenate((matrix[pivots], extra[None, ...]))
    inverse = span.base_inverse(square[:, :, 0, 0])
    ell = np.zeros((len(basis), n, 2), dtype=np.int64)
    ell[:, 0, 0] = inverse[:, -1]
    for k in range(1, n):
        rhs = sum(span.matmul(square[:, :, j], ell[:, k - j]) for j in range(1, k + 1))
        ell[:, k] = -inverse @ rhs % 23
    for k in range(n):
        assert not np.any(
            sum(span.matmul(matrix[:, :, j], ell[:, k - j]) for j in range(k + 1)) % 23
        )

    def project(poly):
        return sum((chart.smul(ell[index[m]], c) for m, c in poly.items()), chart.zero()) % 23

    def partial(poly, i):
        return {
            tuple(x - int(j == i) for j, x in enumerate(m)): c * m[i] % 23
            for m, c in poly.items()
            if m[i]
        }

    jacobian = {}
    for perm in permutations(range(4)):
        product = {(0, 0, 0, 0): chart.ONE}
        for i, j in enumerate(perm):
            product = chart.pmul(product, partial(forms[i], j))
        inversions = sum(x > y for i, x in enumerate(perm) for y in perm[i + 1 :])
        jacobian = chart.padd((jacobian, 1), (product, (-1) ** inversions))
    scale = 96 * chart.sinv(chart.smul(lead, project(jacobian))) % 23
    sources = [tuple(m) for m in archive["minus_monomials"]]
    target = sources.index((3, 0, 3, 0, 0))
    canonical = {}
    for i, m in enumerate(sources):
        if m[4] == 0:
            product = {chart.madd(e, m[:4]): c for e, c in dh.items()}
            canonical[i] = chart.smul(scale, project(product))
    assert len(canonical) == 40
    assert canonical[target][0].tolist() == [5, 5]
    # Independent obstruction rows in the packet are normalized differently.
    # All 40 canonically normalized source jets must differ by the same unit.
    obstruction = archive["obstruction"][:n]
    unit = chart.smul(canonical[target], chart.sinv(obstruction[:, target]))
    for i, jet in canonical.items():
        assert np.array_equal(jet, chart.smul(unit, obstruction[:, i]))
    return canonical[target]


def local_source(chart, coordinates, normalized, s):
    """Lift all three simple roots of q=h=0, evaluating in F_529[r]/r^4."""
    n, one = chart.N, chart.ONE
    mul, inv, const = chart.smul, chart.sinv, chart.const
    r = chart.zero()
    r[1, 0] = 1
    bases = (chart.QBASE, [chart.AB] + chart.QBASE, chart.ABASE + chart.BBASE)

    def power(x, k):
        out = one
        for _ in range(k):
            out = mul(out, x)
        return out

    def compose(coefficients, tau):
        out = chart.zero()
        for c in coefficients[::-1]:
            out = mul(out, tau)
            out[0] = (out[0] + c) % 23
        return out

    def evaluate(tau, alpha):
        coefficients = [compose(c, tau) for c in coordinates.swapaxes(0, 1)]
        tangent = np.arange(1, len(coordinates))[:, None, None] * coordinates[1:]
        derivatives = [compose(c % 23, tau) for c in tangent.swapaxes(0, 1)]
        all_values = []
        for base, cvals, dvals in zip(
            bases,
            (coefficients[:5], [one] + coefficients[5:10], coefficients[10:]),
            (derivatives[:5], [chart.zero()] + derivatives[5:10], derivatives[10:]),
        ):
            vals = [chart.zero() for _ in range(5)]
            for m, coeff, dcoeff in zip(base, cvals, dvals):
                _, i, j, k = m
                alpha_degree = j + k
                for col, (constant, rp, sp, ap) in enumerate(
                    (
                        (1, i, k, alpha_degree),
                        (i, i - 1, k, alpha_degree),
                        (k, i, k - 1, alpha_degree),
                        (alpha_degree, i, k, alpha_degree - 1),
                        (1, i, k, alpha_degree),
                    )
                ):
                    if constant:
                        monomial = mul(power(r, rp), power(alpha, ap))
                        term = mul(dcoeff if col == 4 else coeff, monomial)
                        vals[col] = (vals[col] + constant * s**sp * term) % 23
            all_values.append(vals)
        return all_values

    B = (4 + 6 * s + 17 * s * s + 22 * s**3) % 23
    assert B
    roots = []
    for a in range(23):
        for b in range(23):
            z = np.array([a, b])
            if np.array_equal((B * chart.fmul(z, chart.fmul(z, z))) % 23, [22, 0]):
                roots.append(z)
    assert len(roots) == 3
    total = np.zeros((3, n, 2), dtype=np.int64)
    g = 1 + s**6 + 4 * s + 12 * s**5
    gs = 6 * s**5 + 4 + 60 * s**4
    for root in roots:
        tau, alpha = chart.zero(), chart.zero()
        alpha[0] = root
        for k in range(1, n):
            _, (q, _, _, qa, qt), (h, _, _, ha, ht) = evaluate(tau, alpha)
            jac = (mul(qt, ha) - mul(ht, qa)) % 23
            detinv = inv(jac)[0]
            tau[k] = chart.fmul(detinv, chart.fmul(qa[0], h[k]) - chart.fmul(ha[0], q[k]))
            alpha[k] = chart.fmul(detinv, chart.fmul(ht[0], q[k]) - chart.fmul(qt[0], h[k]))
        (beta, br, bs, ba, _), (q, qr, qs, qa, qt), (h, hr, hs, ha, ht) = evaluate(tau, alpha)
        assert not np.any(q) and not np.any(h)
        beta, br, bs, ba = [-v % 23 for v in (beta, br, bs, ba)]
        p0 = (one + power(r, 6) + 4 * r + 8 * power(r, 5)) % 23
        Pa = (
            p0
            + mul(4 * s * r, power(alpha, 2))
            + mul(8 * s * r, power(alpha, 4))
            + g * power(alpha, 6)
        ) % 23
        assert not np.any((power(beta, 3) + Pa) % 23)
        jac = (mul(qt, ha) - mul(ht, qa)) % 23
        er, es = (mul(qa, hr) - mul(ha, qr)) % 23, (mul(qa, hs) - mul(ha, qs)) % 23
        ez = (mul(qr, hs) - mul(qs, hr)) % 23
        ey = (mul(er, bs) - mul(es, br) + mul(ba, ez)) % 23
        entries = []
        for i, (p, pr, ps) in enumerate(
            (
                (4 * s * r, const(4 * s), 4 * r),
                (8 * s * r, const(8 * s), 8 * r),
                (const(g), const(0), const(gs)),
            ),
            start=1,
        ):
            term = 3 * mul(p, mul(beta, ey)) - mul(power(beta, 2), mul(er, ps) - mul(es, pr))
            entries.append((mul(alpha, term) - 2 * i * mul(p, mul(power(beta, 2), ez))) % 23)
        last = entries[2]
        middle = (entries[1] + mul(power(alpha, 2), last)) % 23
        first = (entries[0] + mul(power(alpha, 2), middle)) % 23
        scalar = mul(compose(normalized, tau), inv(jac))
        for i, value in enumerate((first, middle, last)):
            total[i] += mul(scalar, value)
    # The current has factor 4/3 and the local trace has factor 1/3.
    return 4 * pow(9, -1, 23) * total % 23


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("packet", type=Path)
    parser.add_argument("--s", type=int, nargs="+", default=[1, 2, 3, 4])
    args = parser.parse_args()
    sys.path[:0] = [str(args.packet / sub) for sub in ("algebra", "contact-span")]
    import intrinsic_norm_incidence as chart
    import contact_series_span as span

    chart.N = 4
    chart.ONE = chart.const(1)
    coordinates = np.load(args.packet / "algebra/intrinsic_branch_four_jets.npz")["coordinates"]
    eta = socle_coefficients(args.packet, chart, span)
    print("eta_A / d(tau), in ascending powers of tau:", eta.tolist())
    for s in args.s:
        values = local_source(chart, coordinates, eta, s % 23)
        print("s =", s, "source rows (1, t², t⁴), each in ascending powers of r:", values.tolist())


if __name__ == "__main__":
    main()
