"""Test actual reduced residue numerators modulo the 49 resultant sections, using coupled jets."""

import sys
from pathlib import Path
from itertools import product
import importlib.util
from contextlib import redirect_stdout
from io import StringIO
from tempfile import TemporaryDirectory

import numpy as np

from deformed_trace import socle_coefficients
from fiber_degree import PRIME, multiply, monomials, resultant_jets, determinant


def power(a, n):
    result = np.zeros_like(a)
    result[0, 0] = 1
    for _ in range(n):
        result = multiply(result, a)
    return result


def pp_add(*terms):
    length = max(len(p) for k, p in terms)
    out = np.zeros((length, *terms[0][1].shape[1:]), dtype=np.int64)
    for k, p in terms:
        out[: len(p)] += k * p
    return out % PRIME


def pp_mul(p, q):
    out = np.zeros((len(p) + len(q) - 1, *p.shape[1:]), dtype=np.int64)
    for i, x in enumerate(p):
        for j, y in enumerate(q):
            out[i + j] += multiply(x, y)
    return out % PRIME


def shift(p, k):
    return np.pad(p, ((k, 0), (0, 0), (0, 0)))


def diff_alpha(p):
    return np.array([j * x % PRIME for j, x in enumerate(p) if j])


def evaluate(coeffs, monomials, r, s):
    n = len(coeffs)
    ans = np.zeros((3, 1 + max(i[2] + i[3] for i in monomials), n, 2), dtype=np.int64)
    for k, (_, b, c, d) in enumerate(monomials):
        mult = r**b * s**d
        dr = b * r ** (b - 1) * s**d if b else 0
        ds = d * r**b * s ** (d - 1) if d else 0
        for slot, z in enumerate((mult, dr, ds)):
            ans[slot, c + d] += (z % PRIME) * coeffs[:, k]
    return ans % PRIME


def fixed_contacts(coords, intrinsic, r, s):
    len(coords)
    one = power(coords[:, 0], 0)
    qcoeff = np.concatenate((one[:, None, :], coords[:, 5:10, :]), axis=1)
    q = evaluate(qcoeff, [(1, 1, 0, 0)] + intrinsic.QBASE, r, s)
    h = evaluate(coords[:, 10:26, :], intrinsic.ABASE + intrinsic.BBASE, r, s)
    beta = evaluate(-coords[:, :5, :], intrinsic.QBASE, r, s)
    return q, h, beta


def linear_reduction(q, h):
    u, v = q[0, 0], q[0, 2]
    even = multiply(v, h[0, 0]) - multiply(u, h[0, 2])
    odd = multiply(v, h[0, 1]) - multiply(u, h[0, 3])
    g = (multiply(v, multiply(even, even)) + multiply(u, multiply(odd, odd))) % PRIME
    return u, v, even, odd, g


def kernel_polynomials(q, h, beta, r, s):
    a, hr, hs = h
    z, qr, qs = q
    b, br, bs = beta
    qa, ha, ba = (diff_alpha(i) for i in (z, a, b))
    er = pp_add((1, pp_mul(qa, hr)), (-1, pp_mul(ha, qr)))
    es = pp_add((1, pp_mul(qa, hs)), (-1, pp_mul(ha, qs)))
    ez = pp_add((1, pp_mul(qr, hs)), (-1, pp_mul(qs, hr)))
    ey = pp_add((1, pp_mul(er, bs)), (-1, pp_mul(es, br)), (1, pp_mul(ba, ez)))
    b2 = pp_mul(b, b)
    S = []
    for i, (p, pr, ps) in enumerate(
        (
            (4 * r * s, 4 * s, 4 * r),
            (8 * r * s, 8 * s, 8 * r),
            (1 + s**6 + 4 * s + 12 * s**5, 0, 6 * s**5 + 4 + 60 * s**4),
        ),
        1,
    ):
        A = pp_add((3 * p, pp_mul(b, ey)), (-1, pp_mul(b2, pp_add((ps, er), (-pr, es)))))
        S.append(pp_add((1, shift(A, 1)), (-2 * i * p, pp_mul(b2, ez))))
    return [
        pp_add((1, S[0]), (1, shift(S[1], 2)), (1, shift(S[2], 4))),
        pp_add((1, S[1]), (1, shift(S[2], 2))),
        S[2],
    ]


def clear_kernel(poly, u, v, even, odd, den_degree):
    """Return V^den_degree times V*(O*rem_even−E*rem_odd), using alpha²=−U/V."""
    n = len(v)
    even_remainder = np.zeros((n, 2), dtype=np.int64)
    odd_remainder = even_remainder.copy()
    for j, coef in enumerate(poly):
        b = j // 2
        term = multiply(coef, multiply(power(-u, b), power(v, den_degree - b)))
        if j % 2:
            odd_remainder += term
        else:
            even_remainder += term
    return multiply(v, multiply(odd, even_remainder) - multiply(even, odd_remainder))


def gaussian_pivots(matrix):
    a = matrix.copy()
    ids = list(range(len(a)))
    cols = []
    for col in range(a.shape[1]):
        j = len(cols)
        eligible = np.flatnonzero(np.any(a[j:, col] != 0, axis=1))
        if not len(eligible):
            continue
        k = j + int(eligible[0])
        ids[j], ids[k] = ids[k], ids[j]
        a[[j, k]] = a[[k, j]]
        x, y = a[j, col]
        inv = pow(int((x * x - 15 * y * y) % 23), -1, 23)
        b, c = inv * x % 23, -inv * y % 23
        p, w = a[j, :, 0].copy(), a[j, :, 1].copy()
        A = (b * p + 15 * c * w) % 23
        B = (b * w + c * p) % 23
        d, e = a[j + 1 :, col, 0, None].copy(), a[j + 1 :, col, 1, None].copy()
        a[j + 1 :, :, 0] = (a[j + 1 :, :, 0] - d * A - 15 * e * B) % 23
        a[j + 1 :, :, 1] = (a[j + 1 :, :, 1] - d * B - e * A) % 23
        cols.append(col)
    return sorted(ids[: len(cols)]), cols


def run(archive):
    norm = archive / "algebra/intrinsic_norm_incidence.py"
    spec = importlib.util.spec_from_file_location("upstream_norm", norm)
    intrinsic = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(intrinsic)
    source = archive / "contact-span/contact_span_coefficients.npz"
    raw = np.load(source, allow_pickle=False)["branch"]
    with TemporaryDirectory() as scratch:
        intrinsic.HERE = Path(scratch)
        intrinsic.INPUT = source
        intrinsic.N = len(raw)
        intrinsic.ONE = intrinsic.const(1)
        with redirect_stdout(StringIO()):
            intrinsic.main()
        coords = np.load(intrinsic.HERE / "intrinsic_branch_four_jets.npz", allow_pickle=False)[
            "coordinates"
        ]
    assert coords.shape == (67, 26, 2)
    # Raw Q = delta*Qintr and H = f1_u*Hintr + (linear form)*Qintr.
    base = monomials(2, 0)
    u = base.index((0, 0, 0, 0, 1))
    ab = base.index((1, 1, 0, 0, 0))
    a, p, b, f2_ab = raw[:, u], raw[:, ab], raw[:, len(base) + u], raw[:, len(base) + ab]
    delta = (multiply(a, f2_ab) - multiply(b, p)) % PRIME
    gauge = multiply(power(delta, 3), power(a, 2))
    gs = []
    for r, s in product(range(7), repeat=2):
        q, h, _ = fixed_contacts(coords, intrinsic, r, s)
        g = linear_reduction(q, h)[-1]
        assert np.array_equal(multiply(g, gauge), resultant_jets(raw, r, s))
        gs.append(g)
    sys.path.insert(0, str(archive / "contact-span"))
    import contact_series_span as span

    intrinsic.N = 4
    intrinsic.ONE = intrinsic.const(1)
    eta = np.zeros((6, 2), dtype=np.int64)
    eta[:4] = socle_coefficients(archive, intrinsic, span)
    for r, s in ((1, 2), (4, 6)):
        q, h, beta = fixed_contacts(coords, intrinsic, r, s)
        u, v, even, odd, g = linear_reduction(q, h)
        S = kernel_polynomials(q, h, beta, r, s)
        den_degree = max((len(poly) - 1) // 2 for poly in S)
        print("sample", r, s, "degree", den_degree)
        section_cols = [multiply(power(v, den_degree), b) for b in gs]
        src = [clear_kernel(poly, u, v, even, odd, den_degree) for poly in S]
        matrix = np.stack(section_cols + src, axis=1)
        rows, cols = gaussian_pivots(matrix)
        det = determinant(matrix[rows][:, cols])
        first_row = 4 if (r, s) == (1, 2) else 3
        expected = [8, 17] if (r, s) == (1, 2) else [6, 5]
        assert rows == list(range(first_row, first_row + 52)) and cols == list(range(52))
        assert det == expected
        orders = [int(np.flatnonzero(np.any(x != 0, axis=1))[0]) - den_degree for x in src]
        print(
            "52 by 52 determinant:", det, "first row:", first_row, "numerator valuations:", orders
        )
        den = multiply(power(v, 6), g)
        assert not den[:6].any() and den[6].any()
        # With pole orders at most three, the displayed residue needs only
        # the canonical input coefficients through degree two.
        assert min(orders) >= -3
        unit = den[6:12]
        a0, b0 = unit[0]
        inv0 = np.array([a0, -b0]) * pow(int((a0 * a0 - 15 * b0 * b0) % PRIME), -1, PRIME) % PRIME
        inverse = np.zeros_like(unit)
        inverse[0] = inv0
        for degree in range(1, 6):
            residual = multiply(unit, inverse)[degree]
            inverse[degree] = -multiply(inv0[None, :], residual[None, :])[0] % PRIME
        assert np.array_equal(multiply(unit, inverse), power(unit, 0))
        scaled = [
            (4 * pow(3, -1, PRIME) * multiply(multiply(eta, num[:6]), inverse)[5] % PRIME).tolist()
            for num in src
        ]
        expected_residues = (
            [[14, 1], [22, 22], [0, 0]] if (r, s) == (1, 2) else [[16, 8], [19, 8], [11, 22]]
        )
        assert scaled == expected_residues
        print("three scaled seed residues for input A:", scaled)


if __name__ == "__main__":
    run(Path(sys.argv[1]))
