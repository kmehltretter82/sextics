"""Verify the virtual intersection and rational transverse local lengths.

Run with Python, NumPy and python-flint. Projective q is in an affine patch,
and the nonzero proportionality scalar is absorbed in an etale choice of the
sign of the affine cubic h. All coefficient vectors use the Segre chart.
"""

from itertools import product
from math import comb, gcd
import numpy as np
from flint import fmpz_mat, nmod_mat

P = 101


def porteous():
    value = 0
    for i in range(6):
        for j in (4, 5):
            k = 23 - i - j
            for z in (0, 1):
                m = k - 14 - z
                if 0 <= m <= 4:
                    x = i + j - 4 + m
                    if x + z == 5 and x >= 4:
                        value += (
                            comb(18 + i, i)
                            * 3**j
                            * 2**k
                            * (-1) ** (j - 4 + m)
                            * comb(4, m)
                            * 3**z
                            * (1, 2)[x - 4]
                        )
    assert value == -2210955264
    print("Virtual Porteous integral", value, flush=True)


def mul(f, g):
    result = {}
    for (a, b), v in f.items():
        for (c, d), w in g.items():
            key = (a + c, b + d)
            result[key] = result.get(key, 0) + v * w
    return {k: v for k, v in result.items() if v}


def power(f, k):
    out = {(0, 0): 1}
    for _ in range(k):
        out = mul(out, f)
    return out


def sigma(f):
    return {e: (-1) ** sum(e) * v for e, v in f.items()}


def add(f, g):
    return {k: f.get(k, 0) + g.get(k, 0) for k in f.keys() | g.keys()}


def fixed_Q(i, marked):
    L = {(0, 0): 1, (1, 0): 1, (0, 1): 2, (1, 1): 3}
    Lb = sigma(L)
    q = mul(L, Lb)
    dqs = [{(0, 0): 1, (2, 0): 1, (0, 2): 2, (1, 1): 3, (2, 2): 5}]
    if marked:
        dqs.append(add(L, Lb))
    h = mul(power(L, i), power(Lb, 3 - i))
    hb = sigma(h)
    assert mul(h, hb) == power(q, 3)
    even = [e for e in product(range(7), repeat=2) if sum(e) % 2 == 0]
    hm = list(product(range(4), repeat=2))

    def vector(f):
        return [f.get(k, 0) for k in even]

    columns = [vector(add(mul(hb, {k: 1}), mul(h, sigma({k: 1})))) for k in hm]
    for dq in dqs:
        columns.append([-3 * x for x in vector(mul(power(q, 2), dq))])
    J = np.array(columns, dtype=object).T
    if marked:
        # At p=(1,-2/5), multiply the evaluation equation by 125.
        J = np.vstack([J, [(-2) ** b * 5 ** (3 - b) for _, b in hm] + [0] * len(dqs)])
        pcols = columns[:16] + [
            [-3 * x for x in vector(mul(power(q, 2), {k: 1}))]
            for k in product(range(3), repeat=2)
            if sum(k) % 2 == 0
        ]
        PJ = np.vstack(
            [
                np.array(pcols, dtype=object).T,
                [(-2) ** b * 5 ** (3 - b) for _, b in hm] + [0] * 5,
            ]
        )
        # Remove the common weighted scaling from this full tangent.
        assert 21 - fmpz_mat(PJ.tolist()).rank() - 1 == {1: 7, 2: 8, 3: 3}[i]
    quadratic = [[] for _ in range(len(J))]
    for j, k in enumerate(hm):
        for m, other in enumerate(hm):
            key = (k[0] + other[0], k[1] + other[1])
            if sum(key) % 2 == 0:
                quadratic[even.index(key)].append((j, m, (-1) ** sum(other)))
    cubic = [[] for _ in range(len(J))]
    for a, b in product(range(len(dqs)), repeat=2):
        pp = mul(dqs[a], dqs[b])
        for j, v in enumerate(vector(mul(q, pp))):
            if v:
                quadratic[j].append((16 + a, 16 + b, -3 * v))
        for c, dq in enumerate(dqs):
            for j, v in enumerate(vector(mul(pp, dq))):
                if v:
                    cubic[j].append((a, b, c, -v))
    return J, quadratic, cubic


def columns(matrix):
    A, rank = matrix.rref()
    return [
        next(j for j, value in enumerate(row) if value) for row in A.tolist()[:rank]
    ]


def hom(degree, n):
    if n == 1:
        yield (degree,)
    else:
        for j in range(degree + 1):
            for tail in hom(degree - j, n - 1):
                yield (j,) + tail


def residual_series(J, quadratic, cubic, order, rational):
    m, n = J.shape
    modJ = nmod_mat(J.tolist(), P)
    piv = columns(modJ)
    rows = columns(modJ.transpose())
    free = [k for k in range(n) if k not in piv]
    other = [j for j in range(m) if j not in rows]
    assert fmpz_mat(J.tolist()).rank() == len(piv)
    if not free:
        return [], {}
    M = fmpz_mat(J[np.ix_(rows, piv)].tolist())
    adj, denominator = M.inv().numer_denom()
    D = int(denominator)
    assert D % P
    adj = np.array([[int(c) for c in row] for row in adj.tolist()], dtype=object)
    # The d-th solution coefficients are multiplied by D^(2d-1). Thus the
    # degree-d quadratic convolution uses the previous scaled coefficients
    # directly, while the cubic convolution carries one extra factor of D.
    zero = (0,) * len(free)
    by_degree = [[zero]]
    H, pair, F = {}, {}, {}
    for d in range(1, order + 1):
        terms = list(hom(d, len(free)))
        by_degree.append(terms)
        for e in terms:
            rhs = np.zeros(m, dtype=object)
            sq = np.zeros((n - 16, n - 16), dtype=object)
            cube = np.zeros((n - 16, n - 16, n - 16), dtype=object)
            for x in product(*[range(k + 1) for k in e]):
                if x == zero or x == e:
                    continue
                y = tuple(a - b for a, b in zip(e, x))
                X, Y = H[x], H[y]
                for j, row in enumerate(quadratic):
                    rhs[j] += sum(c * X[k] * Y[other] for k, other, c in row)
                sq += np.outer(X[16:], Y[16:])
                cube += pair[x][:, :, None] * Y[None, None, 16:]
            pair[e] = sq if rational else sq % P
            for j, row in enumerate(cubic):
                rhs[j] += D * sum(v * cube[a, b, c] for a, b, c, v in row)
            v = np.zeros(n, dtype=object)
            if d == 1:
                v[free[e.index(1)]] = D
                v[piv] = -adj @ J[rows, free[e.index(1)]]
            else:
                v[piv] = -adj @ rhs[rows]
            ff = J @ v + D * rhs
            H[e] = v if rational else v % P
            ff = ff if rational else ff % P
            assert all(c == 0 for c in ff[rows])
            F[e] = ff[other]
    return by_degree, F


def macaulay_blocks(by, F, order):
    terms = sum(by[: order + 1], [])
    index = {e: j for j, e in enumerate(terms)}
    num = len(next(iter(F.values())))
    for ds in reversed(range(order - 1)):
        rows = []
        for shift in by[ds]:
            block = [[0] * len(terms) for _ in range(num)]
            for z in sum(by[2 : order - ds + 1], []):
                j = index[tuple(a + b for a, b in zip(shift, z))]
                for k, coeff in enumerate(F[z]):
                    block[k][j] = coeff
            rows.extend(row for row in block if any(row))
        if rows:
            yield rows


def rational_length(by, F, d):
    n = sum(map(len, by[: d + 1]))
    blocks = []
    for batch in macaulay_blocks(by, F, d):
        for row in batch:
            factor = gcd(*row)
            blocks.append([c // factor for c in row])
    if not blocks:
        return n
    mod = nmod_mat([[c % P for c in row] for row in blocks], P)
    selected = columns(mod.transpose())
    S = fmpz_mat([blocks[j] for j in selected])
    U, nullity = S.nullspace()
    assert nullity == n - len(selected)
    K = fmpz_mat([[U[k, j] for j in range(nullity)] for k in range(n)])
    assert (fmpz_mat(blocks) * K).is_zero()
    return nullity


def modular_length(by, F, d):
    n = sum(map(len, by[: d + 1]))
    independent = []
    for batch in macaulay_blocks(by, F, d):
        red, rank = nmod_mat(independent + batch, P).rref()
        independent = red.tolist()[:rank]
    return n - len(independent)


def length(i, marked, low_order, expected):
    J, quadratic, cubic = fixed_Q(i, marked)
    if low_order == 0:
        assert fmpz_mat(J.tolist()).rank() == J.shape[1]
    else:
        by, F = residual_series(J, quadratic, cubic, low_order, True)
        dims = [1] + [rational_length(by, F, d) for d in range(1, low_order + 1)]
        assert dims[-1] == expected
        by, F = residual_series(J, quadratic, cubic, low_order + 1, False)
        # Integral unit pivots at 101 imply char-zero quotient dimension <=
        # modular quotient dimension, while the preceding truncation is a
        # quotient of this truncation. Equality proves Nakayama stabilization.
        assert modular_length(by, F, low_order + 1) == expected
        hilbert = [1] + [b - a for a, b in zip(dims, dims[1:])]
        print(
            "Marked" if marked else "Unmarked",
            "power",
            i,
            "Hilbert",
            hilbert,
            "exact transverse length",
            expected,
            flush=True,
        )
        return
    print(
        "Marked" if marked else "Unmarked",
        "power",
        i,
        "full tangent rank;",
        "exact transverse length",
        expected,
        flush=True,
    )


if __name__ == "__main__":
    porteous()
    length(3, False, 0, 1)
    length(2, False, 5, 49)  # conjugation gives power one
    length(1, True, 5, 49)
    length(2, True, 6, 98)
    length(3, True, 2, 3)
