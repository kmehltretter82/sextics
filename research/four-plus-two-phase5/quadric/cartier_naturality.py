"""Match ordinary Cartier on a smooth contact curve to scaled pole-two Cartier.

Arguments: upstream contact-span directory; directory containing target_cartier.py.
The latter's run(11) returns (mons, rest, V).
"""

from itertools import permutations
from functools import reduce
from pathlib import Path
import sys
import numpy as np
import smooth_span

p = 11
N = 55


def zero():
    return np.zeros(N, dtype=np.int64)


def one():
    v = zero()
    v[0] = 1
    return v


def mul(a, b):
    return np.convolve(a, b)[:N] % p


def inverse(a):
    out = zero()
    out[0] = pow(int(a[0]), -1, p)
    for j in range(1, N):
        out[j] = -out[0] * (a[1 : j + 1] @ out[:j][::-1]) % p
    assert np.array_equal(mul(a, out), one())
    return out


def add(*entries):
    out = {}
    for form, sign in entries:
        for m, c in form.items():
            out[m] = (out.get(m, 0) + sign * c) % p
    return {m: c for m, c in out.items() if np.any(c)}


def times(a, b):
    f = {}
    for m, c in a.items():
        for n, d in b.items():
            k = tuple(i + j for i, j in zip(m, n))
            f[k] = (f.get(k, 0) + mul(c, d)) % p
    return {m: c for m, c in f.items() if np.any(c)}


def partial(f, j):
    return {
        tuple(i - int(t == j) for t, i in enumerate(m)): m[j] * c % p
        for m, c in f.items()
        if m[j]
    }


def gens(jets):
    st = smooth_span.upstream
    forms = []
    off = 0
    for d, pa in st.blocks:
        f = {}
        for m in st.monomials(d, pa):
            f[m] = jets[:N, off]
            off += 1
        forms.append(f)
    return forms


def canonical(jets):
    st = smooth_span.upstream
    f1, f2, g1, g2, A, B = gens(jets)
    lead = f1[(0, 0, 0, 0, 1)]
    minusQ = {m[:4]: -mul(c, inverse(lead)) % p for m, c in f1.items() if m[4] == 0}
    Qp = [{(0, 0, 0, 0): one()}, minusQ, times(minusQ, minusQ)]

    def substitute(f):
        g = {}
        for m, c in f.items():
            g = add((g, 1), (times({m[:4]: c}, Qp[m[4]]), 1))
        return g

    forms = [substitute(f) for f in (f2, add((A, 1), (B, 1)), g1, g2)]
    djet = jets[1 : N + 1] * np.arange(1, N + 1, dtype=np.int64)[:, None] % p
    *_, Ad, Bd = gens(djet)
    dh = substitute(add((Ad, 1), (Bd, 1)))
    base = smooth_span.mon(9)
    index = {m: i for i, m in enumerate(base)}
    coefficients = []
    edges = []
    row = 0
    for f, d in zip(forms, [2, 3, 4, 4]):
        offset = len(coefficients)
        coefficients.extend(f.values())
        for m in smooth_span.mon(9 - d):
            for k, e in enumerate(f):
                edges.append(
                    (row, index[tuple(i + j for i, j in zip(m, e))], offset + k)
                )
            row += 1
    rowids, colids, entries = np.array(edges).T
    weights = np.array(coefficients, dtype=np.float64)
    M0 = np.zeros((row, len(base)), dtype=np.int64)
    M0[rowids, colids] = weights[entries, 0].astype(np.int64)
    _, piv = st.base_rref(M0.T)
    assert len(piv) == len(base) - 1
    _, cols = st.base_rref(M0[piv])
    (free,) = set(range(len(base))) - set(cols)
    fix = np.zeros(len(base), dtype=np.int64)
    fix[free] = 1
    inv = st.base_inverse(np.vstack((M0[piv], fix)))
    v = np.zeros((N, len(base)), dtype=np.int64)
    v[0] = inv[:, -1]
    # All products are of nonnegative residue representatives. Matrix
    # multiplication accumulates exact integers within binary64 precision.
    assert 100 * N * len(edges) < 2**52
    vf = np.zeros(v.shape, dtype=np.float64)
    vf[0] = v[0]
    for k in range(1, N):
        gram = weights[:, 1 : k + 1] @ vf[k - 1 :: -1]
        acc = np.bincount(rowids, weights=gram[entries, colids], minlength=row).astype(
            np.int64
        )
        v[k] = -(inv[:, :-1] @ acc[piv]) % p
        vf[k] = v[k]
        # The update must annihilate every original relation.
        assert not np.any((M0 @ v[k] + acc) % p)

    def project(f):
        return sum((mul(c, v[:, index[m]]) for m, c in f.items()), zero()) % p

    jac = {}
    for perm in permutations(range(4)):
        sign = (-1) ** sum(x > y for j, x in enumerate(perm) for y in perm[j + 1 :])
        jac = add(
            (jac, 1),
            (reduce(times, (partial(f, i) for f, i in zip(forms, perm))), sign),
        )
    scale = 96 * inverse(mul(lead, project(jac))) % p
    return np.stack(
        [
            mul(scale, project(times(dh, times({m[:4]: one()}, Qp[m[4]]))))
            for m in st.minus_basis
        ],
        axis=1,
    )


def main():
    ell, jets = smooth_span.main()
    sys.path.insert(0, str(Path(sys.argv[2]).resolve()))
    import target_cartier as operator

    mons, rest, V = operator.run(p)
    physical = [
        smooth_span.upstream.minus_basis.index(m[:4] + (m[-1] // 2 - 1,)) for m in mons
    ]
    indices = np.array(physical)[rest]
    E = ell[:, indices]
    C = canonical(jets)
    s = smooth_span.upstream.minus_basis.index((3, 0, 3, 0, 0))
    canon_g = mul(C[:, s], inverse(ell[:N, s]))
    g = canon_g
    fit = np.zeros(N, dtype=np.int64)
    fit[0] = 7
    for j in range((N + 1) // 11):
        k = 11 * j + 10
        cut = 1 if j == 0 else k - 10
        rhs = sum(fit[i] * (E[j - i] @ V) for i in range(j + 1)) - sum(
            fit[i] * E[k - i] for i in range(cut)
        )
        matrix = E[: k - cut + 1][::-1].T
        reduced, pivots = operator.rref(np.column_stack((matrix, rhs)), p)
        assert pivots == list(range(len(matrix[0])))
        fit[cut : k + 1] = reduced[: len(pivots), -1]
    assert np.array_equal(fit, g)
    for i in range(len(ell[0])):
        assert np.array_equal(mul(g, ell[:N, i]), C[:, i])
    for j in range((N + 1) // 11):
        assert np.array_equal(C[11 * j + 10, indices], C[j, indices] @ V % p)
    print("canonical unit through order 10", g[:11].tolist())
    print("target canonical form through order 10", C[:11, s].tolist())
    print(
        "ordinary/scaled Cartier identities verified on 48 quotient generators at output orders 0 through 4."
    )


def full():
    global N
    N = 583
    from smooth_sweep import extend

    params = extend(N + 1)
    jets = np.zeros((N + 1, len(params)), dtype=np.int64)
    for j, v in enumerate(params):
        vals = list(v)
        jets[: len(vals), j] = vals
    C = canonical(jets)
    sys.path.insert(0, str(Path(sys.argv[2]).resolve()))
    import target_cartier as operator

    mons, rest, V = operator.run(p)
    indices = [
        smooth_span.upstream.minus_basis.index(mons[i][:4] + (mons[i][-1] // 2 - 1,))
        for i in rest
    ]
    low = C[:48, indices]
    high = C[10::11, indices]
    recovered = smooth_span.upstream.base_inverse(low) @ high[:48] % 11
    assert np.array_equal(recovered, V)
    assert np.array_equal(C[48:53, indices] @ recovered % 11, high[48:53])
    from fiber_point import det

    print(
        "all 2304 matrix entries recovered from canonical curve jets through order 527; initial-row determinant",
        det(low),
    )
    print(
        "five holdout output rows 48..52 (input degrees through 582) verified on all 48 quotient generators."
    )


if __name__ == "__main__":
    if "--full" in sys.argv:
        full()
    else:
        main()
