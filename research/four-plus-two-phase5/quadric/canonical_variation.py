"""Differentiate canonical AJ residues for the b^6 and d^6 core variations.

Usage: python canonical_variation.py UPSTREAM_CONTACT_SPAN_DIRECTORY
"""

import numpy as np
import smooth_span as span
import cartier_naturality as residue

p = 11


def branches(N):
    st = span.upstream
    y0, v = span.seed()
    D = st.differential(y0)
    _, cols = st.base_rref(D)
    inv = st.base_inverse(D[:, cols])
    z = np.zeros((N, 90), dtype=np.int64)
    z[0], z[1] = y0, v
    for n in range(2, N):
        rhs = sum(span.quadratic(z[k], z[n - k]) for k in range(1, n))
        z[n, cols] = -inv @ rhs % p
        assert not np.any(
            (sum(span.quadratic(z[k], z[n - k]) for k in range(n + 1))) % p
        )
    return z, cols, inv


def variation(z, cols, inv, exp):
    st = span.upstream
    S = np.zeros(70, dtype=np.int64)
    S[st.plus_basis.index(exp)] = 1
    D = np.stack([st.differential(c) for c in z])
    dz = np.zeros_like(z)
    dz[0, cols] = inv @ S % p
    for n in range(1, len(z)):
        rhs = sum(D[k] @ dz[n - k] for k in range(1, n + 1))
        dz[n, cols] = -inv @ rhs % p
    for n in range(len(z)):
        L = sum(D[k] @ dz[n - k] for k in range(n + 1))
        assert not np.any((L - (S if n == 0 else 0)) % p)
    return S, dz


def compute(z, dz, M, n):
    assert M % p == 0 and n + 1 <= M
    N = M + n
    residue.N = N
    pert = z[: N + 1].copy()
    pert[M:] = (pert[M:] + dz[: N + 1 - M]) % p
    # Independently verify the substituted family equation, coefficient by coefficient.
    diff = span.upstream.differential(z[0]) @ dz[0] % p
    for k in range(N + 1):
        e = sum(
            span.quadratic(pert[i], pert[k - i]) - span.quadratic(z[i], z[k - i])
            for i in range(k + 1)
        )
        assert not np.any((e - (diff if k == M else 0)) % p)
    ref = residue.canonical(z[: N + 1])
    changed = residue.canonical(pert)
    delta = (changed - ref) % p
    assert not np.any(delta[:M])
    return delta[M:]


def modular_difference(cols, exp, n=11):
    # Work with two genuine deformations over Z/121. Their T-coordinate
    # and all the other free affine coordinates are held fixed.
    from unittest.mock import patch

    st = span.upstream
    y0, d = span.seed()
    squared = p * p
    S = np.zeros(70, dtype=np.int64)
    S[st.plus_basis.index(exp)] = 1

    def quad(v, w):
        return (
            np.bincount(
                st.term_row,
                weights=st.term_sign * v[st.term_i] * w[st.term_j],
                minlength=70,
            ).astype(np.int64)
            % squared
        )

    def diff(v):
        D = np.zeros((70, 90), dtype=np.int64)
        np.add.at(D, (st.term_row, st.term_i), st.term_sign * v[st.term_j])
        np.add.at(D, (st.term_row, st.term_j), st.term_sign * v[st.term_i])
        return D % squared

    # The four-variable kappa = 5/7, rather than its residue-class lift 7.
    kappa = 5 * pow(7, -1, squared) % squared
    Hvec = np.zeros(70, dtype=np.int64)
    for m, c in [
        ((6, 0, 0, 0, 0), 1),
        ((0, 6, 0, 0, 0), 1),
        ((0, 0, 6, 0, 0), 1),
        ((0, 0, 0, 6, 0), 1),
        ((0, 0, 0, 0, 3), 1),
        ((5, 1, 0, 0, 0), kappa),
        ((1, 5, 0, 0, 0), 2 * kappa),
        ((0, 0, 5, 1, 0), kappa),
        ((0, 0, 1, 5, 0), 3 * kappa),
        ((3, 1, 1, 1, 0), kappa),
        ((1, 1, 3, 1, 0), 2 * kappa),
    ]:
        Hvec[st.plus_basis.index(m)] = c % squared
    oldinv = st.base_inverse

    def inverse121(matrix):
        iv = oldinv(matrix % p)
        return iv @ (2 * np.eye(len(matrix), dtype=np.int64) - matrix @ iv) % squared

    answers = []
    for step in (0, p):
        target = (Hvec + step * S) % squared
        z = np.zeros((n + 1, 90), dtype=np.int64)
        z[0] = y0
        err = (target - quad(y0, y0)) % squared
        assert not np.any(err % p)
        z[0, cols] = (
            z[0, cols] + p * (oldinv(diff(y0)[:, cols]) @ (err // p))
        ) % squared
        assert np.array_equal(quad(z[0], z[0]), target)
        D = diff(z[0])
        iv = inverse121(D[:, cols])
        z[1] = d
        z[1, cols] = (z[1, cols] - iv @ (D @ d)) % squared
        for degree in range(2, n + 1):
            z[degree, cols] = (
                -iv @ sum(quad(z[i], z[degree - i]) for i in range(1, degree)) % squared
            )
        for degree in range(1, n + 1):
            assert not np.any(
                sum(quad(z[i], z[degree - i]) for i in range(degree + 1)) % squared
            )
        with (
            patch.object(residue, "p", squared),
            patch.object(residue, "N", n),
            patch.object(st, "base_inverse", inverse121),
        ):
            answers.append(residue.canonical(z))
    delta = (answers[1] - answers[0]) % squared
    assert not np.any(delta % p)
    return delta // p


def main():
    n = 28
    z, cols, inv = branches(44 + n + 1)
    mono = span.upstream.minus_basis
    target = mono.index((3, 0, 3, 0, 0))
    for name, coord in [("b^6", (0, 6, 0, 0, 0)), ("d^6", (0, 0, 0, 6, 0))]:
        S, dz = variation(z, cols, inv, coord)
        first = compute(z, dz, 33, n)
        second = compute(z, dz, 44, n)
        assert np.array_equal(first, second)
        assert np.array_equal(first[:11], modular_difference(cols, coord))
        print(name, "target derivative through order 10", first[:11, target].tolist())
        print("D_lambda z(0), number of nonzero coordinates", np.count_nonzero(dz[0]))
        print("derivative row ranks through order 27", span.rank(first)[0])
    return z


def complete_basis(selection):
    import json
    from pathlib import Path

    n = 96
    M = 110
    z, cols, inv = branches(M + n + 1)
    residue.N = n
    C = residue.canonical(z[: n + 1])
    basis = span.upstream.minus_basis
    J = span.upstream.fixed_jacobian_directions(z[0])
    _, piv = span.upstream.base_rref(J.T)
    kept = [i for i in range(60) if i not in piv]
    assert len(kept) == 48
    sheets = {}
    for name, exp in [("b6", (0, 6, 0, 0, 0)), ("d6", (0, 0, 0, 6, 0))]:
        _, dz = variation(z, cols, inv, exp)
        sheets[name] = compute(z, dz, M, n)
    tests = []
    for record in json.loads(Path(selection).read_text())["certificates"]:
        u = record["u_exponent"]
        for term in record["selected_test_numerators"]:
            R = tuple(term["monomial"]) + (u,)
            tests.append(sheets[term["direction"]][:, basis.index(R)])
    X = np.column_stack([C[:, kept]] + tests)
    assert X.shape == (96, 96) and span.det(X) == 5
    print(
        "48 holomorphic and 48 specified pole-three test pullbacks: first 96 jet rows, determinant 5 mod 11"
    )
    return X


if __name__ == "__main__":
    import sys

    main()
    if len(sys.argv) > 2:
        complete_basis(sys.argv[2])
