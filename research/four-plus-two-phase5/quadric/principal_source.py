"""Source residue at the smooth seed for the formal compact-pairing principal part."""

from math import comb
from itertools import product
import numpy as np
import canonical_variation as cv
import cartier_naturality as sr
from principal_solution import main as primitive

p = 11


def cst(v):
    return int(v) * sr.one() % p


def spow(v, k):
    if not k:
        return sr.one()
    u = spow(sr.mul(v, v), k // 2)
    return sr.mul(v, u) if k % 2 else u


def collect(*terms):
    n = max(len(A) for _, A in terms)
    out = np.zeros((n, sr.N), dtype=np.int64)
    for s, A in terms:
        out[: len(A)] += s * np.array(A, dtype=np.int64)
    return out % p


def prod(A, B):
    out = np.zeros((len(A) + len(B) - 1, sr.N), dtype=np.int64)
    for i, a in enumerate(A):
        for j, b in enumerate(B):
            out[i + j] += sr.mul(a, b)
    return out % p


def shift(A, k):
    return np.concatenate((np.zeros((k, sr.N), dtype=np.int64), A))


def eval_ser(f, z):
    out = sr.zero()
    for c in f[::-1]:
        out = (sr.mul(out, z) + c) % p
    return out


def ssolve(mat, b):
    inv = cv.span.upstream.base_inverse(mat[0])
    out = np.zeros_like(b)
    for k in range(sr.N):
        acc = np.einsum("nmq,nq->m", mat[1 : k + 1], out[:k][::-1])
        out[k] = inv @ (b[k] - acc) % p
    for k in range(sr.N):
        assert not np.any(
            (np.einsum("nmq,nq->m", mat[: k + 1], out[: k + 1][::-1]) - b[k]) % p
        )
    return out


def intrinsic(z):
    f1, f2, _, _, A, B = sr.gens(z)
    ab = (1, 1, 0, 0, 0)
    um = (0, 0, 0, 0, 1)
    det = sr.mul(f1[um], f2[ab]) - sr.mul(f1[ab], f2[um])
    iv = sr.inverse(det)

    def physical(first, second, a, b):
        return {
            e[:4]: sr.mul(iv, (sr.mul(a, first[e]) - sr.mul(b, second[e])) % p)
            for e in first
            if e[4] == 0
        }

    Q = physical(f1, f2, f2[ab], f1[ab])
    q = physical(f2, f1, f1[um], f2[um])
    beta = {m: -w % p for m, w in Q.items()}
    polys = sr.add((A, 1), (B, 1))
    h = {}
    for m, coeff in polys.items():
        g = {m[:4]: coeff}
        h = sr.add((h, 1), (sr.times(g, beta) if m[4] else g, 1))
    # Canonical cubic representative has neither a*b*linear nor weighted u terms.
    linear = list(product(range(2), repeat=4))
    linear = [a for a in linear if sum(a) == 1]
    U = np.zeros((sr.N, 4, 4), dtype=np.int64)
    rhs = np.zeros((sr.N, 4), dtype=np.int64)
    for i, m in enumerate(linear):
        pivot = tuple(a + b for a, b in zip(m, ab[:4]))
        rhs[:, i] = h.get(pivot, sr.zero())
        for j, L in enumerate(linear):
            U[:, i, j] = q.get(tuple(a - b for a, b in zip(pivot, L)), sr.zero())
    b = ssolve(U, rhs)
    h = sr.add((h, 1), (sr.times(q, dict(zip(linear, b.T))), -1))
    assert all(not (m[0] and m[1]) for m in h)
    return q, h, beta


def eval_partials(f, r, s):
    n = 1 + max(sum(m[2:]) for m in f)
    out = np.zeros((3, n, sr.N), dtype=np.int64)
    for (_, i, j, k), v in f.items():
        out[0, j + k] += v * (pow(r, i, p) * pow(s, k, p) % p)
        if i:
            out[1, j + k] += v * (i * pow(r, i - 1, p) * pow(s, k, p) % p)
        if k:
            out[2, j + k] += v * (k * pow(r, i, p) * pow(s, k - 1, p) % p)
    return out % p


def kernels(q, h, beta, r, s):
    (z, zr, zs), (a, ar, ass), (b, br, bs) = [
        eval_partials(f, r, s) for f in [q, h, beta]
    ]

    def da(v):
        return np.array([i * c % p for i, c in enumerate(v) if i])

    er = collect((1, prod(da(z), ar)), (-1, prod(da(a), zr)))
    es = collect((1, prod(da(z), ass)), (-1, prod(da(a), zs)))
    ez = collect((1, prod(zr, ass)), (-1, prod(zs, ar)))
    ey = collect((1, prod(er, bs)), (-1, prod(es, br)), (1, prod(da(b), ez)))
    entries = []
    for i, (v, pr, ps) in enumerate(
        [
            (7 * r * s, 7 * s, 7 * r),
            (14 * r * s, 14 * s, 14 * r),
            (1 + s**6 + 7 * (s + 3 * s**5), 0, 6 * s**5 + 7 * (1 + 15 * s**4)),
        ],
        1,
    ):
        a1 = collect(
            (3 * v, prod(b, ey)), (-1, prod(prod(b, b), collect((ps, er), (-pr, es))))
        )
        entries.append(collect((1, shift(a1, 1)), (-2 * i * v, prod(prod(b, b), ez))))
    tri = [
        collect((1, entries[0]), (1, shift(entries[1], 2)), (1, shift(entries[2], 4))),
        collect((1, entries[1]), (1, shift(entries[2], 2))),
        entries[2],
    ]
    return z, a, tri


def polyrem(f, a):
    rem = f.copy()
    iv = sr.inverse(a[-1])
    for j in range(len(rem) - len(a), -1, -1):
        factor = sr.mul(rem[j + len(a) - 1], iv)
        for k in range(len(a)):
            rem[j + k] = (rem[j + k] - sr.mul(factor, a[k])) % p
    return rem[: len(a) - 1]


def bezout(q, h, s):
    V = [q.get((0, 0, 2 - j, j), sr.zero()) for j in range(3)]
    b = [h.get((0, 0, 3 - j, j), sr.zero()) for j in range(4)]
    mat = np.zeros((sr.N, 5, 5), dtype=np.int64)
    for i in range(3):
        for j in range(3):
            mat[:, i + j, i] = V[j]
    for i in range(2):
        for j in range(4):
            mat[:, i + j, 3 + i] = b[j]
    rhs = np.zeros((sr.N, 5), dtype=np.int64)
    rhs[0, 0] = 1
    sol = ssolve(mat, rhs).T
    A, B = sol[:3], sol[3:]
    U = sr.mul(eval_ser(A, cst(s)), eval_ser(V, cst(s)))
    C = sr.mul(eval_ser(B, cst(s)), eval_ser(b, cst(s)))
    assert np.array_equal((U + C) % p, sr.one())
    return (
        sum(comb(19, i) * sr.mul(spow(U, i), spow(C, 19 - i)) for i in range(6, 20)) % p
    )


def sample(q, h, beta, r, s, coeff):
    z, a, tri = kernels(q, h, beta, r, s)
    u, v = z[0], z[2]
    even_remainder = (sr.mul(v, a[0]) - sr.mul(u, a[2])) % p
    odd_remainder = (sr.mul(v, a[1]) - sr.mul(u, a[3])) % p
    g = (sr.mul(v, spow(even_remainder, 2)) + sr.mul(u, spow(odd_remainder, 2))) % p
    vb = sr.mul(v, a[3])
    if not (vb[0] and g[0]):
        return None
    weight = bezout(q, h, s)
    out = []
    for source in tri:
        even, odd = polyrem(source, z)
        Kq = sr.mul(
            sr.mul(v, sr.inverse(g)),
            sr.mul(odd_remainder, even) - sr.mul(even_remainder, odd),
        )
        # At infinity the residue correction equals the leading remainder coefficient.
        infinity = sr.mul(polyrem(source, prod(z, a))[-1], sr.inverse(vb))
        Kh = (Kq + infinity) % p
        # Independent inverse in the cubic algebra gives the same Kh.
        c0, c1, c2 = [sr.mul(y, sr.inverse(a[3])) for y in a[:3]]
        uv = sr.mul(u, v)
        vv = sr.mul(v, v)
        adj = [
            (
                spow(u, 2)
                + sr.mul(uv, spow(c2, 2) - 2 * c1)
                + sr.mul(vv, spow(c1, 2) - sr.mul(c0, c2))
            )
            % p,
            (sr.mul(vv, sr.mul(c1, c2) - c0)) % p,
            (sr.mul(vv, c1) - uv) % p,
        ]
        check = sr.mul(sr.mul(a[3], sr.inverse(g)), polyrem(prod(source, adj), a)[2])
        assert np.array_equal(check, Kh)
        K = (sr.mul(weight, Kq) + sr.mul(sr.one() - weight, Kh)) % p
        out.append(int((4 * pow(3, -1, p) * (np.arange(1, 97) * coeff @ K[1:])) % p))
    return out


def main():
    c, _ = primitive()
    sr.N = 97
    z, _, _ = cv.branches(sr.N)
    q, h, beta = intrinsic(z)
    out = {}
    for r, s in [(4, 6), (1, 2), (2, 3), (6, 5), (1, 8), (3, 9)]:
        val = sample(q, h, beta, r, s, c)
        if val is not None:
            out[(r, s)] = val
    assert out == {
        (4, 6): [3, 3, 6],
        (1, 2): [8, 0, 10],
        (2, 3): [6, 8, 5],
        (6, 5): [9, 2, 5],
        (1, 8): [3, 7, 5],
        (3, 9): [9, 9, 1],
    }
    print(
        "formal principal-part seed source (includes 4/3, excludes global degree):",
        out,
    )


if __name__ == "__main__":
    main()
