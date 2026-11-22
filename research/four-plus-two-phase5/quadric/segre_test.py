"""Two-quadrics determinantal model: check the selected smoothed norm contact."""

from pathlib import Path
from itertools import product
import sys
import numpy as np

P = 23


def pack(a, b=0, c=0, d=0):
    return (a % P) + P * ((b % P) + P * ((c % P) + P * (d % P)))


def coords(z):
    return tuple(z // P**i % P for i in range(4))


def add(a, b):
    return pack(*(x + y for x, y in zip(coords(a), coords(b))))


def neg(a):
    return pack(*(-x for x in coords(a)))


def mul(a, b):
    x, y, z, w = coords(a)
    X, Y, Z, W = coords(b)
    u, t = z * Z + 15 * w * W, z * W + w * Z
    return pack(
        x * X + 15 * y * Y + 14 * u + 15 * 13 * t,
        x * Y + y * X + 14 * t + 13 * u,
        x * Z + 15 * y * W + z * X + 15 * w * Y,
        x * W + y * Z + z * Y + w * X,
    )


def fpow(x, k):
    y = 1
    while k:
        if k % 2:
            y = mul(y, x)
        k //= 2
        x = mul(x, x)
    return y


def inv(x):
    assert x
    return fpow(x, P**4 - 2)


N = 3


class S:
    def __init__(self, v):
        self.v = tuple(v) if not isinstance(v, int) else (v % P, 0, 0)

    def __add__(self, other):
        if not isinstance(other, S):
            other = S(other)
        return S([add(a, b) for a, b in zip(self.v, other.v)])

    __radd__ = __add__

    def __neg__(self):
        return S([neg(a) for a in self.v])

    def __sub__(self, b):
        return self + (-b)

    def __rsub__(self, b):
        return b + (-self)

    def __mul__(self, other):
        if not isinstance(other, S):
            other = S(other)
        accum = [0] * N
        for i in range(N):
            for j in range(N - i):
                accum[i + j] = add(accum[i + j], mul(self.v[i], other.v[j]))
        return S(accum)

    __rmul__ = __mul__

    def inverse(self):
        z = inv(self.v[0])
        out = S([z, 0, 0])
        delta = S(1) - self * out
        return out * (S(1) + delta + delta * delta)

    def __truediv__(self, other):
        if not isinstance(other, S):
            other = S(other)
        return self * other.inverse()

    def __pow__(self, k):
        if not k:
            return S(1)
        return self * self ** (k - 1)


def sqrt(x, root):
    v = S([root, 0, 0])
    assert (v * v).v[0] == x.v[0]
    for _ in range(2):
        v = (v + x / v) / 2
    assert (v * v - x).v == (0, 0, 0)
    return v


ZERO = S(0)


def clean(f):
    return {m: c for m, c in f.items() if any(c.v)}


def fadd(*args):
    out = {}
    for f in args:
        for m, c in f.items():
            out[m] = out.get(m, ZERO) + c
    return clean(out)


def scale(f, c):
    return clean({m: b * c for m, b in f.items()})


def fmul(f, g):
    out = {}
    for (i, j), c in f.items():
        for (k, l), d in g.items():
            m = i + k, j + l
            out[m] = out.get(m, ZERO) + c * d
    return clean(out)


def power(f, n):
    if not n:
        return {(0, 0): S(1)}
    return fmul(f, power(f, n - 1))


def form(terms, x):
    acc = {}
    for e, coef in terms:
        v = {(0, 0): coef}
        for i, k in enumerate(e):
            if k:
                v = fmul(v, power(x[i], k))
        acc = fadd(acc, v)
    return acc


QBASE = [(2, 0, 0, 0), (0, 2, 0, 0), (0, 0, 2, 0), (0, 0, 1, 1), (0, 0, 0, 2)]
AE = [
    (3, 0, 0, 0),
    (0, 3, 0, 0),
    (1, 0, 2, 0),
    (1, 0, 1, 1),
    (1, 0, 0, 2),
    (0, 1, 2, 0),
    (0, 1, 1, 1),
    (0, 1, 0, 2),
]
BO = [
    (2, 0, 1, 0),
    (2, 0, 0, 1),
    (0, 2, 1, 0),
    (0, 2, 0, 1),
    (0, 0, 3, 0),
    (0, 0, 2, 1),
    (0, 0, 1, 2),
    (0, 0, 0, 3),
]


def source(x):
    terms = [
        ((6, 0, 0, 0), 1),
        ((0, 6, 0, 0), 1),
        ((0, 0, 6, 0), 1),
        ((0, 0, 0, 6), 1),
        ((5, 1, 0, 0), 4),
        ((1, 5, 0, 0), 8),
        ((0, 0, 5, 1), 4),
        ((0, 0, 1, 5), 12),
        ((3, 1, 1, 1), 4),
        ((1, 1, 3, 1), 8),
    ]
    return form([(m, S(c)) for m, c in terms], x)


def seed_chart(data):
    def get(i, offset=0):
        return S([pack(*map(int, x)) for x in data[offset : offset + N, i]])

    q0, q1 = get(5), get(6)
    A, B, C = [get(i, 1) for i in (7, 8, 9)]
    disc = B * B - 4 * A * C
    root = sqrt(disc, pack(0, 0, 1))
    tau = S([0, 1, 0])
    theta = (1 + sqrt(1 - 4 * q0 * q1, 1)) / 2
    d = {(1, 0): S(1) / (2 * root), (0, 1): S(-1) / (2 * root)}
    c = scale(fadd({(1, 0): S(1), (0, 1): S(1)}, scale(d, -2 * B)), S(1) / (4 * A))
    b = scale({(1, 1): -tau / (4 * A), (0, 0): -q0}, S(1) / (theta - q0 * q1 / theta))
    a = fadd({(0, 0): S(1)}, scale(b, -q1 / theta))
    x = [a, b, c, d]
    q = form(
        list(zip(QBASE, [get(i) for i in range(5, 10)])) + [((1, 1, 0, 0), S(1))], x
    )
    assert q == {}
    G = form(zip(QBASE, [get(i) for i in range(5)]), x)
    F = fadd(source(x), scale(power(G, 3), S(-1)))
    HA = form(zip(AE, [get(i) for i in range(10, 18)]), x)
    HB = form(zip(BO, [get(i) for i in range(18, 26)]), x)
    h = fadd(HA, HB)
    k = fadd(HA, scale(HB, -1))
    assert not fadd(F, scale(fmul(h, k), -1))
    return x, G, F, h, k


def diff(p, coord):
    out = {}
    for key, v in p.items():
        i = key[coord]
        if i:
            out[(key[0] - (coord == 0), key[1] - (coord == 1))] = v * i
    return clean(out)


def logarithmic(x, y, i):
    return fadd(fmul(y, diff(x, i)), scale(fmul(x, diff(y, i)), -1))


def kernel_map(f, p, q):
    return fadd(
        fmul(f, fadd(diff(q, 0), scale(diff(p, 1), -1))),
        fmul(p, diff(f, 1)),
        scale(fmul(q, diff(f, 0)), -1),
    )


def matrix(f):
    directions = [
        (j, (a, b))
        for j, (da, db) in enumerate(((4, 6), (6, 4)))
        for a, b in product(range(da + 1), range(db + 1))
        if (a + b) % 2 == 0
    ]
    basis = [(a, b) for a, b in product(range(11), repeat=2) if (a + b) % 2 == 1]
    assert (len(directions), len(basis)) == (36, 60)
    index = {k: i for i, k in enumerate(basis)}
    mat = np.zeros((len(basis), len(directions)), dtype=np.int64)
    for j, (c, key) in enumerate(directions):
        omega = [{}, {}]
        omega[c] = {key: S(1)}
        v = kernel_map(f, *omega)
        for k, t in v.items():
            mat[index[k], j] = t.v[0]
    return mat, basis, directions


def rank_pivots(A):
    m, n = A.shape
    B = A.copy()
    rows = list(range(m))
    columns = []
    rowids = []
    for j in range(n):
        i = len(columns)
        k = next((t for t in range(i, m) if B[t, j]), None)
        if k is None:
            continue
        B[[i, k]] = B[[k, i]]
        rows[i], rows[k] = rows[k], rows[i]
        rowids.append(rows[i])
        columns.append(j)
        mult = inv(int(B[i, j]))
        for col in range(j, n):
            B[i, col] = mul(mult, int(B[i, col]))
        for row in range(i + 1, m):
            mlt = neg(int(B[row, j]))
            if mlt:
                for col in range(j, n):
                    B[row, col] = add(int(B[row, col]), mul(mlt, int(B[i, col])))
    return rowids, columns


def main():
    arr = np.load(sys.argv[1])["coordinates"]
    x, G, F, h, k = seed_chart(arr)
    p, q = [logarithmic(h, k, i) for i in range(2)]
    assert not kernel_map(F, p, q)
    mat, basis, directions = matrix(F)
    rows, cols = rank_pivots(mat)
    assert len(cols) == 35
    print("selected seed rank", len(cols), "rows", rows, "columns", cols)
    print(
        "kernel degrees",
        max(i for i, j in p),
        max(j for i, j in p),
        max(i for i, j in q),
        max(j for i, j in q),
    )
    # A different smooth quadric, preserving the exact H, with Q1 = 0.
    x = [{(0, 0): S(1)}, {(1, 1): S(-1)}, {(1, 0): S(1)}, {(0, 1): S(1)}]
    ctrl = source(x)
    ctrl_matrix, _, _ = matrix(ctrl)
    rr, cc = rank_pivots(ctrl_matrix)
    assert len(cc) == 36
    print("non-contact control rank", len(cc))


if __name__ == "__main__":
    main()
