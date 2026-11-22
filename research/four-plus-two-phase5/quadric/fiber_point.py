"""Independent direct verification of a rational contact of the actual coupled core modulo 11."""

from itertools import product
import numpy as np

p = 11
k = 7
xbase = [2, 6, 1, 7, 7, 4, 8, 4, 3, 0]
hbase = np.array(
    [4, 4, 9, 4, 6, 0, 6, 9, 9, 0, 10, 0, 7, 3, 9, 6], dtype=np.int64
).reshape(4, 4)


class Dual:
    def __init__(self, val, d=0):
        self.val = val % p
        self.d = d % p

    def __add__(self, b):
        if not isinstance(b, Dual):
            b = Dual(b)
        return Dual(self.val + b.val, self.d + b.d)

    __radd__ = __add__

    def __neg__(self):
        return Dual(-self.val, -self.d)

    def __sub__(self, b):
        return self + -b

    def __rsub__(self, b):
        return -self + b

    def __mul__(self, b):
        if not isinstance(b, Dual):
            b = Dual(b)
        return Dual(self.val * b.val, self.val * b.d + self.d * b.val)

    __rmul__ = __mul__

    def __truediv__(self, b):
        if not isinstance(b, Dual):
            b = Dual(b)
        inv = pow(b.val, -1, p)
        return self * Dual(inv, -inv * inv * b.d)

    def __pow__(self, n):
        if n == 0:
            return Dual(1)
        return Dual(pow(self.val, n, p), n * pow(self.val, n - 1, p) * self.d)


def chart(x, w, z):
    m, n, u, v, e = x[:5]
    K = e * (v - u) ** 2 * (m + n) / (n - m) ** 2
    return 1 + K * w * z, n + m * K * w * z, w + z, v * w + u * z


def potential(x, w, z):
    a, b, c, d = chart(x, w, z)
    H = (
        a**6
        + b**6
        + c**6
        + d**6
        + k
        * (
            a**5 * b
            + 2 * a * b**5
            + c**5 * d
            + 3 * c * d**5
            + a * b * c * d * (a**2 + 2 * c**2)
        )
    )
    P, Q, R, T, S = x[5:]
    g = P * a * a + Q * b * b + R * c * c + T * c * d + S * d * d
    return H - g**3


def cubic(w, z):
    return sum(int(hbase[i, j]) * w**i * z**j for i, j in product(range(4), repeat=2))


def rank(a):
    a = np.array(a, dtype=np.int64) % p
    row = 0
    cols = []
    rowids = list(range(len(a)))
    pivots = []
    for j in range(a.shape[1]):
        match = next((i for i in range(row, len(a)) if a[i, j]), None)
        if match is None:
            continue
        a[[row, match]] = a[[match, row]]
        rowids[row], rowids[match] = rowids[match], rowids[row]
        pivots.append(rowids[row])
        cols.append(j)
        a[row] = a[row] * pow(int(a[row, j]), -1, p) % p
        for i in range(row + 1, len(a)):
            a[i] = (a[i] - a[i, j] * a[row]) % p
        row += 1
    return len(cols), pivots, cols


def det(rows):
    a = np.array(rows, dtype=np.int64) % p
    g = 1
    for row in range(len(a)):
        i = next((i for i in range(row, len(a)) if a[i, row]), None)
        if i is None:
            return 0
        if i != row:
            a[[row, i]] = a[[i, row]]
            g = -g
        b = int(a[row, row])
        g = g * b % p
        for i in range(row + 1, len(a)):
            a[i] = (a[i] - a[i, row] * pow(b, -1, p) * a[row]) % p
    return g


def check_curve_smooth():
    """Test the four affine Jacobian ideals over an algebraic closure."""
    import sympy as sp

    w, z = sp.symbols("w z")
    for reverse_w, reverse_z in product((False, True), repeat=2):
        h = sum(
            int(hbase[i, j])
            * w ** (3 - i if reverse_w else i)
            * z ** (3 - j if reverse_z else j)
            for i, j in product(range(4), repeat=2)
        )
        basis = sp.groebner([h, sp.diff(h, w), sp.diff(h, z)], w, z, modulus=11)
        assert [q.as_expr() for q in basis.polys] == [1]
    print(
        "bidegree-(3,3) curve is smooth over the algebraic closure on all four charts."
    )


def main():
    check_curve_smooth()
    fixed = [Dual(z) for z in xbase]
    each = [[Dual(z, int(i == j)) for i, z in enumerate(xbase)] for j in range(10)]
    norm = []
    for w, z in product(range(7), repeat=2):
        val = potential(fixed, w, z)
        h = cubic(w, z) % p
        b = cubic(-w, -z) % p
        assert val.val == h * b % p
        dg = [potential(x, w, z).d for x in each]
        dh = [
            -(w**i) * z**j * (b + (-1) ** (i + j) * h)
            for i, j in product(range(4), repeat=2)
        ]
        norm.append(dg + dh + [0, 0])
    assert rank(norm)[0] == 25
    print("direct stable contact on 7x7 grid, norm Jacobian rank 25.")
    found = False
    for w, z in product(range(p), repeat=2):
        a, b, c, d = chart(fixed, w, z)
        if not a.val or not c % p or cubic(w, z) % p:
            continue
        r = b.val * pow(a.val, -1, p) % p
        s = d.val * pow(c, -1, p) % p
        if r in xbase[:2] or s in xbase[2:4]:
            continue
        tangent = []
        for j in range(12):
            x = each[j] if j < 10 else fixed
            a, b, c, d = chart(x, Dual(w, int(j == 10)), Dual(z, int(j == 11)))
            tangent.append(((b - r * a).d, (d - s * c).d))
        hrow = (
            [0] * 10
            + [w**i * z**j for i, j in product(range(4), repeat=2)]
            + [cubic(Dual(w, 1), Dual(z)).d, cubic(Dual(w), Dual(z, 1)).d]
        )
        f = (
            np.array(
                norm
                + [
                    [t[q] for t in tangent[:10]]
                    + [0] * 16
                    + [t[q] for t in tangent[10:]]
                    for q in range(2)
                ]
                + [hrow],
                dtype=np.int64,
            )
            % p
        )
        rk, rows, _ = rank(f)
        if rk == 28:
            if (r, s, w, z) == (4, 2, 2, 10):
                assert det(f[rows]) == 8
                found = True
            print(
                "r,s",
                r,
                s,
                "w,z",
                w,
                z,
                "alpha",
                c.val * pow(a.val, -1, p) % p,
                "fiber rank",
                rk,
                "det",
                det(f[rows]),
                "rows",
                rows,
            )
    assert found


if __name__ == "__main__":
    main()
