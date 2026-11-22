"""Find a function with its only projective pole at the chosen actual mod-11 contact point."""

from itertools import product
import numpy as np
import sympy as sp

P = 11
b, c, d = sp.symbols("b c d")
q = 4 + b + 4 * b * b + 5 * c * c + 10 * c * d + 7 * d * d
h = (
    5
    + 3 * b**3
    + c * c
    + 8 * c * d
    + 9 * d * d
    + 5 * b * c * c
    + 6 * b * c * d
    + 4 * b * d * d
    + 7 * c
    + 9 * d
    + 5 * b * b * c
    + 2 * b * b * d
    + 4 * c**3
    + 5 * c * c * d
    + 6 * d**3
)
f = sp.Poly(sp.resultant(q, h, d), b, c, modulus=P)
f0 = sp.Poly(f.as_expr().subs(b, 4), c, modulus=P).monic()
linear_remainder = sp.Poly(343 * sp.rem(h, q, d), b, c, d, modulus=P).as_expr()
A, B = [
    sp.Poly(sp.expand(linear_remainder).coeff(d, i), b, c, modulus=P) for i in (1, 0)
]
seed = sp.Poly(c - 6, c, modulus=P)
others = f0.exquo(seed)
assert f0.gcd(f0.diff()).degree() == 0
assert f0.gcd(sp.Poly(A.as_expr().subs(b, 4), c, modulus=P)).degree() == 0


class Series:
    def __init__(self, factor, n):
        self.factor = factor
        self.n = n
        self.e = factor.degree()
        self.mod = np.array([int(factor.nth(i)) % P for i in range(self.e)])
        self.zero = np.zeros((n, self.e), dtype=np.int64)
        self.one = self.zero.copy()
        self.one[0, 0] = 1

    def mul(self, x, y):
        out = self.zero.copy()
        for j in range(self.n):
            z = (
                sum(
                    (np.convolve(x[i], y[j - i]) for i in range(j + 1)),
                    np.zeros(2 * self.e - 1, dtype=np.int64),
                )
                % P
            )
            for k in range(len(z) - 1, self.e - 1, -1):
                z[k - self.e : k] -= z[k] * self.mod
                z[k] = 0
                z %= P
            out[j] = z[: self.e]
        return out

    def power(self, x, n):
        ans = self.one.copy()
        for _ in range(n):
            ans = self.mul(ans, x)
        return ans

    def inv(self, x):
        pol = sp.Poly(
            sum(int(v) * c**j for j, v in enumerate(x[0])), c, modulus=P
        ).invert(self.factor)
        ans = self.zero.copy()
        ans[0] = [int(pol.nth(j)) % P for j in range(self.e)]
        for j in range(1, self.n):
            error = self.mul(x, ans)[j]
            e = self.zero.copy()
            e[0] = error
            ans[j] = -self.mul(e, ans)[0] % P
        assert np.array_equal(self.mul(x, ans), self.one)
        return ans

    def eval(self, F, args):
        z = self.zero.copy()
        for powers, v in F.terms():
            term = int(v) * self.one
            for x, n in zip(args, powers):
                term = self.mul(term, self.power(x, n))
            z = (z + term) % P
        return z


def lifts(factor, n):
    ring = Series(factor, n)
    xc = ring.zero.copy()
    xb = 4 * ring.one
    if n > 1:
        xb[1, 0] = 1
    ct = sp.Poly(c, c, modulus=P) % factor
    xc[0] = [int(ct.nth(j)) % P for j in range(ring.e)]
    inv = ring.inv(ring.eval(f.diff(c), [xb, xc]))
    for k in range(1, n):
        dr = ring.eval(f, [xb, xc])[k]
        e = ring.zero.copy()
        e[0] = dr
        xc[k] = -ring.mul(e, inv)[0] % P
    assert not np.any(ring.eval(f, [xb, xc]))
    xd = -ring.mul(ring.eval(B, [xb, xc]), ring.inv(ring.eval(A, [xb, xc]))) % P
    for eq in (q, h):
        assert not np.any(ring.eval(sp.Poly(eq, b, c, d, modulus=P), [xb, xc, xd]))
    return ring, [xb, xc, xd]


def linear(M, rhs):
    a = np.column_stack([M, rhs]) % P
    piv = []
    pivot_row = 0
    for j in range(M.shape[1]):
        k = next((k for k in range(pivot_row, len(a)) if a[k, j]), None)
        if k is None:
            continue
        a[[pivot_row, k]] = a[[k, pivot_row]]
        a[pivot_row] = a[pivot_row] * pow(int(a[pivot_row, j]), -1, P) % P
        for k in range(len(a)):
            if k != pivot_row:
                a[k] = (a[k] - a[k, j] * a[pivot_row]) % P
        piv.append(j)
        pivot_row += 1
        if pivot_row == len(a):
            break
    if np.any(a[pivot_row:, -1]):
        return None
    x = np.zeros(M.shape[1], dtype=np.int64)
    x[piv] = a[:pivot_row, -1]
    assert not np.any((M @ x - rhs) % P)
    return x


def run():
    print("other five points in one finite algebra:", others.as_expr())
    # a=b=0,c!=0 at projective infinity on b=4a would force incompatible roots.
    qinf = sp.Poly(5 + 10 * d + 7 * d * d, d, modulus=P)
    hinf = sp.Poly(4 + 5 * d + 6 * d**3, d, modulus=P)
    assert qinf.gcd(hinf).degree() == 0 and sp.Poly(qinf, d).LC() != 0
    for n in range(1, 7):
        ring, args = lifts(others, n)
        powers = [[ring.power(x, i) for i in range(n + 1)] for x in args]
        mons = sorted(
            (u for u in product(range(n + 1), repeat=3) if sum(u) <= n),
            key=lambda u: (sum(u), u),
        )
        cols = []
        for i, j, k in mons:
            comp = ring.mul(ring.mul(powers[0][i], powers[1][j]), powers[2][k]).ravel()
            cols.append(np.append(comp, pow(4, i, P) * pow(6, j, P) % P))
        M = np.array(cols).T
        rhs = np.zeros(n * 5 + 1, dtype=np.int64)
        rhs[-1] = 1
        x = linear(M, rhs)
        print("pole order", n, "feasible", x is not None)
        if x is not None:
            num = sp.Poly(
                sum(int(v) * b**i * c**j * d**k for v, (i, j, k) in zip(x, mons)),
                b,
                c,
                d,
                modulus=P,
            )
            print("normalized numerator:", num.as_expr())
            assert not np.any(ring.eval(num, args))
            assert int(num.as_expr().subs({b: 4, c: 6, d: 1})) % P == 1
            assert n == 5
            break
    else:
        raise AssertionError("Riemann–Roch control failed")


if __name__ == "__main__":
    run()
