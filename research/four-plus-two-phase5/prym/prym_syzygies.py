"""Exact arithmetic check of all quadratic/cubic principal-part syzygies at 11."""

from itertools import product
import numpy as np
from flint import nmod_mat, nmod_poly
import prym_seed as seed

P = 11
cut = (1, 3, 5)


def shifts(f, deg):
    return [seed.mul(f, {k: 1}) for k in product(range(deg + 1), repeat=2)]


def selected(fs, deg):
    mons = list(product(range(deg + 1), repeat=2))
    M = nmod_mat([[f.get(k, 0) for k in mons] for f in fs], P)
    R, r = M.transpose().rref()
    return [fs[next(j for j, v in enumerate(row) if v)] for row in R.tolist()[:r]]


def barred(f):
    return {k: v * (-1) ** sum(k) % P for k, v in f.items()}


def line(axis, i):
    return {(0, 0): -i % P, ((1, 0) if axis == 0 else (0, 1)): 1}


def trace_matrix(M, p):
    return nmod_mat((M % p).tolist(), p).rank()


class Axis:
    """The local etale algebra of six points, with their ruling parameter t."""

    def __init__(self, H, axis, val):
        self.axis = axis

        def xy(k):
            return (k[axis], k[1 - axis])

        g = [0] * 7
        for k, c in H.items():
            j, l = xy(k)
            g[l] = (g[l] + c * val**j) % P
        self.g = nmod_poly(g, P)
        assert self.g.degree() == 6 and self.g.gcd(self.g.derivative()).degree() == 0
        self.zero = nmod_poly([], P)
        self.one = nmod_poly([1], P)
        self.x = [nmod_poly([val], P), self.one, self.zero]
        self.y = [nmod_poly([0, 1], P), self.zero, self.zero]
        self.iy = self.inverse(self.g.derivative())
        # Solve H(val+t,y(t))=0 through t^2 in the six-point algebra.
        for i in (1, 2):
            self.y[i] = -self.ev(H)[i] * self.iy % self.g
        assert all(not c for c in self.ev(H))
        self.partial = {}
        for k, c in H.items():
            u = list(k)
            v = u[1 - axis]
            if v:
                u[1 - axis] -= 1
                self.partial[tuple(u)] = v * c % P
        self.omega = self.series_inverse(self.ev(self.partial))
        if axis == 1:
            self.omega = [-c for c in self.omega]
        x = self.y[0]
        self.powers = [self.one]
        for _ in range(20):
            self.powers.append(self.powers[-1] * x % self.g)
        self.ilc = pow(int(self.g.leading_coefficient()), -1, P)

    def inverse(self, f):
        d, s, _ = f.xgcd(self.g)
        assert d == self.one
        return s % self.g

    def times(self, A, B):
        return [
            sum((A[j] * B[i - j] for j in range(i + 1)), self.zero) % self.g
            for i in range(3)
        ]

    def series_inverse(self, A):
        a = self.inverse(A[0])
        return [a, -A[1] * a**2 % self.g, (A[1] ** 2 * a**3 - A[2] * a**2) % self.g]

    def ev(self, f):
        n = max(max(k) for k in f)
        E = [self.one, self.zero, self.zero]
        xp, zp = [E], [E]
        x, z = (self.x, self.y) if self.axis == 0 else (self.y, self.x)
        for _ in range(n):
            xp.append(self.times(xp[-1], x))
            zp.append(self.times(zp[-1], z))
        t = [self.zero] * 3
        for (i, j), v in f.items():
            prod = self.times(xp[i], zp[j])
            t = [t[k] + v * prod[k] for k in range(3)]
        return t

    def trace(self, f, j):
        # For a separable polynomial of degree six, Tr(a/g') is the
        # coefficient of y^5 in a mod g divided by the leading coefficient.
        a = f * self.powers[j] * self.g.derivative() % self.g
        return int(a[5]) * self.ilc % P


def main():
    seed.p = P
    H, q, h, hb, _ = seed.data()
    all_lines = [(axis, i, line(axis, i)) for axis in (0, 1) for i in cut]
    B = seed.mul(*(f for _, _, f in all_lines))
    E = selected(shifts(q, 3) + shifts(hb, 2), 5)
    assert len(E) == 24
    ell = [{key: 1} for key in product(range(2), repeat=2)]
    b = [{key: 1} for key in product(range(3), repeat=2) if sum(key) % 2 == 0]
    quadr = np.zeros((5, 24, 24), dtype=np.int64)
    cubic = np.zeros((4, 24, 24, 24), dtype=np.int64)
    zero_rows = []
    evals = []
    for axis, i, _ in all_lines:
        ring = Axis(H, axis, i)
        Breg = seed.mul(*(f for j, k, f in all_lines if (j, k) != (axis, i)))
        denom = ring.series_inverse(ring.ev(seed.mul(q, Breg)))
        RU = [ring.times(ring.ev(f), denom) for f in E]
        R = np.array(
            [[[int(a[j]) for j in range(6)] for a in v] for v in RU], dtype=np.int64
        )
        evals.append(R[:, 0, :].T)
        # The known holomorphic section has numerator qB and represents s=1.
        unit = ring.times(ring.ev(seed.mul(q, B)), denom)
        assert unit == [ring.zero, ring.one, ring.zero]
        invB = ring.inverse(ring.ev(barred(B))[0])
        quadrows = np.zeros((5, 24, 6), dtype=np.int64)
        qunit = np.zeros((5, 6), dtype=np.int64)
        for k, bb in enumerate(b):
            ll = ring.ev(bb)[0] * ring.omega[0] % ring.g
            s1 = ll * ring.ev(q)[0] % ring.g
            qunit[k] = [ring.trace(s1, t) for t in range(6)]
            for j, f in enumerate(E):
                a = ll * invB * ring.ev(barred(f))[0] % ring.g
                quadrows[k, j] = [ring.trace(a, t) for t in range(6)]
        quadr += np.einsum("ia,kja->kij", R[:, 0, :], quadrows, optimize=True)
        cunits = np.zeros((4, 6), dtype=np.int64)
        K = np.zeros((4, 3, 6, 6, 6), dtype=np.int64)
        for k, ll in enumerate(ell):
            s3 = ring.times(ring.ev(seed.mul(h, ll)), ring.omega)
            cunits[k] = [ring.trace(3 * s3[0], j) for j in range(6)]
            for order in range(3):
                vals = [ring.trace(s3[order], j) for j in range(16)]
                for a, bb, c in product(range(6), repeat=3):
                    K[k, order, a, bb, c] = vals[a + bb + c]
        for a, b1, c, d in product(range(3), repeat=4):
            if a + b1 + c + d == 2:
                cubic += np.einsum(
                    "kxyz,ix,jy,lz->kijl",
                    K[:, a],
                    R[:, b1],
                    R[:, c],
                    R[:, d],
                    optimize=True,
                )
        zero_rows.append(np.vstack([qunit, cunits]))
    # A degree-two expression has a generally nonsymmetric coefficient
    # matrix; polarizing gives all coefficients. The cubic tensor is symmetric.
    assert not np.any((quadr + quadr.transpose(0, 2, 1)) % P)
    assert not np.any(cubic % P)
    relation = np.hstack(zero_rows) % P
    evalmatrix = np.vstack(evals) % P
    assert trace_matrix(relation, P) == 9 and trace_matrix(evalmatrix, P) == 23
    assert not np.any(relation @ evalmatrix % P)
    print("evaluation matrix shape/rank", evalmatrix.shape, trace_matrix(evalmatrix, P))
    print(
        "nine independent syzygies on the 36-dimensional principal parts",
        relation.shape,
        trace_matrix(relation, P),
    )
    print(
        "identically zero polarized quadratic / cubic coefficients",
        5 * (24 * 25 // 2),
        4 * (24 * 25 * 26 // 6),
    )


if __name__ == "__main__":
    main()
