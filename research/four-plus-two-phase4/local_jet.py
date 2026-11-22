"""Two exact boundary equations for each fiber of the chosen nodal branch."""

from pathlib import Path
from argparse import ArgumentParser
import numpy as np

p = 23
v2 = 15


def encode(a, b=0):
    return int(a) % p + p * (int(b) % p)


def mulfield(a, b):
    a0, a1 = a % p, a // p
    b0, b1 = b % p, b // p
    return encode(a0 * b0 + v2 * a1 * b1, a0 * b1 + a1 * b0)


def addfield(a, b):
    return encode(a % p + b % p, a // p + b // p)


def negfield(a):
    return encode(-int(a % p), -int(a // p))


mul = [[mulfield(a, b) for b in range(p * p)] for a in range(p * p)]
add = [[addfield(a, b) for b in range(p * p)] for a in range(p * p)]
neg = [negfield(i) for i in range(p * p)]
inv = [0] + [next(b for b in range(1, p * p) if mul[a][b] == 1) for a in range(1, p * p)]


# truncated power series mod (r^3,dS^2); order is (r^0 dS^0,r^0 dS^1,...).
class J:
    def __init__(self, u):
        self.f = (encode(u), 0, 0, 0, 0, 0) if isinstance(u, int) else tuple(u)

    def __add__(self, other):
        if not isinstance(other, J):
            other = J(other)
        return J([add[a][b] for a, b in zip(self.f, other.f)])

    __radd__ = __add__

    def __neg__(self):
        return J([neg[a] for a in self.f])

    def __sub__(self, other):
        return self + -other

    def __rsub__(self, other):
        return -self + other

    def __mul__(self, other):
        if not isinstance(other, J):
            other = J(other)
        out = [0] * 6
        for n in range(3):
            for s_power in range(2):
                ix = 2 * n + s_power
                for i in range(n + 1):
                    for k in range(s_power + 1):
                        out[ix] = add[out[ix]][
                            mul[self.f[2 * i + k]][other.f[2 * (n - i) + s_power - k]]
                        ]
        return J(out)

    __rmul__ = __mul__

    def inv(self):
        h = self.f[0]
        if not h:
            raise ZeroDivisionError
        c = J([inv[h], 0, 0, 0, 0, 0])
        d = 1 - c * self
        return c * (1 + d * (1 + d * (1 + d)))

    def __truediv__(self, other):
        if not isinstance(other, J):
            other = J(other)
        return self * other.inv()

    def __rtruediv__(self, other):
        return other * self.inv()

    def __pow__(self, n):
        ans = J(1)
        while n:
            if n % 2:
                ans = ans * self
            self = self * self
            n //= 2
        return ans


def field(a):
    return J([a, 0, 0, 0, 0, 0])


def intrinsic_data(jets_file):
    jets = np.load(jets_file)["coordinates"]
    keys = [(0, 0, 0), (2, 0, 0), (0, 0, 2), (0, 1, 2), (0, 2, 2)]
    aval = [(0, 0, 0), (3, 0, 0), (0, 0, 2), (0, 1, 2), (0, 2, 2), (1, 0, 2), (1, 1, 2), (1, 2, 2)]
    bval = [(0, 0, 1), (0, 1, 1), (2, 0, 1), (2, 1, 1), (0, 0, 3), (0, 1, 3), (0, 2, 3), (0, 3, 3)]

    def data(keys, i, sign=1):
        terms = []
        for j, key in enumerate(keys):
            for n in range(5):
                coef = encode(*(sign * jets[n, i + j]))
                if coef:
                    terms.append((n, *key, coef))
        return terms

    return [data(keys, 5) + [(0, 1, 0, 0, 1)], data(aval, 10) + data(bval, 18), data(keys, 0, -1)]


def poly_and_deriv(terms, coords):
    powers = [[x**i for i in range(5)] for x in coords]
    vals = [J(0) for _ in range(5)]
    for n, r, s, a, coef in terms:
        exp = (n, r, s, a)
        for i in range(-1, 4):
            if i >= 0 and not exp[i]:
                continue
            val = field(coef) * (exp[i] if i >= 0 else 1)
            for j in range(4):
                e = exp[j] - (i == j)
                if e:
                    val = val * powers[j][e]
            vals[i + 1] += val
    return vals


def sources(dat, S):
    R = J([0, 0, 1, 0, 0, 0])
    s = J([S, 1, 0, 0, 0, 0])
    B = encode(4 + 6 * S + 17 * S**2 + 22 * S**3)
    if not B:
        raise ZeroDivisionError
    roots = [i for i in range(1, p * p) if mul[B][mul[i][mul[i][i]]] == 22]
    assert len(roots) == 3
    pk = [1 + R**6 + 4 * R + 8 * R**5, 4 * R * s, 8 * R * s, 1 + s**6 + 4 * s + 12 * s**5]
    pr = [6 * R**5 + 4 + 40 * R**4, 4 * s, 8 * s, J(0)]
    ps = [J(0), 4 * R, 8 * R, 6 * s**5 + 4 + 60 * s**4]
    tr = [[J(0) for _ in range(4)] for _ in range(2)]
    for root in roots:
        tau = J(0)
        alpha = field(root)
        for step in range(3):
            coord = (tau, R, s, alpha)
            q, qt, qr, qs, qa = poly_and_deriv(dat[0], coord)
            h, ht, hr, hs, ha = poly_and_deriv(dat[1], coord)
            j = qt * ha - ht * qa
            tau, alpha = tau - (q * ha - h * qa) / j, alpha - (qt * h - ht * q) / j
        coord = (tau, R, s, alpha)
        q, qt, qr, qs, qa = poly_and_deriv(dat[0], coord)
        h, ht, hr, hs, ha = poly_and_deriv(dat[1], coord)
        assert not any(q.f + h.f)
        beta, bt, br, bs, ba = poly_and_deriv(dat[2], coord)
        j = qt * ha - ht * qa
        Er, Es, Ez = qa * hr - ha * qr, qa * hs - ha * qs, qr * hs - qs * hr
        Ey = Er * bs - Es * br + ba * Ez
        Svec = [
            alpha * (3 * pi * beta * Ey - beta**2 * (Er * p_s - Es * p_r))
            - 2 * i * pi * beta**2 * Ez
            for i, (pi, p_r, p_s) in enumerate(zip(pk, pr, ps))
        ]
        sk = [
            Svec[1] + alpha**2 * Svec[2] + alpha**4 * Svec[3],
            Svec[2] + alpha**2 * Svec[3],
            Svec[3],
            J(0),
        ]
        for n in range(2):
            for i in range(4):
                tr[n][i] += 4 * tau**n * sk[i] / (9 * j)
    return tr, pk, pr, ps, R, s


def obstruction(k, pk, pr, ps, R, s):
    f, P, Q, g = pk
    fr = pr[0]
    gs = ps[3]
    A = 3 * (k[2] / 2 - k[1]) / (16 * s)  # r*a since kappa=4.
    den = P - 4 * s * f / fr - 4 * R * g / gs
    L = (
        R * k[1] - 4 * P * A / 3 - 4 * s * R * k[0] / fr - 4 * R * (R * k[3] - 4 * g * A) / gs
    ) / den
    B = 3 * (R * k[0] - f * L) / (2 * fr)
    C = 3 * (R * k[3] - 4 * g * A - g * L) / (2 * gs)
    assert B.f[0] == 0 and L.f[0] == 0
    return [add[add[L.f[0]][A.f[0]]][C.f[1]], add[add[add[L.f[2]][A.f[2]]][B.f[4]]][C.f[3]]]


def obstruction_by_cramer(k, pk, pr, ps, R, s):
    """Independent jet calculation of r J by the adjugate of the moment matrix."""
    from itertools import permutations

    def determinant(matrix):
        value = J(0)
        for perm in permutations(range(4)):
            term = J(1)
            for i, col in enumerate(perm):
                term = term * matrix[i][col]
            if sum(a > b for i, a in enumerate(perm) for b in perm[i + 1 :]) % 2:
                term = -term
            value += term
        return value

    M = [
        [pi, 4 * i * pi / 3, 2 * p_r / 3, 2 * p_s / 3]
        for i, (pi, p_r, p_s) in enumerate(zip(pk, pr, ps))
    ]
    denom = determinant(M)
    ns = [
        determinant([[k[i] if j == col else M[i][j] for j in range(4)] for i in range(4)])
        for col in range(4)
    ]

    def ds(x):
        return J([x.f[1], 0, x.f[3], 0, x.f[5], 0])

    def dr(x):
        return J([x.f[2], x.f[3], mul[2][x.f[4]], mul[2][x.f[5]], 0, 0])

    num = denom * (ns[0] + ns[1] + dr(ns[2]) + ds(ns[3])) - dr(denom) * ns[2] - ds(denom) * ns[3]
    assert not denom.f[0] and not ns[2].f[0] and not num.f[0]

    # Since denom=r*D, r*J=(num/r)/D². Only the first two coefficients are needed.
    def over_r(x):
        return J([x.f[2], x.f[3], x.f[4], x.f[5], 0, 0])

    result = over_r(num) / (over_r(denom) ** 2)
    return [result.f[0], result.f[2]]


def run(path):
    data = intrinsic_data(path)
    points = []
    rows = []
    for S in range(1, p):
        try:
            tr, pk, pr, ps, R, s = sources(data, S)
            # Four x_j*dH/dx_i outside the base-current kernel on this chart,
            # in order (i,j)=(a,c),(a,d),(b,c),(b,d).
            eac = [J(0)] + [(6 - 2 * i) * pk[i] - R * pr[i] for i in range(3)]
            assert not any(((6 - 6) * pk[3] - R * pr[3]).f)
            ebc = [J(0)] + pr[:3]
            eq = [eac, [s * a for a in eac], ebc, [s * a for a in ebc]]
            tar = [J(0), J(0), J(1), J(0)]
            cc = [obstruction(k, pk, pr, ps, R, s) for k in eq + tr + [tar]]
            if S <= 4:
                for k, cert in zip(eq + tr + [tar], cc):
                    assert obstruction_by_cramer(k, pk, pr, ps, R, s) == cert
            for row in zip(*cc):
                rows.append(list(row))
            points.append(S)
        except ZeroDivisionError:
            pass

    def rank(a):
        a = [row[:] for row in a]
        piv = []
        pivot_row = 0
        rowids = list(range(len(a)))
        inds = []
        det = 1
        for j in range(len(a[0])):
            k = next((k for k in range(pivot_row, len(a)) if a[k][j]), None)
            if k is None:
                continue
            a[pivot_row], a[k] = a[k], a[pivot_row]
            rowids[pivot_row], rowids[k] = rowids[k], rowids[pivot_row]
            lead = a[pivot_row][j]
            det = mul[det][lead]
            a[pivot_row] = [mul[v][inv[lead]] for v in a[pivot_row]]
            for k in range(pivot_row + 1, len(a)):
                if a[k][j]:
                    lead = neg[a[k][j]]
                    a[k] = [add[u][mul[lead][v]] for u, v in zip(a[k], a[pivot_row])]
            pivot_row += 1
            piv.append(j)
            inds.append(rowids[pivot_row - 1])
        return piv, inds, det

    for n in (4, 5, 6, 7):
        piv, indices, det = rank([row[:n] for row in rows])
        print(
            "rank,columns,rows,minor for first",
            n,
            ":",
            len(piv),
            piv,
            [(points[i // 2], i % 2) for i in indices],
            (det % p, det // p),
        )
    pivot_rows = [(1, 0), (1, 1), (2, 0), (2, 1), (3, 0), (3, 1), (4, 0)]
    minor = [rows[2 * points.index(S) + n] for S, n in pivot_rows]
    import itertools

    det = 0
    for sigma in itertools.permutations(range(7)):
        term = 1
        for i, j in enumerate(sigma):
            term = mul[term][minor[i][j]]
        if sum(a > b for i, a in enumerate(sigma) for b in sigma[i + 1 :]) % 2:
            term = neg[term]
        det = add[det][term]
    assert det == encode(21, 11)
    for row in minor:
        print([(a % p, a // p) for a in row])
    print("Leibniz determinant:", det % p, det // p)
    return minor


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("jets_file", type=Path)
    run(parser.parse_args().jets_file)
