"""One-pole trace with moving hyperplane lambda and resolvent variable t, both indeterminate."""

from itertools import combinations
import numpy as np
import sympy as sp
from flint import nmod_mpoly_ctx
from one_pole_example import P, b, c, d, q, h, f
from one_pole_trace import fixed_space, rref


def det(matrix, ctx):
    # Determinants of successive first-row minors indexed by column subsets.
    current = {(): ctx.constant(1)}
    for row in range(len(matrix)):
        nexts = {}
        for S in combinations(range(len(matrix)), row + 1):
            answer = ctx.constant(0)
            for i, j in enumerate(S):
                answer += (
                    (-1) ** (row + i) * matrix[row][j] * current[S[:i] + S[i + 1 :]]
                )
            nexts[S] = answer
        current = nexts
    return next(iter(current.values()))


def main():
    n = 16
    mons, coeffs, jets = fixed_space(n)
    ctx = nmod_mpoly_ctx.get(["lam", "t"], P)
    lam, t = ctx.gens()
    zero = ctx.constant(0)
    G = sp.groebner([q, h], d, c, b, order="grevlex", modulus=P)
    ph = [(0, 0), (1, 0), (2, 0), (3, 0), (0, 1), (1, 1)]  # c,d exponents

    def vector(expr):
        rem = sp.Poly(G.reduce(expr)[1], d, c, b, modulus=P)
        v = [ctx.constant(0) for _ in ph]
        for (k, j, i), val in rem.terms():
            v[ph.index((j, k))] += int(val) * lam**i
        return v

    pols = [
        sum(int(s) * b**i * c**j * d**k for s, (i, j, k) in zip(col, mons))
        for col in coeffs.T
    ]
    Vs, Ws = [
        list(map(list, zip(*(vector(sp.expand(mult * p)) for p in pols))))
        for mult in (1, c)
    ]
    evalu = np.array([[int(u(1, 0)) for u in row] for row in Vs])
    _, piv = rref(evalu)
    assert len(piv) == 6
    V = [[row[j] for j in piv] for row in Vs]
    W = [[row[j] for j in piv] for row in Ws]
    M = [[t * v - w for v, w in zip(rv, rw)] for rv, rw in zip(V, W)]
    ev = [
        sum((int(x) * (lam - 4) ** k for k, x in enumerate(jets[:, j])), zero)
        for j in piv
    ]
    bordered = [row + [ctx.constant(int(i == 0))] for i, row in enumerate(M)] + [
        [-v for v in ev] + [zero]
    ]
    den = det(M, ctx)
    num = det(bordered, ctx)
    g = det(V, ctx)
    want = ctx.from_dict({e: int(co) % P for e, co in f.terms()})
    # f uses monomials (b,c); replace them by (lambda,t).
    assert g and den == g * want and num == g * want.derivative(1)
    identity = [[ctx.constant(int(i == j)) for j in range(6)] for i in range(6)]

    def matmul(a, b):
        return [
            [sum((a[i][k] * b[k][j] for k in range(6)), zero) for j in range(6)]
            for i in range(6)
        ]

    C = list(map(list, zip(*(vector(c ** (j + 1) * d**k) for j, k in ph))))
    cp = [identity]
    for _ in range(6):
        cp.append(matmul(cp[-1], C))

    def evaluate_at_C(poly):
        terms = [(int(qc) * lam**a, cp[b]) for (a, b), qc in poly.to_dict().items()]
        return [
            [sum((q * A[i][j] for q, A in terms), zero) for j in range(6)]
            for i in range(6)
        ]

    derivative_c = evaluate_at_C(want.derivative(1))
    derivative_lam = evaluate_at_C(want.derivative(0))
    A = matmul(derivative_c, M)
    rhs = [rw[0] for rw in derivative_lam]
    lam_bordered = [rw + [e] for rw, e in zip(A, rhs)] + [[-v for v in ev] + [zero]]
    lam_den = det(A, ctx)
    lam_num = det(lam_bordered, ctx)
    assert lam_den and lam_num * want == lam_den * want.derivative(0)
    assert (num.derivative(0) * den - num * den.derivative(0)) * lam_den**2 == (
        lam_num.derivative(1) * lam_den - lam_num * lam_den.derivative(1)
    ) * den**2
    print(
        "both dlam and dt coefficients of the full logarithmic trace equal d log F; exterior derivative zero, checked over F11(lam,t)"
    )
    print(
        "interpolation columns",
        piv,
        "geometric polynomial degrees",
        want.degrees(),
        "number of terms",
        len(want),
    )
    print(
        "moving interpolation factor degrees",
        g.degrees(),
        "bordered numerator/denominator degrees",
        num.degrees(),
        den.degrees(),
    )
    print(
        "exact polynomial identities det(M)=det(V)*F(lambda,t), bordered=det(V)*dF/dt"
    )


if __name__ == "__main__":
    main()
