"""Full field-valued trace from one fixed point on the literal smooth F11 contact curve."""

from itertools import product
import numpy as np
import sympy as sp
from one_pole_example import P, b, c, d, q, h, f, A, B, others, seed, lifts


def rref(matrix):
    a = matrix.copy() % P
    piv = []
    r = 0
    for k in range(a.shape[1]):
        i = next((i for i in range(r, len(a)) if a[i, k]), None)
        if i is None:
            continue
        a[[r, i]] = a[[i, r]]
        a[r] = a[r] * pow(int(a[r, k]), -1, P) % P
        z = a[:, k].copy()
        z[r] = 0
        a = (a - z[:, None] * a[r]) % P
        piv.append(k)
        r += 1
        if r == len(a):
            break
    return a, piv


def kernel(a):
    R, pi = rref(a)
    free = [i for i in range(a.shape[1]) if i not in pi]
    out = np.zeros((a.shape[1], len(free)), dtype=np.int64)
    out[free] = np.eye(len(free), dtype=np.int64)
    out[pi] = -R[: len(pi), free] % P
    assert not np.any(a @ out % P)
    return out


def monvals(ring, args, mons):
    caps = np.max(mons, axis=0)
    powers = []
    for x, cap in zip(args, caps):
        row = [ring.one]
        for i in range(cap):
            row.append(ring.mul(row[-1], x))
        powers.append(row)
    return np.array(
        [
            ring.mul(ring.mul(powers[0][i], powers[1][j]), powers[2][k]).ravel()
            for i, j, k in mons
        ]
    ).T


def fixed_space(n):
    # Homogenizing these degree <=n-1 numerators to degree n forces vanishing on a=0.
    gb = sp.groebner([q, h], d, c, b, order="grevlex", modulus=P)
    lm = [tuple(reversed(g.LM(order=gb.order).exponents)) for g in gb.polys]
    mons = [
        e
        for e in product(range(n), repeat=3)
        if sum(e) < n and not any(all(i >= j for i, j in zip(e, m)) for m in lm)
    ]
    r0, aa = lifts(others, n)
    table = monvals(r0, aa, mons)
    basis = kernel(table)
    rp, pa = lifts(seed, n)
    local = monvals(rp, pa, mons) @ basis % P
    return mons, basis, local


def run(n=16):
    mons, space, local = fixed_space(n)
    t = sp.Symbol("t")
    print(
        "pole order",
        n,
        "affine monomials",
        len(mons),
        "one-point functions vanishing at infinity",
        space.shape[1],
        flush=True,
    )
    for lam in (1, 2, 3, 6, 7):
        F = sp.Poly(f.as_expr().subs(b, lam), c, modulus=P).monic()
        assert F.gcd(F.diff()).degree() == 0
        dval = (
            -sp.Poly(B.as_expr().subs(b, lam), c, modulus=P)
            * sp.Poly(A.as_expr().subs(b, lam), c, modulus=P).invert(F)
        ) % F
        C = sp.Poly(c, c, modulus=P)
        res = []
        cp = [C**0]
        dp = [C**0]
        for _ in range(n):
            cp.append(cp[-1] * C % F)
            dp.append(dp[-1] * dval % F)
        for i, j, k in mons:
            co = (pow(lam, i, P) * cp[j] * dp[k]) % F
            res.append([int(co.nth(j)) % P for j in range(6)])
        T = np.array(res).T @ space * pow((lam - 4) ** n, -1, P) % P
        _, piv = rref(T)
        assert len(piv) == 6
        T = T[:, piv]
        Ct = np.zeros((6, 6), dtype=np.int64)
        for j in range(5):
            Ct[j + 1, j] = 1
        Ct[:, -1] = [-int(F.nth(j)) % P for j in range(6)]
        mat = t * sp.Matrix(T) - sp.Matrix(Ct @ T % P)
        # Compute det and the one required adjugate column.  The t-dependent
        # interpolation denominators follow from this system, not from a divisor-root list.
        det = sp.Poly(mat.det(method="domain-ge"), t, modulus=P)
        rhs = []
        for j in range(6):
            rhs.append(
                sp.Poly(
                    (-1) ** j * mat.minor_submatrix(0, j).det(method="domain-ge"),
                    t,
                    modulus=P,
                )
            )
        for i in range(6):
            assert sum(
                (sp.Poly(mat[i, j], t, modulus=P) * rhs[j] for j in range(6)),
                sp.Poly(0, t, modulus=P),
            ) == (det if i == 0 else 0)
        # Negative seed residue, since g has a zero at every pole of b-lambda.
        coeff = [-((-1) ** j) * pow(4 - lam, -j - 1, P) % P for j in range(n)]
        hvals = np.array(coeff[::-1]) @ local[:, piv] % P
        numer = sum((int(z) * h for z, h in zip(hvals, rhs)), sp.Poly(0, t, modulus=P))
        want = sp.Poly(F.as_expr().subs(c, t), t, modulus=P)
        assert det.monic() == want and numer * want == det * want.diff()
        print(
            "b=",
            lam,
            "fixed-seed residue equals full trace, denominator",
            det.monic().as_expr(),
            flush=True,
        )


if __name__ == "__main__":
    run()
