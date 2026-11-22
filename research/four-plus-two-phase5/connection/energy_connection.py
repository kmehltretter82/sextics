"""Relative plane-curve reduction for the actual affine four-plus-two polynomial."""

from functools import lru_cache

import sympy as sp
from sympy.polys.rings import ring

R, r, energy = ring("r,E", sp.QQ)
k = R(5) / 7


# Polynomials in c,d are sparse dictionaries, with coefficients in Q[r,E].
def add(*pieces):
    result = {}
    for coeff, poly in pieces:
        for exp, value in poly.items():
            v = result.get(exp, R.zero) + coeff * value
            if v:
                result[exp] = v
            else:
                result.pop(exp, None)
    return result


def mul(a, b):
    return add(
        *[
            (v, {(i + u, j + w): z for (u, w), z in b.items()})
            for (i, j), v in a.items()
        ]
    )


def diff(polynomial, direction):
    derivative = {}
    for (i, j), value in polynomial.items():
        if direction == 0 and i:
            derivative[i - 1, j] = i * value
        elif direction == 1 and j:
            derivative[i, j - 1] = j * value
    return derivative


def mon(i, j):
    return {(i, j): R.one}


g = {(6, 0): R.one, (0, 6): R.one, (5, 1): k, (1, 5): 3 * k}
f = 1 + k * r + 2 * k * r**5 + r**6
H = add((1, g), (1, {(0, 0): f, (1, 1): k * r, (3, 1): 2 * k * r}))
hgrad = [diff(H, i) for i in range(2)]
Hminus = add((1, H), (-1, {(0, 0): energy}))
hr = {ij: v.diff(r) for ij, v in H.items() if v.diff(r)}


@lru_cache(maxsize=None)
def leading_basis(n):
    """Select an exact basis of the homogeneous graded piece, keeping gradient columns first."""
    labels = []
    cols = []
    for direction in range(2):
        leading = diff(g, direction)
        for i in range(n - 4):
            j = n - 5 - i
            prod = mul(leading, mon(i, j))
            cols.append([prod.get((h, n - h), R.zero).as_expr() for h in range(n + 1)])
            labels.append((direction, i, j))
    for i in range(n + 1):
        cols.append([int(t == i) for t in range(n + 1)])
        labels.append((2, i, n - i))
    matrix = sp.Matrix(cols).T
    piv = matrix.rref()[1]
    inv = matrix[:, list(piv)].inv().tolist()
    return [labels[j] for j in piv], [[R.from_expr(x) for x in a] for a in inv]


def reduce_form(polynomial, petrov=False):
    """Reduce modulo the gradient ideal or modulo (h-E)L + {h,K}.

    Verify the resulting polynomial identity in Q[r,E,c,d] before returning.
    Only the gradient reduction needs its two lifting polynomials downstream.
    """
    work = polynomial.copy()
    rem = {}
    V = [{}, {}]
    K, L = {}, {}
    for n in reversed(range(max(map(sum, polynomial), default=-1) + 1)):
        top = [work.get((i, n - i), R.zero) for i in range(n + 1)]
        if not any(top):
            continue
        labels, inv = leading_basis(n)
        A, B, basis = {}, {}, {}
        for (group, i, j), row in zip(labels, inv):
            value = sum((c * v for c, v in zip(row, top) if c and v), R.zero)
            if value:
                (A, B, basis)[group][i, j] = value
        rem = add((1, rem), (1, basis))
        work = add((1, work), (-1, basis))
        if not (A or B):
            continue
        if petrov:
            co = R.one / (n - 4)
            lowerK = add((co, mul(mon(0, 1), A)), (-co, mul(mon(1, 0), B)))
            lowerL = add((6 * co, diff(A, 0)), (6 * co, diff(B, 1)))
            K = add((1, K), (1, lowerK))
            L = add((1, L), (1, lowerL))
            work = add(
                (1, work),
                (-1, mul(hgrad[0], diff(lowerK, 1))),
                (1, mul(hgrad[1], diff(lowerK, 0))),
                (-1, mul(Hminus, lowerL)),
            )
        else:
            V[0] = add((1, V[0]), (1, A))
            V[1] = add((1, V[1]), (1, B))
            work = add((1, work), (-1, mul(hgrad[0], A)), (-1, mul(hgrad[1], B)))
        assert not any(sum(m) >= n for m in work)
    assert not work
    if petrov:
        assert polynomial == add(
            (1, rem),
            (1, mul(Hminus, L)),
            (1, mul(hgrad[0], diff(K, 1))),
            (-1, mul(hgrad[1], diff(K, 0))),
        )
        return rem
    assert polynomial == add(
        (1, rem), (1, mul(hgrad[0], V[0])), (1, mul(hgrad[1], V[1]))
    )
    return rem, V


def matrices():
    """Build multiplication, divergence, and parameter-derivative matrices by columns."""
    basis = [
        (a, b) for d in (1, 3, 5, 7) for g0, a, b in leading_basis(d)[0] if g0 == 2
    ]
    assert [sum(i + j == degree for i, j in basis) for degree in (1, 3, 5, 7)] == [
        2,
        4,
        4,
        2,
    ]
    assert (3, 0) in basis
    T, B, C = [sp.zeros(12) for _ in range(3)]
    for j, ij in enumerate(basis):
        p = mon(*ij)
        jrem, V = reduce_form(mul(H, p))
        brem = reduce_form(add((1, diff(V[0], 0)), (1, diff(V[1], 1))), True)
        crem = reduce_form(mul(hr, p), True)
        for matrix, values in zip((T, B, C), (jrem, brem, crem)):
            assert set(values) <= set(basis)
            for m, v in values.items():
                matrix[basis.index(m), j] = v.as_expr()
    assert all(energy.as_expr() not in entry.free_symbols for entry in T)
    return basis, T, B, C


def critical_check(rr, matrix):
    """Verify the critical-value polynomial from the 12 sign pairs at fixed r."""
    s, E = sp.symbols("s E")
    kq = sp.Rational(5, 7)
    gs = 1 + s**6 + kq * s + 3 * kq * s**5
    top = sp.Poly(gs, s, domain=sp.QQ)
    assert top.degree() == 6 and top.gcd(top.diff()).degree() == 0
    A = 3 * gs - s * gs.diff(s)
    b = 2 * s * gs.diff(s) - 3 * gs
    root = sp.Poly(A**2 + 4 * kq * rr * s * b, s, domain=sp.QQ)
    assert root.degree() == 12 and root.gcd(root.diff()).degree() == 0
    root = root.monic()
    bq = sp.Poly(b, s, domain=sp.QQ)
    assert root.eval(0) and root.gcd(bq).degree() == 0
    u = (sp.Poly(A, s, domain=sp.QQ) * (2 * bq).invert(root)) % root
    assert root.gcd(u).degree() == 0

    def q(expr):
        return sp.Poly(expr, s, domain=sp.QQ)

    assert (3 * q(gs) * u**2 + 4 * kq * rr * q(s) * u + kq * rr * q(s)) % root == 0
    assert (q(gs.diff(s)) * u**2 + 2 * kq * rr * u + kq * rr) % root == 0
    fr = 1 + rr**6 + kq * rr + 2 * kq * rr**5
    vals = (q(fr) + q(gs) * u**3 + kq * rr * q(s) * (u + 2 * u**2)) % root
    param_value = q(fr) + (q(A) ** 3 * q(gs - s * gs.diff(s))) * (8 * bq**3).invert(
        root
    )
    assert (vals - param_value) % root == 0
    cp = sp.Poly(sp.resultant(root.as_expr(), E - vals.as_expr(), s), E, domain=sp.QQ)
    assert cp == sp.Poly(matrix.charpoly(E).as_expr(), E, domain=sp.QQ)
    assert cp.gcd(cp.diff()).degree() == 0 and cp.eval(fr) != 0
    return cp


def check():
    basis, T, B, C = matrices()
    x = r.as_expr()
    E = energy.as_expr()
    delta = sp.eye(12) * E - T
    identity = sp.eye(12)
    # Evaluate analytic directional derivatives by the inverse matrix rule.
    for rr, ee in [(1, 2), (2, -1), (-1, 3)]:

        def at(mat):
            return mat.subs({x: rr, E: ee})

        inv = at(delta).inv()
        GE = at(B - identity) * inv
        dGE_r = at(B.diff(x)) * inv + GE * at(T.diff(x)) * inv
        dGE_E = at(B.diff(E)) * inv - GE * inv
        Cr = at(C)
        CE = at(C.diff(E))
        Gr = -GE * Cr - CE
        dGr_E = -dGE_E * Cr - GE * CE - at(C.diff(E, 2))
        assert dGE_r - dGr_E + Gr * GE - GE * Gr == sp.zeros(12)
        # Verify the pole-cleared equation for nonconstant sample X,Z and Y=CX+delta Z.
        X = sp.Matrix([x**i + (i + 1) * E for i in range(12)])
        Z = sp.Matrix([E**2 - (i + 1) * x for i in range(12)])
        twist = -sp.Rational(5, 3 * ee)
        Y = C * X + delta * Z
        lhs = (
            at(X.diff(x)) + Gr * at(X) + at(Y.diff(E)) + (GE + twist * identity) * at(Y)
        )
        rhs = at(X.diff(x)) + Cr * (at(X.diff(E)) + twist * at(X))
        rhs += at(delta) * (at(Z.diff(E)) + twist * at(Z)) + at(B) * at(Z)
        assert lhs == rhs
        cp = critical_check(rr, at(T))
        print(
            "r,E",
            rr,
            ee,
            "simple noncentral sign pairs",
            cp.degree(),
            "flat curvature",
            True,
            flush=True,
        )
    print(
        "odd basis in (c,d)",
        basis,
        "target index",
        basis.index((3, 0)),
        "nonzeros T,B,C",
        *[sum(bool(a) for a in M) for M in (T, B, C)],
        flush=True,
    )


if __name__ == "__main__":
    check()
