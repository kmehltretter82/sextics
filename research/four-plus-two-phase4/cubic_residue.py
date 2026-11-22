#!/usr/bin/env python3
"""Compute a moving-point residue representative regular at the generic seed.

For q = U+V*alpha^2 and h cubic with leading coefficient b, write g=Res(q,h).
The positive Jacobian is j = (Dq)*h_alpha-(Dh)*q_alpha.  This computes a
representative whose residue at each simple common root is N/j.
"""

import sympy as sp

a = sp.Symbol("alpha")


def cubic_numerator(N, U, V, h):
    """Return L_h and g so that sum_{h(x)=0} N(x)/(q(x)*h'(x)) = L_h/g."""
    hp = sp.Poly(h, a, domain="EX")
    b = hp.LC()
    assert hp.degree() == 3
    c0, c1, c2 = [hp.nth(i) / b for i in range(3)]
    monic = a**3 + c2 * a**2 + c1 * a + c0
    A0 = U**2 + U * V * (c2**2 - 2 * c1) + V**2 * (c1**2 - c0 * c2)
    A1 = V**2 * (c1 * c2 - c0)
    A2 = V**2 * c1 - U * V
    adj = A0 + A1 * a + A2 * a**2
    norm = U**3 + U**2 * V * (c2**2 - 2 * c1) + U * V**2 * (c1**2 - 2 * c0 * c2) + V**3 * c0**2
    reduced = sp.rem(sp.Poly(N * adj, a, domain="EX"), sp.Poly(monic, a, domain="EX"))
    return sp.factor(b * reduced.nth(2)), sp.factor(b**2 * norm)


def conic_numerator(N, U, V, h):
    """The conic-based representative from the Phase 4 trace calculation."""
    h = sp.Poly(h, a)
    even, odd = V * h.nth(0) - U * h.nth(2), V * h.nth(1) - U * h.nth(3)
    red = sp.rem(sp.Poly(N, a, domain="EX"), sp.Poly(a**2 + U / V, a, domain="EX"))
    return sp.factor(V * (odd * red.nth(0) - even * red.nth(1)))


def bezout_patch(V, b, A, C, m, n):
    """Unnormalized patch weights; their sum is (A*V+C*b)^(m+n-1)."""
    M = m + n - 1
    wq = sum(sp.binomial(M, i) * (A * V) ** i * (C * b) ** (M - i) for i in range(m, M + 1))
    wh = sum(sp.binomial(M, i) * (A * V) ** i * (C * b) ** (M - i) for i in range(m))
    return wq, wh


def check_patch():
    # Every binomial term is divisible by a required clearing power.
    for m, n in ((1, 1), (3, 8), (6, 14)):
        M = m + n - 1
        assert all(i >= m or M - i >= n for i in range(M + 1))
    # Genuine moving s resultant, whose zeros are fixed values of tau.
    tau, s = sp.symbols("tau s")
    V, b = s - tau, s**2 + s + 1
    A, C = -(s + tau + 1), sp.Integer(1)
    rho = tau**2 + tau + 1
    assert sp.expand(A * V + C * b - rho) == 0
    assert sp.resultant(V, b, s) == rho
    wq, wh = bezout_patch(V, b, A, C, 2, 3)
    assert sp.expand(wq + wh - rho**4) == 0
    assert sp.denom(sp.cancel(wq / V**2)) == 1 and sp.denom(sp.cancel(wh / b**3)) == 1
    # At this one-parameter sample both local presentations acquire poles.
    # The patched numerator is polynomial and has the same residues over g=0.
    V, b = tau, 1 - tau
    h = b * a**3 + 2 * a**2 + 3 * a + 1
    N = a**7 + 3 * a + 1
    Lh, g = cubic_numerator(N, sp.Integer(1), V, h)
    Lq = conic_numerator(N, sp.Integer(1), V, h)
    assert sp.cancel(Lq).as_numer_denom()[1].subs(tau, 0) == 0
    assert sp.cancel(Lh).as_numer_denom()[1].subs(tau, 1) == 0
    wq, wh = bezout_patch(V, b, 1, 1, 3, 8)
    patched = sp.cancel(wq * Lq + wh * Lh)
    assert sp.denom(patched) == 1
    for old in [Lq, Lh]:
        residual = sp.cancel((patched - old) / g)
        assert sp.gcd(sp.denom(residual), g, tau) == 1


def run():
    U, V, b, c0, c1, c2 = sp.symbols("U V b c0 c1 c2", nonzero=True)
    h = b * (a**3 + c2 * a**2 + c1 * a + c0)
    q = U + V * a**2
    # In the power basis the adjugate's first column is the scalar inverse
    # numerator of q.  This verifies all three displayed coefficients.
    X = sp.Matrix([[0, 0, -c0], [1, 0, -c1], [0, 1, -c2]])
    Q = U * sp.eye(3) + V * X**2
    adj = Q.adjugate()
    A0, A1, A2 = adj[:, 0]
    assert sp.expand(A0 - (U**2 + U * V * (c2**2 - 2 * c1) + V**2 * (c1**2 - c0 * c2))) == 0
    assert sp.expand(A1 - V**2 * (c1 * c2 - c0)) == 0
    assert sp.expand(A2 - (V**2 * c1 - U * V)) == 0
    assert (adj - A0 * sp.eye(3) - A1 * X - A2 * X**2).applyfunc(sp.expand) == sp.zeros(3)
    _, g = cubic_numerator(1, U, V, h)
    assert sp.expand(g - b**2 * Q.det()) == 0
    assert sp.expand(g - sp.resultant(q, h, a)) == 0
    # Lagrange interpolation: [alpha^2](N/q mod (h/b))/b equals the sum
    # over h-roots.  At V=0 the representative is regular provided U*b != 0.
    n0, n1, n2 = sp.symbols("n0 n1 n2")
    N = n0 + n1 * a + n2 * a**2
    L, g = cubic_numerator(N, U, V, h)
    assert sp.cancel((L / g).subs(V, 0) - n2 / (b * U)) == 0

    # A moving simple point: h has root 2+tau. At tau=0 it meets q=0.
    tau = sp.Symbol("tau")
    h = 7 * (a - 1) * (a - 2 - tau) * (a - 3)
    U, V = tau - 4, sp.Integer(1)
    q = U + V * a**2
    N = a**10 + 2 * a + 3 + tau * a**4
    L, g = cubic_numerator(N, U, V, h)

    def at(z, root):
        return z.subs(a, root)

    roots = (1, 2 + tau, 3)
    hprime = sp.diff(h, a)
    direct = sum(at(N, r) / (at(q, r) * at(hprime, r)) for r in roots)
    assert sp.cancel(L / g - direct) == 0
    j = sp.diff(q, tau) * hprime - sp.diff(h, tau) * sp.diff(q, a)
    assert sp.cancel(sp.residue(L / g, tau, 0) - (N / j).subs({tau: 0, a: 2})) == 0
    # For the two conic roots the residues equal -L_q/g. The coefficient at
    # alpha^-1 at infinity is the sum of all finite residues (the polynomial
    # part has no such coefficient).  Use polynomial division to compute it.
    W = conic_numerator(N, U, V, h)
    denominator = sp.Poly(q * h, a)
    remainder = sp.rem(sp.Poly(N, a, domain="EX"), sp.Poly(q * h, a, domain="EX"))
    infinity_coeff = remainder.nth(denominator.degree() - 1) / denominator.LC()
    assert sp.cancel((L - W) / g - infinity_coeff) == 0
    # All three points meet q simultaneously when q is a multiple of tau:
    # the incidence resultant vanishes cubically, but the desired pole is simple.
    h0 = h.subs(tau, 0)
    L, g = cubic_numerator(N, 2 * tau, 3 * tau, h0)
    expected = sum(
        at(N.subs(tau, 0), r) / ((2 + 3 * r**2) * at(sp.diff(h0, a), r)) for r in (1, 2, 3)
    )
    assert sp.cancel(sp.residue(L / g, tau, 0) - expected) == 0
    check_patch()
    print("Cubic adjugate, resultant, three-point formula, V=0, moving Jacobian sign,")
    print("residue at infinity, common three-sheet fiber, and Bezout patch all verified.")


if __name__ == "__main__":
    run()
