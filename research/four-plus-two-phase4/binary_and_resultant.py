#!/usr/bin/env python3
"""Check the binary normalized trace and the point-resultant residue weights.

Run with Python 3 and SymPy >= 1.12. All symbols represent characteristic-zero
indeterminates. The last check works in the localization v o != 0; it also
gives a residue numerator which is defined wherever v is invertible.
"""

import sympy as sp


def zero(expr):
    assert sp.cancel(expr) == 0, expr


def check_binary():
    C, R, D, S, Cr, Rr, Ds, Ss, a, t = sp.symbols("C R D S C_r R_r D_s S_s alpha t")
    k = sp.I * R / S
    f, g = C**3 + R**2, D**3 + S**2
    fr, gs = 3 * C**2 * Cr + 2 * R * Rr, 3 * D**2 * Ds + 2 * S * Ss
    jx, jy = 2 * C * Rr - 3 * R * Cr, 2 * D * Ss - 3 * S * Ds
    scale = -4 * sp.I * jx * jy / 9

    # In K[alpha]/(alpha^3-k) the denominators used below are rational in K.
    def reduce_a(expr):
        n, d = sp.fraction(sp.cancel(expr))
        d = sp.rem(d, a**3 - k, a)
        assert not d.has(a) and d != 0
        return sp.factor(sp.rem(n, a**3 - k, a) / d)

    # b^2 eta_r, b^2 eta_s for the closed form d(alpha^2 D/C)/(b/C)^2.
    ar, ass = a * Rr / (3 * R), -a * Ss / (3 * S)
    U, V = a**2 * D * jx / (3 * R), -(a**2) * C * jy / (3 * S)
    Z = U * ass - V * ar
    T = reduce_a((6 * g * a**5 * Z + U * gs * a**6 - V * fr) / (f + g * a**6))
    top = reduce_a(sp.Rational(4, 3) * (a * g * T - a * U * gs - 6 * g * Z))
    zero(top - scale)
    # The coefficients before trace include alpha^4 and alpha^2.
    coeffs = [reduce_a(top * a**m) for m in (4, 2, 0)]

    # Compute the trace as the diagonal sum of the *regular representation*.
    def trace(expr):
        diagonal = sum(sp.expand(reduce_a(expr * a**j)).coeff(a, j) for j in range(3))
        return sp.factor(diagonal / 3)

    for observed, expected in zip([trace(k0) for k0 in coeffs], (0, 0, scale)):
        zero(observed - expected)

    # Trace the three flux components in the same cubic algebra *before*
    # applying the differential operator; this checks its sign independently.
    # (t^2-alpha^2)^(-1) = (t^4+t^2 alpha^2+alpha^4)/(t^6-k^2).
    inv = (t**4 + t**2 * a**2 + a**4) / (t**6 - k**2)
    zero(reduce_a(inv * (t**2 - a**2)) - 1)
    flux = [trace(expr * inv) for expr in (2 * a * V, -2 * a * U, -2 * t * Z)]

    def dr(expr):
        return sp.diff(expr, C) * Cr + sp.diff(expr, R) * Rr

    def ds(expr):
        return sp.diff(expr, D) * Ds + sp.diff(expr, S) * Ss

    P = f + g * t**6
    defect = P * (dr(flux[0]) + ds(flux[1]) + sp.diff(flux[2], t))
    defect -= sp.Rational(2, 3) * (fr * flux[0] + gs * t**6 * flux[1] + 6 * g * t**5 * flux[2])
    zero(defect + scale * t**4)
    return "(0, 0, -4*i*(2*C*R_r-3*R*C_r)*(2*D*S_s-3*S*D_s)/9)"


def check_point_resultant():
    # The coefficients have arbitrary first derivatives in each base/curve
    # direction. We substitute point relations only after differentiating.
    u, v, e0, e2, o0, o2, a = sp.symbols("u v e0 e2 o0 o2 alpha")
    variables = (u, v, e0, e2, o0, o2)
    derivatives = sp.symbols("u_x v_x e0_x e2_x o0_x o2_x")

    def dx(p):
        return sum(sp.diff(p, c) * dc for c, dc in zip(variables, derivatives))

    q = u + v * a**2
    h = e0 + e2 * a**2 + a * (o0 + o2 * a**2)
    e, o = v * e0 - u * e2, v * o0 - u * o2
    F = v * e**2 + u * o**2
    zero(sp.resultant(q, h, a) - F)

    def at_root(p):
        return p.subs({u: -v * a**2, e0: -e2 * a**2 - a * (o0 + o2 * a**2)})

    jac = dx(q) * sp.diff(h, a) - dx(h) * sp.diff(q, a)
    zero(at_root(dx(F) - v * o * jac))

    # For any N polynomial (remainder n0+n1*a), the residue of
    # eta_R * v*(o*n0-e*n1)/F at a simple zero equals ell_R*N/j.
    n0, n1 = sp.symbols("n0 n1")
    zero(at_root(v * (o * n0 - e * n1) - v * o * (n0 + a * n1)))
    # A second derivation: sum residues at the two q-roots of N/(q*h).
    # It equals (1/v)*coeff_alpha (N/h mod q), with the opposite sign.
    num = sp.rem((n0 + n1 * a) * v**2 * (e - o * a), a**2 + u / v, a)
    zero(sp.expand(num).coeff(a) / v - v * (e * n1 - o * n0))
    zero(sp.rem(v**2 * (e - o * a) * h - F, a**2 + u / v, a))
    return "F=Res(q,h); dF=v*o*(dq*h_alpha-dh*q_alpha); residue weight=v*(o*n0-e*n1)"


def check_seed():
    # Restrict to E_b at r=0. Q2=ab, so q=r and h=1+alpha^3 B.
    # The tangent of the parameter curve on the normalized branch need not
    # coincide with the signed-minor tangent at the singular central point.
    a, B, C, Br, Cr, u, v = sp.symbols("alpha B C B_s C_s u v")
    k = -1 / B
    beta = -C * a**2
    Er, EZ = -3 * a**2 * B, a**3 * Br
    EY = Er * (-Cr * a**2) + sp.diff(beta, a) * EZ
    j = (u + v * a**2) * 3 * a**2 * B
    g = C**3 - B**2
    wronskian = 2 * C * Br - 3 * B * Cr

    def rem(expr):
        return sp.rem(sp.expand(expr), a**3 - k, a)

    # p1=p2=0, as are their s-partials. The r-partials do not
    # contribute because Es=0. Thus only the S3 term remains.
    S3 = (
        a * (3 * g * beta * EY - beta**2 * Er * (3 * C**2 * Cr - 2 * B * Br)) - 6 * g * beta**2 * EZ
    )
    zero(rem((u + v * a**2) * S3 - j * a**2 * C * wronskian))
    X = sp.Matrix([[0, 0, k * k], [1, 0, 0], [0, 1, 0]])  # multiplication by alpha^2
    T = (u * sp.eye(3) + v * X).inv()
    scalar = B**2 * u**3 + v**3
    for power, expected in zip((3, 2, 1), (u**2, -u * v, v**2)):
        zero((X**power * T).trace() * scalar / 3 - expected)
    return "(4/3)*ell*C*(2*C*B_s-3*B*C_s)/(B^2*u^3+v^3)*(u^2,-u*v,v^2)"


if __name__ == "__main__":
    print("Binary trace:", check_binary())
    print("Coupled point pushforward:", check_point_resultant())
    print("Coupled three-sheet seed trace at r=0:", check_seed())
