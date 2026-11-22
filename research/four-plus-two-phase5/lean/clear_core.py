"""Clear a Q(r,s,t) current to homogeneous core polynomials and verify the identity.

Import clear_current(qr, qs, qt) and pass SymPy expressions in r,s,t.
The default core is the protected four-variable sextic. Run as a script for
an exact nontrivial test against a degenerate core admitting a known current.
"""

from functools import reduce

import sympy as sp

a, b, c, d, r, s, t = sp.symbols("a b c d r s t")
x = (a, b, c, d)
k = sp.Rational(5, 7)
core = (
    a**6
    + k * a**5 * b
    + k * a**3 * b * c * d
    + 2 * k * a * b**5
    + 2 * k * a * b * c**3 * d
    + b**6
    + c**6
    + k * c**5 * d
    + 3 * k * c * d**5
    + d**6
)


def clear_current(qr, qs, qt, H=core, domain=sp.QQ):
    """Return (denominator, numerators) after checking both chart and core."""
    q = tuple(map(sp.sympify, (qr, qs, qt)))
    chart = (r, s, t)
    P = H.subs(dict(zip(x, (1, r, t, s * t))), simultaneous=True)
    chart_defect = (
        sp.Rational(2, 3) * sum(v * P.diff(y) for v, y in zip(q, chart))
        - P * sum(v.diff(y) for v, y in zip(q, chart))
        - t**4
    )

    def is_zero(expr, variables):
        num = sp.together(expr).as_numer_denom()[0]
        return sp.Poly(num, *variables, domain=domain).is_zero

    if not is_zero(chart_defect, chart):
        raise ValueError("The chart derivative does not equal the prescribed t^4.")

    # With p_a=0 the b,c,d components are the Piola pushforward.
    qr, qs, qt = q
    inverses = {r: b / a, s: d / c, t: c / a}
    densities = (0, qr, qt, s * qt + t * qs)
    physical = [
        sp.cancel((2 * a * e / (3 * t)).subs(inverses, simultaneous=True))
        for e in densities
    ]
    # A common polynomial denominator may have degree zero, while the
    # CoreCertificate interface asks for a positive degree.  Multiply by a.
    denominators = [sp.Poly(sp.denom(z), *x, domain=domain) for z in physical]
    D = reduce(sp.Poly.lcm, denominators).as_expr()
    if sp.Poly(D, *x).total_degree() == 0:
        D *= a
    N = tuple(sp.cancel(D * z) for z in physical)
    Dp = sp.Poly(D, *x, domain=domain)
    Np = [sp.Poly(z, *x, domain=domain) for z in N]
    if not Dp.is_homogeneous or not all(
        p.is_zero or (p.is_homogeneous and p.total_degree() == Dp.total_degree() + 1)
        for p in Np
    ):
        raise ArithmeticError("Unexpected degrees in the physical lift")
    contraction = sum(h * H.diff(y) for h, y in zip(N, x))
    divergence = D * sum(h.diff(y) for h, y in zip(N, x)) - sum(
        h * D.diff(y) for h, y in zip(N, x)
    )
    if not is_zero(
        2 * a**3 * c**3 * D**2 - 2 * D * contraction + 3 * H * divergence, x
    ):
        raise ArithmeticError("The cleared polynomial certificate identity failed")
    return D, N


if __name__ == "__main__":
    # At P=t^6 the current (0,0,1/(5t)) has derivative t^4.
    print(clear_current(0, 0, 1 / (5 * t), H=c**6))
    quadratic_field = sp.QQ.algebraic_field(sp.sqrt(2))
    print(clear_current(sp.sqrt(2), 0, 1 / (5 * t), H=c**6, domain=quadratic_field))
    try:
        clear_current(0, 0, 1 / (5 * t))
    except ValueError as e:
        print("Negative control on the exact protected core:", e)
    else:
        raise AssertionError(
            "The degenerate test current cannot solve the coupled core"
        )
