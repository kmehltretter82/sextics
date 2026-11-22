"""Canny support projection on a genuine coupled contact-curve fiber with an excess line."""

import sympy as sp


def main():
    t, d, e = sp.symbols("t d e")
    r = 4
    q = 4 + r + 4 * r * r + 5 * t * t + 10 * t * d + 7 * d * d
    h = (
        5
        + 3 * r**3
        + t * t
        + 8 * t * d
        + 9 * d * d
        + 5 * r * t * t
        + 6 * r * t * d
        + 4 * r * d * d
        + 7 * t
        + 9 * d
        + 5 * r * r * t
        + 2 * r * r * d
        + 4 * t**3
        + 5 * t * t * d
        + 6 * d**3
    )
    core = (
        1
        + r**6
        + t**6
        + d**6
        + 7 * (r + 2 * r**5 + t**5 * d + 3 * t * d**5 + r * t * d * (1 + 2 * t * t))
    )
    q1 = 4 + 8 * r * r + 4 * t * t + 3 * t * d
    sig = h.subs({t: -t, d: -d}, simultaneous=True)
    assert (
        sp.groebner([q], d, t, modulus=11).reduce(sp.expand(core - q1**3 - h * sig))[1]
        == 0
    )
    # Fix a=1,b=4. Multiplying both equations by d creates the whole line d=0.
    # The true point (c,d)=(6,1) is isolated and persists after this modification.
    ev = lambda f: int(f.subs({t: 6, d: 1})) % 11
    assert ev(q) == ev(h) == 0
    jac = ev(q.diff(t) * h.diff(d) - q.diff(d) * h.diff(t))
    assert jac
    f = sp.Poly(sp.resultant(q, h, d), t, modulus=11)
    assert sp.resultant(d * q, d * h, d) == 0
    res = sp.Poly(sp.resultant(d * q + e * t**3, d * h + e * d**4, d), e, t, modulus=11)
    power = min(exp[0] for exp, c in res.terms())
    low = sp.Poly.from_dict(
        {(j,): c for (i, j), c in res.terms() if i == power}, (t,), modulus=11
    )
    fac = low.exquo(f)
    assert power == 1 and f.eval(6) == 0 and low.eval(6) == 0
    print("coupled norm checked; contact Jacobian at (6,1):", jac)
    print("direct sextic in t:", f.as_expr())
    print("excess ideal has identically zero ordinary resultant")
    print("first resultant perturbation order:", power, "excess factor:", fac.as_expr())


if __name__ == "__main__":
    main()
