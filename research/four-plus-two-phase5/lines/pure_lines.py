"""Certify isolation of 24 sign-invariant pure sixth-power lines in the full line Grassmannian."""

import sympy as sp
from sympy.polys.agca.extensions import FiniteExtension
from sympy.polys.rings import ring


def check(infinite_block, prime=11):
    s = sp.symbols("s")
    F = sp.GF(prime)
    k0 = 5 * pow(7, -1, prime) % prime
    binary = sp.Poly(
        1 + s**6 + k0 * s + (3 if infinite_block == "X" else 2) * k0 * s**5,
        s,
        modulus=prime,
    )
    K = FiniteExtension(binary)
    R, v, T = ring("v,T", K)
    ground = R.ground_new
    root = K.generator
    k = K.convert(k0)
    D = K.convert(binary.diff().as_expr())
    C = (2 if infinite_block == "X" else 1) * k * root
    assert binary.gcd(binary.diff()).degree() == 0
    assert D * D.inverse() == C * C.inverse() == K.one
    sigma = ground(root)
    u = w = z = R.zero
    Q = R.one

    def defect():
        # X=1,Y=T. The first block has its unit coordinate at infinity; the
        # second ratio is a root of the squarefree complementary binary sextic.
        P0, P1 = u + v * T, R.one
        S0, S1 = T, w + (sigma + z) * T
        a, b, c, d = (P0, P1, S0, S1) if infinite_block == "X" else (S0, S1, P0, P1)
        H = (
            a**6
            + b**6
            + c**6
            + d**6
            + ground(k)
            * (
                a**5 * b
                + a**3 * b * c * d
                + 2 * a * b**5
                + 2 * a * b * c**3 * d
                + c**5 * d
                + 3 * c * d**5
            )
        )
        return H - Q**3

    # These six rows have invertible linear term in u,w,z and the three
    # coefficients of the quadratic cube root. Thus they solve uniquely
    # as series in v; the remaining row (T^3) tests the single formal branch.
    for m in range(1, 4):
        weight = v**m
        err = defect()
        u -= ground(err.get((m, 4), K.zero) / C) * weight
        w -= ground(err.get((m, 5), K.zero) / D) * weight
        z -= ground(err.get((m, 6), K.zero) / D) * weight
        err = defect()
        for j in range(3):
            Q += ground(err.get((m, j), K.zero) / K.convert(3)) * weight * T**j
        err = defect()
        assert all(err.get((m, j), K.zero) == K.zero for j in (0, 1, 2, 4, 5, 6))
        obstruction = err.get((m, 3), K.zero)
        if m < 3:
            assert obstruction == K.zero
    obstruction = sp.Poly(obstruction.as_expr(), s, modulus=prime)
    gcd = obstruction.gcd(binary).monic()
    print(
        infinite_block,
        "infinity; binary",
        binary.as_expr(),
        "cubic obstruction",
        obstruction.as_expr(),
        "gcd",
        gcd.as_expr(),
    )
    assert gcd.degree() == 0


if __name__ == "__main__":
    check("X")
    check("Y")
