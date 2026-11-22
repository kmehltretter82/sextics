"""Single-border determinant for normalized current traces; exact binary check."""

import sympy as sp


def coefficient_row(poly, variable, degree):
    coeff = sp.Poly(poly, variable)
    return [coeff.nth(i) for i in range(degree + 1)]


def trace_border(relation, variable, denominator, numerators):
    """Macaulay border for the finite complete univariate algebra k[x]/(relation)."""
    n = sp.degree(relation, variable)
    degree = (
        n - 1 + max(sp.degree(poly, variable) for poly in [denominator, *numerators])
    )
    parameters = sp.symbols(f"X:{len(numerators)}")
    test = denominator + sum(x * poly for x, poly in zip(parameters, numerators))
    relations = [relation * variable**i for i in range(degree - n + 1)]
    border = [test * variable**i for i in range(n)]
    matrix = sp.Matrix(
        [coefficient_row(poly, variable, degree) for poly in relations + border]
    )
    zero = dict.fromkeys(parameters, 0)
    constant = matrix.subs(zero)
    assert constant.det() != 0
    # Jacobi's determinant identity computes all logarithmic derivatives at once.
    inverse = constant.inv()
    return [sp.factor((inverse * matrix.diff(x)).trace() / n) for x in parameters]


def main():
    alpha, rho, scale, z = sp.symbols("alpha rho scale z")
    # rho=i R/S and scale=−4 i J_X J_Y/9 in the geometric binary family.
    # The arbitrary nonconstant denominator checks that it cancels in the norm.
    denominator = 1 - z * alpha
    source = [scale * alpha**4, scale * alpha**2, scale]
    result = trace_border(
        alpha**3 - rho, alpha, denominator, [denominator * h for h in source]
    )
    assert result == [0, 0, scale]
    print("Full binary normalized source:", result)
    # An independent coupled cubic chart from the actual nodal kernel: this
    # checks all three alpha roots at one fixed parameter, not the ratio fiber.
    B, C, dC, dB, g, dg, U, V = sp.symbols("B C dC dB g dg U V")
    beta = -C * alpha**2
    ha = 3 * B * alpha**2
    Er, EZ = -ha, dB * alpha**3
    EY = Er * (-dC * alpha**2) - 2 * C * alpha * EZ
    S3 = alpha * (3 * g * beta * EY - beta**2 * Er * dg) - 6 * g * beta**2 * EZ
    # Actual nodal seed at (r,s)=(0,1), reduced mod 23; v²=15.
    v = sp.symbols("v")
    vals = {B: 3, C: 3, dC: 4, dB: 14, g: 18, dg: 1, U: 8 + 21 * v, V: 9}
    f = (B * alpha**3 + 1).subs(vals)
    jac = (ha * (U + V * alpha**2)).subs(vals)
    nums = [sp.Rational(4, 3) * S3.subs(vals) * alpha**j for j in (4, 2, 0)]

    def reduced(poly):
        return sp.rem(poly, f, alpha)

    out = trace_border(f, alpha, reduced(jac), [reduced(poly) for poly in nums])
    companion = sp.Matrix([[0, 0, -sp.Rational(1, 3)], [1, 0, 0], [0, 1, 0]])

    def evaluate(poly):
        return sum(
            (
                coef * companion ** degree[0]
                for degree, coef in sp.Poly(poly, alpha).terms()
            ),
            start=sp.zeros(3),
        )

    expected = [
        sp.factor((evaluate(jac).inv() * evaluate(v)).trace() / 3) for v in nums
    ]
    assert out == expected

    def field_pair(expr):
        num, den = sp.fraction(sp.cancel(expr))
        relation = sp.Poly(v * v - 15, v, modulus=23)

        def convert(poly):
            q = sp.Poly(poly, v)
            return sp.Poly.from_dict(
                {mon: int(c.p) * pow(int(c.q), -1, 23) % 23 for mon, c in q.terms()},
                v,
                modulus=23,
            )

        value = (convert(num) * convert(den).invert(relation)).rem(relation)
        return tuple(int(value.nth(j)) % 23 for j in range(2))

    pairs = [field_pair(entry) for entry in out]
    print("Actual coupled three-sheet control for dτ at (0,1), F23[v]/(v²−15):", pairs)
    assert pairs == [(5, 21), (22, 2), (21, 3)]


if __name__ == "__main__":
    main()
