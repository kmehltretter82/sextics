"""An exact first-order test for the regular base currents on the coupled core."""

import sympy as sp

r, s, t = sp.symbols("r s t")
kappa = sp.Rational(5, 7)
p = (
    1 + r**6 + kappa * r + 2 * kappa * r**5,
    kappa * r * s,
    2 * kappa * r * s,
    1 + s**6 + kappa * s + 3 * kappa * s**5,
)


# Unknowns (L,a,b,c) for q_t=ta, q_r=b, q_s=c and L=-(a+b_r+c_s).
def matrix():
    return sp.Matrix(
        [
            [
                pi,
                sp.Rational(4 * i, 3) * pi,
                sp.Rational(2, 3) * sp.diff(pi, r),
                sp.Rational(2, 3) * sp.diff(pi, s),
            ]
            for i, pi in enumerate(p)
        ]
    )


def evaluate(source, point):
    """Return q and the scalar obstruction to an exact base current at one good point."""
    sub = {r: sp.Integer(point[0]), s: sp.Integer(point[1])}
    source = sp.Matrix(source)
    M = matrix()

    def value(f):
        return f.subs(sub)

    inv = value(M).inv()
    x = inv * value(source)
    dr = inv * (value(source.diff(r)) - value(M.diff(r)) * x)
    ds = inv * (value(source.diff(s)) - value(M.diff(s)) * x)
    return x, x[0] + x[1] + dr[2] + ds[3]


def unpack(H):
    """Convert a homogeneous source in (a,b,c,d) to its ratio-chart density."""
    a, b, c, d = sp.symbols("a b c d")
    chart = sp.Poly(sp.expand(t * H.subs({a: 1, b: r, c: t, d: s * t}, simultaneous=True)), t)
    out = sp.Matrix([chart.nth(2 * i) for i in range(4)])
    assert sp.expand(sum(out[i] * t ** (2 * i) for i in range(4)) - chart.as_expr()) == 0
    return out


def polynomial_exact():
    a, b, c, d = sp.symbols("a b c d")
    x = [a, b, c, d]
    H = (
        a**6
        + b**6
        + c**6
        + d**6
        + kappa
        * (a**5 * b + 2 * a * b**5 + c**5 * d + 3 * c * d**5 + a * b * c * d * (a * a + 2 * c * c))
    )
    # Off-block linear fields are trace free, so their source is x_j dH/dx_i.
    return [unpack(x[j] * sp.diff(H, x[i])) for i in range(4) for j in range(4) if i // 2 != j // 2]


def monomial_sources():
    return [
        ([0] * (int((k + 1) / 2)) + [r**i * s**j] + [0] * (3 - int((k + 1) / 2)))
        for k in (1, 3, 5)
        for i in range(7 - k)
        for j in range(k + 1)
    ]


def get_sample_rows(sources, prime=1009):
    """Exact finite-field evaluation of J on polynomial sources and their derivatives."""
    X = sp.Matrix.hstack(*(sp.Matrix(u) for u in sources))
    M = matrix()
    exprs = [M, M.diff(r), M.diff(s), X, X.diff(r), X.diff(s)]

    # Clearing the fixed scalar denominators gives integer polynomial evaluators.
    def func(mat):
        mult = sp.ilcm(1, *[int(c.q) for term in mat for c in sp.Poly(term, r, s).coeffs()])
        ff = sp.lambdify((r, s), list(mult * mat), modules="math")
        return lambda R, S: sp.Matrix(
            mat.rows, mat.cols, [v * pow(mult, -1, prime) % prime for v in ff(R, S)]
        )

    fs = list(map(func, exprs))
    from random import Random

    rng = Random(324769)
    points, rows = [], []
    while len(rows) < 55:
        point = tuple(rng.randrange(1, prime) for _ in range(2))
        Mv, Mr, Ms, Xv, Xr, Xs = [fn(*point) for fn in fs]
        if Mv.det() % prime == 0:
            continue
        iv = Mv.inv_mod(prime)
        sol = (iv * Xv) % prime
        dr = iv * (Xr - Mr * sol)
        ds = iv * (Xs - Ms * sol)
        vals = (sol.row(0) + sol.row(1) + dr.row(2) + ds.row(3)) % prime
        rows.append(list(vals))
        points.append(point)
    return sp.Matrix(rows), points


def mod_rref(mat, p):
    mat = [[int(a) % p for a in row] for row in mat.tolist()]
    pivot_row = 0
    pivots = []
    for j in range(len(mat[0])):
        k = next((k for k in range(pivot_row, len(mat)) if mat[k][j]), None)
        if k is None:
            continue
        mat[pivot_row], mat[k] = mat[k], mat[pivot_row]
        mat[pivot_row] = [a * pow(mat[pivot_row][j], -1, p) % p for a in mat[pivot_row]]
        for k in range(len(mat)):
            if k != pivot_row and mat[k][j]:
                lead = mat[k][j]
                mat[k] = [(a - lead * b) % p for a, b in zip(mat[k], mat[pivot_row])]
        pivot_row += 1
        pivots.append(j)
    return pivots, mat


def check_single_component_fields():
    """Each single-component rational ODE has a polynomial unknown of degree at most four."""
    P = sum(coefficient * t ** (2 * i) for i, coefficient in enumerate(p))
    expected = (
        sp.Rational(14707736000, 194481),
        sp.Rational(22495616000, 194481),
        sp.Rational(-20432613056, 151263),
    )
    for variable, witness in zip((r, s, t), expected):
        specialize = {x: 1 for x in (r, s, t) if x != variable}
        fiber = sp.Poly(P.subs(specialize), variable)
        assert fiber.degree() == 6 and fiber.gcd(fiber.diff()).degree() == 0
        columns = [
            sp.Poly(
                -fiber.as_expr() * sp.diff(variable**j, variable)
                + sp.Rational(2, 3) * fiber.diff().as_expr() * variable**j,
                variable,
            )
            for j in range(5)
        ]
        columns.append(sp.Poly((t**4).subs(specialize), variable))
        minor = sp.Matrix([[col.nth(i) for col in columns] for i in range(6)]).det()
        assert minor == witness
        print("Single component", variable, "augmented determinant:", minor)


def main():
    check_single_component_fields()
    target = sp.Matrix([0, 0, 1, 0])
    M = matrix()
    print("M(1,1) determinant:", M.subs({r: 1, s: 1}).det())
    sol, obstruction = evaluate(target, (1, 1))
    print("target at (1,1):", list(sol), "obstruction", obstruction)
    b = r / (1 + s)
    c = s / (1 + r)
    a = (sp.diff(b, r) + sp.diff(c, s)) / 3 - sp.Rational(2, 9) * c * sp.diff(p[3], s) / p[3]
    L = -a - sp.diff(b, r) - sp.diff(c, s)
    z = M * sp.Matrix([L, a, b, c])
    assert sp.cancel(z[3]) == 0
    xs, J = evaluate(z, (1, 1))
    assert xs == sp.Matrix([e.subs({r: 1, s: 1}) for e in (L, a, b, c)])
    assert J == 0
    poly = polynomial_exact()
    mon = monomial_sources()
    rows, points = get_sample_rows(poly + [target] + mon)
    # Signs and order follow (derivative index, multiplier index):
    # H_c, r H_c, H_d, r H_d after passing to the ratio chart.
    for idx, (a0, c0) in enumerate(
        [
            (sp.Rational(3, 2), -3 * s / 2),
            (3 * r / 2, -3 * r * s / 2),
            (0, sp.Rational(3, 2)),
            (0, 3 * r / 2),
        ]
    ):
        expected = M * sp.Matrix([0, a0, 0, c0])
        assert (expected - poly[idx + 4]).applyfunc(sp.expand) == sp.zeros(4, 1)
        assert sp.simplify(a0 + sp.diff(c0, s)) == 0
    assert (
        sp.Matrix.hstack(
            *[
                sp.Matrix(
                    [
                        sp.Poly(b, r, s).coeff_monomial((i, j))
                        for b in u
                        for i in range(7)
                        for j in range(7)
                    ]
                )
                for u in poly
            ]
        ).rank()
        == 8
    )

    for width in (8, 9):
        cols, _ = mod_rref(rows[:, :width], 1009)
        assert cols == ([0, 1, 2, 3] if width == 8 else [0, 1, 2, 3, 8])
        print("E and target" if width == 9 else "E", len(cols), cols)
    mono = rows[:, 9:]
    cols, _ = mod_rref(mono, 1009)
    assert len(cols) == 36
    rr, _ = mod_rref(mono[:, cols].T, 1009)
    # Taking independent rows and columns gives a square invertible minor.
    minor = mono.extract(rr, cols)
    from sympy.polys.matrices import DomainMatrix
    from sympy import GF

    det = int(DomainMatrix.from_Matrix(minor).convert_to(GF(1009)).det()) % 1009
    print(
        "rank(W) witness mod 1009:",
        det,
        "on points",
        [points[i] for i in rr],
        "with omitted zero-based monomial indices",
        [i for i in range(40) if i not in cols],
    )
    assert det
    assert obstruction == sp.Rational(22379593, 11895240)


if __name__ == "__main__":
    main()
