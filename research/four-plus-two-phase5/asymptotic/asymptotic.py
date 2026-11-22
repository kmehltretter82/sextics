"""Verify contact order, dominance and the target's base coefficient on the literal sextic."""

from itertools import product
import sympy as sp


x, y, z, q, r, e, t = sp.symbols("x y z q r e t")
coords = (x, y, z)
variables = (x, y, z, q, r, e)
direction = (1, q, r)
# Use 7h to keep the evaluators integral. The actual perturbation is h+e*y^3.
H = (
    7
    + 5 * x
    + 5 * x * y * z
    + 10 * x**5
    + 10 * x * y**3 * z
    + 7 * x**6
    + 7 * y**6
    + 5 * y**5 * z
    + 15 * y * z**5
    + 7 * z**6
    + 7 * e * y**3
)
along = sp.Poly(
    H.subs(
        dict(zip(coords, [v + t * d for v, d in zip(coords, direction)])),
        simultaneous=True,
    ),
    t,
)
contact = [along.nth(i) for i in range(3)]
cubic = [along.nth(i) for i in range(3, 7)]


def evaluate(expr):
    return sp.lambdify(variables, expr, cse=True)


value = evaluate(contact)
jac = evaluate([[f.diff(v) for v in variables] for f in contact])
hess = evaluate(
    [[[f.diff(v, w) for w in variables] for v in variables] for f in contact]
)
cval = evaluate(cubic)
cjac = evaluate([[f.diff(v) for v in variables] for f in cubic])
abc = sp.symbols("c0:4")
disc_expr = sp.discriminant(sum(a * t**i for i, a in enumerate(abc)), t)
disc = sp.lambdify(abc, disc_expr)
disc_jac = sp.lambdify(abc, [disc_expr.diff(c) for c in abc])


def determinant(A):
    (a, b, c), (d, e, f), (g, h, i) = A
    return a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)


def run(p=11):
    def inv(k):
        return pow(int(k) % p, -1, p)

    def solve3(A, v):
        di = inv(determinant(A))
        return [
            determinant(
                [
                    [v[k] if i == j else a for i, a in enumerate(row)]
                    for k, row in enumerate(A)
                ]
            )
            * di
            % p
            for j in range(3)
        ]

    def dot(a, b):
        return sum(u * v for u, v in zip(a, b)) % p

    for x0, y0, z0 in product(range(p), repeat=3):
        if value(x0, y0, z0, 0, 0, 0)[0] % p:
            continue
        for q0, r0 in product(range(p), repeat=2):
            point = x0, y0, z0, q0, r0, 0
            if any(i % p for i in value(*point)):
                continue
            D = jac(*point)
            dep = (2, 3, 4)  # solve for z,q,r over x,y,e
            A = [[row[i] % p for i in dep] for row in D]
            if determinant(A) % p == 0:
                continue

            def tangent(j):
                v = [int(k == j) for k in range(6)]
                v[2:5] = solve3(A, [-row[j] for row in D])
                return v

            X, Y, E = map(tangent, (0, 1, 5))
            dd = hess(*point)

            def mixed(v):
                return solve3(A, [-dot(v, [dot(row, E) for row in ha]) for ha in dd])

            X1, Y1 = mixed(X), mixed(Y)
            J = (X[3] * Y[4] - Y[3] * X[4]) % p
            Jdot = (X1[1] * Y[4] + X[3] * Y1[2] - Y1[1] * X[4] - Y[3] * X1[2]) % p
            C = [v * inv(7) % p for v in cval(*point)]
            Cdot = [dot(v, E) * inv(7) % p for v in cjac(*point)]
            delta = disc(*C) % p
            ddelta = dot(disc_jac(*C), Cdot)
            if not all((J, C[0], C[3], delta)):
                continue
            beta = 3 * inv(2) * (Jdot - J * ddelta * inv(6 * delta)) % p
            if beta == 0:
                continue
            # Expand the physical three-dimensional volume along the congruence.
            dxp, dyp = X[:3], Y[:3]
            dxv, dyv = [0, *X[3:5]], [0, *Y[3:5]]
            v = [1, q0, r0]
            assert determinant([dxp, dyp, v]) % p == 0
            assert (determinant([dxv, dyp, v]) + determinant([dxp, dyv, v])) % p == 0
            assert determinant([dxv, dyv, v]) % p == J
            print(
                "prime", p, "point", point, "det(contact Jacobian)", determinant(A) % p
            )
            print("tangents x,y,e", X, Y, E, "mixed x/e,y/e (z,q,r)", X1, Y1)
            print("residual cubic constant-first", C, "derivative", Cdot)
            print(
                "J2, J2dot, Delta, Deltadot, beta/gamma", J, Jdot, delta, ddelta, beta
            )
            return
    raise ArithmeticError("No nondegenerate sample")


if __name__ == "__main__":
    run()
