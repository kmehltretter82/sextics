"""Solve bounded instances of FixedSupportEquation over Q, with exact checks."""

import argparse
import ast
import json
from itertools import combinations_with_replacement
from pathlib import Path

import sympy as sp
from sympy.polys.matrices import DomainMatrix

from clear_core import a, b, c, core, d


coords = x, y, z = sp.symbols("x y z")


def poly(expr):
    return sp.Poly(expr, *coords, domain=sp.QQ)


h = poly(core.xreplace({a: 1, b: x, c: y, d: z}))


def expression(text):
    """Parse integer arithmetic in x,y,z,h without executing user input."""
    names = dict(zip(("x", "y", "z", "h"), (*coords, h.as_expr())))

    def visit(node):
        if isinstance(node, ast.Constant) and type(node.value) is int:
            return sp.Integer(node.value)
        if isinstance(node, ast.Name) and node.id in names:
            return names[node.id]
        if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.USub, ast.UAdd)):
            value = visit(node.operand)
            return -value if isinstance(node.op, ast.USub) else value
        ops = {
            ast.Add: sp.Add,
            ast.Sub: lambda u, v: u - v,
            ast.Mult: sp.Mul,
            ast.Div: lambda u, v: u / v,
            ast.Pow: sp.Pow,
        }
        if isinstance(node, ast.BinOp) and type(node.op) in ops:
            return ops[type(node.op)](visit(node.left), visit(node.right))
        raise ValueError("Use rational arithmetic in x,y,z,h")

    return visit(ast.parse(text.replace("^", "**"), mode="eval").body)


def basis(degree):
    for total in range(degree + 1):
        for indices in combinations_with_replacement(range(3), total):
            exp = tuple(indices.count(i) for i in range(3))
            yield exp, sp.prod(v**e for v, e in zip(coords, exp))


def cleared(B, n, N):
    """Direct polynomial formula, independent of the sparse column construction."""
    action = sum((v * h.diff(i) for i, v in enumerate(N)), poly(0))
    div = sum((v.diff(i) for i, v in enumerate(N)), poly(0))
    normal = sum((v * B.diff(i) for i, v in enumerate(N)), poly(0))
    return 2 * B * action - 3 * h * (B * div - n * normal)


def verify(result):
    """Replay either saved witness using only polynomial arithmetic on the literal core."""
    B, n, degree = poly(expression(result["B"])), result["n"], result["degree"]
    if B.is_zero or min(n, degree) < 0:
        raise ValueError("Require a nonzero B and nonnegative exponent and degree")
    rhs = poly(sp.cancel(3 * B.as_expr() ** (n + 1) * expression(result["source"])))
    if "numerators" in result:
        N = [poly(expression(v)) for v in result["numerators"]]
        if (
            len(N) != 3
            or max(v.total_degree() for v in N) > degree
            or cleared(B, n, N) != rhs
        ):
            raise ArithmeticError("Polynomial re-expansion failed")
        return
    weight = {
        tuple(powers): sp.Rational(value) for *powers, value in result["obstruction"]
    }

    def pair(p):
        terms = p.as_dict()
        return sum(value * terms.get(exp, 0) for exp, value in weight.items())

    for i in range(3):
        for _, u in basis(degree):
            N = [poly(0)] * 3
            N[i] = poly(u)
            if pair(cleared(B, n, N)) != 0:
                raise ArithmeticError("Dual functional failed on a column")
    if pair(rhs) != 1:
        raise ArithmeticError("Dual functional failed on the source")


def solve(denominator, n, degree, source=y**3, dual=False):
    B = poly(denominator)
    if B.is_zero or min(n, degree) < 0:
        raise ValueError("Require a nonzero B and nonnegative exponent and degree")
    source = sp.cancel(source)
    rhs = poly(sp.cancel(3 * B.as_expr() ** (n + 1) * source))
    monomials = list(basis(degree))
    width = 3 * len(monomials)
    rows = {}

    def add_column(j, coefficients, shift=(0, 0, 0), scale=1):
        for powers, value in coefficients:
            exp = tuple(u + v for u, v in zip(powers, shift))
            row = rows.setdefault(exp, {})
            new = row.get(j, 0) + scale * value
            if new:
                row[j] = new
            else:
                row.pop(j, None)

    div_factor = (-3 * h * B).as_dict(native=True).items()
    for i in range(3):
        action = (
            (2 * B * h.diff(i) + 3 * n * h * B.diff(i)).as_dict(native=True).items()
        )
        for k, (exp, _) in enumerate(monomials):
            j = i * len(monomials) + k
            add_column(j, action, exp)
            if exp[i]:
                shift = list(exp)
                shift[i] -= 1
                add_column(j, div_factor, shift, exp[i])
    add_column(width, rhs.as_dict(native=True).items())
    exps = sorted((e for e, row in rows.items() if row), key=lambda e: (sum(e), e))
    data = {i: rows[e] for i, e in enumerate(exps)}
    matrix = DomainMatrix.from_dod(data, (len(exps), width + 1), sp.QQ)
    reduced, pivots = matrix.rref()
    inconsistent = width in pivots
    result = dict(
        B=str(B.as_expr()),
        n=n,
        degree=degree,
        source=str(source),
        rows=len(exps),
        unknowns=width,
        entries=sum(len(row) for row in data.values()),
        rank=len(pivots) - int(inconsistent),
        augmented_rank=len(pivots),
    )
    if inconsistent:
        if dual:
            # Find a coefficient functional lambda with lambda(A)=0, lambda(rhs)=1.
            target = DomainMatrix.from_dod(
                {width: {0: sp.QQ.one}}, (width + 1, 1), sp.QQ
            )
            rr, pp = matrix.transpose().hstack(target).rref()
            last = len(exps)
            if last in pp:
                raise ArithmeticError("Inconsistent dual system")
            dual_rows = rr.to_dod()
            weight = {
                exps[j]: sp.QQ.to_sympy(dual_rows[i][last])
                for i, j in enumerate(pp)
                if last in dual_rows[i]
            }

            result["obstruction"] = [
                [*exp, str(v)] for exp, v in sorted(weight.items())
            ]
            verify(result)
        return result

    rows_rref = reduced.to_dod()
    values = [sp.S.Zero] * width
    for i, j in enumerate(pivots):
        values[j] = sp.QQ.to_sympy(rows_rref[i].get(width, sp.QQ.zero))
    N = [
        poly(
            sum(
                values[i * len(monomials) + j] * u for j, (_, u) in enumerate(monomials)
            )
        )
        for i in range(3)
    ]
    result["numerators"] = [str(v.as_expr()) for v in N]
    verify(result)
    return result


def control():
    B = 1 + x + y
    N = (x * y + z, 1 + z**2, x - y * z)
    q = [v / B**2 for v in N]
    S = sum(
        sp.Rational(2, 3) * v * h.as_expr().diff(u) - h.as_expr() * v.diff(u)
        for v, u in zip(q, coords)
    )
    result = solve(B, 2, 2, S)
    if (
        result.get("numerators") != [str(poly(v).as_expr()) for v in N]
        or result["rank"] != 30
    ):
        raise ArithmeticError("Control did not recover the unique known flux")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--B", default="h")
    parser.add_argument("--n", type=int, default=1)
    parser.add_argument("--degree", type=int, default=3)
    parser.add_argument("--source", default="y^3")
    parser.add_argument(
        "--dual",
        action="store_true",
        help="verify and output an inconsistency certificate",
    )
    parser.add_argument(
        "--control",
        action="store_true",
        help="recover a manufactured coupled-core flux",
    )
    parser.add_argument(
        "--verify",
        type=Path,
        help="re-expand a saved witness or obstruction without elimination",
    )
    args = parser.parse_args()
    if args.verify:
        verify(json.loads(args.verify.read_text()))
        print("Exact polynomial verification passed")
    else:
        result = (
            control()
            if args.control
            else solve(
                expression(args.B),
                args.n,
                args.degree,
                expression(args.source),
                args.dual,
            )
        )
        print(json.dumps(result, indent=2))
