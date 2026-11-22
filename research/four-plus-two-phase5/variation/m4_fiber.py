"""Emit the 39-variable sparse M4GB input from the public norm circuit, mod 23."""

import argparse
from collections import Counter
import json
from pathlib import Path
import sys
import sympy as sp


def reduce_mod(expr):
    variables = sorted(expr.free_symbols, key=str)
    if not variables:
        return sp.Integer(int(expr.p) * pow(int(expr.q), -1, 23) % 23)
    terms = sp.Poly(expr, *variables).terms()
    return sum(
        (int(c.p) * pow(int(c.q), -1, 23) % 23)
        * sp.prod(x**j for x, j in zip(variables, m))
        for m, c in terms
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("circuit", type=Path)
    parser.add_argument("--r", type=int, default=4)
    parser.add_argument("--s", type=int, default=6)
    parser.add_argument(
        "--names", type=Path, help="Optional variable and reverse-substitution map"
    )
    opts = parser.parse_args()
    data = json.loads(opts.circuit.read_text())
    w_names = [f"W{i}" for i in range(len(data["quartic_multiplier_basis"]))]
    symbols = {
        key: sp.Symbol(key)
        for key in data["parameter_order"] + ["alpha"] + w_names + ["r", "s"]
    }

    def parse(x):
        return sp.sympify(x, locals=symbols)

    equations = []
    for block in ("pivot", "free"):
        for value, row in zip(data["R_" + block], data["M_" + block]):
            equation = parse(value) - sum(
                symbols[name] * parse(c) for name, c in zip(w_names, row) if c != "0"
            )
            equations.append(reduce_mod(equation))
    alpha = symbols["alpha"]
    a2, a3, z, ia, i1, i2, i3 = sp.symbols("a2 a3 z ia i1 i2 i3")
    substitutions = {symbols["r"]: opts.r, symbols["s"]: opts.s}
    for key in ("point_q", "point_h"):
        poly = sp.Poly(reduce_mod(parse(data[key]).subs(substitutions)), alpha)
        equations.append(
            sum(poly.nth(j) * factor for j, factor in enumerate([1, alpha, a2, a3]))
        )
    equations.extend(
        [
            a2 - alpha**2,
            a3 - alpha * a2,
            z - symbols["mu0"] * symbols["mu1"],
            ia * alpha - 1,
            i1 * (1 - z) - 1,
            i2 * (1 - 2 * z) - 1,
            i3 * (1 - 3 * z) - 1,
        ]
    )
    variables = [symbols[k] for k in data["parameter_order"] + ["alpha"] + w_names] + [
        a2,
        a3,
        z,
        ia,
        i1,
        i2,
        i3,
    ]
    removed = []
    while True:
        choices = []
        for i, equation in enumerate(equations):
            for x in sorted(equation.free_symbols, key=str):
                if str(x) not in w_names:
                    continue
                coefficient = equation.diff(x)
                if coefficient.is_Integer and int(coefficient) % 23:
                    cost = sp.count_ops(equation) * sum(e.has(x) for e in equations)
                    choices.append((cost, str(x), i, x, int(coefficient) % 23))
        if not choices:
            break
        _, _, index, x, coefficient = min(choices)
        equation = equations.pop(index)
        variables.remove(x)
        value = sp.expand(-equation.subs(x, 0) * pow(coefficient, -1, 23))
        removed.append([str(x), str(value)])
        for i, equation in enumerate(equations):
            if equation.has(x):
                updated = sp.expand(equation.subs(x, value))
                equations[i] = sp.Poly(
                    updated, *sorted(updated.free_symbols, key=str), modulus=23
                ).as_expr()
    equations = [sp.expand(e) for e in equations]
    names = {x: sp.Symbol(f"X{i}") for i, x in enumerate(variables)}
    if opts.names:
        opts.names.write_text(
            json.dumps(
                {
                    "variables_in_order": list(map(str, variables)),
                    "reverse_substitutions": removed[::-1],
                },
                indent=2,
            )
            + "\n"
        )
    print("$fieldsize 23")
    print(f"$vars {len(variables)} X")
    for equation in equations:
        line = str(equation.xreplace(names)).replace("**", "^").replace(" ", "")
        if "(" in line or ")" in line:
            raise ValueError("The M4GB input format requires expanded monomials")
        print(line)
    degrees = Counter(sp.Poly(e, *variables).total_degree() for e in equations)
    term_count = sum(len(sp.Add.make_args(e)) for e in equations)
    print(
        f"{len(variables)} variables, degrees {dict(degrees)}, {term_count} terms",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
