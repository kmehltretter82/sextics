#!/usr/bin/env python3
"""Emit the full incidence equations at (r,s)=(4,6) modulo 23 for Singular."""

import argparse
import json
from pathlib import Path
import sympy as sp


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("circuit", type=Path)
    parser.add_argument(
        "--presolve", action="store_true", help="eliminate variables by exact substitution"
    )
    args = parser.parse_args()
    data = json.loads(args.circuit.read_text())
    names = data["parameter_order"]
    multipliers = ["W" + str(i) for i in range(len(data["quartic_multiplier_basis"]))]
    variables = names + ["alpha"] + multipliers + ["invChart"]
    syms = {x: sp.Symbol(x) for x in variables + ["r", "s"]}
    equations = []
    for block in ("pivot", "free"):
        for value, row in zip(data["R_" + block], data["M_" + block]):
            entry = sp.sympify(value, locals=syms)
            for col, coeff in zip(multipliers, row):
                if coeff != "0":
                    entry -= syms[col] * sp.sympify(coeff, locals=syms)
            equations.append(entry)
    for key in ("point_q", "point_h"):
        expr = sp.sympify(data[key], locals=syms)
        equations.append(expr.subs({syms["r"]: 4, syms["s"]: 6}))
    z = syms["mu0"] * syms["mu1"]
    equations.append(syms["invChart"] * syms["alpha"] * (1 - z) * (1 - 2 * z) * (1 - 3 * z) - 1)
    assert len(variables) == len(equations) == 47
    print("ring base = 23, (" + ",".join(variables) + "), dp;")
    print("ideal incidence = " + ",\n".join(str(e).replace("**", "^") for e in equations) + ";")
    print('print("47 input equations and 47 variables on the open incidence chart");')
    if args.presolve:
        print('LIB "presolve.lib";')
        print("list reduction = elimpart(incidence);")
        print('print("Number of variables removed by substitution"); size(reduction[2]);')
        print("ideal simplified = reduction[1];")
        print("def reducedRing = tolessvars(simplified);")
        print("setring reducedRing;")
        print('print("Remaining variables and equations"); nvars(basering); size(IMAG);')
        print("ideal basis = slimgb(IMAG);")
    else:
        print("ideal basis = slimgb(incidence);")
    print('print("dimension"); dim(basis);')
    print('print("vector space dimension"); vdim(basis);')
    print("quit;")


if __name__ == "__main__":
    main()
