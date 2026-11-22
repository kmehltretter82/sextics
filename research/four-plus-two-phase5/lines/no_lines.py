"""Generate the seven restriction coefficients on every Schubert cell of Gr(2,4)."""

from itertools import combinations
from pathlib import Path
from sys import argv
from sympy import GF
from sympy.polys.rings import ring


def make(prime):
    out = []
    for i, j in combinations(range(4), 2):
        unknown = [
            (row, k)
            for row, piv in enumerate([i, j])
            for k in range(piv + 1, 4)
            if k not in [i, j]
        ]
        names = [f"n{r}{a}" for r, a in unknown]
        R, *v = ring(",".join(names + ["X", "Y"]), GF(prime))
        X, Y = v[-2:]
        co = dict(zip(unknown, v))
        a, b, c, d = [
            (
                X
                if k == i
                else Y
                if k == j
                else co.get((0, k), 0) * X + co.get((1, k), 0) * Y
            )
            for k in range(4)
        ]
        h = 7 * (a**6 + b**6 + c**6 + d**6) + 5 * (
            a**5 * b
            + 2 * a * b**5
            + c**5 * d
            + 3 * c * d**5
            + a * b * c * d * (a**2 + 2 * c**2)
        )
        eq = [[] for _ in range(7)]
        for mon, coef in h.terms():
            term = "*".join(
                [str(int(coef))]
                + [
                    var if e == 1 else f"{var}^{e}"
                    for var, e in zip(names, mon[:-2])
                    if e
                ]
            )
            eq[mon[-1]].append(term)
        out.append(f"ring r{i}{j}={prime},(" + ",".join(names or ["u"]) + "),dp;")
        out.append("ideal I=" + ",".join("+".join(e) or "0" for e in eq) + ";")
        out.append(f'"chart {i}{j}"; ideal J=slimgb(I); J;')
    return "\n".join(out + ["quit;"])


if __name__ == "__main__":
    Path(argv[1]).write_text(make(11))
