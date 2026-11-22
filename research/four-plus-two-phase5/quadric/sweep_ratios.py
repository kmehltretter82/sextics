"""Revert the smooth-contact parameter using public python-flint polynomial APIs."""

from argparse import ArgumentParser
from pathlib import Path
from time import monotonic
import numpy as np
from flint import nmod_poly


P = 11


def poly(coeff):
    return nmod_poly(coeff, P)


def compose(coeff, inner, n):
    """Compose modulo x^n via f(X)=sum X^j f_j(X)^11 for coefficient decimations f_j."""
    assert inner[0] == 0
    powers = {}

    def visit(coeff, n):
        coeff = coeff[:n]
        if n <= 24:
            if not len(coeff):
                return poly([])
            return poly(coeff.tolist()).compose_mod(
                inner.truncate(n), poly([1]).left_shift(n)
            )
        if n not in powers:
            g = inner.truncate(n)
            parts = [poly([1]), g]
            for j in range(2, P):
                parts.append(parts[-1].mul_low(g, n))
            powers[n] = parts
        m = (n + P - 1) // P
        out = poly([])
        for j, gp in enumerate(powers[n]):
            b = visit(coeff[j::P], m)
            if not b:
                continue
            lifted = [0] * (P * (len(b) - 1) + 1)
            lifted[::P] = list(b)
            if j == 0:
                out += poly(lifted).truncate(n)
            else:
                out += gp.mul_low(poly(lifted), n)
        return out

    return visit(coeff, n)


def revert(coeff, n):
    assert coeff[0] == 0 and coeff[1]
    deriv = np.array([(i * int(coeff[i])) % P for i in range(1, n)], dtype=np.uint8)
    x = poly([0, 1])
    v = poly([0, pow(int(coeff[1]), -1, P)])
    k = 2
    while k < n:
        k = min(2 * k, n)
        diff = compose(coeff, v, k) - x
        slope = compose(deriv, v, k).inverse_series_trunc(k)
        v -= diff.mul_low(slope, k)
    return v


def array(f, n):
    a = np.zeros(n, dtype=np.uint8)
    a[: len(f)] = list(f)
    return a


def main():
    ap = ArgumentParser(description=__doc__)
    ap.add_argument("source", type=Path)
    ap.add_argument("out", type=Path)
    ap.add_argument("-n", type=int)
    ap.add_argument("--sheets", type=int, nargs="+", default=[0])
    args = ap.parse_args()
    src = np.load(args.source)
    out = {}
    for j in args.sheets:
        n = args.n or len(src[f"r{j}"])
        rc = src[f"r{j}"][:n].copy()
        r0 = int(rc[0])
        rc[0] = 0
        start = monotonic()
        zeta = poly(rc.tolist())
        v = revert(rc, n)
        z = compose(src[f"Z{j}"], v, n)
        t0 = poly(src[f"a{j}"][:n].tolist()).inverse_series_trunc(n)
        t = compose(array(t0, n), v, n)
        assert t.mul_low(t, n) == z
        # Both directions are checked in the truncated power series ring.
        assert compose(rc, v, n) == poly([0, 1])
        assert compose(array(z, n), zeta, n) == poly(src[f"Z{j}"][:n].tolist())
        out[f"Z{j}"] = array(z, n)
        out[f"t{j}"] = array(t, n)
        out[f"r0_{j}"] = r0
        print(
            "sheet",
            j,
            "r0",
            r0,
            "Z(r0+x)",
            out[f"Z{j}"][:10].tolist(),
            round(monotonic() - start, 2),
            "seconds",
            flush=True,
        )
    np.savez_compressed(args.out, **out)


if __name__ == "__main__":
    main()
