"""Refine each root over C using the original packet's coefficient circuit."""

import json, sys
from pathlib import Path
import mpmath as mp
import sympy as sp


def compile_circuit(packet):
    d = json.load(open(packet))
    syms = sp.symbols(" ".join(d["variable_order"]))
    r, s = sp.symbols("r s")
    g = sp.symbols("g0:19")
    env = {k.name: k for k in syms + (r, s)}
    par = lambda x: sp.sympify(x, locals=env)
    M = sp.Matrix([[par(x) for x in row] for row in d["M_pivot"] + d["M_free"]])
    R = sp.Matrix([par(x) for x in d["R_pivot"] + d["R_free"]])
    f = (R - M * sp.Matrix(g)).col_join(
        sp.Matrix([par(d["point_q"]), par(d["point_h"])])
    )
    z = syms + g
    entries = [(i, j, v) for (i, j), v in f.jacobian(z).todok().items()]

    def fn(v):
        return sp.lambdify(z + (r, s), v, modules="mpmath", cse=True, docstring_limit=0)

    F = fn(list(f))
    J = fn([t[2] for t in entries])
    beta = fn(par(d["point_beta"]))

    def val(z, r, s):
        return mp.matrix(F(*z, r, s))

    def jac(z, r, s):
        m = mp.matrix(46, 46)
        for (i, j, _), v in zip(entries, J(*z, r, s)):
            m[i, j] = v
        return m

    return val, jac, beta


def mpv(p):
    return mp.mpc(*map(str, p))


def fmt(v):
    return [mp.nstr(v.real, 57), mp.nstr(v.imag, 57)]


def certify(end, F, J, beta, r, s):
    z = list(map(mpv, end))
    # Our full solver orders the quartic before alpha. The packet orders alpha first.
    z = mp.matrix(z[:26] + [z[-1]] + z[26:-1])
    initial = mp.mnorm(F(z, r, s), "inf")
    for t in range(4):
        z -= mp.lu_solve(J(z, r, s), F(z, r, s))
    f = F(z, r, s)
    j = J(z, r, s)
    # Eliminate the 19 quartic multipliers in exactly the packet's prescribed pivot order.
    P = j[:19, 27:]
    L = j[19:, 27:]
    invP = P**-1
    original = f[19:, :] - L * (invP * f[:19, :])
    S = j[19:, :27] - L * invP * j[:19, :27]
    jac = mp.det(S)
    d = mp.det(P)
    be = beta(*z, r, s)
    al = z[26]
    k = mp.mpf(5) / 7
    Pval = (
        1
        + r**6
        + k * r
        + 2 * k * r**5
        + k * r * s * al**2
        + 2 * k * r * s * al**4
        + (1 + s**6 + k * s + 3 * k * s**5) * al**6
    )
    assert max(abs(v) for v in original) < mp.mpf("1e-55") and abs(
        be**3 + Pval
    ) < mp.mpf("1e-55")
    return dict(
        initial_full_residual=mp.nstr(initial, 9),
        direct27_residual=mp.nstr(max(map(abs, original)), 9),
        norm_residual=mp.nstr(max(map(abs, f)), 9),
        alpha=fmt(al),
        beta=fmt(be),
        cover_residual=mp.nstr(abs(be**3 + Pval), 9),
        twentyseven_det=fmt(jac),
        pivot_det_abs=mp.nstr(abs(d), 15),
        coordinates=[fmt(v) for v in z[:27]],
    )


def main():
    packet, continuation, out, *new = sys.argv[1:]
    doc = json.load(open(continuation))
    mp.mp.dps = 72
    F, J, beta = compile_circuit(packet)
    print("compiled literal packet circuit", flush=True)
    inputs = [b["end"] for b in doc["branches"]]
    for name in new:
        record = json.load(open(name))
        candidates = (
            [record["endpoint"]]
            if "endpoint" in record
            else [
                item["z"]
                for item in (record["reached"] if isinstance(record, dict) else record)
            ]
        )
        for z in candidates:
            if all(
                max(abs(complex(*a) - complex(*b)) for a, b in zip(z, w)) > 1e-7
                for w in inputs
            ):
                inputs.append(z)
    results = []
    for i, z in enumerate(inputs):
        item = certify(z, F, J, beta, mp.mpf("0.1"), mp.mpf(str(doc["s"])))
        results.append(item)
        print(
            i,
            item["alpha"],
            item["direct27_residual"],
            item["twentyseven_det"],
            flush=True,
        )
    Path(out).write_text(json.dumps(results, indent=2) + "\n")


if __name__ == "__main__":
    main()
