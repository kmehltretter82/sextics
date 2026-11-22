"""Blow up a nodal solution and follow its neighboring sheets at fixed ratio s."""

from pathlib import Path
import itertools, json
import numpy as np
from contact import nodal_embedding

K = 5 / 7


def basis(d, p=0):
    return sorted(
        e
        for e in itertools.product(range(d + 1), repeat=4)
        if sum(e) == d and (e[2] + e[3]) % 2 == p
    )


q = [(2, 0, 0, 0), (0, 2, 0, 0), (0, 0, 2, 0), (0, 0, 1, 1), (0, 0, 0, 2)]
ab = (1, 1, 0, 0)
aa = [
    (3, 0, 0, 0),
    (0, 3, 0, 0),
    (1, 0, 2, 0),
    (1, 0, 1, 1),
    (1, 0, 0, 2),
    (0, 1, 2, 0),
    (0, 1, 1, 1),
    (0, 1, 0, 2),
]
bb = [
    (2, 0, 1, 0),
    (2, 0, 0, 1),
    (0, 2, 1, 0),
    (0, 2, 0, 1),
    (0, 0, 3, 0),
    (0, 0, 2, 1),
    (0, 0, 1, 2),
    (0, 0, 0, 3),
]
out = basis(6)
gb = basis(4)
qp = [ab] + q
Hterms = [
    ((6, 0, 0, 0), 1),
    ((5, 1, 0, 0), K),
    ((3, 1, 1, 1), K),
    ((1, 5, 0, 0), 2 * K),
    ((1, 1, 3, 1), 2 * K),
    ((0, 6, 0, 0), 1),
    ((0, 0, 6, 0), 1),
    ((0, 0, 5, 1), K),
    ((0, 0, 1, 5), 3 * K),
    ((0, 0, 0, 6), 1),
]
H = np.zeros(44)
for e, t in Hterms:
    H[out.index(e)] = t


def tensor(*bases):
    a = np.zeros((44,) + tuple(map(len, bases)))
    for idx in itertools.product(*[range(len(b)) for b in bases]):
        row = out.index(tuple(map(sum, zip(*[b[i] for b, i in zip(bases, idx)]))))
        a[(row,) + idx] = 1
    return a


tq, ta, tb, tg = tensor(q, q, q), tensor(aa, aa), tensor(bb, bb), tensor(qp, gb)
eins = np.einsum


def pack(z):
    return z[:5], np.r_[1, z[5:10]], z[10:18], z[18:26], z[26:45]


def norm(z):
    Q, M, A, B, G = pack(z)
    return (
        H
        - eins("xijk,i,j,k->x", tq, Q, Q, Q)
        - eins("xij,i,j->x", ta, A, A)
        + eins("xij,i,j->x", tb, B, B)
        - eins("xij,i,j->x", tg, M, G)
    )


def dnorm(z):
    Q, M, A, B, G = pack(z)
    return np.column_stack(
        (
            -3 * eins("xijk,i,j->xk", tq, Q, Q),
            -eins("xij,j->xi", tg[:, 1:, :], G),
            -2 * eins("xij,i->xj", ta, A),
            2 * eins("xij,i->xj", tb, B),
            -eins("xij,i->xj", tg, M),
        )
    )


def norm2(z, v):
    Q, _, _, _, _ = pack(z)
    V, M, A, B, G = pack(v)
    return (
        -3 * eins("xijk,i,j,k->x", tq, Q, V, V)
        - eins("xij,i,j->x", ta, A, A)
        + eins("xij,i,j->x", tb, B, B)
        - eins("xij,i,j->x", tg[:, 1:, :], M[1:], G)
    )


def vals(b, r, s, a, dr=0, da=0):
    return np.array(
        [
            (j if dr else 1)
            * (k + l if da else 1)
            * r ** (j - dr)
            * s**l
            * a ** (k + l - da)
            if j >= dr and k + l >= da
            else 0
            for _, j, k, l in b
        ]
    )


def point(z, r, s, dr=0, da=0):
    Q, M, A, B, G = pack(z)
    a = z[-1]
    return np.array(
        [
            M @ vals(qp, r, s, a, dr, da),
            A @ vals(aa, r, s, a, dr, da) + B @ vals(bb, r, s, a, dr, da),
        ]
    )


def rhs(z, r, s):
    return np.r_[norm(z), point(z, r, s)]


def jac(z, r, s):
    J = np.zeros((46, 46), dtype=z.dtype)
    J[:44, :45] = dnorm(z)
    J[44, 5:10] = vals(q, r, s, z[-1])
    J[45, 10:18] = vals(aa, r, s, z[-1])
    J[45, 18:26] = vals(bb, r, s, z[-1])
    J[44:, -1] = point(z, r, s, da=1)
    return J


def nodalseed(C, B):
    z = np.zeros(45)
    z[2:5] = C
    z[10:12] = 1
    z[22:26] = B
    for e, v in [
        ((4, 0, 0, 0), K),
        ((0, 4, 0, 0), 2 * K),
        ((2, 2, 0, 0), -2),
        ((2, 0, 1, 1), K),
        ((0, 0, 3, 1), 2 * K),
    ]:
        z[26 + gb.index(e)] = v
    return z


def slopes(z):
    # The two-dimensional right kernel, one-dimensional left cokernel,
    # and quadratic form on the kernel give the two node tangents.
    L, S, VH = np.linalg.svd(dnorm(z), full_matrices=True)
    assert S[-1] < 1e-12 and S[-2] > 1e-6
    e, f = VH[-2], VH[-1]
    ell = L[:, -1]
    F = ell @ norm2(z, e)
    G = ell @ norm2(z, f)
    J = ell @ norm2(z, e + f) - F - G
    tangents = [e + k * f for k in np.roots([G, J, F])]
    # Sort by a chart-coordinate ratio, invariant under changing the SVD kernel basis.
    return S, sorted(
        [(v[2] / v[0], v) for v in tangents], key=lambda v: v[0].real, reverse=True
    )


def newton(z, r, s):
    for i in range(15):
        e = rhs(z, r, s)
        J = jac(z, r, s)
        dz = np.linalg.solve(J, -e)
        if max(abs(e)) < 4e-13 and max(abs(dz)) < 7e-10:
            break
        z = z + dz
    if max(abs(rhs(z, r, s))) > 1e-10 or max(abs(dz)) > 1e-8:
        raise ArithmeticError((r, i, max(abs(rhs(z, r, s))), max(abs(dz))))
    return z, i, float(max(abs(rhs(z, r, s)))), float(np.linalg.cond(J))


def data(z0, r0, s, sl, rgrid, alpha_side):
    B0 = pack(z0)[3]
    BS = sum(B0[4 + i] * s**i for i in range(4))
    a0 = np.cbrt(-1 / BS) * np.exp(2j * np.pi * alpha_side / 3)
    z = np.r_[z0.astype(complex), a0]
    J = jac(z, r0, s)
    zdir = sl.astype(complex)
    vel = -zdir / (J[44, :45] @ zdir)
    alpha_vel = -(J[45, :45] @ vel + point(z, r0, s, dr=1)[1]) / J[45, 45]
    vel = np.r_[vel, alpha_vel]
    assert max(abs(J @ vel + np.r_[np.zeros(44), point(z, r0, s, dr=1)])) < 1e-12
    records = []
    for r in rgrid:
        dest, steps, res, cond = newton(z + (r - r0) * vel, r, s)
        Q, M, A, B, G = pack(dest)
        if (
            abs((1 - M[1] * M[2]) * (1 - 3 * M[1] * M[2]) * (1 - 2 * M[1] * M[2]) ** 3)
            < 1e-5
        ):
            raise ArithmeticError("chart pole")
        nerr = max(abs(norm(dest)))
        records.append(
            dict(
                r=float(r),
                steps=int(steps),
                residual=res,
                cond=cond,
                alpha=[dest[-1].real, dest[-1].imag],
                norm_residual=float(nerr),
            )
        )
        z = dest
        r0 = r
        J = jac(z, r, s)
        vel = np.linalg.solve(J, -np.r_[np.zeros(44), point(z, r, s, dr=1)])
    return records, z


def main():
    import argparse

    p = argparse.ArgumentParser()
    p.add_argument("rur")
    p.add_argument("--s", type=float, default=1.0)
    p.add_argument(
        "--verifier",
        type=Path,
        default=Path(__file__).resolve().parent.parent / "verify_rur.py",
    )
    p.add_argument("--output", type=Path, default=Path("continuation.json"))
    a = p.parse_args()
    m, C, B, err = nodal_embedding(a.rur, a.verifier)
    z0 = nodalseed(C, B)
    sv, branches = slopes(z0)
    print(
        "nodal binary root",
        m,
        "C",
        C,
        "B",
        B,
        "identity error",
        err,
        "nodal norm",
        max(abs(norm(z0))),
        "two smallest singular values",
        sv[-2:],
        flush=True,
    )
    out = dict(
        m=m,
        s=a.s,
        C=C.tolist(),
        B=B.tolist(),
        base_error=err,
        start_norm=float(max(abs(norm(z0)))),
        nodal_singular_values=sv[-2:].tolist(),
        branches=[],
    )
    rgrid = np.r_[0.002, 0.004, 0.008, 0.012, np.arange(0.016, 0.1001, 0.004)]
    for j, (k, sl) in enumerate(branches):
        for side in range(3):
            chain, final = data(z0, 0.0, a.s, sl, rgrid, side)
            print(
                j,
                side,
                "slope",
                k,
                "alpha",
                chain[-1]["alpha"],
                "max res",
                max(w["residual"] for w in chain),
                "final condition",
                chain[-1]["cond"],
                flush=True,
            )
            out["branches"].append(
                dict(
                    tangent=j,
                    alpha_choice=side,
                    lambda2_over_lambda0=[k.real, k.imag],
                    states=chain,
                    end=[(z.real, z.imag) for z in final],
                )
            )
    a.output.write_text(json.dumps(out, indent=2) + "\n")


if __name__ == "__main__":
    main()
