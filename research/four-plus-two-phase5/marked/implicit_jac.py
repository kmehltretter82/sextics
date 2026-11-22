"""Check the seven-parameter implicit equations to first order modulo 121."""

import numpy as np
from marked_axis import RES, geometry, seed, solve3, norm
from marked_det import inv_det


def residuals(point, mod):
    x, y, v, j, Z, q00, q20, T, lam, mu, vz, ez = point
    q02, q11 = 4 + mu, 7 + lam
    q22 = 1 - q00 - q20 - q02 - q11
    F, _ = geometry([x, y, v, j, Z, q02, q11, q20, q22], 4, 2, mod)
    inv0 = pow(F[0, 0], -1, mod)
    g = lambda i, j: inv0 * F.get((j, i), 0) % mod  # row order z,w
    half = pow(2, -1, mod)
    L = -g(0, 6)
    A = (T * T + L + g(0, 2)) * half % mod
    ew = (A * A + T * T * L + 2 * L - g(0, 4)) * pow(2 * T * (A + 1), -1, mod) % mod
    vw, yw = (T - ew) % mod, (A - T * ew) % mod
    yz = (vz * vz + g(2, 0)) * half % mod
    h = {(0, i): u for i, u in enumerate([1, vw, yw, ew])}
    h.update({(1, 0): vz, (2, 0): yz, (3, 0): ez})
    for i in (1, 2, 3):
        if i % 2:
            jj = [1, 3, 5]
            mat = [[1, 0, 0], [yw, -vw, 1], [0, -ew, yw]]
        else:
            jj = [2, 4, 6]
            mat = [[-vw, 1, 0], [-ew, yw, -vw], [0, 0, -ew]]
        rhs = [(g(i, k) - norm(h, i, k, mod)) * half % mod for k in jj]
        for k, value in enumerate(solve3(mat, rhs, mod), 1):
            h[i, k] = value

    def row(i, k):
        return (norm(h, i, k, mod) - g(i, k)) % mod

    vertical = [
        (ew * ew - L) % mod,
        row(4, 0),
        row(6, 0),
        sum(h.values()) % mod,
        row(5, 1),
    ]
    return np.array(vertical + [row(i, k) for i, k in RES if (i, k) != (5, 1)])


def check():
    m, p = 121, 11
    params, Q, _ = seed()
    point = params + [Q[0, 0], Q[2, 0], 6, 0, 0, 10, 10]
    assert np.all(residuals(point, p) == 0)
    base = residuals(point, m)
    cols = []
    for i in range(12):
        moved = point.copy()
        moved[i] += p
        diff = (residuals(moved, m) - base) % m
        assert np.all(diff % p == 0)
        cols.append(diff // p)
    J = np.array(cols).T
    inverse, vdet = inv_det(J[:5, 7:], p, p)
    reduced = (J[5:, :7] - J[5:, 7:] @ inverse @ J[:5, :7]) % p
    _, rdet = inv_det(reduced, p, p)
    _, total = inv_det(J, p, p)
    assert total == (-vdet * rdet) % p  # 7*5 odd
    print(
        "vertical Jacobian",
        vdet,
        "reduced seven-parameter Jacobian",
        rdet,
        "full Jacobian",
        total,
    )


if __name__ == "__main__":
    check()
