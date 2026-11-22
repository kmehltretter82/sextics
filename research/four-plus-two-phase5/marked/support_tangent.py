"""The physical tangent plane of the 11-adic smooth seed, computed independently from its rational chart."""

import numpy as np
from marked_det import inv_det
from marked_axis import reconstruct, seed, QMON

p = 11
b, Q, h = seed()
vals = b + [Q[k] for k in QMON] + [6, 9]
r, s, alpha = 4, 2, 6
src = vals + [r, s]


def F(u):
    return np.array(reconstruct(u[:11], u[11], u[12], p * p)[0])


base = F(src)
cols = []
for i in range(13):
    moved = src.copy()
    moved[i] += p
    diff = (F(moved) - base) % (p * p)
    assert not np.any(diff % p)
    cols.append(diff // p)
J = np.array(cols).T
inv, det = inv_det(J[:, :11], p, p)
der = -inv @ J[:, 11:] % p
assert not np.any((J[:, :11] @ der + J[:, 11:]) % p)
dt = der[4] * pow(2 * alpha, -1, p) % p
dd = (s * dt + np.array([0, alpha])) % p
print("full determinant", det, "Z_r,Z_s", der[4], "t_r,t_s", dt, "d_r,d_s", dd)
assert dt[1] or dd[1]
# Affine coordinates a=1,b=r,c=t,d=st; find a normal vector.
normal = np.cross([1, dt[0], dd[0]], [0, dt[1], dd[1]]) % p
print(
    "normal vector in (b,c,d)",
    normal,
    "tangent plane constant",
    normal @ np.array([r, alpha, s * alpha]) % p,
)
