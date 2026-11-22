"""Anti-sign degree-three-cover Frobenius traces via F_(p^2) point counts."""

import numpy as np
from itertools import product
from time import perf_counter
import sys


def run(p):
    k = 5 * pow(7, -1, p) % p
    q = p * p
    nonsquare = next(w for w in range(2, p) if pow(w, (p - 1) // 2, p) == p - 1)
    el = np.arange(q)
    re = el % p
    im = el // p
    add = (re[:, None] + re) % p + p * ((im[:, None] + im) % p)
    mul = (re[:, None] * re + nonsquare * im[:, None] * im) % p + p * (
        (re[:, None] * im + im[:, None] * re) % p
    )

    def power(x, n):
        a = 1
        while n:
            if n % 2:
                a = mul[a, x]
            x = mul[x, x]
            n //= 2
        return a

    f = add[add[1, power(el, 6)], add[mul[k, el], mul[2 * k % p, power(el, 5)]]]
    g = add[add[1, power(el, 6)], add[mul[k, el], mul[3 * k % p, power(el, 5)]]]
    f = np.r_[f, 1]
    g = np.r_[g, 1]
    ab = np.r_[el, 0]
    a2 = np.r_[np.ones(q, dtype=int), 0]
    zeta = next(x for x in range(2, q) if power(x, 3) == 1)
    # sum of the numbers of cube roots for three cubic twists
    cub = power(el, (q - 1) // 3)
    cube_weight = np.array(
        [
            [
                1 if j == 0 else (3 if cub[j] == power(zeta, (-e) % 3) else 0)
                for j in range(q)
            ]
            for e in range(3)
        ]
    )
    quad = power(el, (q - 1) // 2)
    cats = np.zeros(q, dtype=np.int64)
    for z in range(1, q):
        unit = 1 if quad[z] == 1 else -1
        z2 = mul[z, z]
        z3 = mul[z, z2]
        # Coordinates X=(1,r) or (0,1); Y=(1,s) or (0,1). t^2=z.
        v = add[f[:, None], mul[g, z3]]
        coeff = mul[k, mul[ab[:, None], ab]]
        mixed = mul[coeff, add[mul[a2, z][:, None], mul[a2, mul[2, z2]]]]
        hv = add[v, mixed]
        cats -= unit * np.bincount(hv.ravel(), minlength=q)
    out = cube_weight @ cats
    print(
        "p",
        p,
        "trace signs",
        out.tolist(),
        "tau(E) rational coefficients",
        int((out[0] - out[1]) // 3),
        int((out[2] - out[1]) // 3),
        "trace mod p^3",
        out % p**3,
        flush=True,
    )


if __name__ == "__main__":
    for arg in sys.argv[1:]:
        run(int(arg))
