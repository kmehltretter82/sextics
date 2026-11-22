"""Recover the selected cubic from two axis sextics, using only conic data and seed signs."""

from itertools import product
import numpy as np
import segre_test as st

S, ZERO = st.S, st.ZERO


def coefficient(f, i, j):
    return f.get((i, j), ZERO)


def axis_factor(a, u, xseed, vseed):
    x = st.sqrt(-a[3], xseed)
    v = S([vseed, 0, 0])
    for _ in range(2):
        fun = (v * v + a[1]) ** 2 - 8 * u * u * x * v - 4 * u * u * a[2]
        der = 4 * v * (v * v + a[1]) - 8 * u * u * x
        v = v - fun / der
    w = (v * v + a[1]) / (2 * u)
    assert not any((w * w - 2 * v * x - a[2]).v)
    return [u, v, w, x]


def sparse(h, k):
    return {(i, j): h[i][j] for i in range(k) for j in range(4)}


def conjugate(f):
    return {m: (-c if sum(m) % 2 else c) for m, c in f.items()}


def solve(f, seed):
    u = st.sqrt(coefficient(f, 0, 0), coefficient(seed, 0, 0).v[0])
    z = axis_factor(
        [coefficient(f, 0, 2 * i) for i in range(4)],
        u,
        coefficient(seed, 0, 3).v[0],
        coefficient(seed, 0, 1).v[0],
    )
    w = axis_factor(
        [coefficient(f, 2 * i, 0) for i in range(4)],
        u,
        coefficient(seed, 3, 0).v[0],
        coefficient(seed, 1, 0).v[0],
    )
    _, v, y, x = z
    unit = u * x - v * y
    assert unit.v[0] and u.v[0] and x.v[0]
    h = [z]
    for k in range(1, 4):
        old = sparse(h, k)
        prod = st.fmul(old, conjugate(old))

        def half(j):
            return (coefficient(f, k, j) - coefficient(prod, k, j)) / 2

        c0 = w[k]
        if k % 2:
            c1 = (half(1) + v * c0) / u
            base = half(3) - y * c1 + x * c0
            # -v c2 + u c3 = base and -x c2 + y c3 = half(5).
            c2 = (y * base - u * half(5)) / unit
            c3 = (x * base - v * half(5)) / unit
        else:
            c3 = -half(6) / x
            base2 = half(2) - y * c0
            base4 = half(4) + v * c3
            c1 = (y * base2 - u * base4) / unit
            c2 = (x * base2 - v * base4) / unit
        h.append([c0, c1, c2, c3])
    h = st.clean(sparse(h, 4))
    rem = st.fadd(f, st.scale(st.fmul(h, conjugate(h)), -1))
    positions = [(i, j) for i in (4, 5, 6) for j in range(1, 7) if (i + j) % 2 == 0]
    assert len(positions) == 9 and set(rem) <= set(positions)
    return h, [coefficient(rem, *m) for m in positions]


def main():
    import sys

    data = np.load(sys.argv[1])["coordinates"]
    x, C, f, old, _ = st.seed_chart(data)
    h, rem = solve(f, old)
    assert not st.fadd(h, st.scale(old, -1))
    assert not any(any(v.v) for v in rem)
    print(
        "recovered 16 cubic coefficients through tau^2 using only two quadric jets and three axis seed signs/roots."
    )
    pert = st.fadd(C, st.scale(st.power(x[0], 2), S([0, 1, 0])))
    fpert = st.fadd(st.source(x), st.scale(st.power(pert, 3), -1))
    _, rpert = solve(fpert, old)
    print(
        "coupled Q1 perturbation first nonzero degree: 2; second-order residuals:",
        [st.coords(z.v[2]) for z in rpert],
    )
    assert all(not any(z.v[:2]) for z in rpert) and any(z.v[2] for z in rpert)


if __name__ == "__main__":
    main()
