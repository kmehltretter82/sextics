"""Verify all previously isolated ruling-coordinate fibers for a cubic conic restriction."""

from itertools import product
from flint import fmpz_mod_poly_ctx, fq_default_ctx, fq_default_poly_ctx
from osculate import osculate, data, derivative


def rev(f, i, j):
    return {(6 - k if i else k, 6 - l if j else l): v for (k, l), v in f.items()}


def fixed_fiber(f, K, x):
    R = fq_default_poly_ctx(K)
    y = R.gen()

    def atx(ff):
        poly = R.zero()
        for (k, l), c in ff.items():
            poly += R([K(c) * x**k]) * y**l
        return poly

    F = atx(f)
    assert F
    count = 0
    for g, mult in F.factor()[1]:

        def ev(h):
            return atx(h) % g

        U, V = derivative(f, 0), derivative(f, 1)
        fw, fz = ev(U), ev(V)
        if not fw or not fz:
            continue
        w2, wz, zz = [
            ev(v) for v in (derivative(U, 0), derivative(U, 1), derivative(V, 1))
        ]
        E = 2 * fw * fz**2 % g
        A = (w2 * fz**2 - 2 * wz * fw * fz + zz * fw**2) % g
        n0 = y * E % g
        n1 = (-y * A - 2 * fw**2 * fz) % g
        d0 = E
        d1 = -A

        def linpow(a, b):
            rows = [[R.one()]]
            for _ in range(6):
                prev = rows[-1]
                rows.append(
                    [
                        (a * i + b * j) % g
                        for i, j in zip(prev + [R.zero()], [R.zero()] + prev)
                    ]
                )
            return rows

        xx = linpow(R([x]), R.one())
        nn = linpow(n0, n1)
        dd = linpow(d0, d1)
        Fline = [R.zero() for _ in range(13)]
        for (k, l), c in f.items():
            for i, a in enumerate(xx[k]):
                for j, b in enumerate(nn[l]):
                    for m, d in enumerate(dd[6 - l]):
                        Fline[i + j + m] += c * a * b * d
        Fline = [v % g for v in Fline]
        assert not any(Fline[:3])
        Fline = Fline[3:]
        v = [R.one()]
        for m in range(4):
            out = [R.zero() for _ in range(len(Fline) + len(v) - 2)]
            for i, a in enumerate(Fline):
                for j, b in enumerate(v):
                    if i + j:
                        out[i + j - 1] += (3 * j + (1 - 3 * m) * i) * a * b
            v = [a % g for a in out]
        assert any(v), ("possible cube conic", x, g)
        count += g.degree()
    return count


def main():
    factors = osculate(data).factor()[1]
    checks = 0
    matrices = fmpz_mod_poly_ctx(11)
    for f, mult in factors:
        assert mult == 1
        K = fq_default_ctx(modulus=matrices([int(c) for c in f]), var="a")
        x = K.gen()
        for r in (False, True):
            checks += f.degree() * fixed_fiber(rev(data, False, r), K, x)
    K = fq_default_ctx(11, 1, "a")
    for r in (False, True):
        checks += fixed_fiber(rev(data, True, r), K, K(0))
    print(
        "smooth noncritical geometric point checks including chart overlaps",
        checks,
        flush=True,
    )


if __name__ == "__main__":
    main()
