"""Cube-conic filter by osculating the actual (6,6) branch on its quadric."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from itertools import product
from math import comb
from flint import nmod_mpoly_ctx, nmod_poly
import prym_seed as seed

P = 11
R = nmod_mpoly_ctx.get(("Y", "w"), P, "lex")
Y, w = R.gens()
zero = R.constant(0)
one = R.constant(1)
data = seed.data()[0]


def derivative(f, i):
    out = {}
    for k, c in f.items():
        v = k[i]
        if v:
            t = list(k)
            t[i] -= 1
            out[tuple(t)] = c * v % P
    return out


def upoly(F):
    d = F.to_dict()
    return nmod_poly([int(d.get((0, i), 0)) for i in range(int(F.degrees()[1]) + 1)], P)


def content(F):
    d = F.to_dict()
    if not d:
        return nmod_poly([], P)
    n, m = map(int, F.degrees())
    res = nmod_poly([], P)
    for y in range(n + 1):
        poly = nmod_poly([int(d.get((y, i), 0)) for i in range(m + 1)], P)
        res = res.gcd(poly)
        if res.degree() == 0:
            return res
    return res


def embeds(f):
    return R.from_dict({(j, i): int(c) for (i, j), c in f.items()})


EXCEPTION = []


def sdiv(forms, remember=True):
    g = nmod_poly([], P)
    for f in forms:
        g = g.gcd(content(f))
        if g.degree() == 0:
            break
    if g.degree() > 0:
        if remember:
            EXCEPTION.append(g)
        g = R.from_dict({(0, i): int(c) for i, c in enumerate(g) if c})
        return [f // g for f in forms]
    return forms


def primitive(g):
    return sdiv([g])[0]


def ruling(f):
    g = embeds(f)
    # numerator of 27*g^(8/3)*(g^(1/3))'''.
    v = one
    for j in range(3):
        v = 3 * g * v.derivative(0) + (1 - 3 * j) * g.derivative(0) * v
    d = v.to_dict()
    common = nmod_poly([], P)
    for j in range(int(v.degrees()[0]) + 1):
        q = nmod_poly(
            [int(d.get((j, k), 0)) for k in range(int(v.degrees()[1]) + 1)], P
        )
        common = common.gcd(q)
    return common


def osculate(f, extra=0):
    # Change variable Y=a(w)*z to work with a monic sextic P(Y) over F_11[w].
    a = sum((c * w**i for (i, j), c in f.items() if j == 6), zero)
    g = embeds(f)
    EXCEPTION.append(upoly(a))
    A = [a**i for i in range(7)]

    def transform(f, n=6):
        return sum((int(c) * w**i * Y**j * A[n - j] for (i, j), c in f.items()), zero)

    F = transform(f) // a  # monic in Y
    red = lambda b: b % F

    def mul(*xs):
        g = one
        for x in xs:
            g = red(g * x)
        return g

    U, W = [derivative(f, i) for i in (0, 1)]
    fx, fz, fxx, fxz, fzz = [
        red(transform(v))
        for v in [U, W, derivative(U, 0), derivative(U, 1), derivative(W, 1)]
    ]
    D = mul(2 * fx, fz, fz)
    C = mul(fxx, fz, fz) - mul(2 * fxz, fx, fz) + mul(fzz, fx, fx)
    n0, n1, d0, d1 = sdiv(
        [mul(Y, D), -mul(Y, C) - mul(2 * a, fx, fx, fz), mul(a, D), -mul(a, C)]
    )

    def lpower(x, y):
        v = [[one]]
        for _ in range(6):
            prev = v[-1]
            v.append(
                [mul(x, i) + mul(y, j) for i, j in zip(prev + [zero], [zero] + prev)]
            )
        return v

    xp = lpower(w, one)
    yp = lpower(n0, n1)
    zp = lpower(d0, d1)
    # F(w+t, (n0+n1*t)/(d0+d1*t)) times the sixth denominator power.
    tt = [zero for k in range(13)]
    for j in range(7):
        ypz = [zero] * 7
        for k, a0 in enumerate(yp[j]):
            for l, b0 in enumerate(zp[6 - j]):
                ypz[k + l] += mul(a0, b0)
        for i in range(7):
            c = f.get((i, j), 0)
            if not c:
                continue
            for k, a0 in enumerate(xp[i]):
                for l, b0 in enumerate(ypz):
                    tt[k + l] += mul(c * a0, b0)
    assert not any(tt[:3])
    t = sdiv(tt[3:])
    print("residual coefficient degrees", t[0].degrees(), t[-1].degrees(), flush=True)
    # Shift chart w through extra if desired; here use two successive
    # Taylor coefficients of the required fourth derivative of the cube root.
    v = [one]
    for k in range(4):
        out = [zero] * (len(t) + len(v) - 2)
        for i, aa in enumerate(t):
            for j, bb in enumerate(v):
                if i + j:
                    out[i + j - 1] += mul((3 * j + (1 - 3 * k) * i), aa, bb)
        v = out
    gcs = None
    for j in (0, 1, 2):
        b = primitive(v[j]) if j == 0 else v[j]
        print("constraint", j, "degrees", b.degrees(), flush=True)
        if gcs is None:
            res = upoly(F.resultant(b, 0))
            print("first resultant degree", res.degree(), flush=True)
            assert res
            gcs = res // res.gcd(res.derivative()) if res.derivative() else res
            # Squarefree radical in characteristic p, including p-th powers.
            gcs = onefactor(res)
        else:
            am = R.from_dict({(0, i): int(c) for i, c in enumerate(gcs) if c})
            bb = b % am
            FF = F % am
            res = upoly(FF.resultant(bb, 0))
            gcs = gcs.gcd(res)
        print("candidate degree", gcs.degree(), flush=True)
        if gcs.degree() < 25:
            break
    # Discard exceptional x where a vanishes or a point of C is singular;
    # their entire finite fibers require direct verification separately.
    for aa in EXCEPTION:
        aa = onefactor(aa)
        gcs = gcs * aa // gcs.gcd(aa)
    print(
        "candidate degree including content and denominators", gcs.degree(), flush=True
    )
    return gcs


def onefactor(f):
    out = nmod_poly([1], P)
    for k, e in f.factor_squarefree()[1]:
        out *= k
    return out


if __name__ == "__main__":
    for f in (data, {(j, i): v for (i, j), v in data.items()}):
        print(
            "ruling candidates",
            ruling(f),
            "reverse",
            ruling({(6 - i, j): v for (i, j), v in f.items()}),
            flush=True,
        )
    osculate(data)
