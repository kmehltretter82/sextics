"""Exclude all projective cube conics whose intersections lie at singular branch points."""

from itertools import combinations, permutations
from math import comb
from flint import fq_default_ctx, fq_default_poly_ctx
from osculate import data, ruling, seed
from conic_check import rev

K = fq_default_ctx(11, 4, "a")
R = fq_default_poly_ctx(K)
t = R.gen()


def determinant(rows):
    n = len(rows)
    out = K(0)
    for perm in permutations(range(n)):
        val = K(
            (-1) ** sum(perm[i] > perm[j] for i in range(n) for j in range(i + 1, n))
        )
        for i, j in enumerate(perm):
            val *= rows[i][j]
        out += val
    return out


def node_points():
    points = [(True, K(0), K(0))]
    for z, _ in (t**4 + 1).roots():
        w = (z**3 + 2 * z) / K(5)
        points.append((False, w, z))
    assert len(points) == 5
    return points


def local_conics(point):
    inv, w, z = point
    f = rev(data, inv, False)

    def H(a, b):
        val = K(0)
        for (i, j), v in f.items():
            if i >= a and j >= b:
                val += K(v * comb(i, a) * comb(j, b)) * w ** (i - a) * z ** (j - b)
        return val

    assert all(H(a, k - a) == 0 for k in (0, 1) for a in range(k + 1))
    quad = R([H(2, 0), H(1, 1), H(0, 2)])
    assert H(1, 1) ** 2 - 4 * H(2, 0) * H(0, 2) != 0
    cubic = R([H(3, 0), H(2, 1), H(1, 2), H(0, 3)])
    roots = quad.roots()
    assert len(roots) == quad.degree()
    for slope, mul in roots:
        if not slope:
            continue
        k = -cubic(slope) / (quad.derivative()(slope) * slope)
        row = [-z * (1 + k * w) + slope * w, k * z - slope, 1 + k * w, -k]
        yield [row[1], row[0], row[3], row[2]] if inv else row


def triple_conics(points):
    for triple in combinations(points, 3):
        rows = [
            [K(0), K(1), K(0), z] if inv else [K(1), w, z, w * z]
            for inv, w, z in triple
        ]
        co = [
            K((-1) ** j) * determinant([v[:j] + v[j + 1 :] for v in rows])
            for j in range(4)
        ]
        if co[0] * co[3] != co[1] * co[2]:
            yield co


def cobound(coeff):
    a, b, c, d = coeff
    assert a * d - b * c
    num = R([-a, -b])
    den = R([c, d])
    val = R.zero()
    for (i, j), v in data.items():
        val += v * t**i * num**j * den ** (6 - j)
    assert val
    return [12 - val.degree()] + [mul for pol, mul in val.factor()[1]]


def main():
    for f in (data, {(j, i): v for (i, j), v in data.items()}):
        for inv in (False, True):
            assert ruling(rev(f, inv, False)).degree() == 0
    ranks = seed.macaulay(
        list(seed.homogeneous_derivatives(data, (6, 6))), (19, 19), False
    )
    assert ranks[1:3] == (400, 395)
    pts = node_points()
    conics = [c for p in pts for c in local_conics(p)] + list(triple_conics(pts))
    tests = 0
    for co in conics:
        es = cobound(co)
        assert any(e % 3 for e in es), ("cube conic", co, es)
        tests += 1
    print(
        "Geometric ruling restrictions excluded in both directions; gradient matrix",
        ranks[:3],
    )
    print("Ordinary nodes", len(pts), "osculating and triple conic checks", tests)


if __name__ == "__main__":
    main()
