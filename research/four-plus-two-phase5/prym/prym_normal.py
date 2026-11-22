"""Check the four intrinsic normal equations of the Prym contact locus."""

from itertools import product
from flint import nmod_mat
import prym_seed as seed

P = 11
base = (2, 6, 1, 7, 7)


def restricted_core(point, modulus):
    """The physical sextic on the marked Segre chart for an arbitrary Q."""
    m, n, u, v, e = point
    seed.p = modulus
    K = e * (v - u) ** 2 * (m + n) * pow((n - m) ** 2, -1, modulus) % modulus
    a = {(0, 0): 1, (1, 1): K}
    b = {(0, 0): n, (1, 1): m * K % modulus}
    c = {(1, 0): 1, (0, 1): 1}
    d = {(1, 0): v, (0, 1): u}
    k = 5 * pow(7, -1, modulus) % modulus
    terms = [
        (1, [a] * 6),
        (1, [b] * 6),
        (1, [c] * 6),
        (1, [d] * 6),
        (k, [a] * 5 + [b]),
        (2 * k, [a] + [b] * 5),
        (k, [c] * 5 + [d]),
        (3 * k, [c] + [d] * 5),
        (k, [a] * 3 + [b, c, d]),
        (2 * k, [a, b] + [c] * 3 + [d]),
    ]
    h = seed.add(*[(i, seed.mul(*f)) for i, f in terms])
    return h, (a, b, c, d)


def flat(f, degree, parity=None):
    return [
        f.get(k, 0)
        for k in product(range(degree + 1), repeat=2)
        if parity is None or sum(k) % 2 == parity
    ]


def normalized_kernel(M):
    S, dim = M.nullspace()
    assert dim == 1
    vector = [int(S[i, 0]) % P for i in range(S.nrows())]
    first = next(v for v in vector if v)
    return [v * pow(first, -1, P) % P for v in vector]


def main():
    seed.p = P
    H, q, h, hb, _ = seed.data()
    HH, (a, b, c, d) = restricted_core(base, P)
    assert HH == H
    basis = [seed.mul(x, y) for x in (a, b) for y in (c, d)]
    q2 = seed.mul(q, q)
    socle = nmod_mat(
        seed.rows_for([(q2, (4, 4)), (h, (3, 3)), (hb, (3, 3))], (8, 8)), P
    )
    pairing = normalized_kernel(socle)
    assert socle.rank() == 80
    assert all(
        not coeff or (i // 9 + i % 9) % 2 == 1 for i, coeff in enumerate(pairing)
    )

    def value(f):
        return sum(v * u for v, u in zip(pairing, flat(f, 8))) % P

    # Differentiate H on the quadric chart by exact dual differences modulo p².
    # Keeping intrinsic q,h fixed is a local trivialization; the socle pairing
    # annihilates all changes in these coordinates.
    H0, _ = restricted_core(base, P**2)
    dH = []
    for i in range(5):
        test = list(base)
        test[i] += P
        H1, _ = restricted_core(test, P**2)
        keys = H0.keys() | H1.keys()
        delta = {f: (H1.get(f, 0) - H0.get(f, 0)) % P**2 for f in keys}
        assert all(v % P == 0 for v in delta.values())
        dH.append({f: v // P for f, v in delta.items() if v})
    seed.p = P
    normal = nmod_mat([[value(seed.mul(m, t)) for t in dH] for m in basis], P)
    print("four mixed bilinear normal rows in (m,n,u,v,e)", normal.tolist())
    tangent = normalized_kernel(normal)
    print(
        "projected parameter tangent",
        tangent,
        "first four columns minor",
        int(nmod_mat([row[:4] for row in normal.tolist()], P).det()),
    )
    assert normal.rank() == 4
    # Independently quotient the unrestricted 25-condition linearized norm map
    # by all 21 fixed-Q directions and recover the same line in the five Q variables.
    unknowns = [
        seed.mul(q2, {k: 3}) for k in product(range(3), repeat=2) if sum(k) % 2 == 0
    ]
    for k in product(range(4), repeat=2):
        s = seed.mul(hb, {k: 1})
        t = seed.mul(h, {k: (-1) ** sum(k)})
        unknowns.append(seed.add((1, s), (1, t)))
    assert len(unknowns) == 21
    fixed = nmod_mat([flat(x, 6, 0) for x in unknowns], P)
    full = nmod_mat([flat(x, 6, 0) for x in unknowns + dH], P)
    assert fixed.rank() == 21 and full.rank() == 25
    assert all(value(seed.mul(m, x)) == 0 for x in unknowns for m in basis)
    assert all(int(v) % P == 0 for v in normal * nmod_mat([[k] for k in tangent], P))
    combo = seed.add(*[(v, f) for v, f in zip(tangent, dH)])
    assert nmod_mat(fixed.tolist() + [flat(combo, 6, 0)], P).rank() == 21
    print("independent fixed and full ranks", fixed.rank(), full.rank())


if __name__ == "__main__":
    main()
