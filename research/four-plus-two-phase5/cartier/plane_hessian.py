"""Check absence of triple tangent planes on the fixed rational sextic."""

from itertools import permutations
from flint import nmod_mat


def add_monomial(p, e, k):
    p[e] = p.get(e, 0) + k
    if p[e] == 0:
        del p[e]


def mul(f, g):
    h = {}
    for e, k in f.items():
        for u, coefficient in g.items():
            add_monomial(h, tuple(a + b for a, b in zip(e, u)), k * coefficient)
    return h


def derivative(f, i):
    return {e[:i] + (e[i] - 1,) + e[i + 1 :]: v * e[i] for e, v in f.items() if e[i]}


def determinant(matrix):
    h = {}
    for perm in permutations(range(len(matrix))):
        sign = (-1) ** sum(
            perm[i] > perm[j] for i in range(len(perm)) for j in range(i + 1, len(perm))
        )
        f = {(0, 0, 0, 0): sign}
        for i, j in enumerate(perm):
            f = mul(f, matrix[i][j])
        for e, v in f.items():
            add_monomial(h, e, v)
    return h


def monomials(d):
    return [
        (a, b, c, d - a - b - c)
        for a in range(d + 1)
        for b in range(d - a + 1)
        for c in range(d - a - b + 1)
    ]


def generators():
    # 7H, in the physical order (a,b,c,d).
    h = {
        (6, 0, 0, 0): 7,
        (0, 6, 0, 0): 7,
        (0, 0, 6, 0): 7,
        (0, 0, 0, 6): 7,
        (5, 1, 0, 0): 5,
        (1, 5, 0, 0): 10,
        (0, 0, 5, 1): 5,
        (0, 0, 1, 5): 15,
        (3, 1, 1, 1): 5,
        (1, 1, 3, 1): 10,
    }
    hes = [[derivative(derivative(h, i), j) for j in range(4)] for i in range(4)]
    minors = []
    for i in range(4):
        for j in range(i + 1):
            mat = [[hes[u][v] for v in range(4) if v != j] for u in range(4) if u != i]
            minors.append(determinant(mat))
    return [h] + minors


def run(p=13, d=23):
    gens = generators()
    base = monomials(d)
    ranks = []
    for parity in (0, 1):
        cols = [e for e in base if sum(e[2:]) % 2 == parity]
        pos = {e: i for i, e in enumerate(cols)}
        matrices = []
        for f in gens:
            poly = {e: v % p for e, v in f.items() if v % p}
            if not poly:
                continue
            ed = next(iter(poly))
            deg = sum(ed)
            char = sum(ed[2:]) % 2
            assert all(sum(e) == deg and sum(e[2:]) % 2 == char for e in poly)
            for u in monomials(d - deg):
                if (sum(u[2:]) + char) % 2 != parity:
                    continue
                row = [0] * len(cols)
                for e, v in poly.items():
                    row[pos[tuple(a + b for a, b in zip(e, u))]] = v
                matrices.append(row)
        M = nmod_mat(matrices, p)
        rank = M.rank()
        print(
            "modulus",
            p,
            "degree",
            d,
            "parity",
            parity,
            "matrix",
            len(matrices),
            len(cols),
            "rank",
            rank,
            flush=True,
        )
        if rank == len(cols):
            reduced, rr = M.transpose().rref()
            chosen = [
                next(j for j, k in enumerate(row) if k) for row in reduced.tolist()[:rr]
            ]
            det = int(nmod_mat([matrices[j] for j in chosen], p).det())
            print("chosen full minor", det, flush=True)
            assert det

        ranks.append((rank, len(cols)))
    return ranks


if __name__ == "__main__":
    from sys import argv

    rr = run(*map(int, argv[1:]))
    assert all(a == b for a, b in rr)
