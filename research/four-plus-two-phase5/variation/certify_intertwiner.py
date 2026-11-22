"""A three-direction scalar intertwiner certificate in the cyclic family at the literal H."""

from flint import nmod_mat
from certify_span import jacobian_columns, plus

P = 23
EIGENVECTOR = [14, 14, 21, 3, 3, 10, 21, 19, 9, 0, 18, 0, 0, 5, 3, 1]


def quotient_basis(n):
    raw, jacobian, _ = jacobian_columns(n)
    relations, rank = nmod_mat(jacobian, P).rref()
    rows = relations.tolist()[:rank]
    pivot = [next(i for i, v in enumerate(row) if v) for row in rows]
    free = [i for i in range(len(raw)) if i not in pivot]
    quotient = nmod_mat(len(free), len(raw), P)
    for i, k in enumerate(free):
        quotient[i, k] = 1
        for j, h in enumerate(pivot):
            quotient[i, h] = -relations[j, k]
    return [raw[i] for i in free], {v: i for i, v in enumerate(raw)}, quotient


def run():
    source, _, _ = quotient_basis(6)
    target, index, reduce = quotient_basis(12)
    assert (len(source), len(target)) == (32, 16)

    def mult(term):
        return nmod_mat(
            [
                [reduce[i, index[plus(v, term)]] for v in source]
                for i in range(len(target))
            ],
            P,
        )

    # L_b and L_d have common source D_6^- and common target D_12^-.
    b = mult((0, 6, 0, 0))
    d = mult((0, 0, 0, 6))
    a = mult((6, 0, 0, 0))
    stacked = nmod_mat(b.tolist() + d.tolist(), P)
    assert int(stacked.det()) == 20
    combined = a * stacked.inv()
    c = nmod_mat([row[:16] for row in combined.tolist()], P)
    h = nmod_mat([row[16:] for row in combined.tolist()], P)
    eye = nmod_mat([[int(i == j) for j in range(16)] for i in range(16)], P)
    v = nmod_mat([[a] for a in EIGENVECTOR], P)
    assert c * v == v and (c - eye).rank() == 15
    columns = []
    for _ in range(16):
        columns.append([v[i, 0] for i in range(16)])
        v = h * v
    krylov = nmod_mat(columns, P).transpose()
    assert int(krylov.det()) == 10
    print("GF(23): det[L_b;L_d]=20, rank(C-I)=15, det[v,Dv,...,D^15v]=10")
    print("Simultaneous centralizer dimension 1; three-map intertwiner scalar over Q.")


if __name__ == "__main__":
    run()
