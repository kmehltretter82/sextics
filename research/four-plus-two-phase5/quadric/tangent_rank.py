"""Check transversality of the coupled norm-contact and projection at the seed."""

from pathlib import Path
import sys
import numpy as np

sys.path.insert(0, str(Path(sys.argv[1]).resolve().parent))
import intrinsic_norm_incidence as n
import segre_test as st

n.N = 2
n.ONE = n.const(1)
data = np.load(sys.argv[1])
points = data["coordinates"][: n.N]
mult = n.mon(4, 0)
Y = [e for e in n.mon(6, 0) if not (e[0] and e[1])]
Q, Q2, AA, BB = [
    dict(zip(basis, points[:, start:end, :].transpose(1, 0, 2)))
    for basis, start, end in [
        (n.QBASE, 0, 5),
        (n.QBASE, 5, 10),
        (n.ABASE, 10, 18),
        (n.BBASE, 18, 26),
    ]
]
Q2[n.AB] = n.ONE
G2 = dict(zip(mult, data["g2"][: n.N].transpose(1, 0, 2)))
J = []
for base, poly, num in [
    (n.QBASE, n.pmul(Q, Q), -3),
    (n.QBASE, G2, -1),
    (n.ABASE, AA, -2),
    (n.BBASE, BB, 2),
]:
    for e in base:
        prod = {n.madd(e, m): (v * num) % n.P for m, v in poly.items()}
        _, rem = n.divide(Q2, prod, mult)
        J.append([rem.get(m, n.zero()) for m in Y])
M = np.array(J).transpose(1, 0, 2, 3)


def encode(C):
    return C[..., 0] + n.P * C[..., 1]


A = encode(M)
M0 = A[:, :, 0]
M1 = A[:, :, 1]
rows, cols = st.rank_pivots(M0)
print("seed differential rank", len(cols), "row pivots", rows, "col pivots", cols)
rc, cc = st.rank_pivots(M0[:, 10:])
print("seed cubic differential rank", len(cc))


def det(mat):
    mat = mat.copy()
    acc = 1
    for j in range(len(mat)):
        i = next((i for i in range(j, len(mat)) if mat[i, j]), None)
        if i is None:
            return 0
        if i != j:
            mat[[i, j]] = mat[[j, i]]
            acc = st.neg(acc)
        pivot = int(mat[j, j])
        acc = st.mul(acc, pivot)
        iv = st.inv(pivot)
        for m in range(j + 1, len(mat)):
            u = st.neg(st.mul(iv, int(mat[m, j])))
            if u:
                for k in range(j + 1, len(mat)):
                    mat[m, k] = st.add(int(mat[m, k]), st.mul(u, int(mat[j, k])))
    return acc


def detjet(m0, m1):
    der = 0
    for j in range(len(m0)):
        m = m0.copy()
        m[:, j] = m1[:, j]
        der = st.add(der, det(m))
    return det(m0), der


assert len(cols) == 24 and len(cc) == 16
for omit in range(26):
    use = [i for i in range(26) if i != omit]
    d, t = detjet(M0[:, use], M1[:, use])
    if d or t:
        assert omit == 0 and d == 0 and st.coords(t) == (5, 19, 0, 0)
        print(
            "omitted coordinate",
            omit,
            "determinant constant and linear coefficients",
            st.coords(d)[:2],
            st.coords(t)[:2],
        )
        break
else:
    raise AssertionError("missing transversality minor")
