#!/usr/bin/env python3
from itertools import product

MOD = 1000000007
def inv(a, mod=MOD): return pow(a % mod, mod - 2, mod)
inv7 = inv(7)
def mons(deg): return [m for m in product(range(deg + 1), repeat=4) if sum(m) == deg]

H = {
    (6,0,0,0): 1,
    (5,1,0,0): (5 * inv7) % MOD,
    (3,1,1,1): (5 * inv7) % MOD,
    (1,5,0,0): (10 * inv7) % MOD,
    (1,1,3,1): (10 * inv7) % MOD,
    (0,6,0,0): 1,
    (0,0,6,0): 1,
    (0,0,5,1): (5 * inv7) % MOD,
    (0,0,1,5): (15 * inv7) % MOD,
    (0,0,0,6): 1,
}
A_mon = (3,0,3,0)

def pderiv(poly, i, mod=MOD):
    out = {}
    for m, c in poly.items():
        if m[i] > 0:
            m2 = list(m); m2[i] -= 1; t = tuple(m2)
            out[t] = (out.get(t, 0) + c * m[i]) % mod
    return {k: v for k, v in out.items() if v}

def poly_mul(p1, p2, mod=MOD):
    out = {}
    for m1, c1 in p1.items():
        for m2, c2 in p2.items():
            mm = (m1[0]+m2[0], m1[1]+m2[1], m1[2]+m2[2], m1[3]+m2[3])
            out[mm] = (out.get(mm, 0) + c1 * c2) % mod
    return {k: v for k, v in out.items() if v}

dH = [pderiv(H, i) for i in range(4)]
m6, m7, m12 = mons(6), mons(7), mons(12)
idx6 = {m: i for i, m in enumerate(m6)}
idx7 = {m: i for i, m in enumerate(m7)}
idx12 = {m: i for i, m in enumerate(m12)}

def rref_nullspace(rows, ncols, mod=MOD):
    pivots = {}
    for r in rows:
        cur = {k: v % mod for k, v in r.items() if v % mod != 0}
        while cur:
            c = min(cur.keys())
            v = cur[c]
            if c in pivots:
                piv = pivots[c]
                for k, pv in piv.items():
                    nv = (cur.get(k, 0) - v * pv) % mod
                    if nv: cur[k] = nv
                    else: cur.pop(k, None)
            else:
                iv = inv(v, mod)
                pivots[c] = {k: (val * iv) % mod for k, val in cur.items()}
                break
    # Keep un-backsubstituted copy or backsubstitute carefully!
    # Wait: in previous script, `higher = [c2 for c2 in row.keys() if c2 > c and c2 in pivots]`
    # needed `sorted(higher)` or iterating in increasing order can introduce NEW pivot cols > c2!
    # Backsubstitute from LARGEST pivot col down, processing higher pivot cols in DECREASING or INCREASING?
    # Since pivots[c2] for c2 > c ONLY has keys >= c2, and pivots[c2] was ALREADY backsubstituted so it has NO pivot keys > c2!
    # Wait: why did mod_j12 fail? Because mod_j12 used piv_j12 AFTER backsubstitution, which is fine, let's check E1!
    for c in sorted(pivots.keys(), reverse=True):
        row = pivots[c]
        for c2 in sorted([k for k in row.keys() if k > c and k in pivots]):
            factor = row.pop(c2)
            # Since pivots[c2] already has NO pivot keys > c2, this CANNOT add any pivot keys!
            for k, pv in pivots[c2].items():
                if k == c2: continue
                nv = (row.get(k, 0) - factor * pv) % mod
                if nv: row[k] = nv
                else: row.pop(k, None)
    free_cols = [j for j in range(ncols) if j not in pivots]
    basis = []
    for f in free_cols:
        vec = [0] * ncols
        vec[f] = 1
        for c, row in pivots.items():
            if f in row:
                vec[c] = (-row[f]) % mod
        basis.append(vec)
    return len(pivots), pivots, basis

ncols = 4 * len(m7) + len(m6)
e1_rows = [{} for _ in range(len(m12))]
for i in range(4):
    base = i * len(m7)
    for j, m in enumerate(m7):
        for mh, ch in dH[i].items():
            prod_m = (m[0]+mh[0], m[1]+mh[1], m[2]+mh[2], m[3]+mh[3])
            r = e1_rows[idx12[prod_m]]
            col = base + j
            r[col] = (r.get(col, 0) + ch) % MOD
base_D = 4 * len(m7)
for j, m in enumerate(m6):
    prod_m = (m[0]+A_mon[0], m[1]+A_mon[1], m[2]+A_mon[2], m[3]+A_mon[3])
    r = e1_rows[idx12[prod_m]]
    col = base_D + j
    r[col] = (r.get(col, 0) - 1) % MOD

rk_e1, piv_e1, null_e1 = rref_nullspace(e1_rows, ncols)
# Check null_e1 actually satisfies e1_rows!
for v in null_e1[:5]:
    for r in e1_rows:
        assert sum(val * v[col] for col, val in r.items()) % MOD == 0
print("E1 nullspace verified!")

koszul_vecs = [v for v in null_e1 if all(v[base_D + j] == 0 for j in range(len(m6)))]
dc_vecs = [v for v in null_e1 if any(v[base_D + j] != 0 for j in range(len(m6)))]

j12_rows = []
for i in range(4):
    for m in m7:
        d = {}
        for mh, ch in dH[i].items():
            mm = (m[0]+mh[0], m[1]+mh[1], m[2]+mh[2], m[3]+mh[3])
            d[idx12[mm]] = (d.get(idx12[mm], 0) + ch) % MOD
        j12_rows.append(d)
rk_j12, piv_j12, _ = rref_nullspace(j12_rows, len(m12))

def mod_j12(poly12):
    cur = {idx12[m]: c % MOD for m, c in poly12.items() if c % MOD != 0}
    for c in sorted(piv_j12.keys()):
        if c in cur:
            v = cur[c]
            for k, pv in piv_j12[c].items():
                nv = (cur.get(k, 0) - v * pv) % MOD
                if nv: cur[k] = nv
                else: cur.pop(k, None)
    return cur

def unpack(v):
    D_p = {m6[j]: v[base_D + j] for j in range(len(m6)) if v[base_D + j]}
    ps = [{m7[j]: v[i * len(m7) + j] for j in range(len(m7)) if v[i * len(m7) + j]} for i in range(4)]
    return D_p, ps

def bracket(D_a, p_b):
    div_p = {}
    for i in range(4):
        for m, c in pderiv(p_b[i], i).items():
            div_p[m] = (div_p.get(m, 0) + c) % MOD
    out = poly_mul(D_a, div_p)
    for i in range(4):
        for m, c in poly_mul(p_b[i], pderiv(D_a, i)).items():
            out[m] = (out.get(m, 0) - c) % MOD
    return {k: v for k, v in out.items() if v}

D0, _ = unpack(dc_vecs[0])
_, K0 = unpack(koszul_vecs[0])
assert mod_j12(bracket(D0, K0)) == {}
print("Verified: bracket(D_c, K) == 0 mod J_12!")

# Check R_12 quadrics by parity!
r12_coords = [j for j in range(len(m12)) if j not in piv_j12]
print("R_12 monomial parities:")
par_r12 = {}
for c in r12_coords:
    m = m12[c]
    p = ((m[0]+m[1])%2, (m[2]+m[3])%2)
    par_r12.setdefault(p, []).append(c)
for p, cs in par_r12.items():
    print(f"  R_12 parity {p}: {len(cs)} monomials")

unpacked_dc = [unpack(v) for v in dc_vecs]
nD = len(unpacked_dc)
sym_pairs = [(a, b) for a in range(nD) for b in range(a, nD)]
q_r12_rows = [{} for _ in range(len(r12_coords))]
for col, (a, b) in enumerate(sym_pairs):
    Da, pa = unpacked_dc[a]
    Db, pb = unpacked_dc[b]
    if a == b:
        poly = bracket(Da, pa)
    else:
        poly = bracket(Da, pb)
        for m, c in bracket(Db, pa).items():
            poly[m] = (poly.get(m, 0) + c) % MOD
    red = mod_j12(poly)
    for r_idx, coord in enumerate(r12_coords):
        if coord in red:
            q_r12_rows[r_idx][col] = red[coord]

rk_qr12, _, _ = rref_nullspace(q_r12_rows, len(sym_pairs))
print(f"Full 49-var R_12 quadrics: 35 eqs, Sym^2={len(sym_pairs)}, rank={rk_qr12}")

# Separate into (0,0) and (1,1) pure sectors AND check degree-3 / degree-4 Macaulay if possible,
# AND check FULL E2 (455 equations in S_12, eliminating K)!
# WAIT: how can we get MORE quadrics on D_c alone without K?!
# Look: for ANY fixed D_c, K enters E2 via L_{D_c}(K) = D_c * div(K) - K(D_c).
# Wait: K = sum_{i<j} phi_{ij} (dH_j e_i - dH_i e_j) where deg(phi_{ij}) = 2!
# What is D_c * div(K) - K(D_c)?
#   = sum_{i<j} [ D_c (d_i phi_{ij} dH_j - d_j phi_{ij} dH_i) - phi_{ij} (d_i D_c dH_j - d_j D_c dH_i) ]!!
# Notice every term is of the form M_k * dH_k where deg(M_k) = 7!
by_parity = {(0,0): [], (1,1): []}
for v in dc_vecs:
    D_p, ps = unpack(v)
    m = next(iter(D_p.keys()))
    by_parity[((m[0]+m[1])%2, (m[2]+m[3])%2)].append((D_p, ps))

for par, unp in by_parity.items():
    nd_p = len(unp)
    sp = [(a, b) for a in range(nd_p) for b in range(a, nd_p)]
    rows_p = [{} for _ in range(len(r12_coords))]
    for col, (a, b) in enumerate(sp):
        Da, pa = unp[a]
        Db, pb = unp[b]
        poly = bracket(Da, pa) if a == b else dict(bracket(Da, pb))
        if a != b:
            for m, c in bracket(Db, pa).items():
                poly[m] = (poly.get(m, 0) + c) % MOD
        red = mod_j12(poly)
        for r_idx, coord in enumerate(r12_coords):
            if coord in red:
                rows_p[r_idx][col] = red[coord]
    rk_p, _, _ = rref_nullspace(rows_p, len(sp))
    print(f"Pure sector {par}: nd={nd_p}, Sym^2={len(sp)}, R_12 quadric rank = {rk_p}")
