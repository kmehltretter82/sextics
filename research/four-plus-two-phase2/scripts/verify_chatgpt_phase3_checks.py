#!/usr/bin/env python3
"""Independently reproduce ChatGPT's Phase 3 exact checks:
1. Mod-11 contact enumeration & 84x60 first-order deformation obstruction + witnesses l_a, l_c.
2. Bad prime p = 13 singular point [0:0:1:9].
3. Twisted Frobenius traces (t0, t1, t2) and tau_zeta at p = 19, 31, 37, 43, 61.
"""


def modinv(x: int, p: int) -> int:
    return pow(x % p, p - 2, p)


def rref_rank_mod(rows: list[list[int]], p: int) -> int:
    if not rows:
        return 0
    m, n = len(rows), len(rows[0])
    a = [[v % p for v in r] for r in rows]
    rk = 0
    for col in range(n):
        piv = next((r for r in range(rk, m) if a[r][col] != 0), None)
        if piv is None:
            continue
        a[rk], a[piv] = a[piv], a[rk]
        inv = modinv(a[rk][col], p)
        a[rk] = [(v * inv) % p for v in a[rk]]
        for r in range(m):
            if r != rk and a[r][col] != 0:
                f = a[r][col]
                a[r] = [(a[r][j] - f * a[rk][j]) % p for j in range(n)]
        rk += 1
        if rk == m:
            break
    return rk


def exps4(deg: int) -> list[tuple[int, int, int, int]]:
    return [
        (a, b, c, deg - a - b - c)
        for a in range(deg + 1)
        for b in range(deg + 1 - a)
        for c in range(deg + 1 - a - b)
    ]


def pmul(p1: dict, p2: dict, mod: int) -> dict:
    out: dict[tuple[int, int, int, int], int] = {}
    for e1, v1 in p1.items():
        for e2, v2 in p2.items():
            e = (e1[0] + e2[0], e1[1] + e2[1], e1[2] + e2[2], e1[3] + e2[3])
            nv = (out.get(e, 0) + v1 * v2) % mod
            if nv:
                out[e] = nv
            else:
                out.pop(e, None)
    return out


def check_mod11_obstruction():
    p = 11
    Ca = {(2, 0, 0, 0): 3, (1, 1, 0, 0): 9, (0, 2, 0, 0): 7}
    Ba = {(3, 0, 0, 0): 2, (2, 1, 0, 0): 4, (0, 3, 0, 0): 1}
    Cc = {(0, 0, 2, 0): 5, (0, 0, 1, 1): 2, (0, 0, 0, 2): 8}
    Bc = {(0, 0, 3, 0): 5, (0, 0, 0, 3): 7}

    # Check C_a^3 - B_a^2 == f_a mod 11 and C_c^3 - B_c^2 == f_c mod 11 (kappa = 5/7 = 7 mod 11)
    Ca3 = pmul(pmul(Ca, Ca, p), Ca, p)
    Ba2 = pmul(Ba, Ba, p)
    fa_chk = {e: (Ca3.get(e, 0) - Ba2.get(e, 0)) % p for e in set(Ca3) | set(Ba2)}
    fa_chk = {k: v for k, v in fa_chk.items() if v}
    assert fa_chk == {(6, 0, 0, 0): 1, (5, 1, 0, 0): 7, (1, 5, 0, 0): 3, (0, 6, 0, 0): 1}

    Cc3 = pmul(pmul(Cc, Cc, p), Cc, p)
    Bc2 = pmul(Bc, Bc, p)
    fc_chk = {e: (Cc3.get(e, 0) - Bc2.get(e, 0)) % p for e in set(Cc3) | set(Bc2)}
    fc_chk = {k: v for k, v in fc_chk.items() if v}
    assert fc_chk == {(0, 0, 6, 0): 1, (0, 0, 5, 1): 7, (0, 0, 1, 5): 10, (0, 0, 0, 6): 1}

    basis6 = exps4(6)
    ridx = {e: i for i, e in enumerate(basis6)}
    assert len(basis6) == 84

    Ca2 = pmul(Ca, Ca, p)
    Cc2 = pmul(Cc, Cc, p)
    cols = []
    for base, mult, deg in [(Ca2, 3, 2), (Cc2, 3, 2), (Ba, -2, 3), (Bc, -2, 3)]:
        for m in exps4(deg):
            pol = pmul(base, {m: mult % p}, p)
            col = [0] * 84
            for e, v in pol.items():
                col[ridx[e]] = v
            cols.append(col)
    assert len(cols) == 60
    M = [[cols[c][r] for c in range(60)] for r in range(84)]
    Kcpl_poly = {(3, 1, 1, 1): 7, (1, 1, 3, 1): 3}
    rhs = [Kcpl_poly.get(e, 0) for e in basis6]
    M_aug = [M[r] + [rhs[r]] for r in range(84)]
    rk_M = rref_rank_mod(M, p)
    rk_aug = rref_rank_mod(M_aug, p)
    assert rk_M == 59 and rk_aug == 60

    # Short witnesses l_a = (0, 1, 5, 7, 0) and l_c = (0, 1, 3, 0, 4)
    la = [0, 1, 5, 7, 0]
    q_a = [(4 - i, i, 0, 0) for i in range(5)]
    for pol in [Ca2, pmul({(1, 0, 0, 0): 1}, Ba, p), pmul({(0, 1, 0, 0): 1}, Ba, p)]:
        assert sum(la[i] * pol.get(q_a[i], 0) for i in range(5)) % p == 0
    assert sum(la[i] * Kcpl_poly.get((q_a[i][0], q_a[i][1], 1, 1), 0) for i in range(5)) % p == 7

    lc = [0, 1, 3, 0, 4]
    q_c = [(0, 0, 4 - i, i) for i in range(5)]
    for pol in [Cc2, pmul({(0, 0, 1, 0): 1}, Bc, p), pmul({(0, 0, 0, 1): 1}, Bc, p)]:
        assert sum(lc[i] * pol.get(q_c[i], 0) for i in range(5)) % p == 0
    assert sum(lc[i] * Kcpl_poly.get((1, 1, q_c[i][2], q_c[i][3]), 0) for i in range(5)) % p == 3
    print("VERIFIED_MOD11_OBSTRUCTION: rank(M)=59, rank([M|K_cpl])=60, short witnesses l_a, l_c OK")


def check_bad_prime_13():
    p = 13
    kap = (5 * modinv(7, p)) % p  # 5 * 2 = 10 mod 13
    # At [a:b:c:d] = [0:0:1:9], partials wrt a, b vanish automatically; check partials wrt c, d of f_c(c,d)
    c, d = 1, 9
    fc = (c**6 + kap * c**5 * d + 3 * kap * c * d**5 + d**6) % p
    dfc_dc = (6 * c**5 + 5 * kap * c**4 * d + 3 * kap * d**5) % p
    dfc_dd = (kap * c**5 + 15 * kap * c * d**4 + 6 * d**5) % p
    assert fc == 0 and dfc_dc == 0 and dfc_dd == 0
    print("VERIFIED_BAD_PRIME_13: [0:0:1:9] is a singular point of H mod 13")


def check_frobenius_traces():
    expected = [
        (19, 7, (133, -95, -38), (76, 19)),
        (31, 25, (-155, 124, 31), (-93, -31)),
        (37, 26, (-407, 481, -74), (-296, -185)),
        (43, 36, (215, -430, 215), (215, 215)),
        (61, 47, (0, 183, -183), (-61, -122)),
    ]
    for p, zeta, exp_t, exp_tau in expected:
        assert (p - 1) % 3 == 0
        assert pow(zeta, 3, p) == 1 and zeta != 1
        kap = (5 * modinv(7, p)) % p
        nu = next(x for x in range(2, p) if pow(x, (p - 1) // 2, p) == p - 1)
        e3 = (p - 1) // 3
        chi = [0] * p
        for x in range(1, p):
            pw = pow(x, e3, p)
            if pw == 1:
                chi[x] = 0
            elif pw == zeta:
                chi[x] = 1
            else:
                chi[x] = 2

        # Notice #X_j = #{[a:b:c:d:u] in P(1,1,1,1,2) : u^3 = -lam_j H(a,b,c,d)}.
        # Summing over affine (a,b,c,d) != 0 and dividing by (p - 1):
        # For a fixed nonzero H_val = H(a,b,c,d), the number of u in F_p with u^3 = -lam_j H_val
        # is 1 if H_val == 0, else 3 if (-lam_j H_val)^e3 == 1, else 0.
        # Since (-1)^e3 == 1 (as e3 is even for odd p = 1 mod 6!), (-H_val)^e3 == H_val^e3!
        # And (lam_j H_val)^e3 == zeta^j * H_val^e3 == 1 iff H_val^e3 == zeta^{(-j) mod 3}!
        fa = [0] * (p * p)
        fc = [0] * (p * p)
        for a in range(p):
            a2 = (a * a) % p
            a5 = (a2 * a2 * a) % p
            a6 = (a5 * a) % p
            for b in range(p):
                b2 = (b * b) % p
                b5 = (b2 * b2 * b) % p
                b6 = (b5 * b) % p
                fa[a * p + b] = (a6 + b6 + kap * a5 * b + 2 * kap * a * b5) % p
                fc[a * p + b] = (a6 + b6 + kap * a5 * b + 3 * kap * a * b5) % p

        nu2 = (nu * nu) % p
        nu3 = (nu2 * nu) % p
        cnt = [0, 0, 0]
        cnt_sig = [0, 0, 0]
        for a in range(p):
            a2 = (a * a) % p
            for b in range(p):
                fab = fa[a * p + b]
                ab = (a * b) % p
                kab = (kap * ab) % p
                for c in range(p):
                    c2_2 = (2 * c * c) % p
                    for d in range(p):
                        if a == 0 and b == 0 and c == 0 and d == 0:
                            continue
                        fcd = fc[c * p + d]
                        cd = (c * d) % p
                        c1 = (kab * cd * a2) % p
                        c2 = (kab * cd * c2_2) % p
                        hv = (fab + fcd + c1 + c2) % p
                        hs = (fab + nu3 * fcd + nu * c1 + nu2 * c2) % p
                        if hv == 0:
                            cnt[0] += 1; cnt[1] += 1; cnt[2] += 1
                        else:
                            cnt[(-chi[hv]) % 3] += 3
                        if hs == 0:
                            cnt_sig[0] += 1; cnt_sig[1] += 1; cnt_sig[2] += 1
                        else:
                            cnt_sig[(-chi[hs]) % 3] += 3
        # Divide affine nonzero counts by (p - 1); vertex [0:0:0:0:1] is not on u^3 + lam H = 0
        X = [c // (p - 1) for c in cnt]
        Xs = [c // (p - 1) for c in cnt_sig]
        t = tuple((Xs[j] - X[j]) // 2 for j in range(3))
        tau = ((t[0] - t[1]) // 3, (t[2] - t[1]) // 3)
        assert t == exp_t and tau == exp_tau, f"Mismatch at p={p}: got {t}, {tau}"
        print(f"VERIFIED_FROBENIUS p={p}: t={t}, tau_zeta={tau[0]} + ({tau[1]})*zeta")


if __name__ == "__main__":
    check_mod11_obstruction()
    check_bad_prime_13()
    check_frobenius_traces()
