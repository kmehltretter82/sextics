#!/usr/bin/env python3
"""Dependency-free exact Q linear-algebra audit for four-plus-two-residual.

Reproduces every exact rational baseline number in four_plus_two_fixed_divisor_audit.py
using Python's standard-library fractions.Fraction (zero external dependencies),
benchmarks the reduced coprime (N, T) formulation (deg 2d and d+6 vs 2d+6),
and audits rational candidate denominators over Q through moderate degree.
"""
from fractions import Fraction
from time import perf_counter


def exponents(deg: int) -> list[tuple[int, int, int, int]]:
    return [
        (a, b, c, deg - a - b - c)
        for a in range(deg + 1)
        for b in range(deg + 1 - a)
        for c in range(deg + 1 - a - b)
    ]


KAPPA = Fraction(5, 7)
H_POLY: dict[tuple[int, int, int, int], Fraction] = {
    (6, 0, 0, 0): Fraction(1),
    (5, 1, 0, 0): KAPPA,
    (1, 5, 0, 0): 2 * KAPPA,
    (0, 6, 0, 0): Fraction(1),
    (0, 0, 6, 0): Fraction(1),
    (0, 0, 5, 1): KAPPA,
    (0, 0, 1, 5): 3 * KAPPA,
    (0, 0, 0, 6): Fraction(1),
    (3, 1, 1, 1): KAPPA,
    (1, 1, 3, 1): 2 * KAPPA,
}
A_POLY: dict[tuple[int, int, int, int], Fraction] = {(3, 0, 3, 0): Fraction(1)}


def poly_add(p: dict, q: dict, scale_q: Fraction = Fraction(1)) -> dict:
    out = dict(p)
    if scale_q == 0:
        return out
    for m, v in q.items():
        nv = out.get(m, Fraction(0)) + scale_q * v
        if nv == 0:
            out.pop(m, None)
        else:
            out[m] = nv
    return out


def poly_mul(p: dict, q: dict) -> dict:
    out: dict[tuple[int, int, int, int], Fraction] = {}
    for m1, v1 in p.items():
        for m2, v2 in q.items():
            m = (m1[0] + m2[0], m1[1] + m2[1], m1[2] + m2[2], m1[3] + m2[3])
            nv = out.get(m, Fraction(0)) + v1 * v2
            if nv == 0:
                out.pop(m, None)
            else:
                out[m] = nv
    return out


def poly_deriv(p: dict, var: int) -> dict:
    out: dict[tuple[int, int, int, int], Fraction] = {}
    for m, v in p.items():
        if m[var] > 0:
            ml = list(m)
            c = ml[var]
            ml[var] -= 1
            out[tuple(ml)] = out.get(tuple(ml), Fraction(0)) + v * c
    return {k: v for k, v in out.items() if v != 0}


GRAD_H = tuple(poly_deriv(H_POLY, i) for i in range(4))

SOURCE_EXPONENTS = [
    e for e in exponents(6) if (e[0] + e[1]) % 2 == 1 and (e[2] + e[3]) % 2 == 1
]


def sparse_rank_qq(cols: list[dict[int, Fraction]]) -> int:
    """Exact sparse Gaussian elimination over Q on column dictionaries {row_idx: Fraction}."""
    pivots: dict[int, dict[int, Fraction]] = {}
    for col in sorted(cols, key=len):
        cur = dict(col)
        while cur:
            lead = min(cur.keys())
            if lead not in pivots:
                inv = Fraction(1, 1) / cur[lead]
                pivots[lead] = {r: v * inv for r, v in cur.items()}
                break
            fac = cur[lead]
            piv = pivots[lead]
            del cur[lead]
            for r, pv in piv.items():
                if r == lead:
                    continue
                nv = cur.get(r, Fraction(0)) - fac * pv
                if nv == 0:
                    cur.pop(r, None)
                else:
                    cur[r] = nv
    return len(pivots)


def det_qq(mat: list[list[Fraction]]) -> Fraction:
    n = len(mat)
    a = [list(r) for r in mat]
    det = Fraction(1)
    for col in range(n):
        piv = next((r for r in range(col, n) if a[r][col] != 0), None)
        if piv is None:
            return Fraction(0)
        if piv != col:
            a[col], a[piv] = a[piv], a[col]
            det = -det
        det *= a[col][col]
        inv = Fraction(1, 1) / a[col][col]
        for r in range(col + 1, n):
            if a[r][col] != 0:
                fac = a[r][col] * inv
                for j in range(col, n):
                    a[r][j] -= fac * a[col][j]
    return det


def audit_baseline() -> tuple[int, list[dict]]:
    x_vars = [{(1, 0, 0, 0): Fraction(1)}, {(0, 1, 0, 0): Fraction(1)}]
    y_vars = [{(0, 0, 1, 0): Fraction(1)}, {(0, 0, 0, 1): Fraction(1)}]
    sources = [poly_mul(y, GRAD_H[i]) for i in (0, 1) for y in y_vars]
    sources += [poly_mul(x, GRAD_H[j]) for j in (2, 3) for x in x_vars]

    row_idx = {e: i for i, e in enumerate(SOURCE_EXPONENTS)}
    cols = [{row_idx[e]: v for e, v in s.items()} for s in sources]
    col_A = {row_idx[e]: v for e, v in A_POLY.items()}

    rank = sparse_rank_qq(cols)
    rank_with_A = sparse_rank_qq(cols + [col_A])

    outer_set = {i for i, e in enumerate(SOURCE_EXPONENTS) if e[0] + e[1] != 3}
    outer_cols = [{r: v for r, v in c.items() if r in outer_set} for c in cols]
    outer_rank = sparse_rank_qq(outer_cols)

    coeff_A_zero = all((3, 0, 3, 0) not in s for s in sources)
    bideg = {
        f"({b},{6 - b})": sum(e[0] + e[1] == b for e in SOURCE_EXPONENTS)
        for b in (1, 3, 5)
    }
    witness_exponents = [
        (0, 1, 0, 5), (0, 1, 1, 4), (0, 1, 3, 2), (0, 1, 4, 1),
        (0, 5, 0, 1), (0, 5, 1, 0), (1, 0, 0, 5), (1, 0, 1, 4),
    ]
    submat = [
        [sources[c].get(e, Fraction(0)) for c in range(8)]
        for e in witness_exponents
    ]
    wdet = det_qq(submat)

    assert len(SOURCE_EXPONENTS) == 40
    assert bideg == {"(1,5)": 12, "(3,3)": 16, "(5,1)": 12}
    assert rank == 8 and rank_with_A == 9 and outer_rank == 8
    assert coeff_A_zero
    assert wdet == Fraction(-4556250000, 117649)

    print("EXACT_QQ_BASELINE:", {
        "dim_W": 40,
        "bidegrees": bideg,
        "dim_E": rank,
        "dim_W_mod_E": 40 - rank,
        "rank_E_union_A": rank_with_A,
        "outer_projection_rank": outer_rank,
        "A_coeff_zero_on_E": coeff_A_zero,
        "witness_8x8_det": f"{wdet.numerator}/{wdet.denominator}",
    }, flush=True)
    return rank, sources


def audit_fixed_divisor_qq(label: str, D: dict, base_rank: int = 8, use_coprime_T: bool = False) -> dict:
    t0 = perf_counter()
    deg = next(iter(D.keys()))
    d = sum(deg)
    dD = [poly_deriv(D, i) for i in range(4)]

    if not use_coprime_T:
        row_exps = exponents(6 + 2 * d)
        r_idx = {e: i for i, e in enumerate(row_exps)}
        D2 = poly_mul(D, D)
        HD3 = {k: 3 * v for k, v in poly_mul(H_POLY, D).items()}
        cNi = []
        for i in range(4):
            t1 = poly_mul(D, GRAD_H[i])
            t2 = poly_mul(H_POLY, dD[i])
            comb = poly_add({k: -2 * v for k, v in t1.items()}, t2, Fraction(-3))
            cNi.append(comb)

        s_cols = []
        for e in SOURCE_EXPONENTS:
            p = poly_mul({e: Fraction(2)}, D2)
            s_cols.append({r_idx[k]: v for k, v in p.items()})

        n_cols = []
        nbasis = exponents(d + 1)
        for i in range(4):
            for mn in nbasis:
                p = poly_mul({mn: Fraction(1)}, cNi[i])
                if mn[i] > 0:
                    dm = list(mn)
                    ev = dm[i]
                    dm[i] -= 1
                    p = poly_add(p, poly_mul({tuple(dm): Fraction(ev)}, HD3))
                n_cols.append({r_idx[k]: v for k, v in p.items()})

        rk_N = sparse_rank_qq(n_cols)
        rk_full = sparse_rank_qq(s_cols + n_cols)
        dim_WD = 40 - rk_full + rk_N
        res = {
            "label": label,
            "mode": "full_2d+6",
            "deg_D": d,
            "shape": (len(row_exps), 40 + len(n_cols)),
            "rank_full": rk_full,
            "rank_N": rk_N,
            "dim_W_D": dim_WD,
            "new_classes_W_D_mod_E": dim_WD - base_rank,
            "sec": round(perf_counter() - t0, 3),
        }
        print("FIXED_DIVISOR_QQ:", res, flush=True)
        return res
    else:
        # Coprime system:
        # Eq1 (deg 2d):     N . grad D - D * T = 0
        # Eq2 (deg d + 6):  2 * S * D - 2 * N . grad H + 3 * H * (div N - T) = 0
        r1 = exponents(2 * d)
        r2 = exponents(d + 6)
        idx1 = {e: i for i, e in enumerate(r1)}
        off2 = len(r1)
        idx2 = {e: off2 + i for i, e in enumerate(r2)}

        s_cols = []
        for e in SOURCE_EXPONENTS:
            p2 = poly_mul({e: Fraction(2)}, D)
            s_cols.append({idx2[k]: v for k, v in p2.items()})

        nt_cols = []
        nbasis = exponents(d + 1)
        for i in range(4):
            for mn in nbasis:
                col: dict[int, Fraction] = {}
                # Eq1: mn * dD[i]
                for k, v in poly_mul({mn: Fraction(1)}, dD[i]).items():
                    col[idx1[k]] = v
                # Eq2: -2 * mn * GRAD_H[i] + 3 * H * d_i(mn)
                p2 = {k: -2 * v for k, v in poly_mul({mn: Fraction(1)}, GRAD_H[i]).items()}
                if mn[i] > 0:
                    dm = list(mn)
                    ev = dm[i]
                    dm[i] -= 1
                    p2 = poly_add(p2, poly_mul({tuple(dm): Fraction(3 * ev)}, H_POLY))
                for k, v in p2.items():
                    col[idx2[k]] = v
                nt_cols.append(col)

        tbasis = exponents(d)
        for mt in tbasis:
            col = {}
            for k, v in poly_mul({mt: Fraction(-1)}, D).items():
                col[idx1[k]] = v
            for k, v in poly_mul({mt: Fraction(-3)}, H_POLY).items():
                col[idx2[k]] = v
            nt_cols.append(col)

        rk_NT = sparse_rank_qq(nt_cols)
        rk_full = sparse_rank_qq(s_cols + nt_cols)
        dim_WD = 40 - rk_full + rk_NT
        res = {
            "label": label,
            "mode": "coprime_NT_(2d,d+6)",
            "deg_D": d,
            "shape": (len(r1) + len(r2), 40 + len(nt_cols)),
            "rank_full": rk_full,
            "rank_NT": rk_NT,
            "dim_W_D": dim_WD,
            "new_classes_W_D_mod_E": dim_WD - base_rank,
            "sec": round(perf_counter() - t0, 3),
        }
        print("FIXED_DIVISOR_QQ:", res, flush=True)
        return res


def main():
    rank, _ = audit_baseline()
    x0 = {(1, 0, 0, 0): Fraction(1)}
    x0x1 = {(1, 1, 0, 0): Fraction(1)}
    qsum = {
        (2, 0, 0, 0): Fraction(1), (0, 2, 0, 0): Fraction(1),
        (0, 0, 2, 0): Fraction(1), (0, 0, 0, 2): Fraction(1),
    }
    fa = {k: v for k, v in H_POLY.items() if k[2] == 0 and k[3] == 0}
    fc = {k: v for k, v in H_POLY.items() if k[0] == 0 and k[1] == 0}

    for lbl, D in [("x0", x0), ("x0*x1", x0x1), ("x0^2+x1^2+y0^2+y1^2", qsum), ("f_a", fa), ("f_c", fc)]:
        r_full = audit_fixed_divisor_qq(lbl, D, rank, use_coprime_T=False)
        r_nt = audit_fixed_divisor_qq(lbl, D, rank, use_coprime_T=True)
        assert r_full["dim_W_D"] == r_nt["dim_W_D"], f"Mismatch on {lbl}"

    # Also test D = H over Q (not coprime to H, so use full_2d+6)
    audit_fixed_divisor_qq("H (core)", H_POLY, rank, use_coprime_T=False)

    # Exact Q verification of additional rational covariants using the 26x faster coprime_NT system
    a2_2c2 = {(2, 0, 0, 0): Fraction(1), (0, 0, 2, 0): Fraction(2)}
    abcd = {(1, 1, 1, 1): Fraction(1)}
    fa_sub_fc = poly_add(fa, fc, Fraction(-1))
    fa_mul_fc = poly_mul(fa, fc)
    beta = {(1, 1, 3, 1): 2 * KAPPA}
    gamma = {(3, 1, 1, 1): KAPPA}
    h2 = poly_add(poly_mul(beta, beta), poly_mul(fc, gamma), Fraction(-3))
    h1 = poly_add(poly_mul(beta, gamma), poly_mul(fc, fa), Fraction(-9))
    h0 = poly_add(poly_mul(gamma, gamma), poly_mul(beta, fa), Fraction(-3))
    H_cub = poly_add(h2, poly_add(h1, h0))

    for lbl, D in [
        ("a^2+2c^2", a2_2c2),
        ("abcd", abcd),
        ("abcd*(a^2+2c^2)", poly_mul(abcd, a2_2c2)),
        ("f_a - f_c", fa_sub_fc),
        ("f_a * f_c (deg 12)", fa_mul_fc),
    ]:
        res = audit_fixed_divisor_qq(lbl, D, rank, use_coprime_T=True)
        assert res["dim_W_D"] == 8 and res["new_classes_W_D_mod_E"] == 0
    print("ALL EXACT Q BASELINE & BENCHMARK CHECKS PASSED!", flush=True)


if __name__ == "__main__":
    main()
