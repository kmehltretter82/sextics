#!/usr/bin/env python3
"""Generate Work/Sol_no_low_degree_or_scalar_core_certificate.lean and its local stub."""
from fractions import Fraction
from functools import reduce
from math import gcd
from pathlib import Path

ROOT = Path("/Users/kmehltretter/codex/collaboration/prove2me_research/hodge-computer-2026-11-20")
WORK = ROOT / "Work"
THM_DIR = ROOT / "local_env/Theorems"


def mons(deg):
    out = []
    for a in range(deg + 1):
        for b in range(deg - a + 1):
            for c in range(deg - a - b + 1):
                out.append((a, b, c, deg - a - b - c))
    return out


core7 = {
    (6, 0, 0, 0): 7,
    (5, 1, 0, 0): 5,
    (3, 1, 1, 1): 5,
    (1, 5, 0, 0): 10,
    (1, 1, 3, 1): 10,
    (0, 6, 0, 0): 7,
    (0, 0, 6, 0): 7,
    (0, 0, 5, 1): 5,
    (0, 0, 1, 5): 15,
    (0, 0, 0, 6): 7,
}


def pderiv_poly(P, i):
    res = {}
    for m, v in P.items():
        if m[i] > 0:
            m2 = list(m)
            m2[i] -= 1
            m2 = tuple(m2)
            res[m2] = res.get(m2, 0) + v * m[i]
    return res


def mul_poly(P, Q):
    res = {}
    for m1, v1 in P.items():
        for m2, v2 in Q.items():
            m = (m1[0] + m2[0], m1[1] + m2[1], m1[2] + m2[2], m1[3] + m2[3])
            res[m] = res.get(m, 0) + v1 * v2
    return {k: v for k, v in res.items() if v != 0}


pd_core7 = [pderiv_poly(core7, i) for i in range(4)]

MOD = 1000000007


def find_sparse_dual(row_mons, cols, target_col):
    C = len(cols)
    row_nnz = sorted((sum(1 for c in range(C) if m in cols[c]), r) for r, m in enumerate(row_mons))
    basis_C = {}
    indep_rows = []
    dep_with_target = None
    for _, r in row_nnz:
        m = row_mons[r]
        v = [cols[c].get(m, 0) % MOD for c in range(C)]
        tv = target_col.get(m, 0) % MOD
        for pcol, (bvec, btv) in basis_C.items():
            if v[pcol]:
                fac = v[pcol]
                for j in range(pcol, C):
                    if bvec[j]:
                        v[j] = (v[j] - fac * bvec[j]) % MOD
                tv = (tv - fac * btv) % MOD
        piv = next((j for j in range(C) if v[j] != 0), None)
        if piv is not None:
            inv = pow(v[piv], MOD - 2, MOD)
            v = [(x * inv) % MOD for x in v]
            tv = (tv * inv) % MOD
            basis_C[piv] = (v, tv)
            indep_rows.append((r, piv))
        else:
            if tv != 0:
                dep_with_target = r
                break
    assert dep_with_target is not None
    S_rows = [r for r, _ in indep_rows] + [dep_with_target]
    piv_cols = [p for _, p in indep_rows]
    K = len(indep_rows)
    mat = [
        [Fraction(cols[piv_cols[j]].get(row_mons[S_rows[i]], 0)) for i in range(K)]
        + [Fraction(-cols[piv_cols[j]].get(row_mons[dep_with_target], 0))]
        for j in range(K)
    ]
    for col in range(K):
        piv = next(r for r in range(col, K) if mat[r][col] != 0)
        mat[col], mat[piv] = mat[piv], mat[col]
        inv = 1 / mat[col][col]
        for j in range(col, K + 1):
            mat[col][j] *= inv
        for r in range(K):
            if r != col and mat[r][col] != 0:
                fac = mat[r][col]
                for j in range(col, K + 1):
                    mat[r][j] -= fac * mat[col][j]
    weights = {row_mons[S_rows[i]]: mat[i][K] for i in range(K) if mat[i][K] != 0}
    weights[row_mons[dep_with_target]] = Fraction(1)
    lcm = lambda a, b: a * b // gcd(a, b)
    den = reduce(lcm, [v.denominator for v in weights.values()], 1)
    w_int = {k: int(v * den) for k, v in weights.items()}
    g = reduce(gcd, [abs(v) for v in w_int.values()])
    w_int = {k: v // g for k, v in w_int.items()}
    for c in range(C):
        assert sum(w_int.get(m, 0) * v for m, v in cols[c].items()) == 0
    t_dot = sum(w_int.get(m, 0) * v for m, v in target_col.items())
    assert t_dot != 0
    if t_dot < 0:
        w_int = {k: -v for k, v in w_int.items()}
        t_dot = -t_dot
    return sorted(w_int.items()), t_dot


# 1. Scalar core certificate
m12_11 = [m for m in mons(12) if (m[0] + m[1]) % 2 == 1 and (m[2] + m[3]) % 2 == 1]
cols_scalar = []
for i in range(4):
    for alpha in mons(7):
        P = {}
        for m, v in mul_poly(pd_core7[i], {alpha: 1}).items():
            if (m[0] + m[1]) % 2 == 1 and (m[2] + m[3]) % 2 == 1:
                P[m] = P.get(m, 0) - 5 * v
        if alpha[i] > 0:
            a2 = list(alpha)
            a2[i] -= 1
            a2 = tuple(a2)
            for m, v in mul_poly(core7, {a2: alpha[i]}).items():
                if (m[0] + m[1]) % 2 == 1 and (m[2] + m[3]) % 2 == 1:
                    P[m] = P.get(m, 0) + 3 * v
        P = {k: v for k, v in P.items() if v != 0}
        if P:
            cols_scalar.append(P)

rhs_scalar = {
    m: v
    for m, v in mul_poly({(3, 0, 3, 0): 2}, core7).items()
    if (m[0] + m[1]) % 2 == 1 and (m[2] + m[3]) % 2 == 1
}
w_sc, t_sc = find_sparse_dual(m12_11, cols_scalar, rhs_scalar)

# 2. d = 4 certificates (35 monomials)
m10 = mons(10)
cols_p5 = [mul_poly(pd_core7[i], {a: 1}) for i in range(4) for a in mons(5)]
mons4 = mons(4)
duals_d4 = []
for idx_mu, mu in enumerate(mons4):
    other_D = [mul_poly({(3, 0, 3, 0): -7}, {mu2: 1}) for j, mu2 in enumerate(mons4) if j != idx_mu]
    target_mu = mul_poly({(3, 0, 3, 0): -7}, {mu: 1})
    p0 = (3 + mu[0] + mu[1]) % 2
    p1 = (3 + mu[2] + mu[3]) % 2
    m10_par = [m for m in m10 if (m[0] + m[1]) % 2 == p0 and (m[2] + m[3]) % 2 == p1]
    cols_par = [c for c in cols_p5 + other_D if any((m[0] + m[1]) % 2 == p0 and (m[2] + m[3]) % 2 == p1 for m in c)]
    w_mu, t_mu = find_sparse_dual(m10_par, cols_par, target_mu)
    duals_d4.append((mu, t_mu, w_mu))


# Exact simulation of Lean's mergeAllAux + compileScalar + compileStrict in Python!
def packKey(tag, a, b, c, d):
    return ((((tag * 32 + a) * 32 + b) * 32 + c) * 32 + d)


def mergeTwoAux(n, xs, ys):
    out = []
    i = j = 0
    while n > 0 and i < len(xs) and j < len(ys):
        n -= 1
        k1, w1 = xs[i]
        k2, w2 = ys[j]
        if k1 == k2:
            if w1 + w2 != 0:
                out.append((k1, w1 + w2))
            i += 1
            j += 1
        elif k1 <= k2:
            out.append((k1, w1))
            i += 1
        else:
            out.append((k2, w2))
            j += 1
    return out + xs[i:] + ys[j:]


def mergeAllAux(fuel, L):
    while fuel > 0:
        if not L:
            return []
        if len(L) == 1:
            return L[0]
        out = []
        i = 0
        while i + 1 < len(L):
            out.append(mergeTwoAux(2048, L[i], L[i + 1]))
            i += 2
        if i < len(L):
            out.append(L[i])
        L = out
        fuel -= 1
    res = []
    for x in reversed(L):
        res = mergeTwoAux(2048, x, res)
    return res


def sim_rowScalar(m, w):
    m0, m1, m2, m3 = m
    terms0 = []
    for ea, eb, ec, ed, c0 in [
        (9, 0, 3, 0, 14), (8, 1, 3, 0, 10), (6, 1, 4, 1, 10), (4, 5, 3, 0, 20), (4, 1, 6, 1, 20),
        (3, 6, 3, 0, 14), (3, 0, 9, 0, 14), (3, 0, 8, 1, 10), (3, 0, 4, 5, 30), (3, 0, 3, 6, 14),
    ]:
        terms0.append([(packKey(5, 0, 0, 0, 0), c0 * w)] if (ea, eb, ec, ed) == m else [])
    for tag, ea, eb, ec, ed, c0 in [
        (0, 5, 0, 0, 0, -210), (0, 4, 1, 0, 0, -125), (0, 2, 1, 1, 1, -75), (0, 0, 5, 0, 0, -50), (0, 0, 1, 3, 1, -50),
        (1, 5, 0, 0, 0, -25), (1, 3, 0, 1, 1, -25), (1, 1, 4, 0, 0, -250), (1, 1, 0, 3, 1, -50), (1, 0, 5, 0, 0, -210),
        (2, 3, 1, 0, 1, -25), (2, 1, 1, 2, 1, -150), (2, 0, 0, 5, 0, -210), (2, 0, 0, 4, 1, -125), (2, 0, 0, 0, 5, -75),
        (3, 3, 1, 1, 0, -25), (3, 1, 1, 3, 0, -50), (3, 0, 0, 5, 0, -25), (3, 0, 0, 1, 4, -375), (3, 0, 0, 0, 5, -210),
    ]:
        if ea <= m0 and eb <= m1 and ec <= m2 and ed <= m3:
            terms0.append([(packKey(tag, m0 - ea, m1 - eb, m2 - ec, m3 - ed), c0 * w)])
        else:
            terms0.append([])
    core_list = [
        (6, 0, 0, 0, 21), (5, 1, 0, 0, 15), (3, 1, 1, 1, 15), (1, 5, 0, 0, 30), (1, 1, 3, 1, 30),
        (0, 6, 0, 0, 21), (0, 0, 6, 0, 21), (0, 0, 5, 1, 15), (0, 0, 1, 5, 45), (0, 0, 0, 6, 21),
    ]
    termsD01 = []
    for tag in (0, 1):
        for ea, eb, ec, ed, c0 in core_list:
            if ea <= m0 and eb <= m1 and ec <= m2 and ed <= m3:
                sh = [m0 - ea, m1 - eb, m2 - ec, m3 - ed]
                mult = sh[tag] + 1
                sh[tag] += 1
                termsD01.append([(packKey(tag, *sh), c0 * mult * w)])
            else:
                termsD01.append([])
    termsD23 = []
    for tag in (2, 3):
        for ea, eb, ec, ed, c0 in core_list:
            if ea <= m0 and eb <= m1 and ec <= m2 and ed <= m3:
                sh = [m0 - ea, m1 - eb, m2 - ec, m3 - ed]
                mult = sh[tag] + 1
                sh[tag] += 1
                termsD23.append([(packKey(tag, *sh), c0 * mult * w)])
            else:
                termsD23.append([])
    return mergeTwoAux(2048, mergeAllAux(6, terms0), mergeTwoAux(2048, mergeAllAux(6, termsD01), mergeAllAux(6, termsD23)))


assert mergeAllAux(10, [sim_rowScalar(m, w) for m, w in w_sc]) == [(packKey(5, 0, 0, 0, 0), t_sc)]
print("EXACT LEAN SIMULATION OF compileScalar PASSED!")


def fmt_weights(ws):
    items = [f"(({a}, {b}, {c}, {d}), ({w} : ℤ))" for (a, b, c, d), w in ws]
    lines = []
    for i in range(0, len(items), 4):
        lines.append("    " + ", ".join(items[i : i + 4]))
    return "[\n" + ",\n".join(lines) + "\n  ]"


idx_map = {mu: idx for idx, (mu, _, _) in enumerate(duals_d4)}
tree_lines = []
for a in range(5):
    prefix_a = "if" if a == 0 else "else if"
    tree_lines.append(f"    {prefix_a} ha{a} : s 0 = {a} then")
    for b in range(5 - a):
        prefix_b = "if" if b == 0 else "else if"
        tree_lines.append(f"      {prefix_b} hb{b} : s 1 = {b} then")
        for c in range(5 - a - b):
            d = 4 - a - b - c
            idx = idx_map[(a, b, c, d)]
            prefix_c = "if" if c == 0 else "else if"
            tree_lines.append(f"        {prefix_c} hc{c} : s 2 = {c} then")
            tree_lines.append(f"          have hd : s 3 = {d} := by omega")
            tree_lines.append(f"          rw [eq_toFS s, ha{a}, hb{b}, hc{c}, hd]; exact hz_{idx}")
        tree_lines.append("        else omega")
    tree_lines.append("      else omega")
tree_lines.append("    else omega")
decision_tree = "\n".join(tree_lines)

cert_defs = []
cert_haves = []
for idx, ((a, b, c, d), k_val, ws) in enumerate(duals_d4):
    cert_defs.append(
        f"def wsD4_{idx} : List (Exps × ℤ) :=\n  {fmt_weights(ws)}\n\n"
        f"lemma hb_d4_{idx} : allBound wsD4_{idx} = true := by decide\n"
        f"lemma hk_d4_{idx} : ({k_val} : ℤ) ≠ 0 := by decide\n"
        f"lemma hc_d4_{idx} : compileStrict wsD4_{idx} = [(packKey 4 {a} {b} {c} {d}, ({k_val} : ℤ))] := by decide\n"
    )
    cert_haves.append(
        f"  have hz_{idx} : coeff (toFS ({a}, {b}, {c}, {d})) D = 0 :=\n"
        f"    zero_of_strict_cert P D h_eq {a} {b} {c} {d} ({k_val} : ℤ) wsD4_{idx}\n"
        f"      hb_d4_{idx} hk_d4_{idx} hc_d4_{idx}"
    )

coprime_text = (WORK / "Sol_coprime_coreCertificate_iff_strict.lean").read_text()
irred_and_coprime_head = (
    coprime_text.split("theorem solution")[0]
    .replace(
        "import Mathlib.Algebra.MvPolynomial.Equiv\n",
        "import Mathlib.Algebra.MvPolynomial.Equiv\nimport Mathlib.Algebra.MvPolynomial.Funext\n",
    )
    .replace(
        "set_option maxHeartbeats 2000000\n",
        "set_option maxHeartbeats 2000000\nset_option maxRecDepth 100000\n",
    )
)

lean_file = (
    irred_and_coprime_head
    + """
namespace FP2Obstructions

open MvPolynomial HodgeComputer.FourPlusTwoResidual
open scoped BigOperators

private abbrev R4 := MvPolynomial (Fin 4) ℂ
private abbrev Exps := ℕ × ℕ × ℕ × ℕ

private theorem coprime_strict
    (D : R4) (hndvd : ¬ (core ∣ D))
    (N : Fin 4 → R4) (hcert : CoreCertificate D N) :
    ∃ p : Fin 4 → R4, StrictCoreCertificate D p := by
  rcases hcert with ⟨hD_ne, ⟨k, hDk, hUk⟩, hc⟩
  have hprime : Prime core :=
    UniqueFactorizationMonoid.irreducible_iff_prime.mp FP2CoreIrreducible.core_irreducible_proof
  let M : R4 := D * (∑ i, pderiv i (N i)) - ∑ i, pderiv i D * N i
  have hdiv_H : core ∣ D * (2 * ((∑ i, pderiv i core * N i) - coreSource * D)) := by
    refine ⟨3 * M, ?_⟩
    dsimp [M]; linear_combination -hc
  have hdiv2 : core ∣ 2 * ((∑ i, pderiv i core * N i) - coreSource * D) :=
    (hprime.dvd_or_dvd hdiv_H).resolve_left hndvd
  obtain ⟨R2, hR2⟩ := hdiv2
  let R : R4 := C (1 / 2 : ℂ) * R2
  have hR : (∑ i, pderiv i core * N i) - coreSource * D = core * R := by
    have hhalf : C (1 / 2 : ℂ) * (2 : R4) = 1 := by
      rw [← map_ofNat C 2, ← map_mul]; norm_num
    calc
      (∑ i, pderiv i core * N i) - coreSource * D
          = C (1 / 2 : ℂ) * (2 * ((∑ i, pderiv i core * N i) - coreSource * D)) := by
            linear_combination -(∑ i, pderiv i core * N i - coreSource * D) * hhalf
      _ = C (1 / 2 : ℂ) * (core * R2) := by rw [hR2]
      _ = core * R := by dsimp [R]; ring
  have hLHS_hom : ((∑ i, pderiv i core * N i) - coreSource * D).IsHomogeneous (6 + (k + 1)) := by
    have h1 : (∑ i, pderiv i core * N i).IsHomogeneous (6 + (k + 1)) := by
      apply IsHomogeneous.sum
      intro i _
      exact (show (6 - 1) + (k + 2) = 6 + (k + 1) by omega) ▸
        (FP2CoprimeStrict.core_hom.pderiv (i := i)).mul (hUk i)
    have h2 : (coreSource * D).IsHomogeneous (6 + (k + 1)) :=
      FP2CoprimeStrict.coreSource_hom.mul hDk
    exact h1.sub h2
  have hR_hom : R.IsHomogeneous (k + 1) :=
    FP2CoprimeStrict.isHomogeneous_of_mul_isHomogeneous
      FP2CoprimeStrict.core_hom FP2CoprimeStrict.core_ne_zero (hR ▸ hLHS_hom)
  have h3M : 3 * M = 2 * D * R := by
    have hc' : core * (3 * M - 2 * D * R) = 0 := by
      dsimp [M]; linear_combination hc + 2 * D * hR
    exact sub_eq_zero.mp ((mul_eq_zero.mp hc').resolve_left FP2CoprimeStrict.core_ne_zero)
  let p : Fin 4 → R4 := fun i => N i - C (1 / 6 : ℂ) * R * X i
  refine ⟨p, hD_ne, ⟨k, hDk, ?_⟩, ?_, ?_⟩
  · intro i
    have hRX : (C (1 / 6 : ℂ) * R * X i).IsHomogeneous (k + 2) := by
      simpa using (hR_hom.C_mul (1 / 6 : ℂ)).mul (isHomogeneous_X ℂ i)
    exact (hUk i).sub hRX
  · have hE_H : (∑ i : Fin 4, pderiv i core * X i) = 6 * core := by
      calc
        (∑ i : Fin 4, pderiv i core * X i) = ∑ i : Fin 4, X i * pderiv i core :=
          Finset.sum_congr rfl (fun i _ => by ring)
        _ = (6 : ℕ) • core := FP2CoprimeStrict.core_hom.sum_X_mul_pderiv
        _ = 6 * core := by simp [nsmul_eq_mul]
    have h_pull : (∑ i : Fin 4, pderiv i core * (C (1 / 6 : ℂ) * R * X i)) =
        C (1 / 6 : ℂ) * R * (∑ i : Fin 4, pderiv i core * X i) := by
      rw [Finset.mul_sum]
      exact Finset.sum_congr rfl (fun i _ => by ring)
    have hsum_p : (∑ i : Fin 4, pderiv i core * p i) =
        (∑ i : Fin 4, pderiv i core * N i) -
          C (1 / 6 : ℂ) * R * (∑ i : Fin 4, pderiv i core * X i) := by
      simp only [p, mul_sub, Finset.sum_sub_distrib, h_pull]
    have hsixth : C (1 / 6 : ℂ) * (6 : R4) = 1 := by
      rw [← map_ofNat C 6, ← map_mul]; norm_num
    rw [hsum_p, hE_H]
    linear_combination hR - core * R * hsixth
  · have hE_R : (∑ i : Fin 4, X i * pderiv i R) = ((k + 1 : ℕ) : R4) * R := by
      simpa [nsmul_eq_mul] using hR_hom.sum_X_mul_pderiv
    have hE_D : (∑ i : Fin 4, pderiv i D * X i) = ((k + 1 : ℕ) : R4) * D := by
      calc
        (∑ i : Fin 4, pderiv i D * X i) = ∑ i : Fin 4, X i * pderiv i D :=
          Finset.sum_congr rfl (fun i _ => by ring)
        _ = (k + 1) • D := hDk.sum_X_mul_pderiv
        _ = ((k + 1 : ℕ) : R4) * D := by simp [nsmul_eq_mul]
    have hdiv_RX : (∑ i : Fin 4, pderiv i (C (1 / 6 : ℂ) * R * X i)) =
        C (1 / 6 : ℂ) * ((∑ i : Fin 4, X i * pderiv i R) + 4 * R) := by
      have h1 : ∀ i : Fin 4, pderiv i (C (1 / 6 : ℂ) * R * X i) =
          C (1 / 6 : ℂ) * (X i * pderiv i R + R) := by
        intro i
        simp only [Derivation.leibniz, pderiv_C, pderiv_X_self, smul_eq_mul, mul_one]
        ring
      simp_rw [h1, ← Finset.mul_sum, Finset.sum_add_distrib, Finset.sum_const, Finset.card_fin,
        nsmul_eq_mul]
      norm_num
    have hdir_D : (∑ i : Fin 4, pderiv i D * (C (1 / 6 : ℂ) * R * X i)) =
        C (1 / 6 : ℂ) * R * (∑ i : Fin 4, pderiv i D * X i) := by
      rw [Finset.mul_sum]
      exact Finset.sum_congr rfl (fun i _ => by ring)
    have hexp : D * (∑ i, pderiv i (p i)) - ∑ i, pderiv i D * p i =
        M - D * (∑ i, pderiv i (C (1 / 6 : ℂ) * R * X i)) +
          ∑ i, pderiv i D * (C (1 / 6 : ℂ) * R * X i) := by
      dsimp [M, p]
      simp only [map_sub, mul_sub, Finset.sum_sub_distrib]
      ring
    have htwo_thirds : C (1 / 6 : ℂ) * (4 : R4) * (3 : R4) = 2 := by
      rw [← map_ofNat C 4, ← map_ofNat C 3, ← map_mul, ← map_mul]
      norm_num
      exact map_ofNat C 2
    have h3_ne : (3 : R4) ≠ 0 := by
      rw [← map_ofNat C 3, Ne, C_eq_zero]; norm_num
    apply mul_left_cancel₀ h3_ne
    rw [mul_zero, hexp, hdiv_RX, hdir_D, hE_R, hE_D]
    linear_combination h3M - D * R * htwo_thirds

def toFS (e : Exps) : Fin 4 →₀ ℕ :=
  Finsupp.equivFunOnFinite.symm ![e.1, e.2.1, e.2.2.1, e.2.2.2]

lemma toFS_apply_0 (a b c d : ℕ) : toFS (a, b, c, d) 0 = a := by
  change (Finsupp.equivFunOnFinite (toFS (a, b, c, d))) 0 = a
  unfold toFS; rw [Equiv.apply_symm_apply]; rfl

lemma toFS_apply_1 (a b c d : ℕ) : toFS (a, b, c, d) 1 = b := by
  change (Finsupp.equivFunOnFinite (toFS (a, b, c, d))) 1 = b
  unfold toFS; rw [Equiv.apply_symm_apply]; rfl

lemma toFS_apply_2 (a b c d : ℕ) : toFS (a, b, c, d) 2 = c := by
  change (Finsupp.equivFunOnFinite (toFS (a, b, c, d))) 2 = c
  unfold toFS; rw [Equiv.apply_symm_apply]; rfl

lemma toFS_apply_3 (a b c d : ℕ) : toFS (a, b, c, d) 3 = d := by
  change (Finsupp.equivFunOnFinite (toFS (a, b, c, d))) 3 = d
  unfold toFS; rw [Equiv.apply_symm_apply]; rfl

lemma eq_toFS (s : Fin 4 →₀ ℕ) : s = toFS (s 0, s 1, s 2, s 3) := by
  ext i
  fin_cases i
  · exact (toFS_apply_0 _ _ _ _).symm
  · exact (toFS_apply_1 _ _ _ _).symm
  · exact (toFS_apply_2 _ _ _ _).symm
  · exact (toFS_apply_3 _ _ _ _).symm

def monP (a b c d : ℕ) : R4 := X 0 ^ a * X 1 ^ b * X 2 ^ c * X 3 ^ d

lemma monP_eq_monomial (a b c d : ℕ) : monP a b c d = monomial (toFS (a, b, c, d)) 1 := by
  rw [monomial_eq, map_one, one_mul, Finsupp.prod_fintype _ _ (fun _ => pow_zero _),
    Fin.prod_univ_four, toFS_apply_0, toFS_apply_1, toFS_apply_2, toFS_apply_3]
  rfl

lemma coeff_monP_mul (a b c d m0 m1 m2 m3 : ℕ) (P : R4) :
    coeff (toFS (m0, m1, m2, m3)) (monP a b c d * P) =
      if a ≤ m0 ∧ b ≤ m1 ∧ c ≤ m2 ∧ d ≤ m3 then
        coeff (toFS (m0 - a, m1 - b, m2 - c, m3 - d)) P
      else 0 := by
  rw [monP_eq_monomial]
  split_ifs with h
  · rcases h with ⟨h0, h1, h2, h3⟩
    have hadd : toFS (m0, m1, m2, m3) =
        toFS (a, b, c, d) + toFS (m0 - a, m1 - b, m2 - c, m3 - d) := by
      ext i; simp only [Finsupp.coe_add, Pi.add_apply]
      fin_cases i
      · change toFS (m0, m1, m2, m3) 0 = toFS (a, b, c, d) 0 + toFS (m0 - a, m1 - b, m2 - c, m3 - d) 0
        rw [toFS_apply_0, toFS_apply_0, toFS_apply_0]; omega
      · change toFS (m0, m1, m2, m3) 1 = toFS (a, b, c, d) 1 + toFS (m0 - a, m1 - b, m2 - c, m3 - d) 1
        rw [toFS_apply_1, toFS_apply_1, toFS_apply_1]; omega
      · change toFS (m0, m1, m2, m3) 2 = toFS (a, b, c, d) 2 + toFS (m0 - a, m1 - b, m2 - c, m3 - d) 2
        rw [toFS_apply_2, toFS_apply_2, toFS_apply_2]; omega
      · change toFS (m0, m1, m2, m3) 3 = toFS (a, b, c, d) 3 + toFS (m0 - a, m1 - b, m2 - c, m3 - d) 3
        rw [toFS_apply_3, toFS_apply_3, toFS_apply_3]; omega
    rw [hadd, coeff_monomial_mul, one_mul]
  · rw [coeff_mul]
    apply Finset.sum_eq_zero
    rintro ⟨u, v⟩ huv
    rw [Finset.mem_antidiagonal] at huv
    rw [coeff_monomial]
    split_ifs with hu
    · exfalso
      apply h
      subst hu
      have k0 := congrArg (fun f : Fin 4 →₀ ℕ => f 0) huv
      have k1 := congrArg (fun f : Fin 4 →₀ ℕ => f 1) huv
      have k2 := congrArg (fun f : Fin 4 →₀ ℕ => f 2) huv
      have k3 := congrArg (fun f : Fin 4 →₀ ℕ => f 3) huv
      simp only [Finsupp.coe_add, Pi.add_apply, toFS_apply_0, toFS_apply_1,
        toFS_apply_2, toFS_apply_3] at k0 k1 k2 k3
      exact ⟨by omega, by omega, by omega, by omega⟩
    · simp

lemma coeff_monP_sc (ea eb ec ed m0 m1 m2 m3 : ℕ) (sc : ℂ) :
    coeff (toFS (m0, m1, m2, m3)) (C sc * monP ea eb ec ed) =
      if ea = m0 ∧ eb = m1 ∧ ec = m2 ∧ ed = m3 then sc else 0 := by
  rw [monP_eq_monomial, C_mul_monomial, mul_one, coeff_monomial]
  by_cases h1 : toFS (ea, eb, ec, ed) = toFS (m0, m1, m2, m3)
  · rw [if_pos h1]
    have k0 := congrArg (fun f : Fin 4 →₀ ℕ => f 0) h1
    have k1 := congrArg (fun f : Fin 4 →₀ ℕ => f 1) h1
    have k2 := congrArg (fun f : Fin 4 →₀ ℕ => f 2) h1
    have k3 := congrArg (fun f : Fin 4 →₀ ℕ => f 3) h1
    simp only [toFS_apply_0, toFS_apply_1, toFS_apply_2, toFS_apply_3] at k0 k1 k2 k3
    rw [if_pos ⟨k0, k1, k2, k3⟩]
  · rw [if_neg h1]
    by_cases h2 : ea = m0 ∧ eb = m1 ∧ ec = m2 ∧ ed = m3
    · exfalso; apply h1; rcases h2 with ⟨rfl, rfl, rfl, rfl⟩; rfl
    · rw [if_neg h2]

lemma coeff_pderiv_0 (m0 m1 m2 m3 : ℕ) (P : R4) :
    coeff (toFS (m0, m1, m2, m3)) (pderiv 0 P) =
      ((m0 + 1 : ℕ) : ℂ) * coeff (toFS (m0 + 1, m1, m2, m3)) P := by
  rw [coeff_pderiv]
  have hs : toFS (m0, m1, m2, m3) + Finsupp.single 0 1 = toFS (m0 + 1, m1, m2, m3) := by
    ext i; fin_cases i <;>
      simp [toFS_apply_0, toFS_apply_1, toFS_apply_2, toFS_apply_3]
  rw [hs, toFS_apply_0]; push_cast; ring

lemma coeff_pderiv_1 (m0 m1 m2 m3 : ℕ) (P : R4) :
    coeff (toFS (m0, m1, m2, m3)) (pderiv 1 P) =
      ((m1 + 1 : ℕ) : ℂ) * coeff (toFS (m0, m1 + 1, m2, m3)) P := by
  rw [coeff_pderiv]
  have hs : toFS (m0, m1, m2, m3) + Finsupp.single 1 1 = toFS (m0, m1 + 1, m2, m3) := by
    ext i; fin_cases i <;>
      simp [toFS_apply_0, toFS_apply_1, toFS_apply_2, toFS_apply_3]
  rw [hs, toFS_apply_1]; push_cast; ring

lemma coeff_pderiv_2 (m0 m1 m2 m3 : ℕ) (P : R4) :
    coeff (toFS (m0, m1, m2, m3)) (pderiv 2 P) =
      ((m2 + 1 : ℕ) : ℂ) * coeff (toFS (m0, m1, m2 + 1, m3)) P := by
  rw [coeff_pderiv]
  have hs : toFS (m0, m1, m2, m3) + Finsupp.single 2 1 = toFS (m0, m1, m2 + 1, m3) := by
    ext i; fin_cases i <;>
      simp [toFS_apply_0, toFS_apply_1, toFS_apply_2, toFS_apply_3]
  rw [hs, toFS_apply_2]; push_cast; ring

lemma coeff_pderiv_3 (m0 m1 m2 m3 : ℕ) (P : R4) :
    coeff (toFS (m0, m1, m2, m3)) (pderiv 3 P) =
      ((m3 + 1 : ℕ) : ℂ) * coeff (toFS (m0, m1, m2, m3 + 1)) P := by
  rw [coeff_pderiv]
  have hs : toFS (m0, m1, m2, m3) + Finsupp.single 3 1 = toFS (m0, m1, m2, m3 + 1) := by
    ext i; fin_cases i <;>
      simp [toFS_apply_0, toFS_apply_1, toFS_apply_2, toFS_apply_3]
  rw [hs, toFS_apply_3]; push_cast; ring

def packKey (tag a b c d : ℕ) : ℕ :=
  ((((tag * 32 + a) * 32 + b) * 32 + c) * 32 + d)

def valFn (N : Fin 4 → R4) (D : R4) (sc : ℂ) (k : ℕ) : ℂ :=
  let d := k % 32
  let k1 := k / 32
  let c := k1 % 32
  let k2 := k1 / 32
  let b := k2 % 32
  let k3 := k2 / 32
  let a := k3 % 32
  let tag := k3 / 32
  match tag with
  | 0 => coeff (toFS (a, b, c, d)) (N 0)
  | 1 => coeff (toFS (a, b, c, d)) (N 1)
  | 2 => coeff (toFS (a, b, c, d)) (N 2)
  | 3 => coeff (toFS (a, b, c, d)) (N 3)
  | 4 => coeff (toFS (a, b, c, d)) D
  | _ => sc

lemma valFn_pack (N : Fin 4 → R4) (D : R4) (sc : ℂ) (tag a b c d : ℕ)
    (ha : a < 32) (hb : b < 32) (hc : c < 32) (hd : d < 32) :
    valFn N D sc (packKey tag a b c d) =
      match tag with
      | 0 => coeff (toFS (a, b, c, d)) (N 0)
      | 1 => coeff (toFS (a, b, c, d)) (N 1)
      | 2 => coeff (toFS (a, b, c, d)) (N 2)
      | 3 => coeff (toFS (a, b, c, d)) (N 3)
      | 4 => coeff (toFS (a, b, c, d)) D
      | _ => sc := by
  unfold valFn packKey
  have h1 : ((((tag * 32 + a) * 32 + b) * 32 + c) * 32 + d) % 32 = d := by omega
  have h2 : ((((tag * 32 + a) * 32 + b) * 32 + c) * 32 + d) / 32 =
      ((tag * 32 + a) * 32 + b) * 32 + c := by omega
  have h3 : (((tag * 32 + a) * 32 + b) * 32 + c) % 32 = c := by omega
  have h4 : (((tag * 32 + a) * 32 + b) * 32 + c) / 32 = (tag * 32 + a) * 32 + b := by omega
  have h5 : ((tag * 32 + a) * 32 + b) % 32 = b := by omega
  have h6 : ((tag * 32 + a) * 32 + b) / 32 = tag * 32 + a := by omega
  have h7 : (tag * 32 + a) % 32 = a := by omega
  have h8 : (tag * 32 + a) / 32 = tag := by omega
  dsimp only
  rw [h1, h2, h3, h4, h5, h6, h7, h8]

def evalPairs (val : ℕ → ℂ) : List (ℕ × ℤ) → ℂ
  | [] => 0
  | (k, w) :: xs => (w : ℂ) * val k + evalPairs val xs

lemma evalPairs_append (val : ℕ → ℂ) :
    ∀ (xs ys : List (ℕ × ℤ)), evalPairs val (xs ++ ys) = evalPairs val xs + evalPairs val ys
  | [], ys => by simp [evalPairs]
  | (k, w) :: xs, ys => by simp [evalPairs, evalPairs_append val xs ys, add_assoc]

def mergeTwoAux : ℕ → List (ℕ × ℤ) → List (ℕ × ℤ) → List (ℕ × ℤ)
  | 0, xs, ys => xs ++ ys
  | _ + 1, [], ys => ys
  | _ + 1, xs, [] => xs
  | n + 1, (k1, w1) :: xs, (k2, w2) :: ys =>
      if k1 == k2 then
        let w := w1 + w2
        if w == 0 then mergeTwoAux n xs ys else (k1, w) :: mergeTwoAux n xs ys
      else if Nat.ble k1 k2 then
        (k1, w1) :: mergeTwoAux n xs ((k2, w2) :: ys)
      else
        (k2, w2) :: mergeTwoAux n ((k1, w1) :: xs) ys

def mergeTwo (bound : ℕ) (xs ys : List (ℕ × ℤ)) : List (ℕ × ℤ) :=
  match xs, ys with
  | [], ys => ys
  | xs, [] => xs
  | xs, ys => mergeTwoAux bound xs ys

lemma evalPairs_mergeTwoAux (val : ℕ → ℂ) :
    ∀ (n : ℕ) (xs ys : List (ℕ × ℤ)),
      evalPairs val (mergeTwoAux n xs ys) = evalPairs val xs + evalPairs val ys := by
  intro n
  induction n with
  | zero => intro xs ys; exact evalPairs_append val xs ys
  | succ n ih =>
      intro xs ys
      cases xs with
      | nil => simp [mergeTwoAux, evalPairs]
      | cons x xs' =>
          cases ys with
          | nil => simp [mergeTwoAux, evalPairs]
          | cons y ys' =>
              rcases x with ⟨k1, w1⟩
              rcases y with ⟨k2, w2⟩
              simp only [mergeTwoAux]
              by_cases h12 : (k1 == k2) = true
              · rw [if_pos h12]
                have hk : k1 = k2 := beq_iff_eq.mp h12
                subst hk
                by_cases hw : (w1 + w2 == 0) = true
                · rw [if_pos hw, ih xs' ys']
                  have hw0 : w1 + w2 = 0 := beq_iff_eq.mp hw
                  have hcast : (w1 : ℂ) + (w2 : ℂ) = 0 := by exact_mod_cast hw0
                  simp only [evalPairs]
                  linear_combination -val k1 * hcast
                · rw [if_neg hw, evalPairs, ih xs' ys']
                  simp only [evalPairs, Int.cast_add]
                  ring
              · rw [if_neg h12]
                by_cases hle : Nat.ble k1 k2 = true
                · rw [if_pos hle, evalPairs, ih xs' ((k2, w2) :: ys')]
                  simp only [evalPairs]; ring
                · rw [if_neg hle, evalPairs, ih ((k1, w1) :: xs') ys']
                  simp only [evalPairs]; ring

lemma evalPairs_mergeTwo (val : ℕ → ℂ) (bound : ℕ) (xs ys : List (ℕ × ℤ)) :
    evalPairs val (mergeTwo bound xs ys) = evalPairs val xs + evalPairs val ys := by
  cases xs with
  | nil => simp [mergeTwo, evalPairs]
  | cons x xs' =>
      cases ys with
      | nil => simp [mergeTwo, evalPairs]
      | cons y ys' => exact evalPairs_mergeTwoAux val bound _ _

def mergeStep (bound : ℕ) : List (List (ℕ × ℤ)) → List (List (ℕ × ℤ))
  | [] => []
  | [x] => [x]
  | x :: y :: rest => mergeTwo bound x y :: mergeStep bound rest

def evalListSum (val : ℕ → ℂ) : List (List (ℕ × ℤ)) → ℂ
  | [] => 0
  | x :: rest => evalPairs val x + evalListSum val rest

lemma evalListSum_mergeStep (val : ℕ → ℂ) (bound : ℕ) :
    ∀ L : List (List (ℕ × ℤ)), evalListSum val (mergeStep bound L) = evalListSum val L
  | [] => rfl
  | [x] => rfl
  | x :: y :: rest => by
      simp only [mergeStep, evalListSum, evalPairs_mergeTwo, evalListSum_mergeStep val bound rest]
      ring

def mergeFold : List (List (ℕ × ℤ)) → List (ℕ × ℤ)
  | [] => []
  | x :: rest => mergeTwo 2048 x (mergeFold rest)

lemma evalPairs_mergeFold (val : ℕ → ℂ) :
    ∀ L : List (List (ℕ × ℤ)), evalPairs val (mergeFold L) = evalListSum val L
  | [] => rfl
  | x :: rest => by
      simp only [mergeFold, evalPairs_mergeTwo, evalPairs_mergeFold val rest, evalListSum]

def mergeAllAux : ℕ → ℕ → List (List (ℕ × ℤ)) → List (ℕ × ℤ)
  | 0, _, L => mergeFold L
  | _ + 1, _, [] => []
  | _ + 1, _, [x] => x
  | fuel + 1, bound, L => mergeAllAux fuel (bound * 2) (mergeStep bound L)

lemma evalPairs_mergeAllAux (val : ℕ → ℂ) :
    ∀ (fuel bound : ℕ) (L : List (List (ℕ × ℤ))),
      evalPairs val (mergeAllAux fuel bound L) = evalListSum val L
  | 0, _, L => evalPairs_mergeFold val L
  | _ + 1, _, [] => rfl
  | _ + 1, _, [x] => by simp [mergeAllAux, evalListSum]
  | fuel + 1, bound, x :: y :: rest => by
      simp only [mergeAllAux,
        evalPairs_mergeAllAux val fuel (bound * 2) (mergeStep bound (x :: y :: rest)),
        evalListSum_mergeStep]

-- Explicit 7 * core and 7 * pderiv i core decompositions
lemma core7_eq : C 7 * core =
    C 7 * monP 6 0 0 0 + C 5 * monP 5 1 0 0 + C 5 * monP 3 1 1 1 + C 10 * monP 1 5 0 0 +
      C 10 * monP 1 1 3 1 + C 7 * monP 0 6 0 0 + C 7 * monP 0 0 6 0 + C 5 * monP 0 0 5 1 +
      C 15 * monP 0 0 1 5 + C 7 * monP 0 0 0 6 := by
  apply MvPolynomial.funext
  intro x
  simp only [core, monP, eval_add, eval_mul, eval_pow, eval_C, eval_X]
  ring

lemma pderiv_X_pow_succ (i j : Fin 4) : ∀ n : ℕ,
    pderiv i (X j ^ (n + 1) : R4) = if j = i then ((n + 1 : ℕ) : R4) * X j ^ n else 0
  | 0 => by
      simp only [zero_add, pow_one, pow_zero, Nat.cast_one, one_mul, pderiv_X]
      by_cases hji : j = i
      · subst hji; simp
      · simp [hji]
  | n + 1 => by
      rw [pow_succ _ (n + 1), pderiv_mul, pderiv_X_pow_succ i j n, pderiv_X]
      by_cases hji : j = i
      · subst hji; simp; ring
      · simp [hji]

lemma pd0_core7_eq : C 7 * pderiv 0 core =
    C 42 * monP 5 0 0 0 + C 25 * monP 4 1 0 0 + C 15 * monP 2 1 1 1 +
      C 10 * monP 0 5 0 0 + C 10 * monP 0 1 3 1 := by
  have h := congrArg (pderiv 0) core7_eq
  rw [pderiv_mul, pderiv_C, zero_mul, zero_add] at h
  rw [h]
  unfold monP
  simp only [pow_zero, mul_one, one_mul, pow_one,
    show ∀ j : Fin 4, (X j : R4) ^ 6 = X j ^ (5 + 1) from fun _ => rfl,
    show ∀ j : Fin 4, (X j : R4) ^ 5 = X j ^ (4 + 1) from fun _ => rfl,
    show ∀ j : Fin 4, (X j : R4) ^ 3 = X j ^ (2 + 1) from fun _ => rfl,
    map_add, pderiv_mul, pderiv_C, pderiv_X_pow_succ, pderiv_X]
  simp
  apply MvPolynomial.funext; intro x
  simp only [map_add, map_mul, map_pow, map_ofNat, eval_X]
  ring

lemma pd1_core7_eq : C 7 * pderiv 1 core =
    C 5 * monP 5 0 0 0 + C 5 * monP 3 0 1 1 + C 50 * monP 1 4 0 0 +
      C 10 * monP 1 0 3 1 + C 42 * monP 0 5 0 0 := by
  have h := congrArg (pderiv 1) core7_eq
  rw [pderiv_mul, pderiv_C, zero_mul, zero_add] at h
  rw [h]
  unfold monP
  simp only [pow_zero, mul_one, one_mul, pow_one,
    show ∀ j : Fin 4, (X j : R4) ^ 6 = X j ^ (5 + 1) from fun _ => rfl,
    show ∀ j : Fin 4, (X j : R4) ^ 5 = X j ^ (4 + 1) from fun _ => rfl,
    show ∀ j : Fin 4, (X j : R4) ^ 3 = X j ^ (2 + 1) from fun _ => rfl,
    map_add, pderiv_mul, pderiv_C, pderiv_X_pow_succ, pderiv_X]
  simp
  apply MvPolynomial.funext; intro x
  simp only [map_add, map_mul, map_pow, map_ofNat, eval_X]
  ring

lemma pd2_core7_eq : C 7 * pderiv 2 core =
    C 5 * monP 3 1 0 1 + C 30 * monP 1 1 2 1 + C 42 * monP 0 0 5 0 +
      C 25 * monP 0 0 4 1 + C 15 * monP 0 0 0 5 := by
  have h := congrArg (pderiv 2) core7_eq
  rw [pderiv_mul, pderiv_C, zero_mul, zero_add] at h
  rw [h]
  unfold monP
  simp only [pow_zero, mul_one, one_mul, pow_one,
    show ∀ j : Fin 4, (X j : R4) ^ 6 = X j ^ (5 + 1) from fun _ => rfl,
    show ∀ j : Fin 4, (X j : R4) ^ 5 = X j ^ (4 + 1) from fun _ => rfl,
    show ∀ j : Fin 4, (X j : R4) ^ 3 = X j ^ (2 + 1) from fun _ => rfl,
    map_add, pderiv_mul, pderiv_C, pderiv_X_pow_succ, pderiv_X]
  simp
  apply MvPolynomial.funext; intro x
  simp only [map_add, map_mul, map_pow, map_ofNat, eval_X]
  ring

lemma pd3_core7_eq : C 7 * pderiv 3 core =
    C 5 * monP 3 1 1 0 + C 10 * monP 1 1 3 0 + C 5 * monP 0 0 5 0 +
      C 75 * monP 0 0 1 4 + C 42 * monP 0 0 0 5 := by
  have h := congrArg (pderiv 3) core7_eq
  rw [pderiv_mul, pderiv_C, zero_mul, zero_add] at h
  rw [h]
  unfold monP
  simp only [pow_zero, mul_one, one_mul, pow_one,
    show ∀ j : Fin 4, (X j : R4) ^ 6 = X j ^ (5 + 1) from fun _ => rfl,
    show ∀ j : Fin 4, (X j : R4) ^ 5 = X j ^ (4 + 1) from fun _ => rfl,
    show ∀ j : Fin 4, (X j : R4) ^ 3 = X j ^ (2 + 1) from fun _ => rfl,
    map_add, pderiv_mul, pderiv_C, pderiv_X_pow_succ, pderiv_X]
  simp
  apply MvPolynomial.funext; intro x
  simp only [map_add, map_mul, map_pow, map_ofNat, eval_X]
  ring

lemma source_core7_eq (sc : ℂ) : C 7 * (2 * C sc * coreSource * core) =
    C sc * (C 14 * monP 9 0 3 0 + C 10 * monP 8 1 3 0 + C 10 * monP 6 1 4 1 +
      C 20 * monP 4 5 3 0 + C 20 * monP 4 1 6 1 + C 14 * monP 3 6 3 0 +
      C 14 * monP 3 0 9 0 + C 10 * monP 3 0 8 1 + C 30 * monP 3 0 4 5 +
      C 14 * monP 3 0 3 6) := by
  have h : C 7 * (2 * C sc * coreSource * core) = 2 * C sc * coreSource * (C 7 * core) := by ring
  rw [h, core7_eq]
  apply MvPolynomial.funext; intro x
  simp only [coreSource, monP, eval_add, eval_mul, eval_pow, eval_C, eval_X, eval_ofNat]
  ring

-- Shift builders for strict and scalar-core rows (no variable shadowing!)
def shMon (tag ea eb ec ed : ℕ) (c0 : ℤ) (m : Exps) (w : ℤ) : List (ℕ × ℤ) :=
  let (m0, m1, m2, m3) := m
  if decide (ea ≤ m0 ∧ eb ≤ m1 ∧ ec ≤ m2 ∧ ed ≤ m3) then
    [(packKey tag (m0 - ea) (m1 - eb) (m2 - ec) (m3 - ed), c0 * w)]
  else []

lemma eval_shMon_v0 (ea eb ec ed : ℕ) (c0 : ℤ) (m0 m1 m2 m3 : ℕ) (w : ℤ)
    (N : Fin 4 → R4) (D : R4) (sc : ℂ)
    (hm : m0 < 31 ∧ m1 < 31 ∧ m2 < 31 ∧ m3 < 31) :
    evalPairs (valFn N D sc) (shMon 0 ea eb ec ed c0 (m0, m1, m2, m3) w) =
      (w : ℂ) * ((c0 : ℂ) * coeff (toFS (m0, m1, m2, m3)) (monP ea eb ec ed * N 0)) := by
  unfold shMon; dsimp only
  rw [coeff_monP_mul]
  by_cases h : ea ≤ m0 ∧ eb ≤ m1 ∧ ec ≤ m2 ∧ ed ≤ m3
  · rw [if_pos (decide_eq_true h), if_pos h, evalPairs, evalPairs, add_zero,
      valFn_pack N D sc 0 _ _ _ _ (by omega) (by omega) (by omega) (by omega)]
    push_cast; ring
  · rw [if_neg (by simpa using h), if_neg h, evalPairs]; ring

lemma eval_shMon_v1 (ea eb ec ed : ℕ) (c0 : ℤ) (m0 m1 m2 m3 : ℕ) (w : ℤ)
    (N : Fin 4 → R4) (D : R4) (sc : ℂ)
    (hm : m0 < 31 ∧ m1 < 31 ∧ m2 < 31 ∧ m3 < 31) :
    evalPairs (valFn N D sc) (shMon 1 ea eb ec ed c0 (m0, m1, m2, m3) w) =
      (w : ℂ) * ((c0 : ℂ) * coeff (toFS (m0, m1, m2, m3)) (monP ea eb ec ed * N 1)) := by
  unfold shMon; dsimp only
  rw [coeff_monP_mul]
  by_cases h : ea ≤ m0 ∧ eb ≤ m1 ∧ ec ≤ m2 ∧ ed ≤ m3
  · rw [if_pos (decide_eq_true h), if_pos h, evalPairs, evalPairs, add_zero,
      valFn_pack N D sc 1 _ _ _ _ (by omega) (by omega) (by omega) (by omega)]
    push_cast; ring
  · rw [if_neg (by simpa using h), if_neg h, evalPairs]; ring

lemma eval_shMon_v2 (ea eb ec ed : ℕ) (c0 : ℤ) (m0 m1 m2 m3 : ℕ) (w : ℤ)
    (N : Fin 4 → R4) (D : R4) (sc : ℂ)
    (hm : m0 < 31 ∧ m1 < 31 ∧ m2 < 31 ∧ m3 < 31) :
    evalPairs (valFn N D sc) (shMon 2 ea eb ec ed c0 (m0, m1, m2, m3) w) =
      (w : ℂ) * ((c0 : ℂ) * coeff (toFS (m0, m1, m2, m3)) (monP ea eb ec ed * N 2)) := by
  unfold shMon; dsimp only
  rw [coeff_monP_mul]
  by_cases h : ea ≤ m0 ∧ eb ≤ m1 ∧ ec ≤ m2 ∧ ed ≤ m3
  · rw [if_pos (decide_eq_true h), if_pos h, evalPairs, evalPairs, add_zero,
      valFn_pack N D sc 2 _ _ _ _ (by omega) (by omega) (by omega) (by omega)]
    push_cast; ring
  · rw [if_neg (by simpa using h), if_neg h, evalPairs]; ring

lemma eval_shMon_v3 (ea eb ec ed : ℕ) (c0 : ℤ) (m0 m1 m2 m3 : ℕ) (w : ℤ)
    (N : Fin 4 → R4) (D : R4) (sc : ℂ)
    (hm : m0 < 31 ∧ m1 < 31 ∧ m2 < 31 ∧ m3 < 31) :
    evalPairs (valFn N D sc) (shMon 3 ea eb ec ed c0 (m0, m1, m2, m3) w) =
      (w : ℂ) * ((c0 : ℂ) * coeff (toFS (m0, m1, m2, m3)) (monP ea eb ec ed * N 3)) := by
  unfold shMon; dsimp only
  rw [coeff_monP_mul]
  by_cases h : ea ≤ m0 ∧ eb ≤ m1 ∧ ec ≤ m2 ∧ ed ≤ m3
  · rw [if_pos (decide_eq_true h), if_pos h, evalPairs, evalPairs, add_zero,
      valFn_pack N D sc 3 _ _ _ _ (by omega) (by omega) (by omega) (by omega)]
    push_cast; ring
  · rw [if_neg (by simpa using h), if_neg h, evalPairs]; ring

lemma eval_shMon_d (ea eb ec ed : ℕ) (c0 : ℤ) (m0 m1 m2 m3 : ℕ) (w : ℤ)
    (N : Fin 4 → R4) (D : R4) (sc : ℂ)
    (hm : m0 < 31 ∧ m1 < 31 ∧ m2 < 31 ∧ m3 < 31) :
    evalPairs (valFn N D sc) (shMon 4 ea eb ec ed c0 (m0, m1, m2, m3) w) =
      (w : ℂ) * ((c0 : ℂ) * coeff (toFS (m0, m1, m2, m3)) (monP ea eb ec ed * D)) := by
  unfold shMon; dsimp only
  rw [coeff_monP_mul]
  by_cases h : ea ≤ m0 ∧ eb ≤ m1 ∧ ec ≤ m2 ∧ ed ≤ m3
  · rw [if_pos (decide_eq_true h), if_pos h, evalPairs, evalPairs, add_zero,
      valFn_pack N D sc 4 _ _ _ _ (by omega) (by omega) (by omega) (by omega)]
    push_cast; ring
  · rw [if_neg (by simpa using h), if_neg h, evalPairs]; ring

def shDer0 (ea eb ec ed : ℕ) (c0 : ℤ) (m : Exps) (w : ℤ) : List (ℕ × ℤ) :=
  let (m0, m1, m2, m3) := m
  if decide (ea ≤ m0 ∧ eb ≤ m1 ∧ ec ≤ m2 ∧ ed ≤ m3) then
    [(packKey 0 (m0 - ea + 1) (m1 - eb) (m2 - ec) (m3 - ed), c0 * ((m0 - ea + 1 : ℕ) : ℤ) * w)]
  else []

def shDer1 (ea eb ec ed : ℕ) (c0 : ℤ) (m : Exps) (w : ℤ) : List (ℕ × ℤ) :=
  let (m0, m1, m2, m3) := m
  if decide (ea ≤ m0 ∧ eb ≤ m1 ∧ ec ≤ m2 ∧ ed ≤ m3) then
    [(packKey 1 (m0 - ea) (m1 - eb + 1) (m2 - ec) (m3 - ed), c0 * ((m1 - eb + 1 : ℕ) : ℤ) * w)]
  else []

def shDer2 (ea eb ec ed : ℕ) (c0 : ℤ) (m : Exps) (w : ℤ) : List (ℕ × ℤ) :=
  let (m0, m1, m2, m3) := m
  if decide (ea ≤ m0 ∧ eb ≤ m1 ∧ ec ≤ m2 ∧ ed ≤ m3) then
    [(packKey 2 (m0 - ea) (m1 - eb) (m2 - ec + 1) (m3 - ed), c0 * ((m2 - ec + 1 : ℕ) : ℤ) * w)]
  else []

def shDer3 (ea eb ec ed : ℕ) (c0 : ℤ) (m : Exps) (w : ℤ) : List (ℕ × ℤ) :=
  let (m0, m1, m2, m3) := m
  if decide (ea ≤ m0 ∧ eb ≤ m1 ∧ ec ≤ m2 ∧ ed ≤ m3) then
    [(packKey 3 (m0 - ea) (m1 - eb) (m2 - ec) (m3 - ed + 1), c0 * ((m3 - ed + 1 : ℕ) : ℤ) * w)]
  else []

lemma eval_shDer0 (ea eb ec ed : ℕ) (c0 : ℤ) (m0 m1 m2 m3 : ℕ) (w : ℤ)
    (N : Fin 4 → R4) (D : R4) (sc : ℂ)
    (hm : m0 < 31 ∧ m1 < 31 ∧ m2 < 31 ∧ m3 < 31) :
    evalPairs (valFn N D sc) (shDer0 ea eb ec ed c0 (m0, m1, m2, m3) w) =
      (w : ℂ) * ((c0 : ℂ) * coeff (toFS (m0, m1, m2, m3)) (monP ea eb ec ed * pderiv 0 (N 0))) := by
  unfold shDer0; dsimp only
  rw [coeff_monP_mul]
  by_cases h : ea ≤ m0 ∧ eb ≤ m1 ∧ ec ≤ m2 ∧ ed ≤ m3
  · rw [if_pos (decide_eq_true h), if_pos h, coeff_pderiv_0, evalPairs, evalPairs, add_zero,
      valFn_pack N D sc 0 _ _ _ _ (by omega) (by omega) (by omega) (by omega)]
    push_cast; ring
  · rw [if_neg (by simpa using h), if_neg h, evalPairs]; ring

lemma eval_shDer1 (ea eb ec ed : ℕ) (c0 : ℤ) (m0 m1 m2 m3 : ℕ) (w : ℤ)
    (N : Fin 4 → R4) (D : R4) (sc : ℂ)
    (hm : m0 < 31 ∧ m1 < 31 ∧ m2 < 31 ∧ m3 < 31) :
    evalPairs (valFn N D sc) (shDer1 ea eb ec ed c0 (m0, m1, m2, m3) w) =
      (w : ℂ) * ((c0 : ℂ) * coeff (toFS (m0, m1, m2, m3)) (monP ea eb ec ed * pderiv 1 (N 1))) := by
  unfold shDer1; dsimp only
  rw [coeff_monP_mul]
  by_cases h : ea ≤ m0 ∧ eb ≤ m1 ∧ ec ≤ m2 ∧ ed ≤ m3
  · rw [if_pos (decide_eq_true h), if_pos h, coeff_pderiv_1, evalPairs, evalPairs, add_zero,
      valFn_pack N D sc 1 _ _ _ _ (by omega) (by omega) (by omega) (by omega)]
    push_cast; ring
  · rw [if_neg (by simpa using h), if_neg h, evalPairs]; ring

lemma eval_shDer2 (ea eb ec ed : ℕ) (c0 : ℤ) (m0 m1 m2 m3 : ℕ) (w : ℤ)
    (N : Fin 4 → R4) (D : R4) (sc : ℂ)
    (hm : m0 < 31 ∧ m1 < 31 ∧ m2 < 31 ∧ m3 < 31) :
    evalPairs (valFn N D sc) (shDer2 ea eb ec ed c0 (m0, m1, m2, m3) w) =
      (w : ℂ) * ((c0 : ℂ) * coeff (toFS (m0, m1, m2, m3)) (monP ea eb ec ed * pderiv 2 (N 2))) := by
  unfold shDer2; dsimp only
  rw [coeff_monP_mul]
  by_cases h : ea ≤ m0 ∧ eb ≤ m1 ∧ ec ≤ m2 ∧ ed ≤ m3
  · rw [if_pos (decide_eq_true h), if_pos h, coeff_pderiv_2, evalPairs, evalPairs, add_zero,
      valFn_pack N D sc 2 _ _ _ _ (by omega) (by omega) (by omega) (by omega)]
    push_cast; ring
  · rw [if_neg (by simpa using h), if_neg h, evalPairs]; ring

lemma eval_shDer3 (ea eb ec ed : ℕ) (c0 : ℤ) (m0 m1 m2 m3 : ℕ) (w : ℤ)
    (N : Fin 4 → R4) (D : R4) (sc : ℂ)
    (hm : m0 < 31 ∧ m1 < 31 ∧ m2 < 31 ∧ m3 < 31) :
    evalPairs (valFn N D sc) (shDer3 ea eb ec ed c0 (m0, m1, m2, m3) w) =
      (w : ℂ) * ((c0 : ℂ) * coeff (toFS (m0, m1, m2, m3)) (monP ea eb ec ed * pderiv 3 (N 3))) := by
  unfold shDer3; dsimp only
  rw [coeff_monP_mul]
  by_cases h : ea ≤ m0 ∧ eb ≤ m1 ∧ ec ≤ m2 ∧ ed ≤ m3
  · rw [if_pos (decide_eq_true h), if_pos h, coeff_pderiv_3, evalPairs, evalPairs, add_zero,
      valFn_pack N D sc 3 _ _ _ _ (by omega) (by omega) (by omega) (by omega)]
    push_cast; ring
  · rw [if_neg (by simpa using h), if_neg h, evalPairs]; ring

def shSc (ea eb ec ed : ℕ) (c0 : ℤ) (m : Exps) (w : ℤ) : List (ℕ × ℤ) :=
  let (m0, m1, m2, m3) := m
  if decide (ea = m0 ∧ eb = m1 ∧ ec = m2 ∧ ed = m3) then
    [(packKey 5 0 0 0 0, c0 * w)]
  else []

lemma eval_shSc (ea eb ec ed : ℕ) (c0 : ℤ) (m0 m1 m2 m3 : ℕ) (w : ℤ)
    (N : Fin 4 → R4) (D : R4) (sc : ℂ) :
    evalPairs (valFn N D sc) (shSc ea eb ec ed c0 (m0, m1, m2, m3) w) =
      (w : ℂ) * ((c0 : ℂ) * coeff (toFS (m0, m1, m2, m3)) (C sc * monP ea eb ec ed)) := by
  unfold shSc; dsimp only
  rw [coeff_monP_sc]
  by_cases h : ea = m0 ∧ eb = m1 ∧ ec = m2 ∧ ed = m3
  · rw [if_pos (decide_eq_true h), if_pos h, evalPairs, evalPairs, add_zero,
      valFn_pack N D sc 5 0 0 0 0 (by omega) (by omega) (by omega) (by omega)]
    push_cast; ring
  · rw [if_neg (by simpa using h), if_neg h, evalPairs]; ring

def rowStrict (m : Exps) (w : ℤ) : List (ℕ × ℤ) :=
  mergeAllAux 6 2 [
    shMon 0 5 0 0 0 42 m w, shMon 0 4 1 0 0 25 m w, shMon 0 2 1 1 1 15 m w,
    shMon 0 0 5 0 0 10 m w, shMon 0 0 1 3 1 10 m w,
    shMon 1 5 0 0 0 5 m w, shMon 1 3 0 1 1 5 m w, shMon 1 1 4 0 0 50 m w,
    shMon 1 1 0 3 1 10 m w, shMon 1 0 5 0 0 42 m w,
    shMon 2 3 1 0 1 5 m w, shMon 2 1 1 2 1 30 m w, shMon 2 0 0 5 0 42 m w,
    shMon 2 0 0 4 1 25 m w, shMon 2 0 0 0 5 15 m w,
    shMon 3 3 1 1 0 5 m w, shMon 3 1 1 3 0 10 m w, shMon 3 0 0 5 0 5 m w,
    shMon 3 0 0 1 4 75 m w, shMon 3 0 0 0 5 42 m w,
    shMon 4 3 0 3 0 (-7) m w
  ]

lemma eval_rowStrict (P : Fin 4 → R4) (D : R4) (m : Exps) (w : ℤ)
    (hm : m.1 < 31 ∧ m.2.1 < 31 ∧ m.2.2.1 < 31 ∧ m.2.2.2 < 31) :
    evalPairs (valFn P D 0) (rowStrict m w) =
      (w : ℂ) * coeff (toFS m) (C 7 * ((∑ i, pderiv i core * P i) - coreSource * D)) := by
  rcases m with ⟨m0, m1, m2, m3⟩
  have hm' : m0 < 31 ∧ m1 < 31 ∧ m2 < 31 ∧ m3 < 31 := hm
  rw [rowStrict, evalPairs_mergeAllAux]
  simp only [evalListSum,
    eval_shMon_v0 _ _ _ _ _ m0 m1 m2 m3 w P D 0 hm',
    eval_shMon_v1 _ _ _ _ _ m0 m1 m2 m3 w P D 0 hm',
    eval_shMon_v2 _ _ _ _ _ m0 m1 m2 m3 w P D 0 hm',
    eval_shMon_v3 _ _ _ _ _ m0 m1 m2 m3 w P D 0 hm',
    eval_shMon_d _ _ _ _ _ m0 m1 m2 m3 w P D 0 hm']
  have hrew : C 7 * ((∑ i, pderiv i core * P i) - coreSource * D) =
      (C 7 * pderiv 0 core) * P 0 + (C 7 * pderiv 1 core) * P 1 +
      (C 7 * pderiv 2 core) * P 2 + (C 7 * pderiv 3 core) * P 3 -
      C 7 * (monP 3 0 3 0 * D) := by
    simp only [Fin.sum_univ_four, coreSource, monP]; ring
  rw [hrew, pd0_core7_eq, pd1_core7_eq, pd2_core7_eq, pd3_core7_eq]
  simp only [add_mul, coeff_add, coeff_sub, mul_assoc, coeff_C_mul]
  push_cast
  ring

def rowScalar (m : Exps) (w : ℤ) : List (ℕ × ℤ) :=
  let terms0 : List (List (ℕ × ℤ)) := [
    shSc 9 0 3 0 14 m w, shSc 8 1 3 0 10 m w, shSc 6 1 4 1 10 m w,
    shSc 4 5 3 0 20 m w, shSc 4 1 6 1 20 m w, shSc 3 6 3 0 14 m w,
    shSc 3 0 9 0 14 m w, shSc 3 0 8 1 10 m w, shSc 3 0 4 5 30 m w,
    shSc 3 0 3 6 14 m w,
    shMon 0 5 0 0 0 (-210) m w, shMon 0 4 1 0 0 (-125) m w, shMon 0 2 1 1 1 (-75) m w,
    shMon 0 0 5 0 0 (-50) m w, shMon 0 0 1 3 1 (-50) m w,
    shMon 1 5 0 0 0 (-25) m w, shMon 1 3 0 1 1 (-25) m w, shMon 1 1 4 0 0 (-250) m w,
    shMon 1 1 0 3 1 (-50) m w, shMon 1 0 5 0 0 (-210) m w,
    shMon 2 3 1 0 1 (-25) m w, shMon 2 1 1 2 1 (-150) m w, shMon 2 0 0 5 0 (-210) m w,
    shMon 2 0 0 4 1 (-125) m w, shMon 2 0 0 0 5 (-75) m w,
    shMon 3 3 1 1 0 (-25) m w, shMon 3 1 1 3 0 (-50) m w, shMon 3 0 0 5 0 (-25) m w,
    shMon 3 0 0 1 4 (-375) m w, shMon 3 0 0 0 5 (-210) m w
  ]
  let termsD01 : List (List (ℕ × ℤ)) := [
    shDer0 6 0 0 0 21 m w, shDer0 5 1 0 0 15 m w, shDer0 3 1 1 1 15 m w,
    shDer0 1 5 0 0 30 m w, shDer0 1 1 3 1 30 m w, shDer0 0 6 0 0 21 m w,
    shDer0 0 0 6 0 21 m w, shDer0 0 0 5 1 15 m w, shDer0 0 0 1 5 45 m w,
    shDer0 0 0 0 6 21 m w,
    shDer1 6 0 0 0 21 m w, shDer1 5 1 0 0 15 m w, shDer1 3 1 1 1 15 m w,
    shDer1 1 5 0 0 30 m w, shDer1 1 1 3 1 30 m w, shDer1 0 6 0 0 21 m w,
    shDer1 0 0 6 0 21 m w, shDer1 0 0 5 1 15 m w, shDer1 0 0 1 5 45 m w,
    shDer1 0 0 0 6 21 m w
  ]
  let termsD23 : List (List (ℕ × ℤ)) := [
    shDer2 6 0 0 0 21 m w, shDer2 5 1 0 0 15 m w, shDer2 3 1 1 1 15 m w,
    shDer2 1 5 0 0 30 m w, shDer2 1 1 3 1 30 m w, shDer2 0 6 0 0 21 m w,
    shDer2 0 0 6 0 21 m w, shDer2 0 0 5 1 15 m w, shDer2 0 0 1 5 45 m w,
    shDer2 0 0 0 6 21 m w,
    shDer3 6 0 0 0 21 m w, shDer3 5 1 0 0 15 m w, shDer3 3 1 1 1 15 m w,
    shDer3 1 5 0 0 30 m w, shDer3 1 1 3 1 30 m w, shDer3 0 6 0 0 21 m w,
    shDer3 0 0 6 0 21 m w, shDer3 0 0 5 1 15 m w, shDer3 0 0 1 5 45 m w,
    shDer3 0 0 0 6 21 m w
  ]
  mergeTwo 128 (mergeAllAux 6 2 terms0) (mergeTwo 128 (mergeAllAux 6 2 termsD01) (mergeAllAux 6 2 termsD23))

lemma eval_rowScalar (N : Fin 4 → R4) (sc : ℂ) (m : Exps) (w : ℤ)
    (hm : m.1 < 31 ∧ m.2.1 < 31 ∧ m.2.2.1 < 31 ∧ m.2.2.2 < 31) :
    evalPairs (valFn N 0 sc) (rowScalar m w) =
      (w : ℂ) * coeff (toFS m)
        (C 7 * (2 * C sc * coreSource * core - 5 * (∑ i, pderiv i core * N i) +
          3 * core * (∑ i, pderiv i (N i)))) := by
  rcases m with ⟨m0, m1, m2, m3⟩
  have hm' : m0 < 31 ∧ m1 < 31 ∧ m2 < 31 ∧ m3 < 31 := hm
  rw [rowScalar, evalPairs_mergeTwo, evalPairs_mergeTwo,
    evalPairs_mergeAllAux, evalPairs_mergeAllAux, evalPairs_mergeAllAux]
  simp only [evalListSum, eval_shSc,
    eval_shMon_v0 _ _ _ _ _ m0 m1 m2 m3 w N 0 sc hm',
    eval_shMon_v1 _ _ _ _ _ m0 m1 m2 m3 w N 0 sc hm',
    eval_shMon_v2 _ _ _ _ _ m0 m1 m2 m3 w N 0 sc hm',
    eval_shMon_v3 _ _ _ _ _ m0 m1 m2 m3 w N 0 sc hm',
    eval_shDer0 _ _ _ _ _ m0 m1 m2 m3 w N 0 sc hm',
    eval_shDer1 _ _ _ _ _ m0 m1 m2 m3 w N 0 sc hm',
    eval_shDer2 _ _ _ _ _ m0 m1 m2 m3 w N 0 sc hm',
    eval_shDer3 _ _ _ _ _ m0 m1 m2 m3 w N 0 sc hm']
  have hrew : C 7 * (2 * C sc * coreSource * core - 5 * (∑ i, pderiv i core * N i) +
      3 * core * (∑ i, pderiv i (N i))) =
      C 7 * (2 * C sc * coreSource * core) -
      C 5 * ((C 7 * pderiv 0 core) * N 0 + (C 7 * pderiv 1 core) * N 1 +
             (C 7 * pderiv 2 core) * N 2 + (C 7 * pderiv 3 core) * N 3) +
      C 3 * ((C 7 * core) * pderiv 0 (N 0) + (C 7 * core) * pderiv 1 (N 1) +
             (C 7 * core) * pderiv 2 (N 2) + (C 7 * core) * pderiv 3 (N 3)) := by
    simp only [Fin.sum_univ_four, ← map_ofNat C 5, ← map_ofNat C 3]; ring
  rw [hrew, source_core7_eq sc, pd0_core7_eq, pd1_core7_eq, pd2_core7_eq, pd3_core7_eq, core7_eq]
  simp only [mul_add, add_mul, coeff_add, coeff_sub, mul_assoc, coeff_C_mul]
  push_cast
  ring

def allBound (ws : List (Exps × ℤ)) : Bool :=
  ws.all fun ⟨(a, b, c, d), _⟩ => decide (a < 31 ∧ b < 31 ∧ c < 31 ∧ d < 31)

def compileStrict (ws : List (Exps × ℤ)) : List (ℕ × ℤ) :=
  mergeAllAux 10 32 (ws.map fun x => rowStrict x.1 x.2)

def compileScalar (ws : List (Exps × ℤ)) : List (ℕ × ℤ) :=
  mergeAllAux 10 32 (ws.map fun x => rowScalar x.1 x.2)

lemma eval_compileStrict_zero (P : Fin 4 → R4) (D : R4)
    (h_eq : (∑ i, pderiv i core * P i) - coreSource * D = 0) :
    ∀ (ws : List (Exps × ℤ)), allBound ws = true →
      evalPairs (valFn P D 0) (compileStrict ws) = 0 := by
  intro ws hb
  rw [compileStrict, evalPairs_mergeAllAux]
  induction ws with
  | nil => rfl
  | cons x rest ih =>
      simp only [allBound, List.all_cons, Bool.and_eq_true] at hb
      rcases x with ⟨m, w⟩
      have hm : m.1 < 31 ∧ m.2.1 < 31 ∧ m.2.2.1 < 31 ∧ m.2.2.2 < 31 := of_decide_eq_true hb.1
      simp only [List.map_cons, evalListSum, eval_rowStrict P D m w hm, ih hb.2,
        h_eq, mul_zero, coeff_zero, add_zero]

lemma zero_of_strict_cert (P : Fin 4 → R4) (D : R4)
    (h_eq : (∑ i, pderiv i core * P i) - coreSource * D = 0)
    (a b c d : ℕ) (k : ℤ) (ws : List (Exps × ℤ))
    (hb : allBound ws = true) (hk : k ≠ 0)
    (hc : compileStrict ws = [(packKey 4 a b c d, k)])
    (habcd : a < 32 ∧ b < 32 ∧ c < 32 ∧ d < 32 := by omega) :
    coeff (toFS (a, b, c, d)) D = 0 := by
  have hz := eval_compileStrict_zero P D h_eq ws hb
  rw [hc, evalPairs, evalPairs, add_zero,
    valFn_pack P D 0 4 a b c d habcd.1 habcd.2.1 habcd.2.2.1 habcd.2.2.2] at hz
  have hkC : (k : ℂ) ≠ 0 := Int.cast_ne_zero.mpr hk
  exact (mul_eq_zero.mp hz).resolve_left hkC

lemma eval_compileScalar_zero (N : Fin 4 → R4) (sc : ℂ)
    (h_eq : 2 * C sc * coreSource * core - 5 * (∑ i, pderiv i core * N i) +
      3 * core * (∑ i, pderiv i (N i)) = 0) :
    ∀ (ws : List (Exps × ℤ)), allBound ws = true →
      evalPairs (valFn N 0 sc) (compileScalar ws) = 0 := by
  intro ws hb
  rw [compileScalar, evalPairs_mergeAllAux]
  induction ws with
  | nil => rfl
  | cons x rest ih =>
      simp only [allBound, List.all_cons, Bool.and_eq_true] at hb
      rcases x with ⟨m, w⟩
      have hm : m.1 < 31 ∧ m.2.1 < 31 ∧ m.2.2.1 < 31 ∧ m.2.2.2 < 31 := of_decide_eq_true hb.1
      simp only [List.map_cons, evalListSum, eval_rowScalar N sc m w hm, ih hb.2,
        h_eq, mul_zero, coeff_zero, add_zero]

"""
    + f"def wsScalar : List (Exps × ℤ) :=\n  {fmt_weights(w_sc)}\n\n"
    + "lemma hb_sc : allBound wsScalar = true := by decide\n"
    + f"lemma hc_sc : compileScalar wsScalar = [(packKey 5 0 0 0 0, ({t_sc} : ℤ))] := by decide\n"
    + f"lemma hk_sc : ({t_sc} : ℤ) ≠ 0 := by decide\n\n"
    + """theorem no_scalar_core_cert (c : ℂ) (N : Fin 4 → R4) :
    ¬ CoreCertificate (C c * core) N := by
  rintro ⟨hD_ne, -, hcore⟩
  have hreduced : (C c * core) *
      (2 * C c * coreSource * core - 5 * (∑ i, pderiv i core * N i) +
        3 * core * (∑ i, pderiv i (N i))) = 0 := by
    have hp : ∀ i : Fin 4, pderiv i (C c * core) = C c * pderiv i core := fun i => by
      rw [pderiv_mul, pderiv_C, zero_mul, zero_add]
    simp only [Fin.sum_univ_four, hp] at hcore ⊢
    linear_combination hcore
  have hlin0 : 2 * C c * coreSource * core - 5 * (∑ i, pderiv i core * N i) +
      3 * core * (∑ i, pderiv i (N i)) = 0 :=
    (mul_eq_zero.mp hreduced).resolve_left hD_ne
  have hz := eval_compileScalar_zero N c hlin0 wsScalar hb_sc
  rw [hc_sc, evalPairs, evalPairs, add_zero,
    valFn_pack N 0 c 5 0 0 0 0 (by omega) (by omega) (by omega) (by omega)] at hz
  have hc0 : c = 0 := (mul_eq_zero.mp hz).resolve_left (Int.cast_ne_zero.mpr hk_sc)
  apply hD_ne
  rw [hc0, map_zero, zero_mul]

"""
    + "\n".join(cert_defs)
    + """
theorem no_deg4_strict (P : Fin 4 → R4) (D : R4)
    (hD_hom : D.IsHomogeneous 4)
    (hdir : ∑ i, pderiv i core * P i = coreSource * D) : D = 0 := by
  have h_eq : (∑ i, pderiv i core * P i) - coreSource * D = 0 := sub_eq_zero.mpr hdir
"""
    + "\n".join(cert_haves)
    + """
  rw [D.as_sum]
  apply Finset.sum_eq_zero
  intro s hs
  have hdeg : s 0 + s 1 + s 2 + s 3 = 4 := by
    have h := hD_hom (mem_support_iff.mp hs)
    simpa [Finsupp.weight_apply, Finsupp.sum_fintype, Fin.sum_univ_four, add_assoc] using h
  have hs0 : coeff s D = 0 := by
"""
    + decision_tree
    + """
  rw [hs0, map_zero]

lemma hom_zero_eq_C {Q : R4} (hQ : Q.IsHomogeneous 0) : Q = C (coeff 0 Q) := by
  conv_lhs => rw [Q.as_sum]
  by_cases h0 : (0 : Fin 4 →₀ ℕ) ∈ Q.support
  · rw [Finset.sum_eq_single 0]
    · simp [monomial_zero']
    · intro s hs hne
      exfalso; apply hne
      have hdeg : s 0 + s 1 + s 2 + s 3 = 0 := by
        have h := hQ (mem_support_iff.mp hs)
        simpa [Finsupp.weight_apply, Finsupp.sum_fintype, Fin.sum_univ_four, add_assoc] using h
      ext i; fin_cases i <;> simp <;> omega
    · intro hnot; exact absurd h0 hnot
  · rw [notMem_support_iff.mp h0, map_zero]
    apply Finset.sum_eq_zero
    intro s hs
    have hdeg : s 0 + s 1 + s 2 + s 3 = 0 := by
      have h := hQ (mem_support_iff.mp hs)
      simpa [Finsupp.weight_apply, Finsupp.sum_fintype, Fin.sum_univ_four, add_assoc] using h
    have hs_eq : s = 0 := by ext i; fin_cases i <;> simp <;> omega
    subst hs_eq
    exact absurd hs h0

lemma not_core_dvd_of_hom_le_four {d : ℕ} {D : R4}
    (hd : d ≤ 4) (hD_ne : D ≠ 0) (hD_hom : D.IsHomogeneous d) :
    ¬ (core ∣ D) := by
  rintro ⟨Q, hQ⟩
  have hQ_ne : Q ≠ 0 := by
    intro h0; apply hD_ne; rw [hQ, h0, mul_zero]
  have hex_comp : ∃ m ∈ Finset.range (Q.totalDegree + 1), homogeneousComponent m Q ≠ 0 := by
    by_contra hall
    apply hQ_ne
    rw [← sum_homogeneousComponent Q]
    apply Finset.sum_eq_zero
    intro m hm
    by_contra hne
    exact hall ⟨m, hm, hne⟩
  obtain ⟨m, hm_mem, hm_ne⟩ := hex_comp
  rw [Finset.mem_range] at hm_mem
  have hmul_comp : homogeneousComponent (6 + m) (core * Q) = core * homogeneousComponent m Q := by
    conv_lhs => rw [← sum_homogeneousComponent Q, Finset.mul_sum, map_sum]
    rw [Finset.sum_eq_single m]
    · exact homogeneousComponent_eq_self
        (FP2CoprimeStrict.core_hom.mul (homogeneousComponent_isHomogeneous m Q))
    · intro j _ hjm
      have hj_mem : core * homogeneousComponent j Q ∈
          homogeneousSubmodule (Fin 4) ℂ (6 + j) :=
        FP2CoprimeStrict.core_hom.mul (homogeneousComponent_isHomogeneous j Q)
      rw [homogeneousComponent_of_mem hj_mem, if_neg (by omega)]
    · intro hnot
      exact absurd (Finset.mem_range.mpr hm_mem) hnot
  have hlhs_zero : homogeneousComponent (6 + m) (core * Q) = 0 := by
    rw [← hQ, homogeneousComponent_of_mem hD_hom, if_neg (by omega)]
  rw [hlhs_zero] at hmul_comp
  exact hm_ne ((mul_eq_zero.mp hmul_comp.symm).resolve_left FP2CoprimeStrict.core_ne_zero)

end FP2Obstructions

theorem solution :
    (∀ (c : ℂ) (N : Fin 4 → MvPolynomial (Fin 4) ℂ),
      ¬ HodgeComputer.FourPlusTwoResidual.CoreCertificate
        (MvPolynomial.C c * HodgeComputer.FourPlusTwoResidual.core) N) ∧
    (∀ (D : MvPolynomial (Fin 4) ℂ) (N : Fin 4 → MvPolynomial (Fin 4) ℂ),
      D.IsHomogeneous 6 →
      HodgeComputer.FourPlusTwoResidual.CoreCertificate D N →
      ¬ (HodgeComputer.FourPlusTwoResidual.core ∣ D) ∧
        ∃ p : Fin 4 → MvPolynomial (Fin 4) ℂ,
          HodgeComputer.FourPlusTwoResidual.StrictCoreCertificate D p) ∧
    (∀ (d : ℕ) (D : MvPolynomial (Fin 4) ℂ) (N : Fin 4 → MvPolynomial (Fin 4) ℂ),
      d ≤ 4 → D.IsHomogeneous d →
      ¬ HodgeComputer.FourPlusTwoResidual.CoreCertificate D N) := by
  refine ⟨FP2Obstructions.no_scalar_core_cert, ?_, ?_⟩
  · intro D N hD6 hcert
    have hndvd : ¬ (HodgeComputer.FourPlusTwoResidual.core ∣ D) := by
      rintro ⟨Q, hQ⟩
      have hQ0 : Q.IsHomogeneous 0 :=
        FP2CoprimeStrict.isHomogeneous_of_mul_isHomogeneous
          FP2CoprimeStrict.core_hom FP2CoprimeStrict.core_ne_zero
          (show (HodgeComputer.FourPlusTwoResidual.core * Q).IsHomogeneous (6 + 0) from hQ ▸ hD6)
      have hQC := FP2Obstructions.hom_zero_eq_C hQ0
      have hD_eq : D = MvPolynomial.C (MvPolynomial.coeff 0 Q) *
          HodgeComputer.FourPlusTwoResidual.core := by
        conv_lhs => rw [hQ, hQC]
        ring
      exact FP2Obstructions.no_scalar_core_cert (MvPolynomial.coeff 0 Q) N (hD_eq ▸ hcert)
    exact ⟨hndvd, FP2Obstructions.coprime_strict D hndvd N hcert⟩
  · intro d D N hd hDd hcert
    have hD_ne : D ≠ 0 := hcert.1
    have hndvd := FP2Obstructions.not_core_dvd_of_hom_le_four hd hD_ne hDd
    obtain ⟨p, -, -, hp_dir, -⟩ := FP2Obstructions.coprime_strict D hndvd N hcert
    let G : MvPolynomial (Fin 4) ℂ := MvPolynomial.X 0 ^ (4 - d)
    have hG_hom : G.IsHomogeneous (4 - d) := by
      dsimp [G]
      have h1 : (MvPolynomial.X 0 ^ (4 - d) : MvPolynomial (Fin 4) ℂ).IsHomogeneous (1 * (4 - d)) :=
        (MvPolynomial.isHomogeneous_X ℂ (0 : Fin 4)).pow (4 - d)
      have heq : 1 * (4 - d) = 4 - d := by omega
      rwa [heq] at h1
    have hGD_hom : (G * D).IsHomogeneous 4 := by
      have h1 : (G * D).IsHomogeneous ((4 - d) + d) := hG_hom.mul hDd
      have heq : (4 - d) + d = 4 := by omega
      rwa [heq] at h1
    have hGD_dir : ∑ i, MvPolynomial.pderiv i HodgeComputer.FourPlusTwoResidual.core * (G * p i) =
        HodgeComputer.FourPlusTwoResidual.coreSource * (G * D) := by
      calc
        ∑ i, MvPolynomial.pderiv i HodgeComputer.FourPlusTwoResidual.core * (G * p i)
            = G * ∑ i, MvPolynomial.pderiv i HodgeComputer.FourPlusTwoResidual.core * p i := by
              rw [Finset.mul_sum]
              exact Finset.sum_congr rfl (fun i _ => by ring)
        _ = G * (HodgeComputer.FourPlusTwoResidual.coreSource * D) := by rw [hp_dir]
        _ = HodgeComputer.FourPlusTwoResidual.coreSource * (G * D) := by ring
    have hGD0 : G * D = 0 := FP2Obstructions.no_deg4_strict (fun i => G * p i) (G * D) hGD_hom hGD_dir
    have hG_ne : G ≠ 0 := by
      intro h0
      have hev := congrArg (MvPolynomial.eval ![1, 0, 0, 0]) h0
      simp [G] at hev
    exact hD_ne ((mul_eq_zero.mp hGD0).resolve_left hG_ne)
"""
)

out_path = WORK / "Sol_no_low_degree_or_scalar_core_certificate.lean"
out_path.write_text(lean_file)
print(f"Wrote {out_path} ({len(lean_file)} bytes)")
