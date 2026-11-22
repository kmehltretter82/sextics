import Mathlib.RingTheory.Derivation.Basic
import Mathlib.Tactic

noncomputable section

namespace FourPlusTwoRatio

variable {K : Type*} [Field K] [CharZero K]

/-- Push the ratio-chart current in coordinates (r,s,t) to the core with
(a,b,c,d) = (a,ar,at,ast).  Its a component is zero. -/
theorem transport (Dr Ds Dt : Derivation ℤ K K)
    (a s t P qr qs qt : K) (ha : a ≠ 0) (ht : t ≠ 0)
    (hra : Dr a = 0) (hsa : Ds a = 0) (hta : Dt a = 0)
    (hrt : Dr t = 0) (hst : Ds t = 0) (htt : Dt t = 1)
    (hss : Ds s = 1)
    (source : -(P * (Dr qr + Ds qs + Dt qt) -
        (2 / 3 : K) * (qr * Dr P + qs * Ds P + qt * Dt P)) = t ^ 4) :
    let H := a ^ 6 * P
    let db := fun z : K => Dr z / a
    let dc := fun z : K => (Dt z - s / t * Ds z) / a
    let dd := fun z : K => Ds z / (a * t)
    let pb := (2 : K) * a * qr / (3 * t)
    let pc := (2 : K) * a * qt / (3 * t)
    let pd := (2 : K) * a * (s * qt + t * qs) / (3 * t)
    pb * db H + pc * dc H + pd * dd H -
      (3 / 2 : K) * H * (db pb + dc pc + dd pd) = a ^ 6 * t ^ 3 := by
  intro H db dc dd pb pc pd
  have h₂ (δ : Derivation ℤ K K) : δ 2 = 0 := δ.map_natCast 2
  have h₃ (δ : Derivation ℤ K K) : δ 3 = 0 := δ.map_natCast 3
  have hdiv : db pb + dc pc + dd pd =
      (2 / (3 * t)) * (Dr qr + Ds qs + Dt qt) := by
    dsimp only [db, dc, dd, pb, pc, pd]
    simp only [Derivation.leibniz_div, Derivation.leibniz, map_add, smul_eq_mul,
      h₂, h₃, hra, hsa, hta, hrt, hst, htt, hss]
    field_simp [ha, ht]
    ring
  have hdir : pb * db H + pc * dc H + pd * dd H =
      (2 * a ^ 6 / (3 * t)) * (qr * Dr P + qs * Ds P + qt * Dt P) := by
    dsimp only [db, dc, dd, pb, pc, pd, H]
    simp only [Derivation.leibniz, Derivation.leibniz_pow, smul_eq_mul, hra, hsa, hta,
      mul_zero]
    field_simp [ha, ht]
    ring
  rw [hdir, hdiv]
  dsimp only [H]
  field_simp [ht] at source ⊢
  linear_combination source

end FourPlusTwoRatio
