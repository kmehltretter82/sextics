import RationalCore
import HodgeComputer.Constructions.ResidueTransport.CoordinateIndependence

/-! A linear polynomial equation on the affine core implies the protected theorem. -/
noncomputable section
namespace FourPlusTwoCertificate
open MvPolynomial
open scoped BigOperators

abbrev AffinePoly := MvPolynomial (Fin 3) ℂ

/-- `H(1,x,y,z)`, with local indices 0,1,2 for x,y,z. -/
def affineCore : AffinePoly :=
  1 + C (5/7 : ℂ) * X 0 + C (5/7 : ℂ) * X 0 * X 1 * X 2 + C (10/7 : ℂ) * X 0^5 +
  C (10/7 : ℂ) * X 0 * X 1^3 * X 2 + X 0^6 + X 1^6 + C (5/7 : ℂ) * X 1^5 * X 2 +
  C (15/7 : ℂ) * X 1 * X 2^5 + X 2^6

/-- Integer-cleared coefficient identity for `q_i = N_i/B^n`. It is linear
in the three polynomial numerators for each fixed `B,n`. -/
def FixedSupportEquation (B : AffinePoly) (n : ℕ) (N : Fin 3 → AffinePoly) : Prop :=
  2*B*(∑ i, N i * pderiv i affineCore) - 3*affineCore*
    (B*(∑ i, pderiv i (N i)) - (n : AffinePoly)*(∑ i, N i * pderiv i B)) =
    3 * X 1^3 * B^(n+1)

private theorem d_div_pow {K : Type*} [Field K] (D : Derivation ℤ K K)
    (x y : K) (n : ℕ) (hy : y ≠ 0) :
    D (x/y^n) = (y*D x - (n : K)*x*D y)/y^(n+1) := by
  cases n with
  | zero => simp [hy]
  | succ m =>
    rw [Derivation.leibniz_div, Derivation.leibniz_pow]
    simp only [smul_eq_mul, nsmul_eq_mul, Nat.succ_sub_one, Nat.cast_succ, pow_succ]
    field_simp [hy] <;> ring

private theorem affine_equation {K : Type*} [Field K] [CharZero K]
    (θ : AffinePoly →+* K) (E : Fin 3 → Derivation ℤ K K)
    (compat : ∀ i p, E i (θ p) = θ (pderiv i p))
    (B : AffinePoly) (n : ℕ) (N : Fin 3 → AffinePoly) (hb : θ B ≠ 0)
    (hraw : FixedSupportEquation B n N) :
    let q : Fin 3 → K := fun i => θ (N i)/(θ B)^n
    2*(∑ i, q i * E i (θ affineCore)) - 3*θ affineCore*(∑ i, E i (q i)) = 3*θ (X 1)^3 := by
  intro q
  have hmul : (θ B)^(n+1)*(∑ i, q i*E i (θ affineCore)) =
      θ B * (∑ i, θ (N i)*θ (pderiv i affineCore)) := by
    rw [Finset.mul_sum, Finset.mul_sum]
    apply Finset.sum_congr rfl
    intro i _
    dsimp only [q]
    rw [compat, pow_succ]
    field_simp [hb]
  have hdiv : (θ B)^(n+1)*(∑ i, E i (q i)) =
      θ B*(∑ i, θ (pderiv i (N i))) - (n : K)*(∑ i, θ (N i)*θ (pderiv i B)) := by
    simp only [Finset.mul_sum, ← Finset.sum_sub_distrib]
    apply Finset.sum_congr rfl
    intro i _
    dsimp only [q]
    rw [d_div_pow _ _ _ _ hb, compat, compat]
    field_simp [hb]
  have h := congrArg θ hraw
  simp only [map_sub, map_mul, map_sum, map_pow, map_ofNat, map_natCast] at h
  apply mul_left_cancel₀ (pow_ne_zero (n+1) hb)
  linear_combination h + 2*hmul - 3*θ affineCore*hdiv

section Field
variable {K : Type*} [Field K] [CharZero K] [Algebra CorePoly K] [IsFractionRing CorePoly K]

private theorem corePartial_map (i : Fin 4) (p : CorePoly) :
    (corePartial i : Derivation ℤ K K) (algebraMap CorePoly K p) =
      algebraMap CorePoly K (pderiv i p) := by
  let : IsScalarTower ℚ CorePoly K := IsScalarTower.of_algebraMap_eq' (Subsingleton.elim _ _)
  exact QuadricContactFraction.extend_map (F := K)
    ((pderiv i : Derivation ℂ CorePoly CorePoly).restrictScalars ℚ) p

omit [CharZero K] [Algebra CorePoly K] [IsFractionRing CorePoly K] in
private theorem chart_identity (φ : CorePoly →+* K) (a : K) (ha : a ≠ 0) (haφ : a=φ (X 0)) :
    let x := fun i : Fin 3 => φ (X i.succ)/a
    let θ := eval₂Hom (φ.comp C) x
    φ core = a^6*θ affineCore ∧ φ coreSource = a^6*(θ (X 1))^3 := by
  intro x θ
  constructor <;>
    simp only [core, affineCore, coreSource, θ, x, map_add, map_mul, map_pow, map_one,
      eval₂Hom_X', eval₂Hom_C, RingHom.comp_apply, ← haφ, show (0 : Fin 3).succ = (1 : Fin 4) from rfl,
      show (1 : Fin 3).succ = (2 : Fin 4) from rfl, show (2 : Fin 3).succ = (3 : Fin 4) from rfl] <;>
    field_simp [ha] <;> ring

private theorem of_affine_equation (B : AffinePoly) (n : ℕ) (N : Fin 3 → AffinePoly) (hB : B ≠ 0)
    (hraw : FixedSupportEquation B n N) :
    ∃ v : Fin 4 → K, CoreFieldEquation v := by
  classical
  let φ := algebraMap CorePoly K
  let a := φ (X 0)
  let x := fun i : Fin 3 => φ (X i.succ)/a
  let θ := eval₂Hom (φ.comp C) x
  let δ := fun i : Fin 4 => (corePartial i : Derivation ℤ K K)
  let E := fun i : Fin 3 => a • δ i.succ
  have ha : a ≠ 0 := IsFractionRing.to_map_ne_zero_of_mem_nonZeroDivisors
    (mem_nonZeroDivisors_iff_ne_zero.mpr (X_ne_zero (R := ℂ) (0 : Fin 4)))
  have hDa (i : Fin 3) : δ i.succ a = 0 := by
    simp [δ, a, φ, corePartial_map, pderiv_X]
  have hDc (i j : Fin 3) : δ i.succ (φ (X j.succ)) = if j=i then 1 else 0 := by
    simp only [δ, φ, corePartial_map, pderiv_X]
    split <;> simp_all
  have hC (i : Fin 3) (c : ℂ) : E i ((φ.comp C) c) = 0 := by
    simp only [E, δ, φ, Derivation.smul_apply, smul_eq_mul, RingHom.comp_apply,
      corePartial_map, pderiv_C, map_zero, mul_zero]
  have hx (i j : Fin 3) : E i (x j) = if j=i then 1 else 0 := by
    simp only [E, x, Derivation.smul_apply, smul_eq_mul, Derivation.leibniz_div, hDa,
      hDc, mul_zero, sub_zero]
    split <;> simp_all
    field_simp [ha]
  have hθ (i p) : E i (θ p) = θ (pderiv i p) := by
    rw [show E i (θ p) = ∑ j, θ (pderiv j p) * E i (x j) from
      HoldoutConstruction.derivation_eval₂_chain (E i) (φ.comp C) x (hC i) p]
    simp [hx]
  have hinj := ResidueCoordinateIndependence.injective_of_coordinate_derivations (φ.comp C) x E hC hx
  have hb : θ B ≠ 0 := by
    intro hz
    exact hB (hinj (by simpa only [map_zero] using hz))
  let q := fun i => θ (N i)/(θ B)^n
  have hq := affine_equation θ E hθ B n N hb hraw
  change 2*(∑ i, q i * E i (θ affineCore)) - 3*θ affineCore*(∑ i, E i (q i)) = 3*θ (X 1)^3 at hq
  obtain ⟨hH, hA⟩ := chart_identity φ a ha rfl
  change φ core = a^6*θ affineCore at hH
  change φ coreSource = a^6*θ (X 1)^3 at hA
  let v : Fin 4 → K := Fin.cases 0 (fun i => (2/3 : K)*a*q i)
  have hdir : (∑ j, v j * δ j (φ core)) = (2/3 : K)*a^6 * ∑ i, q i * E i (θ affineCore) := by
    rw [Fin.sum_univ_succ, Finset.mul_sum]
    simp only [v, Fin.cases_zero, Fin.cases_succ, zero_mul, zero_add]
    apply Finset.sum_congr rfl
    intro i _
    rw [hH, Derivation.leibniz, Derivation.leibniz_pow, hDa]
    simp only [E, Derivation.smul_apply, smul_eq_mul, mul_zero]
    ring
  have h23 (i : Fin 3) : δ i.succ (2/3 : K) = 0 := by
    rw [Derivation.leibniz_div, show δ i.succ (2 : K) = 0 from (δ i.succ).map_natCast 2,
      show δ i.succ (3 : K) = 0 from (δ i.succ).map_natCast 3]
    simp
  have hdiv : (∑ j, δ j (v j)) = (2/3 : K) * ∑ i, E i (q i) := by
    rw [Fin.sum_univ_succ, Finset.mul_sum]
    simp only [v, Fin.cases_zero, Fin.cases_succ, map_zero, zero_add]
    apply Finset.sum_congr rfl
    intro i _
    simp only [Derivation.leibniz, h23, hDa, smul_eq_mul, E, Derivation.smul_apply, mul_zero,
      add_zero]
    ring
  refine ⟨v, ?_⟩
  change (∑ j, v j * δ j (φ core)) - (3/2 : K)*φ core*(∑ j, δ j (v j)) = φ coreSource
  rw [hdir, hdiv, hH, hA]
  linear_combination (a^6/3)*hq

end Field

/-- A three-variable fixed-support polynomial solution discharges the
protected four-plus-two instance. `B` may be any nonzero polynomial. -/
theorem hasCertificate_of_fixed_support (B : AffinePoly) (n : ℕ) (N : Fin 3 → AffinePoly)
    (hB : B ≠ 0) (hraw : FixedSupportEquation B n N) :
    HodgeComputer.HasCertificate (n := HodgeSubmission.n)
      HodgeSubmission.F HodgeSubmission.q HodgeSubmission.A := by
  obtain ⟨v, h⟩ := of_affine_equation (K := FractionRing CorePoly) B n N hB hraw
  exact hasCertificate_of_core_field v h

end FourPlusTwoCertificate
