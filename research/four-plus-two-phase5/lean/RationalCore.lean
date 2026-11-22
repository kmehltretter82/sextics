import UngradedCore
import HodgeComputer.Constructions.FiniteAlgebra.FractionDerivation
import Mathlib.RingTheory.Localization.FractionRing
import Mathlib.RingTheory.Localization.Integer

/-! The exact rational core current implies the protected four-plus-two target. -/
noncomputable section
namespace FourPlusTwoCertificate
open MvPolynomial
open scoped BigOperators

/-- Clear a common denominator in the fraction field of the four literal core
variables. The coordinate derivations are specified by their restrictions to
polynomials. No denominator or homogeneous degree bound is required. -/
theorem hasCertificate_of_rational_current {K : Type*} [Field K] [CharZero K]
    [Algebra CorePoly K] [IsFractionRing CorePoly K]
    (δ : Fin 4 → Derivation ℤ K K)
    (hδ : ∀ i p, δ i (algebraMap CorePoly K p) = algebraMap CorePoly K (pderiv i p))
    (v : Fin 4 → K)
    (hsource : (∑ i, v i * δ i (algebraMap CorePoly K core)) -
      (3 / 2 : K) * algebraMap CorePoly K core * (∑ i, δ i (v i)) = algebraMap CorePoly K coreSource) :
    HodgeComputer.HasCertificate (n := HodgeSubmission.n)
      HodgeSubmission.F HodgeSubmission.q HodgeSubmission.A := by
  classical
  let φ := algebraMap CorePoly K
  obtain ⟨d, hd⟩ := IsLocalization.exist_integer_multiples_of_finite (nonZeroDivisors CorePoly) v
  choose N hN using hd
  let D : CorePoly := d
  have hD : D ≠ 0 := mem_nonZeroDivisors_iff_ne_zero.mp d.property
  have hv (i) : φ (N i) = φ D * v i := by simpa only [Algebra.smul_def] using hN i
  have hv' (i) : φ (pderiv i (N i)) = φ D * δ i (v i) + v i * φ (pderiv i D) := by
    have hi := congrArg (δ i) (hv i)
    simpa only [φ, hδ, Derivation.leibniz, smul_eq_mul] using hi
  have hn : (∑ i, φ (pderiv i core) * φ (N i)) = φ D * ∑ i, v i * δ i (φ core) := by
    rw [Finset.mul_sum]
    apply Finset.sum_congr rfl
    intro i _
    rw [hv i, show δ i (φ core) = φ (pderiv i core) from hδ i core]
    ring
  have hc : (∑ i, φ (pderiv i D) * φ (N i)) = φ D * ∑ i, v i * φ (pderiv i D) := by
    rw [Finset.mul_sum]
    apply Finset.sum_congr rfl
    intro i _
    rw [hv i]
    ring
  have hdiv : φ D * (∑ i, φ (pderiv i (N i))) - (∑ i, φ (pderiv i D) * φ (N i)) =
      (φ D)^2 * ∑ i, δ i (v i) := by
    simp only [hv', Finset.sum_add_distrib, ← Finset.mul_sum, hc]
    ring
  apply hasCertificate_of_ungraded D N hD
  apply IsFractionRing.injective CorePoly K
  change φ _ = φ 0
  simp only [map_add, map_sub, map_mul, map_pow, map_sum, map_ofNat, map_zero, hn, hdiv]
  change (∑ i, v i * δ i (φ core)) - (3 / 2 : K) * φ core * (∑ i, δ i (v i)) = φ coreSource at hsource
  linear_combination -2*(φ D)^2 * hsource

section Canonical
variable {K : Type*} [Field K] [CharZero K] [Algebra CorePoly K] [IsFractionRing CorePoly K]

local instance : IsScalarTower ℚ CorePoly K :=
  IsScalarTower.of_algebraMap_eq' (Subsingleton.elim _ _)

/-- Coordinate partials on any presentation of the fraction field of the core. -/
def corePartial (i : Fin 4) : Derivation ℤ K K :=
  (QuadricContactFraction.extend (F := K)
    ((pderiv i : Derivation ℂ CorePoly CorePoly).restrictScalars ℚ)).restrictScalars ℤ

/-- The rational twisted equation on the four literal core variables. -/
def CoreFieldEquation (v : Fin 4 → K) : Prop :=
  (∑ i, v i * corePartial i (algebraMap CorePoly K core)) -
    (3 / 2 : K) * algebraMap CorePoly K core * (∑ i, corePartial i (v i)) =
    algebraMap CorePoly K coreSource

/-- The only input on the rational current is its equation with the
canonical coordinate partials on the fraction field. -/
theorem hasCertificate_of_fraction_current (v : Fin 4 → K) (h : CoreFieldEquation v) :
    HodgeComputer.HasCertificate (n := HodgeSubmission.n)
      HodgeSubmission.F HodgeSubmission.q HodgeSubmission.A := by
  refine hasCertificate_of_rational_current (K := K) corePartial ?_ v h
  intro i p
  exact QuadricContactFraction.extend_map (F := K)
    ((pderiv i : Derivation ℂ CorePoly CorePoly).restrictScalars ℚ) p

end Canonical

/-- The canonical instance of the result on the literal field of core rational functions. -/
theorem hasCertificate_of_core_field (v : Fin 4 → FractionRing CorePoly) (h : CoreFieldEquation v) :
    HodgeComputer.HasCertificate (n := HodgeSubmission.n)
      HodgeSubmission.F HodgeSubmission.q HodgeSubmission.A :=
  hasCertificate_of_fraction_current v h

end FourPlusTwoCertificate
