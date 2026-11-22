import HodgeComputer.Constructions.RationalFields.Calculus
import Mathlib.Tactic

/-! Specialize a finite polar part for the linear pencil P + ε S. -/

noncomputable section
namespace FourPlusTwoSpecialization

open HodgeTestbed.Conic
open scoped BigOperators

variable {ι K : Type*} [Fintype ι] [Field K] [CharZero K]

/-- The chart source of a rational current. -/
def current (δ : ι → Derivation ℚ K K) (P : K) (q : ι → K) : K :=
  (2 / 3 : K) * fieldAction δ q P - P * fieldDivergence δ q

/-- Running coefficient in the telescoping identity. -/
def d : ℕ → ℚ
  | 0 => 1
  | j+1 => -(2 / (3 * ((j : ℚ) + 1))) * d j - d j

/-- Coefficient of the j-th Laurent current in the specialized primitive. -/
def c : ℕ → ℚ
  | 0 => 1
  | j+1 => -(2 / (3 * ((j : ℚ) + 1))) * d j

lemma d_succ (j : ℕ) : d (j+1) = c (j+1) - d j := rfl

lemma c_succ_mul (j : ℕ) : c (j+1) * ((j : ℚ) + 1) = -(2 / 3 : ℚ) * d j := by
  have h : (j : ℚ) + 1 ≠ 0 := by positivity
  simp only [c]
  field_simp [h]

lemma current_add (δ : ι → Derivation ℚ K K) (P : K) (v w : ι → K) :
    current δ P (v+w) = current δ P v + current δ P w := by
  change current δ P (fun i => v i + w i) = _
  simp only [current, fieldAction_add, fieldDivergence_add]
  ring

/-- A single term of the pole correction, with S = P h. -/
theorem term (δ : ι → Derivation ℚ K K) (P h : K) (v : ι → K) (j : ℕ) :
    current δ P (fun i => algebraMap ℚ K (c (j+1)) * h^(j+1) * v i) =
      algebraMap ℚ K (d (j+1)) * h^(j+1) * current δ P v +
        algebraMap ℚ K (d j) * h^j * current δ (P*h) v := by
  have hc : algebraMap ℚ K (c (j+1)) * ((j : K) + 1) =
      -(2 / 3 : K) * algebraMap ℚ K (d j) := by
    simpa only [map_mul, map_add, map_one, map_neg, map_div₀, map_ofNat, map_natCast] using
      congrArg (algebraMap ℚ K) (c_succ_mul j)
  have hdiff : fieldAction δ v (algebraMap ℚ K (c (j+1)) * h^(j+1)) =
      algebraMap ℚ K (c (j+1)) * (((j : K) + 1) * h^j * fieldAction δ v h) := by
    rw [Derivation.leibniz, Derivation.map_algebraMap, Derivation.leibniz_pow]
    simp [smul_eq_mul, nsmul_eq_mul, mul_assoc]
  simp only [current, fieldAction_scale, fieldDivergence_scale, hdiff, d_succ,
    map_sub, Derivation.leibniz, smul_eq_mul]
  linear_combination -P * h^j * (fieldAction δ v h) * hc

/-- The constant term plus the weighted polar corrections through order m. -/
def corrected (h : K) (q : ℕ → ι → K) (m : ℕ) : ι → K :=
  ∑ j ∈ Finset.range (m+1), fun i => algebraMap ℚ K (c j) * h^j * q j i

/-- Partial summation retains just the next coefficient of the Laurent chain. -/
theorem partial_sum (δ : ι → Derivation ℚ K K) (P h A : K) (q : ℕ → ι → K) (m : ℕ)
    (hzero : current δ P (q 0) + current δ (P*h) (q 1) = A)
    (hpolar : ∀ j, j < m → current δ P (q (j+1)) + current δ (P*h) (q (j+2)) = 0) :
    current δ P (corrected h q m) =
      A - algebraMap ℚ K (d m) * h^m * current δ (P*h) (q (m+1)) := by
  induction m with
  | zero =>
    have hz : corrected h q 0 = q 0 := by ext i; simp [corrected, c]
    rw [hz]
    simpa [d, eq_sub_iff_add_eq] using hzero
  | succ m ih =>
    have prev := ih (fun j hj => hpolar j (by omega))
    have next := hpolar m (by omega)
    rw [show corrected h q (m+1) = corrected h q m +
         (fun i => algebraMap ℚ K (c (m+1)) * h^(m+1) * q (m+1) i) by
           unfold corrected
           rw [Finset.sum_range_succ _ (m+1)], current_add, prev, term]
    linear_combination (algebraMap ℚ K (d (m+1)) * h^(m+1)) * next

/-- A Laurent coefficient chain of any finite pole order specializes to an
actual current with source A at ε=0. -/
theorem specialize (δ : ι → Derivation ℚ K K) (P S A : K) (hP : P ≠ 0)
    (q : ℕ → ι → K) (m : ℕ)
    (hzero : current δ P (q 0) + current δ S (q 1) = A)
    (hpolar : ∀ j, j < m → current δ P (q (j+1)) + current δ S (q (j+2)) = 0)
    (htop : q (m+1) = 0) :
    current δ P (corrected (S/P) q m) = A := by
  have hp : P * (S/P) = S := by field_simp
  have h := partial_sum δ P (S/P) A q m (by simpa only [hp] using hzero)
    (fun j hj => by simpa only [hp] using hpolar j hj)
  rw [h, htop]
  simp [current, fieldDivergence]

end FourPlusTwoSpecialization
