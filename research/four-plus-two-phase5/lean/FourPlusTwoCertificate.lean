import HodgeTarget
import Mathlib.Algebra.MvPolynomial.Equiv
import Mathlib.Algebra.MvPolynomial.Funext
import Mathlib.Algebra.Polynomial.Degree.SmallDegree
import Mathlib.Algebra.Polynomial.Derivative
import Mathlib.Algebra.Polynomial.Eval.Irreducible
import Mathlib.Analysis.Complex.Polynomial.Basic
import Mathlib.RingTheory.MvPolynomial.EulerIdentity
import Mathlib.RingTheory.Polynomial.Eisenstein.Basic
import Mathlib.RingTheory.Polynomial.UniqueFactorization
import Mathlib.Tactic.FinCases
import Mathlib.Tactic.LinearCombination
import Mathlib.Tactic.Ring

/-!
This specializes the suspension construction and irreducibility proof in the
accepted Prove2Me theorem `normalChart_residue_of_coreCertificate`
(5f8e0726-2721-43c9-a74a-c9e0e987f824) to the generated sextics target.
The four-variable and binary sextics come from the project's four-plus-two
core reduction.  The final constructor discharges the protected verifier.
-/

set_option autoImplicit false
set_option maxHeartbeats 4000000
noncomputable section

open MvPolynomial HodgeComputer HodgeSubmission
open scoped BigOperators

namespace FourPlusTwoCertificate

/-- The four-variable core `H(x₀,x₁,y₀,y₁)` of the four-plus-two sextic
`F = H(X₀,X₁,Y₀,Y₁) + B(X₂,Y₂)`. Variables `0,1,2,3` stand for `X₀,X₁,Y₀,Y₁`. -/
def core : MvPolynomial (Fin 4) ℂ :=
  X 0 ^ 6 + C (5 / 7 : ℂ) * X 0 ^ 5 * X 1 + C (5 / 7 : ℂ) * X 0 ^ 3 * X 1 * X 2 * X 3
    + C (10 / 7 : ℂ) * X 0 * X 1 ^ 5 + C (10 / 7 : ℂ) * X 0 * X 1 * X 2 ^ 3 * X 3 + X 1 ^ 6
    + X 2 ^ 6 + C (5 / 7 : ℂ) * X 2 ^ 5 * X 3 + C (15 / 7 : ℂ) * X 2 * X 3 ^ 5 + X 3 ^ 6

/-- The source `x₀³·y₀³` in the core variables. -/
def coreSource : MvPolynomial (Fin 4) ℂ := X 0 ^ 3 * X 2 ^ 3

/-- A strict twisted certificate on the core: a nonconstant homogeneous denominator `D`
and numerators `N` of one degree more such that the rational field `v = N / D` satisfies
`v(H) - (3/2)·H·div(v) = A`, written after clearing `D²` and the factor `2`. -/
def CoreCertificate (D : MvPolynomial (Fin 4) ℂ) (N : Fin 4 → MvPolynomial (Fin 4) ℂ) : Prop :=
  D ≠ 0 ∧ (∃ k : ℕ, D.IsHomogeneous (k + 1) ∧ ∀ i, (N i).IsHomogeneous (k + 2)) ∧
    2 * coreSource * D ^ 2 - 2 * D * (∑ i, pderiv i core * N i)
      + 3 * core * (D * (∑ i, pderiv i (N i)) - ∑ i, pderiv i D * N i) = 0



abbrev CorePoly := MvPolynomial (Fin 4) ℂ
abbrev BinaryPoly := MvPolynomial (Fin 2) ℂ
abbrev SixPoly := MvPolynomial (Fin 6) ℂ

private def F6 : SixPoly := F
private def A6 : SixPoly := A

def binaryB : BinaryPoly :=
  X 0 ^ 6 + C (15 / 7 : ℂ) * X 0 ^ 4 * X 1 ^ 2 +
    C (10 / 7 : ℂ) * X 0 ^ 2 * X 1 ^ 4 + X 1 ^ 6

def liftCore : CorePoly →+* SixPoly := (rename ![0, 1, 3, 4]).toRingHom
def liftBinary : BinaryPoly →+* SixPoly := (rename ![2, 5]).toRingHom

lemma liftCore_X0 : liftCore (X 0) = X 0 := by simp [liftCore]
lemma liftCore_X1 : liftCore (X 1) = X 1 := by simp [liftCore]
lemma liftCore_X2 : liftCore (X 2) = X 3 := by simp [liftCore]
lemma liftCore_X3 : liftCore (X 3) = X 4 := by simp [liftCore]
lemma liftBinary_X0 : liftBinary (X 0) = X 2 := by simp [liftBinary]
lemma liftBinary_X1 : liftBinary (X 1) = X 5 := by simp [liftBinary]

lemma targetF_eq : F6 = liftCore core + liftBinary binaryB := by
  apply MvPolynomial.funext (σ := Fin 6)
  intro x
  simp only [F6, F, core, binaryB, liftCore, liftBinary, AlgHom.toRingHom_eq_coe, RingHom.coe_coe,
    eval_add, eval_mul, eval_pow, eval_C, eval_X, eval_ofNat, eval_rename, Function.comp_apply,
    Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.head_cons, Matrix.cons_val_two,
    Matrix.tail_cons, Matrix.cons_val_three]
  ring

lemma targetA_eq : A6 = liftCore coreSource := by
  simp only [A6, A, coreSource, map_mul, map_pow, liftCore_X0, liftCore_X2]

lemma core_hom : core.IsHomogeneous 6 := by
  have hX := isHomogeneous_X (σ := Fin 4) ℂ
  have hC := isHomogeneous_C (Fin 4) (R := ℂ)
  have t1 : (X 0 ^ 6 : CorePoly).IsHomogeneous 6 := (hX 0).pow 6
  have t2 : (C (5 / 7 : ℂ) * X 0 ^ 5 * X 1 : CorePoly).IsHomogeneous 6 :=
    ((hC _).mul ((hX 0).pow 5)).mul (hX 1)
  have t3 : (C (5 / 7 : ℂ) * X 0 ^ 3 * X 1 * X 2 * X 3 : CorePoly).IsHomogeneous 6 :=
    ((((hC _).mul ((hX 0).pow 3)).mul (hX 1)).mul (hX 2)).mul (hX 3)
  have t4 : (C (10 / 7 : ℂ) * X 0 * X 1 ^ 5 : CorePoly).IsHomogeneous 6 :=
    ((hC _).mul (hX 0)).mul ((hX 1).pow 5)
  have t5 : (C (10 / 7 : ℂ) * X 0 * X 1 * X 2 ^ 3 * X 3 : CorePoly).IsHomogeneous 6 :=
    ((((hC _).mul (hX 0)).mul (hX 1)).mul ((hX 2).pow 3)).mul (hX 3)
  have t6 : (X 1 ^ 6 : CorePoly).IsHomogeneous 6 := (hX 1).pow 6
  have t7 : (X 2 ^ 6 : CorePoly).IsHomogeneous 6 := (hX 2).pow 6
  have t8 : (C (5 / 7 : ℂ) * X 2 ^ 5 * X 3 : CorePoly).IsHomogeneous 6 :=
    ((hC _).mul ((hX 2).pow 5)).mul (hX 3)
  have t9 : (C (15 / 7 : ℂ) * X 2 * X 3 ^ 5 : CorePoly).IsHomogeneous 6 :=
    ((hC _).mul (hX 2)).mul ((hX 3).pow 5)
  have t10 : (X 3 ^ 6 : CorePoly).IsHomogeneous 6 := (hX 3).pow 6
  unfold core
  exact ((((((((t1.add t2).add t3).add t4).add t5).add t6).add t7).add t8).add t9).add t10

lemma binaryB_hom : binaryB.IsHomogeneous 6 := by
  have hX := isHomogeneous_X (σ := Fin 2) ℂ
  have hC := isHomogeneous_C (Fin 2) (R := ℂ)
  have b1 : (X 0 ^ 6 : BinaryPoly).IsHomogeneous 6 := (hX 0).pow 6
  have b2 : (C (15 / 7 : ℂ) * X 0 ^ 4 * X 1 ^ 2 : BinaryPoly).IsHomogeneous 6 :=
    ((hC _).mul ((hX 0).pow 4)).mul ((hX 1).pow 2)
  have b3 : (C (10 / 7 : ℂ) * X 0 ^ 2 * X 1 ^ 4 : BinaryPoly).IsHomogeneous 6 :=
    ((hC _).mul ((hX 0).pow 2)).mul ((hX 1).pow 4)
  have b4 : (X 1 ^ 6 : BinaryPoly).IsHomogeneous 6 := (hX 1).pow 6
  unfold binaryB
  exact ((b1.add b2).add b3).add b4

lemma liftCore_hom {p : CorePoly} {m : ℕ} (hp : p.IsHomogeneous m) :
    (liftCore p).IsHomogeneous m :=
  hp.rename_isHomogeneous (f := ![0, 1, 3, 4])

lemma liftBinary_hom {p : BinaryPoly} {m : ℕ} (hp : p.IsHomogeneous m) :
    (liftBinary p).IsHomogeneous m :=
  hp.rename_isHomogeneous (f := ![2, 5])

lemma pderiv_liftCore (p : CorePoly) (j : Fin 6) :
    pderiv j (liftCore p) =
      if j = 0 then liftCore (pderiv 0 p)
      else if j = 1 then liftCore (pderiv 1 p)
      else if j = 3 then liftCore (pderiv 2 p)
      else if j = 4 then liftCore (pderiv 3 p)
      else 0 := by
  induction p using MvPolynomial.induction_on generalizing j with
  | C z =>
    simp only [liftCore, AlgHom.toRingHom_eq_coe, RingHom.coe_coe, rename_C, pderiv_C, map_zero]
    split_ifs <;> rfl
  | add p q hp hq =>
    simp only [map_add, hp, hq]
    split_ifs <;> ring
  | mul_X p i hp =>
    rw [map_mul]
    fin_cases i <;> fin_cases j <;>
      simp [liftCore_X0, liftCore_X1, liftCore_X2, liftCore_X3, pderiv_X, hp]

lemma pderiv_liftBinary (p : BinaryPoly) (j : Fin 6) :
    pderiv j (liftBinary p) =
      if j = 2 then liftBinary (pderiv 0 p)
      else if j = 5 then liftBinary (pderiv 1 p)
      else 0 := by
  induction p using MvPolynomial.induction_on generalizing j with
  | C z =>
    simp only [liftBinary, AlgHom.toRingHom_eq_coe, RingHom.coe_coe, rename_C, pderiv_C, map_zero]
    split_ifs <;> rfl
  | add p q hp hq =>
    simp only [map_add, hp, hq]
    split_ifs <;> ring
  | mul_X p i hp =>
    rw [map_mul]
    fin_cases i <;> fin_cases j <;>
      simp [liftBinary_X0, liftBinary_X1, pderiv_X, hp]

private def atZero : Ideal (Polynomial ℂ) := Ideal.span {Polynomial.X}

private abbrev Coeff5 := MvPolynomial (Fin 5) ℂ
private def univF : Polynomial Coeff5 := MvPolynomial.finSuccEquiv ℂ 5 F6
private def cF5 : Coeff5 := C (5 / 7 : ℂ) * X 0
private def cF3 : Coeff5 := C (5 / 7 : ℂ) * X 0 * X 2 * X 3
private def cF1 : Coeff5 :=
  C (10 / 7 : ℂ) * X 0 ^ 5 + C (10 / 7 : ℂ) * X 0 * X 2 ^ 3 * X 3
private def cF0 : Coeff5 :=
  X 0 ^ 6 + X 1 ^ 6 + C (15 / 7 : ℂ) * X 1 ^ 4 * X 4 ^ 2 +
    C (10 / 7 : ℂ) * X 1 ^ 2 * X 4 ^ 4 + X 2 ^ 6 +
    C (5 / 7 : ℂ) * X 2 ^ 5 * X 3 + C (15 / 7 : ℂ) * X 2 * X 3 ^ 5 + X 3 ^ 6 + X 4 ^ 6

private lemma univF_formula :
    univF = Polynomial.X ^ 6 +
      (Polynomial.C cF5 * Polynomial.X ^ 5 + Polynomial.C cF3 * Polynomial.X ^ 3 +
        Polynomial.C cF1 * Polynomial.X + Polynomial.C cF0) := by
  have hX0 : MvPolynomial.finSuccEquiv ℂ 5 (X 0) = Polynomial.X := by
    simp [MvPolynomial.finSuccEquiv_apply]
  have hX1 : MvPolynomial.finSuccEquiv ℂ 5 (X 1) = Polynomial.C (X 0) := by
    simp [MvPolynomial.finSuccEquiv_apply]; rfl
  have hX2 : MvPolynomial.finSuccEquiv ℂ 5 (X 2) = Polynomial.C (X 1) := by
    simp [MvPolynomial.finSuccEquiv_apply]; rfl
  have hX3 : MvPolynomial.finSuccEquiv ℂ 5 (X 3) = Polynomial.C (X 2) := by
    simp [MvPolynomial.finSuccEquiv_apply]; rfl
  have hX4 : MvPolynomial.finSuccEquiv ℂ 5 (X 4) = Polynomial.C (X 3) := by
    simp [MvPolynomial.finSuccEquiv_apply]; rfl
  have hX5 : MvPolynomial.finSuccEquiv ℂ 5 (X 5) = Polynomial.C (X 4) := by
    simp [MvPolynomial.finSuccEquiv_apply]; rfl
  have hC : ∀ c : ℂ, MvPolynomial.finSuccEquiv ℂ 5 (C c) = Polynomial.C (C c) := by
    intro c; simp [MvPolynomial.finSuccEquiv_apply]
  unfold univF
  rw [targetF_eq]
  unfold core binaryB
  simp only [map_add, map_mul, map_pow, liftCore_X0, liftCore_X1, liftCore_X2, liftCore_X3,
    liftBinary_X0, liftBinary_X1]
  have hLC : ∀ c : ℂ, liftCore (C c) = C c := fun c => by simp [liftCore]
  have hLB : ∀ c : ℂ, liftBinary (C c) = C c := fun c => by simp [liftBinary]
  simp only [hLC, hLB, map_add, map_mul, map_pow, hX0, hX1, hX2, hX3, hX4, hX5, hC,
    cF5, cF3, cF1, cF0]
  ring

lemma F6_irreducible : Irreducible F6 := by
  have hX0_poly : (Polynomial.X : Polynomial ℂ) = Polynomial.X - Polynomial.C 0 := by simp
  have hatZero : atZero.IsPrime := by
    rw [atZero, hX0_poly]
    exact (Ideal.span_singleton_prime (Polynomial.X_sub_C_ne_zero 0)).mpr
      (Polynomial.prime_X_sub_C 0)
  let sec : Polynomial ℂ := (Polynomial.X + 1) ^ 6 - 1
  let P : Polynomial (Polynomial ℂ) := Polynomial.X ^ 6 + Polynomial.C sec
  have hP_irred : Irreducible P := by
    have hmon : P.Monic := Polynomial.monic_X_pow_add_C sec (by decide)
    have hdP : P.natDegree = 6 := Polynomial.natDegree_X_pow_add_C
    have hmem : sec ∈ atZero := by
      rw [atZero, Ideal.mem_span_singleton, hX0_poly, Polynomial.dvd_iff_isRoot]
      simp [sec, Polynomial.IsRoot]
    have hnot_sq : sec ∉ atZero ^ 2 := by
      intro hsq
      rw [atZero, hX0_poly, Ideal.span_singleton_pow, Ideal.mem_span_singleton] at hsq
      have hd : (Polynomial.X - Polynomial.C (0 : ℂ) : Polynomial ℂ) ∣
          Polynomial.derivative sec := by
        simpa using Polynomial.pow_sub_one_dvd_derivative_of_pow_dvd hsq
      have he := Polynomial.dvd_iff_isRoot.mp hd
      have hsec_d : Polynomial.eval 0 (Polynomial.derivative sec) = 6 := by
        have hcomp : sec = (Polynomial.X ^ 6 - 1 : Polynomial ℂ).comp (Polynomial.X + 1) := by
          simp [sec]
        rw [hcomp, Polynomial.derivative_comp, Polynomial.eval_mul, Polynomial.eval_comp]
        simp
        norm_num
      rw [Polynomial.IsRoot, hsec_d] at he
      norm_num at he
    have heis : P.IsEisensteinAt atZero := by
      refine ⟨?_, ?_, ?_⟩
      · rw [hmon.leadingCoeff, atZero, Ideal.mem_span_singleton, hX0_poly, Polynomial.dvd_iff_isRoot]
        simp [Polynomial.IsRoot]
      · intro m hm; rw [hdP] at hm
        by_cases hm0 : m = 0
        · subst hm0; simpa [P] using hmem
        · have h6m : 6 ≠ m := by omega
          simp [P, Polynomial.coeff_X_pow, Polynomial.coeff_C, hm0, Ne.symm h6m]
      · simpa [P] using hnot_sq
    exact heis.irreducible hatZero hmon.isPrimitive (by rw [hdP]; decide)
  have hmon : univF.Monic := by
    rw [univF_formula]
    have hlt : (Polynomial.C cF5 * Polynomial.X ^ 5 + Polynomial.C cF3 * Polynomial.X ^ 3 +
        Polynomial.C cF1 * Polynomial.X + Polynomial.C cF0).degree < 6 := by compute_degree!
    exact Polynomial.monic_X_pow_add hlt
  let spec5 : Coeff5 →+* Polynomial ℂ :=
    MvPolynomial.eval₂Hom Polynomial.C ![0, 0, Polynomial.X + 1, 0, Polynomial.C Complex.I]
  have hspec : univF.map spec5 = P := by
    have hi2 : (Polynomial.C Complex.I : Polynomial ℂ) ^ 2 = -1 := by
      rw [← Polynomial.C_pow, Complex.I_sq, map_neg, map_one]
    have hi6 : (Polynomial.C Complex.I : Polynomial ℂ) ^ 6 = -1 := by
      calc
        (Polynomial.C Complex.I : Polynomial ℂ) ^ 6 =
            ((Polynomial.C Complex.I : Polynomial ℂ) ^ 2) ^ 3 := by ring
        _ = (-1) ^ 3 := by rw [hi2]
        _ = -1 := by norm_num
    rw [univF_formula]
    simp only [Polynomial.map_add, Polynomial.map_mul, Polynomial.map_pow, Polynomial.map_X,
      Polynomial.map_C]
    have hs5 : spec5 cF5 = 0 := by simp [spec5, cF5]
    have hs3 : spec5 cF3 = 0 := by simp [spec5, cF3]
    have hs1 : spec5 cF1 = 0 := by simp [spec5, cF1]
    have hs0 : spec5 cF0 = sec := by
      have hv0 : spec5 (X 0) = 0 := by simp [spec5]
      have hv1 : spec5 (X 1) = 0 := by simp [spec5]
      have hv2 : spec5 (X 2) = Polynomial.X + 1 := by simp [spec5]
      have hv3 : spec5 (X 3) = 0 := by simp [spec5]
      have hv4 : spec5 (X 4) = Polynomial.C Complex.I := by simp [spec5]
      have hvC : ∀ z : ℂ, spec5 (C z) = Polynomial.C z := fun z => by simp [spec5]
      unfold cF0
      dsimp only [sec]
      simp only [map_add, map_mul, map_pow, hv0, hv1, hv2, hv3, hv4, hvC, hi6]
      ring
    rw [hs5, hs3, hs1, hs0]; simp [P]
  have hu : Irreducible univF :=
    Polynomial.Monic.irreducible_of_irreducible_map spec5 univF hmon (hspec ▸ hP_irred)
  exact (MulEquiv.irreducible_iff (MvPolynomial.finSuccEquiv ℂ 5)).mp hu

lemma not_dvd_liftCore {P : CorePoly} (hP : P ≠ 0) : ¬ (F6 ∣ liftCore P) := by
  have hex : ∃ v : Fin 4 → ℂ, eval v P ≠ 0 := by
    by_contra h
    apply hP
    apply MvPolynomial.funext
    intro v
    rw [map_zero]
    by_contra hv
    exact h ⟨v, hv⟩
  obtain ⟨v, hv⟩ := hex
  let qPoly : Polynomial ℂ := Polynomial.X ^ 6 + Polynomial.C (eval v core)
  have hq_deg : 0 < qPoly.degree := by
    rw [Polynomial.degree_X_pow_add_C (by decide)]; decide
  obtain ⟨z, hz⟩ := Complex.exists_root hq_deg
  let w : Fin 6 → ℂ := ![v 0, v 1, z, v 2, v 3, 0]
  have hev_core : ∀ R : CorePoly, eval w (liftCore R) = eval v R := by
    intro R
    rw [liftCore, AlgHom.toRingHom_eq_coe, RingHom.coe_coe, eval_rename]
    have hw : w ∘ ![0, 1, 3, 4] = v := by ext i; fin_cases i <;> rfl
    rw [hw]
  have hevF : eval w F6 = 0 := by
    rw [targetF_eq, eval_add, hev_core]
    have hevB : eval w (liftBinary binaryB) = z ^ 6 := by
      simp only [liftBinary, AlgHom.toRingHom_eq_coe, RingHom.coe_coe, eval_rename]
      simp [binaryB, w]
    rw [hevB]
    rw [Polynomial.IsRoot] at hz
    simp only [qPoly, Polynomial.eval_add, Polynomial.eval_pow, Polynomial.eval_X,
      Polynomial.eval_C] at hz
    linear_combination hz
  rintro ⟨Q, hQ⟩
  have hev := congrArg (eval w) hQ
  rw [hev_core, eval_mul, hevF, zero_mul] at hev
  exact hv hev

lemma not_dvd_liftBinary : ¬ (F6 ∣ liftBinary binaryB) := by
  let w : Fin 6 → ℂ := ![Complex.I, 0, 1, 0, 0, 0]
  have hi6 : Complex.I ^ 6 = -1 := by
    calc
      Complex.I ^ 6 = (Complex.I ^ 2) ^ 3 := by rw [← pow_mul]
      _ = -1 := by rw [Complex.I_sq]; norm_num
  have hevF : eval w F6 = 0 := by
    rw [targetF_eq, eval_add]
    simp only [liftCore, liftBinary, AlgHom.toRingHom_eq_coe, RingHom.coe_coe, eval_rename]
    simp [core, binaryB, w, hi6]
  have hevB : eval w (liftBinary binaryB) = 1 := by
    simp only [liftBinary, AlgHom.toRingHom_eq_coe, RingHom.coe_coe, eval_rename]
    simp [binaryB, w]
  rintro ⟨Q, hQ⟩
  have hev := congrArg (eval w) hQ
  rw [hevB, eval_mul, hevF, zero_mul] at hev
  exact one_ne_zero hev

/-- A four-variable certificate yields the protected catalogue certificate. -/
theorem ofCoreCertificate (D : CorePoly) (N : Fin 4 → CorePoly) (h : CoreCertificate D N) :
    HodgeComputer.HasCertificate (n := HodgeSubmission.n)
      HodgeSubmission.F HodgeSubmission.q HodgeSubmission.A := by
  rcases h with ⟨hD_ne, ⟨k, hD_hom, hN_hom⟩, hcore⟩
  let H6 : SixPoly := liftCore core
  let B6 : SixPoly := liftBinary binaryB
  let D_ : SixPoly := liftCore D
  let K : SixPoly := liftCore (coreSource * D - ∑ i, pderiv i core * N i)
  let S : SixPoly := 6 * D_ * H6 * B6
  let N6 : Fin 6 → SixPoly :=
    ![6 * H6 * B6 * liftCore (N 0),
      6 * H6 * B6 * liftCore (N 1),
      H6 * K * X 2,
      6 * H6 * B6 * liftCore (N 2),
      6 * H6 * B6 * liftCore (N 3),
      H6 * K * X 5]
  let Q6 : SixPoly := -24 * D_ * H6 * B6 * K
  have hcore_ne : core ≠ 0 := by
    intro hc; have hev := congrArg (eval ![1, 0, 0, 0]) hc; simp [core] at hev
  have hDC_ne : D * core ≠ 0 := mul_ne_zero hD_ne hcore_ne
  have hndvd_DC : ¬ (F6 ∣ D_ * H6) := by
    rw [show D_ * H6 = liftCore (D * core) from (map_mul _ D core).symm]
    exact not_dvd_liftCore hDC_ne
  have hndvd_6 : ¬ (F6 ∣ (6 : SixPoly)) := by
    rintro ⟨Q, hQ⟩
    have hev := congrArg (eval (0 : Fin 6 → ℂ)) hQ
    have hF0 : eval (0 : Fin 6 → ℂ) F6 = 0 := by
      rw [targetF_eq, eval_add]
      simp only [liftCore, liftBinary, AlgHom.toRingHom_eq_coe, RingHom.coe_coe, eval_rename]
      simp [core, binaryB]
    rw [eval_mul, hF0, zero_mul, map_ofNat] at hev
    norm_num at hev
  have hprimeF : Prime F6 :=
    UniqueFactorizationMonoid.irreducible_iff_prime.mp F6_irreducible
  have hndvd_S : ¬ (F6 ∣ S) := by
    have h1 : ¬ (F6 ∣ 6 * (D_ * H6)) := fun hd =>
      (hprimeF.dvd_or_dvd hd).elim hndvd_6 hndvd_DC
    have h2 : ¬ (F6 ∣ 6 * (D_ * H6) * B6) := fun hd =>
      (hprimeF.dvd_or_dvd hd).elim h1 not_dvd_liftBinary
    simpa [S, mul_assoc] using h2
  have hS_ne : S ≠ 0 := fun h0 => hndvd_S (h0 ▸ dvd_zero F6)
  have hrelPrime : IsRelPrime F6 S := F6_irreducible.isRelPrime_iff_not_dvd.mpr hndvd_S
  have h6_hom : (6 : SixPoly).IsHomogeneous 0 := by
    rw [← map_ofNat C 6]; exact isHomogeneous_C _ _
  have hH6_hom : H6.IsHomogeneous 6 := liftCore_hom core_hom
  have hB6_hom : B6.IsHomogeneous 6 := liftBinary_hom binaryB_hom
  have hD_hom6 : D_.IsHomogeneous (k + 1) := liftCore_hom hD_hom
  have hS_deg : S.IsHomogeneous (k + 13) := by
    have hm := ((h6_hom.mul hD_hom6).mul hH6_hom).mul hB6_hom
    exact (show 0 + (k + 1) + 6 + 6 = k + 13 by omega) ▸ hm
  have hS_tot : S.totalDegree = k + 13 := hS_deg.totalDegree hS_ne
  have hK_hom : K.IsHomogeneous (k + 7) := by
    apply liftCore_hom
    have hA_hom : coreSource.IsHomogeneous 6 :=
      ((isHomogeneous_X ℂ (0 : Fin 4)).pow 3).mul ((isHomogeneous_X ℂ (2 : Fin 4)).pow 3)
    have ht1 : (coreSource * D).IsHomogeneous (k + 7) :=
      (show 6 + (k + 1) = k + 7 by omega) ▸ hA_hom.mul hD_hom
    have ht2 : (∑ i, pderiv i core * N i).IsHomogeneous (k + 7) := by
      apply IsHomogeneous.sum
      intro i _
      exact (show (6 - 1) + (k + 2) = k + 7 by omega) ▸ (core_hom.pderiv (i := i)).mul (hN_hom i)
    exact ht1.sub ht2
  have hN6_hom : ∀ i : Fin 6, (N6 i).IsHomogeneous (k + 14) := by
    intro i
    fin_cases i
    · change (6 * H6 * B6 * liftCore (N 0)).IsHomogeneous (k + 14)
      exact (show 0 + 6 + 6 + (k + 2) = k + 14 by omega) ▸
        (((h6_hom.mul hH6_hom).mul hB6_hom).mul (liftCore_hom (hN_hom 0)))
    · change (6 * H6 * B6 * liftCore (N 1)).IsHomogeneous (k + 14)
      exact (show 0 + 6 + 6 + (k + 2) = k + 14 by omega) ▸
        (((h6_hom.mul hH6_hom).mul hB6_hom).mul (liftCore_hom (hN_hom 1)))
    · change (H6 * K * X 2).IsHomogeneous (k + 14)
      exact (show 6 + (k + 7) + 1 = k + 14 by omega) ▸
        ((hH6_hom.mul hK_hom).mul (isHomogeneous_X ℂ (2 : Fin 6)))
    · change (6 * H6 * B6 * liftCore (N 2)).IsHomogeneous (k + 14)
      exact (show 0 + 6 + 6 + (k + 2) = k + 14 by omega) ▸
        (((h6_hom.mul hH6_hom).mul hB6_hom).mul (liftCore_hom (hN_hom 2)))
    · change (6 * H6 * B6 * liftCore (N 3)).IsHomogeneous (k + 14)
      exact (show 0 + 6 + 6 + (k + 2) = k + 14 by omega) ▸
        (((h6_hom.mul hH6_hom).mul hB6_hom).mul (liftCore_hom (hN_hom 3)))
    · change (H6 * K * X 5).IsHomogeneous (k + 14)
      exact (show 6 + (k + 7) + 1 = k + 14 by omega) ▸
        ((hH6_hom.mul hK_hom).mul (isHomogeneous_X ℂ (5 : Fin 6)))
  have hEB : liftBinary (pderiv 0 binaryB) * X 2 +
      liftBinary (pderiv 1 binaryB) * X 5 = 6 * B6 := by
    have hE := congrArg liftBinary binaryB_hom.sum_X_mul_pderiv
    simp only [map_add, map_mul, Fin.sum_univ_two, liftBinary_X0,
      liftBinary_X1, nsmul_eq_mul, map_natCast] at hE
    change _ = 6 * liftBinary binaryB
    linear_combination hE
  have hK_def : liftCore (∑ i, pderiv i core * N i) =
      liftCore coreSource * D_ - K := by
    dsimp [K, D_]; rw [map_sub, map_mul]; ring
  have hdir_core : ∀ p : CorePoly,
      liftCore (pderiv 0 p) * liftCore (N 0) +
      liftCore (pderiv 1 p) * liftCore (N 1) +
      liftCore (pderiv 2 p) * liftCore (N 2) +
      liftCore (pderiv 3 p) * liftCore (N 3) =
      liftCore (∑ i, pderiv i p * N i) := by
    intro p; simp only [Fin.sum_univ_four, map_add, map_mul]
  have hpd6 : ∀ j : Fin 6, pderiv j (6 : SixPoly) = 0 := fun j => by
    rw [← map_ofNat C 6, pderiv_C]
  have hN0 : N6 0 = 6 * H6 * B6 * liftCore (N 0) := rfl
  have hN1 : N6 1 = 6 * H6 * B6 * liftCore (N 1) := rfl
  have hN2 : N6 2 = H6 * K * X 2 := rfl
  have hN3 : N6 3 = 6 * H6 * B6 * liftCore (N 2) := rfl
  have hN4 : N6 4 = 6 * H6 * B6 * liftCore (N 3) := rfl
  have hN5 : N6 5 = H6 * K * X 5 := rfl
  have hpK2 : pderiv 2 K = 0 := by dsimp [K]; rw [pderiv_liftCore]; simp
  have hpK5 : pderiv 5 K = 0 := by dsimp [K]; rw [pderiv_liftCore]; simp
  have hdir6 : (∑ i : Fin 6, pderiv i F6 * N6 i) = A6 * S := by
    rw [targetF_eq, targetA_eq]
    simp only [Fin.sum_univ_six, hN0, hN1, hN2, hN3, hN4, hN5]
    dsimp only [S]
    simp only [map_add, pderiv_liftCore, pderiv_liftBinary, Fin.reduceEq, ite_true, ite_false,
      add_zero, zero_add]
    have hdH := hdir_core core
    linear_combination 6 * H6 * B6 * hdH + 6 * H6 * B6 * hK_def + H6 * K * hEB
  have hdiv6 : S * (∑ i : Fin 6, pderiv i (N6 i)) - ∑ i : Fin 6, pderiv i S * N6 i = F6 * Q6 := by
    have hcore_lift : 3 * H6 *
        (D_ * (liftCore (pderiv 0 (N 0)) +
               liftCore (pderiv 1 (N 1)) +
               liftCore (pderiv 2 (N 2)) +
               liftCore (pderiv 3 (N 3))) -
         (liftCore (pderiv 0 D) * liftCore (N 0) +
          liftCore (pderiv 1 D) * liftCore (N 1) +
          liftCore (pderiv 2 D) * liftCore (N 2) +
          liftCore (pderiv 3 D) * liftCore (N 3))) =
        -2 * D_ * K := by
      have h0 := congrArg liftCore hcore
      simp only [Fin.sum_univ_four, map_zero, map_sub, map_add, map_mul, map_pow, map_ofNat] at h0
      dsimp [H6, D_, K]
      simp only [Fin.sum_univ_four, map_sub, map_add, map_mul]
      linear_combination h0
    rw [targetF_eq]
    simp only [Fin.sum_univ_six, hN0, hN1, hN2, hN3, hN4, hN5]
    dsimp only [S, Q6, H6, B6, D_]
    simp only [pderiv_mul, hpd6, hpK2, hpK5, pderiv_liftCore, pderiv_liftBinary, pderiv_X_self,
      Fin.reduceEq, ite_true, ite_false]
    change (6 * liftCore D * liftCore core * liftBinary binaryB) * _ - _ =
      (liftCore core + liftBinary binaryB) *
        (-24 * liftCore D * liftCore core * liftBinary binaryB * K)
    dsimp only [H6, B6, D_] at hcore_lift hEB
    linear_combination
      12 * liftCore core * liftBinary binaryB ^ 2 * hcore_lift -
      6 * liftCore D * liftCore core ^ 2 * K * hEB
  apply HodgeComputer.Verify.hasCertificate (d := 6) (D := S) (N := N6) (Q := Q6)
  refine ⟨F6_irreducible.ne_zero, ?_, ?_, hS_ne, hS_tot ▸ hS_deg, ?_, hrelPrime, ?_⟩
  · change F6.IsHomogeneous 6
    rw [targetF_eq]
    exact hH6_hom.add hB6_hom
  · change Homogeneous A6 6
    rw [targetA_eq]
    apply Homogeneous.of_isHomogeneous (k := 6)
    apply liftCore_hom
    exact ((isHomogeneous_X ℂ (0 : Fin 4)).pow 3).mul ((isHomogeneous_X ℂ (2 : Fin 4)).pow 3)
  · intro i
    have hi := Homogeneous.of_isHomogeneous (hN6_hom i)
    convert hi using 1
    simp [hS_tot, HodgeSubmission.q]
    omega
  · apply normalized_identity
    · simpa only [HodgeSubmission.q, Nat.cast_one, one_mul, HodgeTestbed.directional,
        F6, A6] using hdir6
    · simpa only [HodgeSubmission.q, pow_one, HodgeTestbed.divergenceNumerator,
        HodgeTestbed.divergence, HodgeTestbed.directional, F6] using hdiv6

end FourPlusTwoCertificate
