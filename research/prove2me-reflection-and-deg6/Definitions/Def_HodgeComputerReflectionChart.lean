import Definitions.Def_DiamondCertificate

set_option autoImplicit false
noncomputable section

namespace HodgeComputer.ReflectionChart
open MvPolynomial
open scoped BigOperators

/-- Five-variable chart polynomial ring in coordinates `(λ, b, c, d, u)` or `(λ, b, c, d, v)`. -/
abbrev ChartPoly := MvPolynomial (Fin 5) ℂ

/-- Reflection cover substitution `u ↦ v²` on the last coordinate (`Fin 5` index `4`). -/
def squareLift : ChartPoly →+* ChartPoly :=
  eval₂Hom C (fun i : Fin 5 => if i = 4 then X i ^ 2 else X i)

/-- Transformed double-pole source numerator `4 v P(λ, b, c, d, v²)` on the reflection chart. -/
def reflectedSource (P : ChartPoly) : ChartPoly :=
  C 4 * X 4 * squareLift P

/-- Pullback of a six-variable polynomial along the affine reflection chart
`(X₀, X₁, X₂, Y₀, Y₁, Y₂) = (1 + v, 1 - v, b, c + λ v, c - λ v, d)`. -/
def affinePullback : AmbientPoly 5 →+* ChartPoly :=
  eval₂Hom C ![1 + X 4, 1 - X 4, X 1, X 2 + X 0 * X 4, X 2 - X 0 * X 4, X 3]

/-- Canonical smooth test point `(0, 0, 0, 0, i)` on the affine reflection chart. -/
def normalPoint : Fin 5 → ℂ := ![0, 0, 0, 0, Complex.I]

/-- Cleared relative-primitive polynomial defect on the weighted cubic threefold `f(λ, b, c, d, u) = 0`:
`f_u² ∑_{j ∈ {b,c,d}} [f_u (D ∂_j C_j - C_j ∂_j D) - f_j (D ∂_u C_j - C_j ∂_u D)] - 2 D² (f_u P_u - f_{uu} P)`. -/
def cubicResidueDefect (f P D : ChartPoly) (Cvec : Fin 3 → ChartPoly) : ChartPoly :=
  (pderiv 4 f) ^ 2 *
    (∑ j : Fin 3,
      (pderiv 4 f * (D * pderiv j.succ.castSucc (Cvec j) - Cvec j * pderiv j.succ.castSucc D) -
        pderiv j.succ.castSucc f * (D * pderiv 4 (Cvec j) - Cvec j * pderiv 4 D))) -
    C 2 * D ^ 2 * (pderiv 4 f * pderiv 4 P - pderiv 4 (pderiv 4 f) * P)

/-- Existence of a cleared relative two-form primitive `β_j = C_j / D` (`j ∈ {b, c, d}`)
on the weighted threefold `f(λ, b, c, d, u) = 0` for the residue density `2(f_u P_u - f_{uu} P) / f_u³`. -/
def RelativePrimitive (f P : ChartPoly) : Prop :=
  ∃ (D : ChartPoly) (Cvec : Fin 3 → ChartPoly) (Q : ChartPoly),
    D ≠ 0 ∧ IsRelPrime f D ∧ cubicResidueDefect f P D Cvec = f * Q

/-- Cleared constant-in-`u` canonical trace coefficient `K₀` for a cubic equation
`f = ℓ u³ + c₂ u² + c₁ u + c₀` and cubic source `P = p₀ + p₁ u + p₂ u² + p₃ u³`. -/
def cubicTraceK0 (ell c2 c1 c0 p0 p2 p3 : ChartPoly) : ChartPoly :=
  C 6 * ell ^ 2 * p0 - C 2 * ell * c1 * p2 + (C 3 * ell * c0 + c2 * c1) * p3

/-- Cleared linear-in-`u` canonical trace coefficient `K₁` for a cubic equation
`f = ℓ u³ + c₂ u² + c₁ u + c₀` and cubic source `P = p₀ + p₁ u + p₂ u² + p₃ u³`. -/
def cubicTraceK1 (ell c2 p1 p2 p3 : ChartPoly) : ChartPoly :=
  C 6 * ell ^ 2 * p1 - C 4 * ell * c2 * p2 + C 2 * c2 ^ 2 * p3

end HodgeComputer.ReflectionChart
