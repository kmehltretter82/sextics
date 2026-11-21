import Definitions.Def_HodgeComputerDiagonalPermutation
import Definitions.Def_HodgeComputerReflectionChart

set_option autoImplicit false
noncomputable section

namespace HodgeComputer.DiagonalPermutation
open MvPolynomial HodgeComputer.ReflectionChart

/-- The `(01)`-transposition-invariant degree-six source
`B₂ = X₀³ Y₀³ + X₁³ Y₁³ - 2 X₂³ Y₂³` (`benediamond`, Discussion #177). -/
def transpositionSource : HodgeComputer.AmbientPoly 5 :=
  X 0 ^ 3 * X 3 ^ 3 + X 1 ^ 3 * X 4 ^ 3 - C 2 * X 2 ^ 3 * X 5 ^ 3

/-- Leading coefficient `ℓ(λ) = 2 + (10/7) λ² + 2 λ⁶` of the `(01)`-reflection cubic. -/
def diagEll : ChartPoly :=
  C 2 + C (10 / 7 : ℂ) * X 0 ^ 2 + C 2 * X 0 ^ 6

/-- Quadratic coefficient `B(λ, b, c, d)` of the `(01)`-reflection cubic. -/
def diagB : ChartPoly :=
  C 30 - C (10 / 11 : ℂ) * X 3 ^ 2 - C (10 / 7 : ℂ) * X 2 ^ 2 +
    C (10 / 7 : ℂ) * X 0 * X 1 * X 3 - C (20 / 11 : ℂ) * X 0 * X 1 * X 2 +
    C (30 / 13 : ℂ) * X 0 ^ 2 * X 2 * X 3 + C (30 / 11 : ℂ) * X 0 ^ 2 * X 1 -
    C (20 / 13 : ℂ) * X 0 ^ 3 * X 3 + C 30 * X 0 ^ 4 * X 2 ^ 2 -
    C (10 / 13 : ℂ) * X 0 ^ 4 * X 1 ^ 2

/-- Linear coefficient `C(λ, b, c, d)` of the `(01)`-reflection cubic. -/
def diagC : ChartPoly :=
  C 30 + C (10 / 13 : ℂ) * X 2 * X 3 ^ 3 + C (10 / 13 : ℂ) * X 2 ^ 3 * X 3 +
    C (30 / 7 : ℂ) * X 1 * X 2 * X 3 + C (30 / 11 : ℂ) * X 1 * X 2 ^ 2 -
    C (20 / 13 : ℂ) * X 0 * X 3 ^ 3 - C (60 / 13 : ℂ) * X 0 * X 2 ^ 2 * X 3 +
    C (30 / 7 : ℂ) * X 0 * X 1 * X 3 - C (60 / 11 : ℂ) * X 0 * X 1 * X 2 +
    C (10 / 7 : ℂ) * X 0 * X 1 ^ 3 * X 3 - C (20 / 11 : ℂ) * X 0 * X 1 ^ 3 * X 2 -
    C (10 / 7 : ℂ) * X 0 ^ 2 + C (30 / 13 : ℂ) * X 0 ^ 2 * X 2 * X 3 +
    C 30 * X 0 ^ 2 * X 2 ^ 4 + C (10 / 11 : ℂ) * X 0 ^ 2 * X 1 +
    C (10 / 11 : ℂ) * X 0 ^ 2 * X 1 ^ 3

/-- Constant coefficient `D(b, c, d)` of the `(01)`-reflection cubic. -/
def diagD : ChartPoly :=
  C 2 + C (10 / 11 : ℂ) * X 3 ^ 2 + X 3 ^ 6 + C (10 / 13 : ℂ) * X 2 * X 3 ^ 3 +
    C (10 / 7 : ℂ) * X 2 ^ 2 + C (10 / 13 : ℂ) * X 2 ^ 3 * X 3 + C 2 * X 2 ^ 6 +
    C (10 / 7 : ℂ) * X 1 * X 2 * X 3 + C (10 / 11 : ℂ) * X 1 * X 2 ^ 2 +
    C (10 / 13 : ℂ) * X 1 ^ 2 * X 2 ^ 4 + C (10 / 7 : ℂ) * X 1 ^ 3 * X 2 * X 3 +
    C (10 / 11 : ℂ) * X 1 ^ 3 * X 2 ^ 2 + X 1 ^ 6

/-- The weighted cubic equation `f_diag(λ, b, c, d, u) = ℓ u³ + B u² + C u + D`
of `diagonal-permutation` on the `(01)`-reflection chart. -/
def diagWeighted : ChartPoly :=
  diagEll * X 4 ^ 3 + diagB * X 4 ^ 2 + diagC * X 4 + diagD

/-- The descended cubic source `P_diag(λ, b, c, d, u)` of `transpositionSource`
on the `(01)`-reflection chart. -/
def diagSource : ChartPoly :=
  C 2 * X 0 ^ 3 * X 4 ^ 3 +
    (C 6 * X 0 * X 2 ^ 2 + C 18 * X 0 ^ 2 * X 2 + C 6 * X 0 ^ 3) * X 4 ^ 2 +
    (C 6 * X 2 ^ 3 + C 18 * X 0 * X 2 ^ 2 + C 6 * X 0 ^ 2 * X 2) * X 4 +
    (C 2 * X 2 ^ 3 - C 2 * X 1 ^ 3 * X 3 ^ 3)

/-- Cleared canonical linear-in-`u` trace source `K₀ + K₁ u` of `diagonal-permutation`,
obtained after subtracting the first-order pole `p₃ / ℓ` and quadratic branch trace modulo `f_diag`. -/
def diagCanonicalSource : ChartPoly :=
  cubicTraceK0 diagEll diagB diagC diagD
    (C 2 * X 2 ^ 3 - C 2 * X 1 ^ 3 * X 3 ^ 3)
    (C 6 * X 0 * X 2 ^ 2 + C 18 * X 0 ^ 2 * X 2 + C 6 * X 0 ^ 3)
    (C 2 * X 0 ^ 3) +
  cubicTraceK1 diagEll diagB
    (C 6 * X 2 ^ 3 + C 18 * X 0 * X 2 ^ 2 + C 6 * X 0 ^ 2 * X 2)
    (C 6 * X 0 * X 2 ^ 2 + C 18 * X 0 ^ 2 * X 2 + C 6 * X 0 ^ 3)
    (C 2 * X 0 ^ 3) * X 4

end HodgeComputer.DiagonalPermutation
