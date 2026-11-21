import Definitions.Def_HodgeComputerSymmetricSexticResidual
import Definitions.Def_HodgeComputerReflectionChart

set_option autoImplicit false
noncomputable section

namespace HodgeComputer.SymmetricSexticResidual
open MvPolynomial HodgeComputer.ReflectionChart

/-- Left quadric factor `1 + 2 b - u` on the reflection chart `(λ, b, c, d, u)`. -/
def leftQuadric : ChartPoly := 1 + C 2 * X 1 - X 4

/-- Right quadric factor `c² + 2 c d - λ² u` on the reflection chart `(λ, b, c, d, u)`. -/
def rightQuadric : ChartPoly := X 2 ^ 2 + C 2 * X 2 * X 3 - X 0 ^ 2 * X 4

/-- Weighted cubic equation `f_sym(λ, b, c, d, u)` of `symmetric-sextic-residual`
on the reflection chart `(1 + v, 1 - v, b, c + λ v, c - λ v, d)`, `u = v²`. -/
def symWeighted : ChartPoly :=
  C 2 + X 1 ^ 6 + C 2 * X 2 ^ 6 + X 3 ^ 6 +
    C 30 * X 4 + C 30 * X 4 ^ 2 + C 2 * X 4 ^ 3 +
    C 30 * X 0 ^ 2 * X 2 ^ 4 * X 4 +
    C 30 * X 0 ^ 4 * X 2 ^ 2 * X 4 ^ 2 + C 2 * X 0 ^ 6 * X 4 ^ 3 +
    C (5 / 7 : ℂ) * leftQuadric ^ 2 * rightQuadric +
    C (10 / 7 : ℂ) * leftQuadric * rightQuadric ^ 2

/-- Descended quadratic source `P_sym(λ, b, c, d, u) = b d (1 - u) (c² - λ² u)`
of `symmetric-sextic-residual` on the reflection chart. -/
def symSource : ChartPoly :=
  X 1 * X 3 * (1 - X 4) * (X 2 ^ 2 - X 0 ^ 2 * X 4)

/-- Leading `u³` coefficient `ℓ_sym(λ)` of `symWeighted`. -/
def symEll : ChartPoly :=
  C 2 + C 2 * X 0 ^ 6 - C (5 / 7 : ℂ) * X 0 ^ 2 - C (10 / 7 : ℂ) * X 0 ^ 4

/-- Quadratic `u²` coefficient `B_sym(λ, b, c, d)` of `symWeighted`. -/
def symB : ChartPoly :=
  let L0 : ChartPoly := 1 + C 2 * X 1
  let R0 : ChartPoly := X 2 ^ 2 + C 2 * X 2 * X 3
  C 30 + C 30 * X 0 ^ 4 * X 2 ^ 2 +
    C (5 / 7 : ℂ) * (R0 + C 2 * X 0 ^ 2 * L0) +
    C (10 / 7 : ℂ) * (X 0 ^ 4 * L0 + C 2 * X 0 ^ 2 * R0)

/-- Linear `u` coefficient `C_sym(λ, b, c, d)` of `symWeighted`. -/
def symC : ChartPoly :=
  let L0 : ChartPoly := 1 + C 2 * X 1
  let R0 : ChartPoly := X 2 ^ 2 + C 2 * X 2 * X 3
  C 30 + C 30 * X 0 ^ 2 * X 2 ^ 4 +
    C (5 / 7 : ℂ) * (-X 0 ^ 2 * L0 ^ 2 - C 2 * L0 * R0) +
    C (10 / 7 : ℂ) * (-C 2 * X 0 ^ 2 * L0 * R0 - R0 ^ 2)

/-- Constant coefficient `D_sym(b, c, d)` of `symWeighted`. -/
def symD : ChartPoly :=
  let L0 : ChartPoly := 1 + C 2 * X 1
  let R0 : ChartPoly := X 2 ^ 2 + C 2 * X 2 * X 3
  C 2 + X 1 ^ 6 + C 2 * X 2 ^ 6 + X 3 ^ 6 +
    C (5 / 7 : ℂ) * L0 ^ 2 * R0 + C (10 / 7 : ℂ) * L0 * R0 ^ 2

/-- Cleared canonical linear-in-`u` trace source `K₀ + K₁ u` of `symmetric-sextic-residual`. -/
def symCanonicalSource : ChartPoly :=
  cubicTraceK0 symEll symB symC symD
    (X 1 * X 3 * X 2 ^ 2) (X 1 * X 3 * X 0 ^ 2) 0 +
  cubicTraceK1 symEll symB
    (-X 1 * X 3 * (X 2 ^ 2 + X 0 ^ 2)) (X 1 * X 3 * X 0 ^ 2) 0 * X 4

end HodgeComputer.SymmetricSexticResidual
