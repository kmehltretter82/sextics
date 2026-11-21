import Definitions.Def_HodgeComputerFourPlusTwoSymmetry

set_option autoImplicit false
noncomputable section

namespace HodgeComputer.FourPlusTwoResidual

open MvPolynomial
open scoped BigOperators

/-- Cleared four-variable core residual `F_H(S, D, N)` for an arbitrary sextic source `S`:
`2 S D² - 2 D (N · ∇H) + 3 H (D div(N) - N · ∇D)`. -/
def coreResidual (S D : MvPolynomial (Fin 4) ℂ)
    (N : Fin 4 → MvPolynomial (Fin 4) ℂ) : MvPolynomial (Fin 4) ℂ :=
  2 * S * D ^ 2 - 2 * D * (∑ i, pderiv i core * N i) +
    3 * core * (D * (∑ i, pderiv i (N i)) - ∑ i, pderiv i D * N i)

/-- Generalized core certificate for an arbitrary source `S`, allowing `deg D = d ≥ 0`
(including `D = 1` with linear numerators `deg Nᵢ = 1`). -/
def CoreCertificateFor (S D : MvPolynomial (Fin 4) ℂ)
    (N : Fin 4 → MvPolynomial (Fin 4) ℂ) : Prop :=
  D ≠ 0 ∧ (∃ d : ℕ, D.IsHomogeneous d ∧ ∀ i, (N i).IsHomogeneous (d + 1)) ∧
    coreResidual S D N = 0

/-- A sextic source `S` is polynomial-exact on the core if `S = L · ∇H - (3/2) H div(L)`
for a homogeneous linear vector field `L` (`deg Lᵢ = 1`), i.e. the `D = 1` baseline. -/
def PolyExactCoreSource (S : MvPolynomial (Fin 4) ℂ) : Prop :=
  ∃ L : Fin 4 → MvPolynomial (Fin 4) ℂ,
    (∀ i, (L i).IsHomogeneous 1) ∧
      2 * S = 2 * (∑ i, pderiv i core * L i) - 3 * core * (∑ i, pderiv i (L i))

/-- Lower-degree linear pair `(N, T)` equivalent to `CoreCertificateFor S D N` when `H ∤ D`,
with equations of degrees `2d` and `d + 6` instead of `2d + 6`. -/
def CoprimeLinearCoreSystem (S D : MvPolynomial (Fin 4) ℂ)
    (N : Fin 4 → MvPolynomial (Fin 4) ℂ) (T : MvPolynomial (Fin 4) ℂ) : Prop :=
  (∑ i, pderiv i D * N i = D * T) ∧
    (2 * S * D - 2 * (∑ i, pderiv i core * N i) +
      3 * core * ((∑ i, pderiv i (N i)) - T) = 0)

end HodgeComputer.FourPlusTwoResidual
