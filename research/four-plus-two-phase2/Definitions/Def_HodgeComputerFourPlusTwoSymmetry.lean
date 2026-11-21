import Definitions.Def_HodgeComputerFourPlusTwoCore

set_option autoImplicit false
noncomputable section

namespace HodgeComputer.FourPlusTwoResidual

open MvPolynomial
open scoped BigOperators

/-- Sign involution `(x₀, x₁, x₂, x₃) ↦ (-x₀, -x₁, x₂, x₃)` on the four-variable core ring. -/
def flipAB : MvPolynomial (Fin 4) ℂ →ₐ[ℂ] MvPolynomial (Fin 4) ℂ :=
  aeval ![-X 0, -X 1, X 2, X 3]

/-- Sign involution `(x₀, x₁, x₂, x₃) ↦ (x₀, x₁, -x₂, -x₃)` on the four-variable core ring. -/
def flipCD : MvPolynomial (Fin 4) ℂ →ₐ[ℂ] MvPolynomial (Fin 4) ℂ :=
  aeval ![X 0, X 1, -X 2, -X 3]

/-- A strict divergence-free core certificate: `D ≠ 0` homogeneous of degree `k + 1` and
`N` homogeneous of degree `k + 2` satisfying both `N · ∇H = A D` and `D div(N) - N · ∇D = 0`. -/
def StrictCoreCertificate (D : MvPolynomial (Fin 4) ℂ)
    (N : Fin 4 → MvPolynomial (Fin 4) ℂ) : Prop :=
  D ≠ 0 ∧ (∃ k : ℕ, D.IsHomogeneous (k + 1) ∧ ∀ i, (N i).IsHomogeneous (k + 2)) ∧
    (∑ i, pderiv i core * N i = coreSource * D) ∧
    (D * (∑ i, pderiv i (N i)) - ∑ i, pderiv i D * N i = 0)

/-- An even-degree `(ℤ/2 × ℤ/2)`-equivariant core certificate: `D` has even positive degree `2m + 2`
and is invariant under both `flipAB` and `flipCD`, while `(N₀, N₁)` is `flipAB`-even and `(N₂, N₃)`
is `flipAB`-odd, satisfying the exact four-variable `CoreCertificate` differential equation. -/
def BipartiteCoreCertificate (D : MvPolynomial (Fin 4) ℂ)
    (N : Fin 4 → MvPolynomial (Fin 4) ℂ) : Prop :=
  D ≠ 0 ∧ (∃ m : ℕ, D.IsHomogeneous (2 * m + 2) ∧ ∀ i, (N i).IsHomogeneous (2 * m + 3)) ∧
    flipAB D = D ∧ flipCD D = D ∧
    flipAB (N 0) = N 0 ∧ flipAB (N 1) = N 1 ∧
    flipAB (N 2) = -N 2 ∧ flipAB (N 3) = -N 3 ∧
    2 * coreSource * D ^ 2 - 2 * D * (∑ i, pderiv i core * N i)
      + 3 * core * (D * (∑ i, pderiv i (N i)) - ∑ i, pderiv i D * N i) = 0

end HodgeComputer.FourPlusTwoResidual
