import Theorems.Thm_HodgeComputer_FourPlusTwoResidual_coreCertificate_iff_bipartite
import Theorems.Thm_HodgeComputer_FourPlusTwoResidual_bipartiteCoreCertificate

set_option autoImplicit false
noncomputable section

theorem solution :
    ∃ (D : MvPolynomial (Fin 4) ℂ) (N : Fin 4 → MvPolynomial (Fin 4) ℂ),
      HodgeComputer.FourPlusTwoResidual.CoreCertificate D N :=
  HodgeComputer.FourPlusTwoResidual.coreCertificate_iff_bipartite.mpr
    HodgeComputer.FourPlusTwoResidual.bipartiteCoreCertificate
