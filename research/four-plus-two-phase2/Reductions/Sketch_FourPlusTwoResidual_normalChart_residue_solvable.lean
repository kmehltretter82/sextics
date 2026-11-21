import Theorems.Thm_HodgeComputer_FourPlusTwoResidual_normalChart_residue_of_coreCertificate
import Theorems.Thm_HodgeComputer_FourPlusTwoResidual_coreCertificate_iff_bipartite
import Theorems.Thm_HodgeComputer_FourPlusTwoResidual_bipartiteCoreCertificate

set_option autoImplicit false
noncomputable section

open HodgeComputer HodgeComputer.FourPlusTwoResidual

theorem solution :
    ∃ (S : AmbientPoly n) (W : Fin 5 → AmbientPoly n) (Q : AmbientPoly n),
      S ≠ 0 ∧ S.IsHomogeneous S.totalDegree ∧ IsRelPrime F S ∧
      (∀ j, (W j).IsHomogeneous (S.totalDegree + 1)) ∧
      NormalChart.residue F A S W = F * Q := by
  obtain ⟨D, N, hDN⟩ :=
    HodgeComputer.FourPlusTwoResidual.coreCertificate_iff_bipartite.mpr
      HodgeComputer.FourPlusTwoResidual.bipartiteCoreCertificate
  exact HodgeComputer.FourPlusTwoResidual.normalChart_residue_of_coreCertificate D N hDN
