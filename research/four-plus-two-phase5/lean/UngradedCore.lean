import FourPlusTwoCertificate
import HodgeComputer.Constructions.HoldoutCover.HomogeneousField

/-! Finite homogeneous extraction for the twisted four-variable current. -/
noncomputable section
namespace FourPlusTwoCertificate
open MvPolynomial HoldoutConstruction
open scoped BigOperators

private def defect (D : CorePoly) (N : Fin 4 → CorePoly) : CorePoly :=
  2 • (coreSource * D^2) - 2 • (D * polynomialDirectional core N) +
    3 • (core * polynomialDivNum D N)

private lemma defect_scale (D S : CorePoly) (N : Fin 4 → CorePoly) :
    defect (S*D) (fun i => S*N i) = S^2 * defect D N := by
  simp only [defect, polynomialDivNum, polynomialDirectional_scale,
    polynomialDivergence_scale, polynomialDirectional_mul, nsmul_eq_mul]
  ring

private lemma defect_change (T U : CorePoly) (M : Fin 4 → CorePoly) :
    defect (T-U) M = defect T M +
      (2 • (coreSource * U * (U-2*T)) + 2 • (U * polynomialDirectional core M) +
        3 • (core * (polynomialDirectional U M - U*polynomialDivergence M))) := by
  have hd : polynomialDirectional (T-U) M = polynomialDirectional T M - polynomialDirectional U M := by
    simp only [polynomialDirectional, map_sub, mul_sub, Finset.sum_sub_distrib]
  simp only [defect, polynomialDivNum, hd, nsmul_eq_mul]
  ring

private lemma scaled_order {n : ℕ} {Q : CorePoly} (h : HasOrder n Q) (c : ℕ) :
    HasOrder n (c • Q) := by
  rw [nsmul_eq_mul]
  simpa using (hasOrder_zero_degree (c : CorePoly)).mul h

private lemma directional_order {n w : ℕ} {P : CorePoly} {M : Fin 4 → CorePoly}
    (hP : HasOrder n P) (hM : ∀ i, HasOrder w (M i)) :
    HasOrder (w+(n-1)) (polynomialDirectional P M) :=
  HasOrder.sum _ _ fun i _ => (hM i).mul (hP.pderiv i)

private lemma inverse_order (D : CorePoly) (d p : ℕ) (hD : HasOrder d D) :
    HasOrder ((p-1)*d) (finiteInverse D (homogeneousComponent d D) p) := by
  let h := homogeneousComponent d D
  have hh : HasOrder d h := isHomogeneous_hasOrder (homogeneousComponent_isHomogeneous d D)
  have hR : HasOrder d (h-D) := hh.sub hD
  apply HasOrder.sum
  intro i hi
  have hip := Finset.mem_range.mp hi
  have hsum : i + (p-1-i) = p-1 := by omega
  simpa only [h, ← Nat.add_mul, hsum] using (hh.pow i).mul (hR.pow (p-1-i))

private lemma least_component (D : CorePoly) (hne : D ≠ 0) :
    ∃ d, HasOrder d D ∧ homogeneousComponent d D ≠ 0 := by
  classical
  have hex : ∃ n, homogeneousComponent n D ≠ 0 := by
    by_contra! hall
    apply hne
    rw [← sum_homogeneousComponent D]
    exact Finset.sum_eq_zero (fun k _ => hall k)
  exact ⟨Nat.find hex, hasOrder_of_components_zero (fun k hk => not_ne_iff.mp (Nat.find_min hex hk)), Nat.find_spec hex⟩

private lemma source_hom : coreSource.IsHomogeneous 6 :=
  ((isHomogeneous_X ℂ (0 : Fin 4)).pow 3).mul ((isHomogeneous_X ℂ (2 : Fin 4)).pow 3)

private lemma defect_component (T : CorePoly) (M : Fin 4 → CorePoly) (e : ℕ)
    (he : 0 < e) (hT : T.IsHomogeneous e) :
    homogeneousComponent (2*e+6) (defect T M) =
      defect T (fun i => homogeneousComponent (e+1) (M i)) := by
  have hA : (coreSource*T^2).IsHomogeneous (2*e+6) := by
    convert source_hom.mul (hT.pow 2) using 1
    omega
  have hdir : homogeneousComponent (2*e+6) (T * polynomialDirectional core M) =
      T * polynomialDirectional core (fun i => homogeneousComponent (e+1) (M i)) := by
    rw [show 2*e+6=e+(5+(e+1)) by omega, component_mul_homogeneous hT,
      component_polynomialDirectional (show core.IsHomogeneous (5+1) from core_hom)]
  have hdiv : homogeneousComponent (2*e+6) (core * polynomialDivNum T M) =
      core * polynomialDivNum T (fun i => homogeneousComponent (e+1) (M i)) := by
    rw [show 2*e+6=6+2*e by omega, component_mul_homogeneous core_hom,
      component_polynomialDivNum he hT]
  simp only [defect, map_add, map_sub, map_nsmul, homogeneousComponent_eq_self hA, hdir, hdiv]

/-- Project the rescaled polynomial identity after its error has moved above
homogeneous degree `2*(p*d)+6`, where `p=d+2`. -/
private theorem project_initial (D : CorePoly) (N : Fin 4 → CorePoly) (d : ℕ)
    (hd : 0 < d) (hD : HasOrder d D) (h0 : homogeneousComponent d D ≠ 0)
    (hc : defect D N = 0) : ∃ T M, CoreCertificate T M := by
  let h := homogeneousComponent d D
  let p := d+2
  let e := p*d
  let u := p*(d+1)
  let w := (p-1)*d
  let T := h^p
  let U := (h-D)^p
  let S := finiteInverse D h p
  let M : Fin 4 → CorePoly := fun i => S*N i
  have hp : 0 < p := by dsimp [p]; omega
  have he : 0 < e := Nat.mul_pos hp hd
  have hu : u = e + (d+2) := by simp only [u, e, p, Nat.mul_add, Nat.mul_one]
  have hw : w+d=e := by
    dsimp only [w, e, p]
    rw [show d+2-1 = d+1 by omega, ← Nat.succ_mul]
  have hh := homogeneousComponent_isHomogeneous d D
  have hT : T.IsHomogeneous e := by simpa only [T, h, e, Nat.mul_comm] using hh.pow p
  have tOrd := isHomogeneous_hasOrder hT
  have uOrd : HasOrder u U := (order_sub_initial hD).pow p
  have mOrd (i : Fin 4) : HasOrder w (M i) := by
    simpa only [Nat.add_zero, M, S, h, w] using (inverse_order D d p hD).mul (hasOrder_zero_degree (N i))
  have divOrd : HasOrder (w-1) (polynomialDivergence M) :=
    HasOrder.sum _ _ (fun i _ => (mOrd i).pderiv i)
  have hErr1 : HasOrder (2*e+7) (coreSource*U*(U-2*T)) := by
    have hUT : HasOrder e (U-2*T) := by
      apply HasOrder.sub (uOrd.mono (by omega))
      simpa only [two_nsmul, two_mul] using tOrd.add tOrd
    apply ((isHomogeneous_hasOrder source_hom).mul uOrd |>.mul hUT).mono
    omega
  have hErr2 : HasOrder (2*e+7) (U*polynomialDirectional core M) := by
    apply (uOrd.mul (directional_order (isHomogeneous_hasOrder core_hom) mOrd)).mono
    omega
  have hErr3 : HasOrder (2*e+1) (polynomialDirectional U M) := by
    apply (directional_order uOrd mOrd).mono
    omega
  have hErr4 : HasOrder (2*e+1) (U*polynomialDivergence M) := by
    apply (uOrd.mul divOrd).mono
    omega
  have hErr5 : HasOrder (2*e+7) (core*(polynomialDirectional U M-U*polynomialDivergence M)) :=
    ((isHomogeneous_hasOrder core_hom).mul (hErr3.sub hErr4)).mono (by omega)
  have hz := (((scaled_order hErr1 2).add (scaled_order hErr2 2)).add (scaled_order hErr5 3)).component_eq_zero
    (k := 2*e+6) (by omega)
  have hscaled : defect (T-U) M = 0 := by
    have hds : D*S = T-U := finiteInverse_identity D h p
    rw [← hds, mul_comm D S, show M = (fun i => S*N i) from rfl, defect_scale, hc, mul_zero]
  rw [defect_change] at hscaled
  have hproj := congrArg (homogeneousComponent (2*e+6)) hscaled
  rw [map_add, hz, add_zero, map_zero, defect_component T M e he hT] at hproj
  refine ⟨T, (fun i => homogeneousComponent (e+1) (M i)), pow_ne_zero p h0, ⟨e-1, ?_, ?_⟩, ?_⟩
  · simpa only [show e-1+1=e by omega] using hT
  · intro i
    simpa only [show e-1+2=e+1 by omega] using homogeneousComponent_isHomogeneous (e+1) (M i)
  · simpa only [defect, nsmul_eq_mul, Nat.cast_ofNat, polynomialDivNum, polynomialDirectional,
      polynomialDivergence, mul_comm, mul_left_comm, mul_assoc] using hproj

/-- A denominator of any degree and four ungraded polynomial numerators
suffice. Multiplication by `X 0` first permits a nonzero constant term in D. -/
theorem coreCertificate_of_ungraded (D : CorePoly) (N : Fin 4 → CorePoly) (hD : D ≠ 0)
    (hraw : 2 * coreSource * D ^ 2 - 2 * D * (∑ i, pderiv i core * N i)
      + 3 * core * (D * (∑ i, pderiv i (N i)) - ∑ i, pderiv i D * N i) = 0) :
    ∃ T M, CoreCertificate T M := by
  have hc : defect D N = 0 := by
    simpa only [defect, nsmul_eq_mul, Nat.cast_ofNat, polynomialDivNum, polynomialDirectional,
      polynomialDivergence, mul_comm, mul_left_comm, mul_assoc] using hraw
  let x : CorePoly := X 0
  obtain ⟨d, ho, hn⟩ := least_component (x*D) (mul_ne_zero (X_ne_zero 0) hD)
  have hOrd1 : HasOrder 1 (x*D) := by
    simpa only [Nat.add_zero] using (isHomogeneous_hasOrder (isHomogeneous_X ℂ (0 : Fin 4))).mul (hasOrder_zero_degree D)
  have hd : 0 < d := by
    by_contra! hv
    exact hn (hOrd1.component_eq_zero (by omega))
  exact project_initial (x*D) (fun i => x*N i) d hd ho hn (by rw [defect_scale, hc, mul_zero])

/-- The protected four-plus-two certificate follows directly from one
ungraded cleared four-variable identity. -/
theorem hasCertificate_of_ungraded (D : CorePoly) (N : Fin 4 → CorePoly) (hD : D ≠ 0)
    (hraw : 2 * coreSource * D ^ 2 - 2 * D * (∑ i, pderiv i core * N i)
      + 3 * core * (D * (∑ i, pderiv i (N i)) - ∑ i, pderiv i D * N i) = 0) :
    HodgeComputer.HasCertificate (n := HodgeSubmission.n)
      HodgeSubmission.F HodgeSubmission.q HodgeSubmission.A := by
  obtain ⟨T, M, h⟩ := coreCertificate_of_ungraded D N hD hraw
  exact ofCoreCertificate T M h

end FourPlusTwoCertificate
