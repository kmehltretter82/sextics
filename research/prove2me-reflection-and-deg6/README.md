# Reflection-Chart Transport, Canonical Trace Reductions, and Degree-Six Core Normalization

Reproducible Lean 4 (`leanprover/lean4:v4.33.1`, Mathlib `0df444a360eaa60ab8c11dca51a86af692955474`) modules and exact modular/rational rank checks accompanying our `prove2.me/collections/hodge-computer` and `hodge.computer` discussion updates on:

- `symmetric-sextic-residual` (`Permutation symmetry · residual sector`, Discussion #117)
- `diagonal-permutation` (`Simultaneous permutations · standard sector`, Discussion #177)
- `four-plus-two-residual` (`Four plus two · suspension remainder`, Discussion #120)

**Scope:** All direct proofs in `Proofs/` are `sorry`-free and depend only on standard axioms `[propext, Classical.choice, Quot.sound]`. The root `HasCertificate` targets for all three frontier instances remain open; the modules in `Reductions/` formalize verified `SKETCH_ACCEPTED` reductions to single affine relative-primitive / canonical-trace targets.

## Contents

- `Definitions/`
  - `Def_HodgeComputerReflectionChart.lean`: shared 5-variable `(01)`-reflection chart, cleared cubic residue defect, `RelativePrimitive`, and cleared degree-3 canonical trace coefficients `(cubicTraceK0, cubicTraceK1)`.
  - `Def_HodgeComputerSymmetricChart.lean`: weighted cubic `symWeighted`, source `symSource`, and linear canonical trace source `symCanonicalSource` for `symmetric-sextic-residual`.
  - `Def_HodgeComputerDiagonalPermutationQuotient.lean`: `(01)`-transposition-invariant source `transpositionSource`, weighted cubic `diagWeighted`, cubic source `diagSource`, and linear canonical trace source `diagCanonicalSource` for `diagonal-permutation`.
- `Proofs/` (all `ACCEPTED` on Prove2Me)
  - `Sol_certificate_of_relative_primitive.lean`: proves `HodgeComputer.ReflectionChart.certificate_of_relative_primitive` (`62cd9479-5060-4477-af7f-35c19b081f1d`).
  - `Sol_relative_primitive_of_canonical_trace.lean`: proves `HodgeComputer.ReflectionChart.relative_primitive_of_canonical_trace` (`50907190-8358-42f0-ac43-d647e3439a5b`).
  - `Sol_symmetric_irreducible.lean`: proves `HodgeComputer.SymmetricSexticResidual.chart_and_target_irreducible` (`d9a76ddf-3d60-4cb0-82af-3f471b6c1e11`).
  - `Sol_diagonal_irreducible.lean`: proves `HodgeComputer.DiagonalPermutation.chart_and_target_irreducible` (`61adb094-0105-41c3-9e23-65eca6c96e25`).
  - `Sol_degreeSix_coreCertificate_iff_strict.lean`: proves `HodgeComputer.FourPlusTwoResidual.degreeSix_coreCertificate_iff_strict` (`e9e09848-70c4-4d05-b1af-3c6a6342d107`), formalizing `benediamond`'s degree-6 strict normalization.
- `Reductions/` (all `SKETCH_ACCEPTED` on Prove2Me)
  - `Sketch_SymmetricSexticResidual_certificate.lean` & `Sketch_SymmetricSexticResidual_relative_primitive_exists.lean`
  - `Sketch_DiagonalPermutation_certificate.lean`, `Sketch_DiagonalPermutation_transposition_quotient_certificate.lean`, & `Sketch_DiagonalPermutation_relative_primitive_exists.lean`
- `scripts/fpt_deg6_check.py`: standalone Python check of the degree-6 core `E1` parity split (`49 = 28_{(0,0)} + 21_{(1,1)}`), `(R_H)_{12}` quadric rank (`35/35`, `16/16` per pure parity), and Koszul bracket span (`420/420` of `J(H)_{12}`).
- `published_manifest.json`: Prove2Me definition, theorem, submission, and mission comment UUIDs.

Prepared for `kmehltretter82` with Antigravity (Gemini) assistance, building on prior contributions by `benediamond`, `fqlx`, `kejace`, and `DiamondJM`.
