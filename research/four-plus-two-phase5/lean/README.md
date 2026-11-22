# From a four-variable current to the protected four-plus-two certificate

`FourPlusTwoCertificate.ofCoreCertificate` takes homogeneous polynomial data
`D, N₀, …, N₃` in `ℂ[a,b,c,d]`, with `D ≠ 0`, positive denominator degree
and `deg Nᵢ = deg D + 1` in the homogeneous sense. Its remaining hypothesis is

```text
2 A D² − 2 D ∑ HᵢNᵢ + 3 H (D ∑ ∂ᵢNᵢ − ∑ DᵢNᵢ) = 0.
```

Here `H` is the protected four-variable core and `A = a³c³`. The theorem
produces exactly `HodgeComputer.HasCertificate HodgeSubmission.F
HodgeSubmission.q HodgeSubmission.A` from those inputs.

`FourPlusTwoRatio.transport` proves the chart identity for the change of
variables `(a,b,c,d) = (a,ar,at,ast)`. Given a rational current
`q = (qᵣ,qₛ,qₜ)` with `(2/3)q(P) − P div(q) = t⁴`, set

```text
(pₐ,pᵦ,p𝒸,p𝒹) = 2a/(3t) (0,qᵣ,qₜ,s qₜ + t qₛ).
```

Then `p(H) − (3/2)H div(p) = a⁶t³ = a³c³`. For a concrete current over the
rationals, `clear_core.clear_current(qr, qs, qt)` substitutes these formulas,
clears denominators, checks homogeneity, and expands the exact core polynomial
identity. For a number field, pass its exact SymPy domain to `domain=`; the script also checks a control over $\mathbb Q(\sqrt{2})$. A rational current for the coupled core remains
the missing input.

## Reproduction

Run `python3 clear_core.py` with SymPy installed. The nonzero control uses
`H=c⁶`, for which `qₜ=1/(5t)` solves the required equation. Its certificate is
`D=c³, N=(0,0,2a³c/15,2a³d/15)`. A second check verifies that this control
current fails for the catalogue core.

On a Linux host with elan, Git, Python 3 and the standard Mathlib build tools,
put the four Lean files (`HodgeTarget.lean`, `FourPlusTwoCertificate.lean`,
`RatioTransport.lean`, `AxiomCheck.lean`) and `upstream-check.sh` in one
directory and run `bash upstream-check.sh`. The script pins upstream sextics
`2388542bd16eb4df849515418e1c2ab291e6dc6a` and loads the toolchain and
Mathlib revision specified by that checkout. `HodgeTarget.lean` was generated
by the protected verifier from the instance record and compared byte-for-byte
with that upstream revision's generator. A solve submission imports the
verifier-generated target.

## Prior work

`DiamondJM` gave the exact core equation and six-variable field,
`benediamond` gave the current and suspension formulas and the intrinsic
trace reduction, and `kmehltretter82` formalized the initial target binding
and suspension identity.
[Core equation](https://github.com/hodge-computer/sextics/discussions/120#discussioncomment-19504157)
[Trace reduction](https://github.com/hodge-computer/sextics/discussions/120#discussioncomment-19448881)
[Prior Lean artifacts](https://github.com/hodge-computer/sextics/discussions/120#discussioncomment-19648129)

The exact
irreducibility proof, coprimality and degree proofs, and polynomial suspension
identities were adapted from the accepted Phase 2 Prove2Me proof
`normalChart_residue_of_coreCertificate`, theorem
`5f8e0726-2721-43c9-a74a-c9e0e987f824`, submission
`83bdfb8e-2661-479f-8cf7-e1dda1aaf525`. The present files bind that
construction directly to the protected sextics verifier and verify the explicit
ratio-chart field transformation. [Discussion 120](https://github.com/hodge-computer/sextics/discussions/120)
