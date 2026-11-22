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

## Arbitrary polynomial core witnesses

A verified four-variable identity with $D\ne0$ and arbitrary polynomial numerators now implies the protected statement directly. Multiply all inputs by $X_0$ if needed to give the lowest denominator component positive degree. If that component is $h$ of degree $d$, set $p=d+2$ and use the finite inverse series for $D$ against $h$. Its product with $D$ equals $h^p$ up to a term of strictly higher order. Taking the homogeneous component of the cleared identity in degree $2dp+6$ gives a homogeneous core witness with denominator $h^p$ and numerators of degree $dp+1$. The proof adapts the Phase 4 packet's expansion and uses the project's existing homogeneous-component and finite-inverse lemmas. The `hasCertificate_of_ungraded` theorem was independently replayed on upstream 4.36.0 with the standard three axioms. [Protected closure from an ungraded identity](UngradedCore.lean)

The same result now accepts a rational vector field in the fraction field of $\mathbf C[a,b,c,d]$. The sole equation is $v(H)-\tfrac32H\operatorname{div}v=a^3c^3$, with the canonical extensions of the coordinate derivatives. Lean chooses a common denominator, differentiates the numerator equations, applies homogeneous extraction, and concludes the protected theorem. This implication also passed independent kernel replay on upstream Lean 4.36.0. [One-equation fraction-field endpoint](RationalCore.lean)

## Specialization along a pencil

Suppose a current for $P+\varepsilon S$ has source with constant term $t^4$ and no negative powers of $\varepsilon$. Let $q_j$ be the coefficient at $\varepsilon^{-j}$, allowing a pole of order $m$. Expanding the identity gives
$$D_P(q_0)+D_S(q_1)=t^4,\qquad
D_P(q_j)+D_S(q_{j+1})=0\quad(1\le j\le m),\qquad q_{m+1}=0.$$
Then the specialized rational current is
$$Q=\sum_{j=0}^m { -2/3 \choose j}(S/P)^j q_j.$$
The factor comes from the identity
$D_P((1+\varepsilon S/P)^{-2/3}q)=(1+\varepsilon S/P)^{-5/3}D_{P+\varepsilon S}(q)$.
The Lean theorem proves the finite coefficient formula by telescoping for any finite set of derivations. It starts from the displayed coefficient equations. A rational function in the pencil parameter has only finitely many negative Laurent coefficients, so an actual generic identity supplies these inputs even when direct evaluation at zero encounters a pole. On the ratio chart the result feeds the existing transport theorem. The proof passed kernel replay locally on 4.33.1 and on the Google Cloud worker with upstream 4.36.0. Its axiom audit lists `propext`, `Classical.choice` and `Quot.sound`. [Finite specialization proof](SpecializeCurrent.lean)

## A fixed affine denominator

Set $(x,y,z)=(b/a,c/a,d/a)$. For any nonzero $B$ in the three-variable polynomial ring, natural $n$, and polynomial triple $N$, `hasCertificate_of_fixed_support` takes the following identity directly to the protected theorem, with $h=H(1,x,y,z)$.

$$
2B\sum_i N_i h_i-3h\left(B\sum_i\partial_iN_i-n\sum_iN_i B_i\right)=3y^3B^{n+1}.
$$

The identity is linear in the unknown coefficients of $N$ for fixed $B,n$. The proof verifies the injection of the affine coordinate ring, its compatibility with partial derivatives, the lift of $N/B^n$ to the four-variable fraction field, and the exact target. It uses coordinate independence from the accepted cycle-three construction. The source passed independent kernel replay on Lean 4.33.1 and upstream Lean 4.36.0, with the three standard axioms. [Fixed-support Lean implication](FixedSupport.lean) · [Geometric criterion for the choice of $B$](supported-residue.md)

## Reproduction

Run `python3 clear_core.py` with SymPy installed. The nonzero control uses
`H=c⁶`, for which `qₜ=1/(5t)` solves the required equation. Its certificate is
`D=c³, N=(0,0,2a³c/15,2a³d/15)`. A second check verifies that this control
current fails for the catalogue core.

On a Linux host with elan, Git, Python 3 and the standard Mathlib build tools,
put the eight Lean files (`HodgeTarget.lean`, `FourPlusTwoCertificate.lean`,
`RatioTransport.lean`, `SpecializeCurrent.lean`, `UngradedCore.lean`, `RationalCore.lean`, `FixedSupport.lean`, `AxiomCheck.lean`) and `upstream-check.sh` in one
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
