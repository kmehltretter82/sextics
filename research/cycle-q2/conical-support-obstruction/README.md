# The conical paired-Euler support alone cannot serve cycle-q2

For the exact protected `cycle-q2` target, there is no homogeneous certificate with denominator

$$h^m,\qquad h=6X_0^6-6Y_0^6+X_0^3X_1Y_0Y_1,$$

for any nonnegative integer `m`. This statement concerns the **physical homogeneous certificate**. It does not exclude other or mixed supports, and the original instance remains open. In particular, it does not by itself exclude a coefficient ansatz in an affine connection problem: reconstruction or homogenization can introduce additional support such as `X2`, normal derivatives, or discriminants.

The computational core is a complete exact exclusion for the literal denominator `D=h²`. **A separate coefficient-adjoint checker verifies the integer obstruction on every degree-19 numerator monomial of normal degree below twelve in all six components, without using a symmetry projection, a generated matrix, or a rank assertion.** The passage from this denominator to every power of its support uses the project's [prescribed-support theorem](https://github.com/hodge-computer/sextics/blob/3f28813a41d183c58785d0844548fa290d128122/docs/01-certificates.md#prescribed-support-and-pole-multiplicity), with all hypotheses checked below. The direct computation is not a finite-degree search extrapolated to arbitrary denominators.

This continues `benediamond`'s [paired-Euler partial current and connection packet](https://github.com/hodge-computer/sextics/pull/223#issuecomment-6712406587), `fqlx`'s [draft #221](https://github.com/hodge-computer/sextics/pull/221), and the [cycle discussion](https://github.com/hodge-computer/sextics/discussions/61). `kejace` already excluded all denominator degrees at most five and degree-six denominators in a single character. The `h^0,h^1` computations below calibrate against that work; the `h²` test supplies the new support exclusion. Prepared for `kmehltretter82` with OpenAI Codex assistance. No Lean formalization or complete certificate is claimed.

## Exact target and geometric hypotheses

Variables are ordered `X0,X1,X2,Y0,Y1,Y2`. Let

$$F=\sum_i(X_i^6+Y_i^6)+\frac12(X_0^3X_1Y_0Y_1+X_1^3X_2Y_1Y_2+X_0X_2^3Y_0Y_2),$$

$$A=X_0^2X_1^2X_2^3Y_0^2Y_1^2Y_2,\qquad q=2.$$

`instance.rec` is an exact copy of the protected target. The checker parses its rational polynomial expressions using Python's AST and compares both polynomials and `q` to its own literal input; no floating-point evaluation is used.

The polynomial `h` is nonzero, homogeneous of degree six, even in total `Y` degree, and independent of `Y2`. Since `F` is monic of degree six in `Y2`, `h` is coprime to `F`: a nonconstant irreducible factor of `h` lies in the coefficient ring, and cannot divide a polynomial with leading coefficient one.

The projective hypersurface `F=0` is smooth over the complex numbers. At a nonzero vector choose a coordinate of maximum modulus `M>0`. Its pure-power derivative has modulus `6 M^5`. In any `X_i` derivative the sum of absolute coupling coefficients is `2`; in any `Y_i` derivative it is `1`. Every remaining monomial has modulus at most `M^5`. Thus the chosen partial derivative cannot vanish, because its coupling part has modulus at most `2 M^5<6 M^5`.

## Why the fixed-denominator search loses no ordinary certificate

Write

$$L_2(v)=2v(F)-F\operatorname{div}v.$$

Fix `D=h^m`. Initially allow an arbitrary quotient in the ordinary condition. Divide each homogeneous numerator by the monic polynomial `F²` in `Y2`:

$$N_i=M_i+F^2K_i,\qquad \deg_{Y2}M_i<12.$$

This is polynomial division over the five-coordinate coefficient ring. It introduces no denominators and preserves total homogeneity, so each `M_i` has degree `6m+7`. The literal identity

$$L_2(F^2K/D)=-F^3\operatorname{div}(K/D)$$

preserves ordinary divisibility with the **same** denominator `D`. More explicitly, after clearing `D²`, the change is a multiple of `F³` with polynomial coefficient `D div K-K(D)`. For the remainder `M`, the cleared error has `Y2` degree at most 17: `D` and its derivatives are `Y2` independent, every derivative of `F` has degree at most five, and each `M_i` has degree at most eleven. Because `F³` is monic of degree eighteen, this ordinary error must be zero. Thus every ordinary certificate with literal denominator `h^m` admits a strict representative with that same denominator and `deg_Y2 M_i<12`. No zero-quotient assumption was imposed on the original field.

## The certificate is verified on the full unsymmetrized numerator space

The standalone `verify_dual_adjoint.py` reads the protected `instance.rec`, constructs `F`, `A`, and `h` directly, then tests the full coefficient adjoints of

$$\sum_i[(4hF_i+4Fh_i)N_i-2Fh\partial_iN_i]=2Ah^3.$$

For each component it accumulates the coefficient of **every** monomial capable of pairing with the finite integer functional, with total degree nineteen and `Y2` degree below twelve. Every such adjoint coefficient vanishes exactly. All other numerator monomials pair to zero automatically because their image has no overlap with the functional's finite support. The source pairs nontrivially. This verifies the complete unsymmetrized search space directly, with no finite-group assumption or matrix construction in the verifier. `ADJOINT_VERIFICATION.json` records its successful readback and source/certificate hashes.

Thus the obstruction itself does not require trusting the symmetry calculation below; the finite group merely made generating a certificate cheaper.

## Symmetry used only to generate the certificate efficiently

The affine diagonal sixth-root group is

$$G=\{w\in(\mathbb Z/6)^6:3w_0+w_1+w_3+w_4=0,\quad3w_1+w_2+w_4+w_5=0,\quad w_0+3w_2+w_3+w_5=0\}.$$

It has order 432, including the six scalar substitutions; its projective quotient has order 72. Both `F` and `h` are invariant. Character averaging over `G` preserves the literal denominator, homogeneity, the strict identity, and the `Y2` degree bound. If `a=(2,2,3,2,2,1)`, the entire surviving numerator space in component `i` consists of monomials `Z^e` satisfying

$$|e|=6m+7,\qquad e_5<12,\qquad e-a-e_i\in W\pmod6,$$

where `W` is the row span modulo six of

$$(3,1,0,1,1,0),\quad(0,3,1,0,1,1),\quad(1,0,3,1,0,1).$$

The checker enumerates `G` and `W` exactly, verifies `|W|=108`, their orthogonality, and `432*108=6^6`; hence `W` is exactly the annihilator of `G`. It enumerates every exponent composition satisfying the displayed conditions. The global parity is already contained in this finite group; no additional parity or selected-support ansatz is imposed.

## Integer linear system and dual certificate

Use `F2=2F` to avoid rational coefficients. For `m>=1`, divide the cleared strict equation by the common nonzero factor `h^(m-1)`. The complete linear equation is

$$2h\sum_iN_i\partial_iF_2-F_2h\sum_i\partial_iN_i+mF_2\sum_iN_i\partial_i h=2Ah^{m+1}.$$

Each allowed numerator monomial gives one column. All polynomial coefficients are exact integers. At `m=2` the matrix has:

- 3,510 numerator columns;
- 3,768 coefficient rows;
- 70,899 nonzero entries.

FLINT's exact integer nullspace routine produced a left relation with 1,817 nonzero integer coefficients, divided by their positive gcd. The independent checker reconstructs every column and verifies, over the integers,

$$\ell(L(e_j))=0\quad\text{for all 3,510 columns},\qquad\ell(2Ah^3)\ne0.$$

The exact nonzero source pairing is recorded in `fixed_cone_m2.json`; as a convenient check it is congruent to 53 modulo 101. **The annihilation identities are checked over the integers, not inferred from ranks modulo a prime.** Consequently this relation excludes solutions with arbitrary complex coefficients. The reported exact rank 2,286 is a generation diagnostic and is not needed to trust the obstruction.

Calibration runs give:

| `m` | columns | rows | generated exact rank | nonzero dual coefficients |
|---|---:|---:|---:|---:|
| 0 | 70 | 87 | 58 | 24 |
| 1 | 728 | 1,293 | 576 | 301 |
| 2 | 3,510 | 3,768 | 2,286 | 1,817 |

The `m=2` generation took about 190 seconds in Python 3.14.6 with python-flint 0.9.0. Independent verification uses only the Python standard library and takes about 0.4 seconds here. No numerical tolerance or search timeout enters the exclusion.

## From h² to the whole physical support

The prescribed-support theorem applies to the smooth sextic, its homogeneous source of degree twelve, `q=2`, and the nonconstant homogeneous even legal equation `h`. If a certificate had all of its physical projective poles confined to `V(h)`, that theorem would supply one with denominator `h²` (a squarefree equation squared also suffices). The fixed-denominator result just proved excludes it. Therefore no certificate has physical projective poles confined to `V(h)`, including any literal denominator `h^m`.

This does not exclude unions `V(h) ∪ V(B)`, denominators with additional support, affine primitives acquiring poles on `X2=0` when homogenized, or unrestricted certificates. The paired-Euler current remains a valid partial current; a completion must acquire additional physical projective support.

## Reproduction

Strongest independent verification, requiring only Python 3 and no symmetry calculation:

```sh
python3 verify_dual_adjoint.py --instance instance.rec \
  --certificate fixed_cone_m2.json --output ADJOINT_VERIFICATION.json
```

A separate reproduction using the original symmetry-reduced matrix generator also passes:

```sh
python3 fixed_cone_search.py --m 2 --verify fixed_cone_m2.json
```

The default `instance.rec` is read beside the script; `--instance /path/to/data/cycle-q2/instance.rec` binds directly to another checkout. The checker reconstructs all monomials, columns, and right-hand coefficients from scratch.

To regenerate the certificate, install python-flint in an isolated environment and run:

```sh
python3 fixed_cone_search.py --m 2
```

The generated relation may have a common integer factor; dividing it and its source pairing by that gcd does not change verification. The supplied JSON is already primitive. No verifier, protected target, or repository source was modified by this research computation.
