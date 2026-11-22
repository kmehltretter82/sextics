# The four-plus-two contact curve and protected Lean target

For the catalogue equation let $(a,b,c,d)=(X_0,X_1,Y_0,Y_1)$ and $A=a^3c^3$. A four-variable homogeneous certificate now implies the protected Lean theorem for the exact six-variable instance. Constructing that certificate is open. On the chart $(a,b,c,d)=(1,r,t,st)$, it suffices to construct a rational current $q$ for the exact coupled polynomial $P$ satisfying

$$\tfrac23 q(P)-P\operatorname{div}(q)=t^4.$$

The preceding packet computes the full normalized trace for the binary control and a local contribution for the coupled family. [Phase 4 calculations](../four-plus-two-phase4/README.md)

## The protected Lean implication

`FourPlusTwoCertificate.ofCoreCertificate` takes the four-variable polynomial identity together with its homogeneity conditions and constructs the six-variable certificate, including the coprimality proof. `FourPlusTwoRatio.transport` verifies the field transformation needed to use the chart equation. Both compiled and passed independent kernel replay on the upstream revision `2388542bd16eb4df849515418e1c2ab291e6dc6a` using Lean 4.36.0 on Linux and on the research toolchain using Lean 4.33.1. The reported axioms are `propext`, `Classical.choice`, and `Quot.sound`. The accompanying exact arithmetic check exercises a nonzero control and rejects its current on this coupled core. [Proofs and reproduction](lean/README.md)

The homogeneous extraction theorem removes the grading assumptions on a cleared identity. Its fractional-field corollary takes a rational vector field satisfying the literal core equation with source $a^3c^3$ all the way to the protected certificate. A separate formal specialization result allows a finite Laurent pole when passing to the coupled member of the pencil. [Lean statements and checks](lean/README.md)

## Smaller equations for the coupled family

On a smooth split quadric $Q_2$, the restriction $f=H-Q_1^3$ has bidegree $(6,6)$. A norm contact gives a factorization $f=h(w,z)h(-w,-z)$ with $h$ of bidegree $(3,3)$. Factoring two axis sextics gives finitely many choices. On the specified open set, three linear solves recover the remaining coefficients of $h$. Nine explicit equations in ten quadric parameters remain. The additional incidence at the prescribed ratios cuts out the moving fiber. At the normalized branch we recover all sixteen cubic coefficients through order two. The intrinsic Jacobian minor omitting the $a^2$ coefficient of $Q_1$ is $(5+19v)\tau+O(\tau^2)$, where $v^2=15$ modulo 23. This proves generic transversality there and also records the rank drop at the nodal seed. [Model, equations and tests](quadric/axis_model.md)

For a fixed moving fiber we also remove the deck and sign choices. Writing $Z=\alpha^2=-U/V$ and $L=Q_1(1,r,\alpha,s\alpha)$, the normalized quadric $Q_1/L$ has value one at the moving point and $L^3=P(r,s,\alpha)$. This leaves nine parameters. The source for the prescribed input is invariant on all twelve presentations of a generic quotient point. An even moving equation gives an explicit formula there. On a selected component the degree reduction is governed by its stabilizer. [Quotient source and trace transfer](quadric/quotient-source.md) · [Four-edge reconstruction and fiber bounds](quadric/finite_chart.md)

On the marked chart, seven equations in seven geometric parameters suffice on an explicit open set. A complete two-parameter slice at the marked seed over the algebraic closure of the field of eleven elements has exactly one simple solution. The finite algebra used in the seven-variable reduction has rank 480. This calculation leaves the degree and trace of the full moving fiber open. [Equations, open conditions and reproduction](marked/README.md)

On a smooth quadric section of the sextic, an ordered projective contact determines a nonzero Prym 3-torsion class with a unique section up to scalar. Conversely an effective twist with divisor disjoint from its conjugate recovers the contact. The characteristic-zero sheet through the 11-adic seed satisfies these conditions. [Geometric argument and arithmetic checks](prym/README.md)

The fixed-component support argument also gives a way to search for the rational current directly. A polynomial containing the swept surface can be defined by a perturbed resultant. With that polynomial fixed, the numerator search uses linear systems and is guaranteed to terminate, though its required size is unknown. [Support construction and coupled control calculation](marked/support.md) · [Residue on the literal affine core](lean/supported-residue.md)

## Exact algebraic starts

Write $g(x)=1+(5/7)x+(15/7)x^5+x^6$. Normalize the binary decomposition $g=C^3-B^2$ with $C=u(1+vx+wx^2)$, $B=b_0(1+b_1x+b_2x^2+b_3x^3)$, $m=u^3$ and $b_0^2=m-1$. We work where $m(m-1)w\ne0$. The first three nonconstant coefficient equations determine $b_1,b_2,b_3$ in succession. Coefficients four through six then give three equations in $m,v,w$. Over $\mathbb Q$, the quotient has dimension 40 as a vector space. Its lexicographic basis expresses $v,w$ and the chart inverse as polynomials in $m$ over one irreducible degree-40 equation. All seven coefficients of the original identity verify in that field.

The given mod-23 choice $(C,B)=(15+18x+16x^2,4+6x+17x^2+22x^3)$ lifts to $m\equiv91534\pmod{23^4}$. A direct calculation checks the algebra through $23^{32}$. Eight of the forty embeddings of the normalized binary field are real. The unique $m$ between 1 and 2 gives $m\approx1.39309512609600217$. Choosing the positive real cube and square roots gives the nodal seed used for complex continuation. Normalizing removes the separate root-of-unity and sign choices in $C,B$. [Generate the exact equations](binary_start.py) · [Check the rational algebra, the 23-adic lift and the real embeddings](verify_rur.py)

## Thirteen complex sheets and a pointwise source match

At $(r,s)=(1/10,1)$ we certify thirteen distinct simple solutions of the coupled equations and enclose all three current coefficients for each of the forty odd polynomial inputs on every sheet. For the input $A$ the unnormalized sum over these thirteen has third coefficient approximately $-0.193016339+1.58558028i$. Three other inputs on this subset span the three source entries at that point. Their unique constant combination giving $(0,0,1)$ there fails when the same sheets are continued to $r=1/10+10^{-14}$ with $s=1$. A uniform contraction bound identifies the local branches at both values, and all three error intervals exclude zero. Other inputs, other sheets and rational coboundaries are still available. [Certified roots, source calculations and continuation](complex/README.md)

A separate three-worker VM run found 115 distinct roots at the same ratios. Independent interval checks certified all 115. Conjugation and sign changes extend the count to 440 distinct moving coordinates and at least 2,640 distinct coefficient lifts across the unquotiented chart. The degree of the chosen component is still open. [Reproducible run and interval checks](complex/README.md#bounded-vm-orbit-and-interval-verification)

## A smooth contact over the 11-adic numbers

A finite-field enumeration found a geometrically smooth contact of genus four with a simple sheet over $(r,s)=(4,2)$. Evaluation of the original flux independently checks its source $(9,9,6)$ modulo 11 for $A$. The lifted source is $(13036912,13328335,93231)$ modulo $11^7$. On the component through this point the first 48 Taylor rows span the full obstruction space, with minor 9 modulo 11. Its relationship to the component studied at the place above 23 is open. [Contact, lift, current and spanning minors](quadric/axis_model.md)

## Degree of the projected swept surface

The surface swept out by the selected family projects to an irreducible surface in ordinary $\mathbb P^3$. The projection from $u^3=-H$ is finite, and the upstream family dominates the ratio plane. Intersect with $d=c$ and put $c=1$. At the nodal parameter there are six points $(a,b)=(\xi,0),(0,\xi)$ with $\xi^3=-B(1,1)$. The conic and cubic have an invertible point Jacobian at each, hence all six lift along the normalized branch. Their Taylor series yield 2,400 rows on the 2,145 plane monomials of degree at most 64. The matrix has full column rank over $\mathbb F_{23}[v]/(v^2-15)$. The branch and point lifts divide only by units, so a nonzero maximal minor survives in characteristic zero.

It follows that the degree of this projected surface is at least 65. Its irreducible defining form restricts nontrivially to $d=c$, since the central elliptic components contain varying $d/c$. Dehomogenizing that restriction would otherwise give a nonzero relation among the 2,145 monomials. The bound applies to this specific surface. On the component through the smooth 11-adic point, a separate test on $d=2c$ gives rank 4,753 on the monomials through degree 96. Its projected surface therefore has degree at least 97. The minimum denominator degree of an unrestricted core certificate remains open. [Branch recurrence](extend_branch.py) · [Independent polynomial check and plane-section rank](swept_slice.py)

## Rational corrections with specified fiber poles

For currents whose fiber poles are confined to $P=0$, $t=0$, and infinity, reduction in $t$ turns the problem into two base derivatives of sections of a rank-three connection. Either base derivative alone is insufficient, with arbitrary rational base denominators. A class-specific Frobenius calculation and a geometric analysis of the possible supporting divisors exclude this whole space for the prescribed source. The moving poles $t^2=\alpha^2$ of the norm-contact current are a separate part of the full trace. [Derivation, restrictions and reproducible tests](connection/gauss-manin.md) · [Equivariant elliptic maps through degree five](connection/elliptic-fibers.md)

## The prescribed class and its fiber poles

At the good primes 11, 17 and 23, the class with numerator $a^3c^3$ has a nonzero Frobenius slope-two component. Two independent coefficient calculations verify the same stable ranks, and direct point counts recover the traces. At 11 the canonical differentials on the smooth contact curve independently recover all 2,304 entries of the ambient Cartier matrix. Another 240 entries on five later coefficient rows check the result. [Class-specific calculation](cartier/target-cartier.md) · [Contact curve check](quadric/cartier_check.md)

For divisors above curves in the ratio base, the part of their degree-one cohomology that can contribute to this sector comes from constant elliptic factors with $j=0$. Their Gysin images have slope $3/2$ at the inert primes. The branch divisor and the chart boundaries do not supply the missing class. Thus any successful rational chart primitive for this target has a finite nonzero fiber pole away from $P=0$, even with arbitrary base denominators. The poles $t^2=\alpha^2$ in the norm-contact construction satisfy this requirement. [Geometry and scope of the obstruction](cartier/vertical-support.md)

The same slope test now excludes all finite unions of linear physical core divisors as additional fiber poles. The degree-23 Hessian calculation controls even singular plane sections. [Linear-support exclusion](cartier/linear-support.md)

## Reducing the unknown boundary residues

In characteristic zero one can move the complementary source residues to a single smooth parameter point by replacing the closed input by a cohomologous meromorphic form. For an exact change of input we give an explicit correction whose moving fiber poles cancel. Producing a globally admissible principal part and evaluating the remaining traces are still necessary. [One-point residue formula and exact-input correction](one-point-residue.md)

## Which closed input to trace

Three cyclic coefficient variations force scalar action on the target sector for a correspondence extending over those directions. A second argument computes the possible actions of the generic total correspondence using its six symmetry graphs. Consequently the literal input $A=a^3c^3$ suffices up to a constant for the generic total family. The normalized rational source and its correction are still unknown. A Lean specialization theorem permits a generic rational witness to have a finite pole at the catalogue point. [Exact matrix checks](variation/README.md) · [Generic total trace and correction](variation/total-trace.md) · [Moving-degree calculation](variation/moving-degree.md) · [Specialization theorem](lean/SpecializeCurrent.lean)

At the marked smooth contact, differentiating in the $b^6,d^6$ directions gives all 48 complementary test forms. The 96 compact tests have a unit jet determinant. Their unique local solution of order at most 96 has order exactly 96 and fails an additional boundary test when the interpolation divisor contains the poles of the slice coordinate. Its six computed seed residues are partial contributions. [Normalized jets, target pairing and boundary obstruction](variation/seed-tests.md)

## Reproduction and credits

Use Python 3 with NumPy, SymPy and python-flint. The binary basis calculation uses Singular 4.3.1. Let `packet` be the extracted public four-plus-two contact archive specified in the Phase 4 README. Run these commands from this directory, with scratch outputs in a directory of your choice.

```sh
python3 binary_start.py --prime 0 > /tmp/four-plus-two-seed.sing
Singular -q /tmp/four-plus-two-seed.sing > /tmp/four-plus-two-seed.out
python3 verify_rur.py /tmp/four-plus-two-seed.out
python3 extend_branch.py "$packet" 399 /tmp/four-plus-two-jets.npz
python3 swept_slice.py "$packet" --jets /tmp/four-plus-two-jets.npz --degree 64
python3 "$packet/algebra/intrinsic_norm_incidence.py"
for check in segre_test axis_reduction tangent_rank; do
  python3 "quadric/$check.py" "$packet/algebra/intrinsic_branch_four_jets.npz"
done
python3 quadric/quotient_source.py
python3 quadric/orbit_check.py
python3 quadric/four_edges.py
python3 quadric/fiber_point.py
python3 quadric/fiber_current.py
python3 quadric/smooth_span.py "$packet/contact-span"
python3 quadric/smooth_sweep.py "$packet/contact-span" 96 1664
python3 complex/certify_roots.py complex/endpoints.json
python3 complex/current.py complex/endpoints.json --transport-test
python3 connection/gauss_manin.py
python3 connection/base_search.py
python3 exact_input.py
python3 cartier/good_reduction.py 11 17 23
python3 cartier/inert_trace.py 11 17 23
python3 cartier/target_cartier.py 11 17 23
python3 cartier/check_coeff.py 11 17 23
python3 cartier/plane_hessian.py
python3 quadric/cartier_naturality.py "$packet/contact-span" cartier --full
python3 variation/certify_span.py
python3 variation/certify_intertwiner.py
python3 variation/weighted_intertwiner.py
python3 variation/trace_determinant.py
python3 variation/moving_intersection.py
python3 connection/fit_trace.py
python3 quadric/canonical_variation.py "$packet/contact-span" variation/span.json
python3 quadric/principal_source.py "$packet/contact-span" variation/span.json
python3 marked/marked_slice.py
python3 marked/marked_plane.py
python3 marked/support_projection.py
python3 marked/support_tangent.py
python3 prym/prym_seed.py
python3 prym/prym_normal.py
python3 prym/prym_syzygies.py
python3 lean/clear_core.py
```

The calculations use the exact contact family, tangent data and current of `benediamond`, the core and suspended field of `DiamondJM`, the earlier work of `kejace`, and the Phase 2–3 work of `kmehltretter82` with Antigravity and ChatGPT. Codex prepared and verified this continuation for `kmehltretter82`. [Original contact archive and formulas](https://github.com/hodge-computer/sextics/discussions/120#discussioncomment-19448881) · [Core equation and suspension](https://github.com/hodge-computer/sextics/discussions/120#discussioncomment-19504157) · [Prior Lean work](https://github.com/hodge-computer/sextics/discussions/120#discussioncomment-19648129)
