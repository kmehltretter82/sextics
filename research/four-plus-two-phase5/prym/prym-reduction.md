# Four normal directions for effective Prym 3-torsion

Use the smooth double cover C→C′, the hyperplane bundle L and a nonzero Prym class τ from the torsion model. Write M=Lτ and K_C=4L. At a genuine norm contact s spans H⁰(M), with divisor D disjoint from σD. All fiberwise formulas are projective in s. Trivializations of the cube and norm of τ can be chosen locally.

## The obstruction dual

Serre duality and the complete intersection (q,h) identify the thirteen-dimensional obstruction dual as
```
H⁰(K_C M⁻¹) = (σs) H⁰(2L) ⊕ s² H⁰(L).
```
Indeed the left side is 2L+σD. Introducing denominator q gives the degree-(4,4) numerator space q·H⁰(O_Q(2,2)) + h·H⁰(O_Q(1,1)). The Koszul relation has degree (5,5), so the sum is direct. The two descriptions correspond by sσs proportional to q and s³ proportional to h.

Let o∈H¹(M) obstruct lifting s across a square-zero equivariant deformation of the quadric and torsion class. The cubic s³ extends as an ambient cubic, so 3s²o=0 in H¹(3L). The even biquadratic sσs extends as an even ambient quadratic, so (σs)o+s(σo)=0 in H¹(2L). These statements follow also from the product rule on Čech representatives. Thus o annihilates the four-dimensional summand s²H⁰(L) and the five-dimensional subspace (σs)H⁰(2L)⁺. Its four remaining coordinates are dual to (σs)H⁰(2L)⁻. The four mixed physical monomials ac,ad,bc,bd form a basis.

A purely polynomial check uses the Artinian complete intersection (Q,q²,h,σh) in four physical variables. Its socle has degree eight and negative σ-character, since the action interchanges the last two generators and has determinant +1 in the coordinates. The degree-six invariant piece of this quotient is the cokernel of the fixed-quadric linearized norm equations. The perfect residue pairing identifies its dual with the four odd quadratics. Compute in Segre coordinates and pair the derivative of H|Q in each chart parameter with these four bilinears. With a common socle scale, the exact mod-11 matrix at the known sheet is
```
       m  n   u  v  e
 ac    8  7   5  5  3
 ad    7  1   6  3  9
 bc    6  5  10  8  2
 bd    9  8   7  9  0
```
The first four columns have determinant 10. The kernel is (1,0,4,2,2). Eliminating all 21 free q,h coefficients in the full 25-condition Jacobian gives the same line. The polynomial calculation remains valid on this special fiber because q,h,σh have no common zero. The smooth-curve Prym interpretation applies on its characteristic-zero lift.

## Nine equations among universal principal parts

Choose a reduced auxiliary divisor A∈|3L| disjoint from σA, locally on the parameter space. There are 36 points geometrically. On the finite Prym cover let E=π_*M(A), of rank 24. On the projective bundle of lines in E the universal meromorphic section s has simple poles in A. Its principal part e takes values in F(1), where F=π_*(M(A)|A) has rank 36 and (1) is the relative projective twist. Vanishing of e is precisely holomorphy of s.

Each b∈H⁰(2L)⁺ gives a quadratic identity on s. The meromorphic differential b sσs has simple poles in A+σA and is σ-invariant. Residues at conjugate points agree. Summing over A alone therefore gives zero, and the residue at each p∈A is linear in the principal part of s at p, with coefficient bσs. This gives five relations whose coefficients are linear in s.

Each l∈H⁰(L) gives a cubic identity from the sum of residues of l s³. In local frames, write s=t⁻¹(r₀+r₁t+r₂t²+⋯) and write the coefficient differential as (k₀+k₁t+k₂t²+⋯)dt. The residue equals
```
r₀ [3k₀(r₁²+r₀r₂) + 3k₁r₀r₁ + k₂r₀²].
```
Division by the leading pole coefficient therefore gives a regular polynomial on the entire space of local jets, including r₀=0. As the residue itself is intrinsic, this quotient transforms as a linear functional on the principal-part line. The four relations thus have coefficients quadratic in s.

Their joint map Ψ on F(1) satisfies Ψe=0 identically. In the family of even quadrics with base hyperplane class x and relative hyperplane class u on P(E), the target, after choosing the cubic and norm trivializations locally, is
```
(H⁰(2L)⁺)∨ ⊗ O(2u−x)  ⊕  (H⁰(L))∨ ⊗ O(3u−x).
```
The factor −x comes from relative adjunction ω=4L⊗O(x), so a section of 4L has residue values in O(−x). The trivializations of τ have finite ambiguity when chosen compatibly.

At a genuine zero, the nine coefficient rows restrict to bσs and 3ls² above. Restriction to A is injective on H⁰(K_C M⁻¹), because K_C M⁻¹(−A) is τ⁻¹ on the fiber and has no sections. Thus Ψ is surjective on an explicit open containing the genuine zero locus in each such auxiliary chart. On that open e is a section of a rank-27 bundle on the 28-dimensional P(E). Its vanishing gives the expected one-dimensional problem, and the moving point imposes one further equation. Extending the kernel or its zero cycle across the discriminant, the auxiliary boundary and the rank-drop locus is required for a complete count. Rational GRR on the smooth open alone loses precisely the data needed for an intersection degree, since a nonzero discriminant hypersurface trivializes a multiple of x there.

Run `python prym_normal.py` and `python prym_syzygies.py` beside `prym_seed.py` with NumPy and python-flint. The second check chooses A cut out by the three w-rulings and three z-rulings of respective values 1,3,5. It checks all intersections are étale, that σA is disjoint there, and that q is invertible there. In the six finite degree-six coordinate algebras it solves H=0 through second order and forms the polarized residue tensors on all 24 sections. All 1,500 quadratic coefficients and all 10,400 cubic coefficients vanish. The arithmetic map on the 36 principal parts has rank nine at the holomorphic section and annihilates the independent rank-23 evaluation matrix.
