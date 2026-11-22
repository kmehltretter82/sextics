# A Prym 3-torsion model for the moving contacts

Work geometrically in characteristic zero and let Q be a smooth even quadric in the four physical coordinates. Write C=Q∩{H=0}, L=O_C(1,1) and σ for the sign involution on the second block. Assume C is smooth. It is a connected bidegree-(6,6) curve of genus 25. The involution on Q has four isolated fixed points, and its tangent action at each is minus the identity. If the invariant sextic vanished there, its local equation in the even linearization would also have zero differential. Thus σ acts freely on C, and C'=C/σ has genus 13. The Prym has dimension 12. [Mumford, *Prym Varieties I*](https://www.dam.brown.edu/people/mumford/alg_geom/papers/1974a--PrymVar-I-NC.pdf)

## Geometric correspondence

A projective contact consists of a biquadratic q, even under σ, and an ordered bicubic h, up to nonzero scalar on each, such that H is a linear combination of q³ and hσh with both coefficients nonzero. The restriction to C makes q³ and hσh proportional. The sections h and σh have disjoint zero sets on C. A common zero on C would also be a zero of q, and the contact identity would annihilate dH on Q there.

The subscheme D={q=h=0} of Q is the complete intersection of the (2,2) and (3,3) divisors. A common component would be contained in the smooth irreducible C, whose bidegree is larger. The intersection has length 12. On C, near the support of D, the section σh is a unit and h is a unit times q³, so D is Cartier on C, div_C(h)=3D and div_C(q)=D+σD. Consequently τ=O_C(D)⊗L⁻¹ has order dividing three, and στ=τ⁻¹. It lies in the connected Prym since τ=(1−σ)(τ²) in multiplicative Picard notation, and 1−σ has connected image in the norm kernel.

For k=1,2,3 the restriction H⁰(Q,O(k,k))→H⁰(C,kL) is an isomorphism by the restriction sequence and Künneth. If τ were trivial, the unique lift of a section with divisor D would be a hyperplane form l. Both lifted projective sections would then be h=l³ and q=lσl. Their powers are proportional on Q and cannot give a nonzero H coefficient in this contact. Thus τ has exact order three.

The effective twist Lτ has exactly one section up to scalar at every such contact. Adjunction gives K_C=4L. The Koszul resolution of D gives h⁰(I_D(4,4))=9+4=13. Restriction identifies this with H⁰(C,K_C−D) since H⁰ and H¹ of O_Q(−2,−2) vanish. Riemann–Roch gives h⁰(O_C(D))=13−12=1.

Conversely take a nonzero Prym 3-torsion class τ, and an effective divisor D of a section s of Lτ such that D and σD are disjoint. Trivialize τ³ to lift s³ uniquely to a projective bicubic h. Use a trivialization of the norm of τ on C' to identify sσs with a section of 2L on C. This gives the standard even linearization, because that trivialization is pulled back from the quotient and the product is symmetric in s and σs. It lifts uniquely to a projective even biquadratic q. Changes of the two trivializations alter only scalars. On C, hσh and q³ have the same divisor and are proportional. After scaling, their difference on Q is cH, since the kernel of restriction in (6,6) has dimension one.

The scalar c is nonzero. Otherwise any common prime component of h and σh on Q would meet the ample curve C and contradict disjointness. Coprimality and hσh proportional to q³ would then force every prime in div_Q(h) to have multiplicity divisible by three. Its third part has class (1,1) on Q, so h is proportional to a cube of a hyperplane section. Restriction gives τ=0, a contradiction. Taking square and cube roots of the two nonzero coefficients over the algebraic closure recovers the affine norm identity. If s and σs share a zero, then c must vanish, because a nonzero c would make that point singular on C. This explains the disjointness condition.

Conjugating h exchanges τ with τ⁻¹. The projective contact is unique for its ordered nonzero torsion class whenever the open condition holds. Scalar lifts contribute six choices geometrically, and the conjugate orientation and moving square root must be tracked separately along the selected component.

## Five base coordinates, a finite cover and a rank condition

Smooth quadrics with smooth section C form an open set in the projective five-space of even quadrics. The nonzero 3-torsion in the relative Prym is finite étale of degree 3²⁴−1=282,429,536,480 there. Étale locally choose a representing line bundle τ. For an auxiliary divisor A∈|3L|, the degree-48 bundle (Lτ)(A)=K_Cτ has h¹=h⁰(τ⁻¹)=0 and h⁰=24. Its restriction to A has dimension 36. Thus effectiveness is exactly rank at most 23 of the 36-by-24 restriction matrix. At all genuine contacts its rank is 23. On this rank stratum the kernel is a line, giving s up to scalar. One then checks disjointness and the separate moving-point equation h(p)=0, where p lies on the fixed ratio line and on Q. Away from its exceptional loci, adjoining such p has degree two over the quadric base.

There is a further intrinsic check at a genuine torsion class. Multiplying by q identifies H⁰(C,2Lτ) with the (3,3) sections vanishing on σD. The complete intersection gives precisely q times the four hyperplane sections, together with σh. Thus h⁰(C,2Lτ)=5.

There are five continuous quadric coordinates. This finite-cover and matrix description does not supply coordinates for the level cover or its monodromy, count isolated moving points, or identify their distribution among components. The good open may miss components generically contained in the singular-branch locus. A moving degree or a normalized trace for an individual component still requires those controls.

## Exact check at the 11-adic moving sheet

The smooth (3,3) contact curve h=0 and the (6,6) sextic C are different curves. Over F₁₁ the fixed-corner values of C are (1,8,0,6). In particular C passes through the fixed corner (w,z)=(∞,0), where its quadratic term in (ξ,z)=(1/w,z) is 5ξ²+5z². The genus-25 étale calculation does not apply to this singular special fiber.

For the actual characteristic-zero Hensel lift at ratios (4,2), the four values modulo 121 are (67,118,11,39), all nonzero in Q₁₁. The 840-by-400 matrix of bidegree-(19,19) multiples of the four bihomogeneous first derivatives of H|Q has 395 unit pivots. Unit elimination modulo 121 leaves a residual matrix divisible by 11 whose quotient has rank five modulo 11. These five additional pivots have valuation one over Q₁₁. All (19,19) coefficients are spanned, proving geometric smoothness of the actual lifted C. The chart and the lift of the full moving norm equations are those of the previously verified 28-by-28 étale marked seed.

The three generators q,h,σh generate the whole (7,7) piece already modulo 11, with an 86-by-64 Macaulay matrix and a full minor 7. The ruling resultant of q,h is squarefree degree 12 and its factors have degrees 1,5 and 6. Thus D is reduced at this specialization, and the lifted section supplies a nonzero rational Prym 3-torsion class over Q₁₁.

There is a direct rank check at the known torsion class. Take B₀=W₀³Z₀³ and A=C∩{B₀=0}. A degree-(10,10) matrix with a full minor 9 verifies q does not vanish on A. Since div_C(q)=D+σD, sections of K_Cτ=3L+D have denominator q and numerators in
```
V = q·H⁰(O_Q(3,3)) + (σh)·H⁰(O_Q(2,2)).
```
This follows from the Koszul resolution of σD. In degree (5,5) the space vanishing on A is W=B₀·H⁰(O_Q(2,2)). The checked dimensions of V,W,V+W modulo 11 are 24,9,32, so the restriction has rank 23. Its one-dimensional kernel is represented by B₀q. The same dimensions hold on the lift by these unit minors, the Koszul formulas, and this kernel vector.

Run `python prym_seed.py` with NumPy and python-flint. The script rebuilds H,q,h from the physical forms and verifies the norm identity in all 49 coefficients over F₁₁ and modulo 121. It then checks the fixed corners, resultants, disjointness, the restriction rank and the two stages of the derivative matrix.
