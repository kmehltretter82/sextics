# The discriminant, monodromy and determinant degrees of the Prym cover

Let S={H=0} be the smooth fixed sextic surface and σ its block involution. This calculation concerns a general pencil of even quadrics on S, and the ambient cover of its smooth members by nonzero Prym 3-torsion labels.

## Two types of nodal degeneration

The quotient S/σ is embedded in P⁵ by the six even quadratic monomials a²,ab,b²,c²,cd,d². The two Veronese relations and the even sextic give a complete intersection of degrees (2,2,3). This follows directly by taking the second Veronese of the invariant homogeneous coordinate ring. The two fixed block lines meet S in twelve simple points. Their images are A₁ quotient singularities, with (−1,−1) action on the tangent coordinates upstairs. Locally one can take the two coordinates from the opposite block, since the pure binary sextics have simple roots. Hyperplane sections of this quotient correspond to even quadrics on S.

At each of the twelve fixed points, the general section through it has one ordinary node upstairs. Indeed the three independent quadratic terms in the opposite block prescribe an arbitrary tangent quadratic form. The invariant local Morse lemma identifies an equivariant smoothing with xy=t and σ(x,y)=(−x,−y). A general pencil avoids intersections of the twelve hyperplanes in parameter space and all their exceptional subloci.

Away from them the other discriminant is the irreducible dual to the quotient surface. If it had codimension at least two, a general pencil would have only the twelve fixed nodes. This is impossible by the Euler-number calculation below. In characteristic zero the conormal biduality theorem gives a unique ordinary tangency at a general point of a hypersurface dual: the conormal map has a single, reduced contact point at a general smooth dual hyperplane, and its nonsingular differential gives a nondegenerate tangent Hessian. The dual is distinct from the twelve linear discriminants by biduality. Its generic hyperplane section has one ordinary node downstairs and two exchanged ordinary nodes upstairs.

The pencil has 24 simple basepoints, the intersection of S with two general quadrics. The blown-up surface has Euler number 108+24=132, since c₂(S)=108 by adjunction. A smooth genus-25 fiber has Euler number −48. Thus 132−2(−48)=228 nodes are counted upstairs. Twelve occur over fixed-point parameters and the remaining 216 occur in 108 pairs. In particular the nonlinear dual hypersurface has degree **108**.

The generic one- and two-node members are irreducible. If the reduced divisor 2H_S were reducible, a component of degree d≤6 would meet its residual in at most two points counted with multiplicity at these nodes. The Hodge index theorem gives intersection at least 2d−d²/6, strictly above two for 2≤d≤6. A degree-one component is a line of self-intersection −4 by adjunction and meets its residual with intersection six. Here H_S denotes the hyperplane class of S.

## The level cover over the pencil

On the 24-dimensional anti-invariant F₃-homology, a fixed-node Dehn twist acts trivially. Its vanishing circle in xy=t is preserved with orientation by σ. At a free pair of nodes let δ,σδ be the two disjoint vanishing cycles, oriented compatibly. The local action on an anti-invariant vector v is, up to the common Picard–Lefschetz sign,
```
v ↦ v + (v,δ)(δ−σδ).
```
The two vanishing cycles are independent and span a primitive integral lattice, because cutting the smooth fiber along both gives the connected normalization with holes. One can choose integral dual arcs closing to loops for each circle independently. Their difference remains nonzero modulo three, and the alternating pairing on the anti-invariant subspace is nondegenerate. Thus this local operator is a rank-one transvection. The mod-3 Prym labels are identified with this anti-invariant space using the fact that two is invertible. [SGA 7, Picard–Lefschetz formula](https://grothendiecksga.com/read/sga7/en/15-3.html)

The geometric monodromy on these labels is Sp₂₄(F₃). To see the spanning assertion modulo three, all the vanishing cycles of the pencil generate H₁ of its fiber: attaching the thimbles and capping the base gives the fundamental group of the total space as the quotient by the vanishing cycles, since the pencil has a section. The total space is a blowup of a smooth hypersurface in P³ and is simply connected. Projecting onto the anti-invariant space kills the fixed-node circles; the free-pair differences therefore span. In the full parameter space the smooth meridians of the free-pair dual hypersurface are conjugate since that hypersurface is irreducible. A nonzero invariant subspace pairs nontrivially with some spanning root, contains that root by the transvection formula, and contains its entire orbit. This proves irreducibility. The ambient P⁵ is simply connected and general-position contraction of loops shows that conjugates of the discriminant meridians generate the image. The twelve fixed-node operators are trivial. The group is consequently irreducible and generated by transvections. Hall’s theorem with r=1, ℓ=3 and empty exceptional set gives the full symplectic group. [Hall, Theorem 3.1](https://arxiv.org/pdf/math/0608718#page=6)

In particular the cover of nonzero labels is irreducible. Its finite projective normalization over P⁵ remains irreducible upon passage to a general pencil by successive Bertini hyperplane sections for the basepoint-free linear system of this finite map. Its normalized restriction has degree N=3²⁴−1=282,429,536,480. Each of the 108 free-pair parameters has 3²³−1 unramified labels and 2·3²² points of ramification index three. The twelve fixed-node points are unramified. Thus the total ramification degree is **108·4·3²²**, and Riemann–Hurwitz gives genus **207·3²²+2 = 6,495,879,339,065**.

## The first compactified GRR correction

Use the balanced 3-root extension of τ with τ³ trivial, and measure rational degrees on the coarse normalized level pencil (equivalently rigidify the constant root gerbe). Compare the determinant of cohomology with that of the untwisted bundle on the same cover. At an index-three point over a free pair, the coarse local chart at either node after base change is xy=s³. The balanced local twisted chart has coordinates uv=s, x=u³,y=v³ with group action (u,v)↦(ζu,ζ⁻¹v). The local character of τ is nontrivial on both nodes: it is the pairing of the 3-torsion character with the vanishing circle, and the two pairings are opposite. At the unramified labels and at the fixed nodes this character is zero.

For weight a=1 or 2, the difference from the trivial character in the degree-one GRR term at each ramified node is
```
(1/3) Σ_{j=1,2} (ζ^(aj)−1)/((1−ζ^j)(1−ζ^(−j))) = −1/3.
```
This is the local inertia contribution on the smooth twisted chart. Equivalently it is (3/2)(B₂(a/3)−B₂(0)); the two ordered branches give the factor two in the boundary formula. Untwisted ambient hyperplane powers have stabilizer weight zero and do not change this difference. Extra root-stack charts at character-zero nodes, if used, cancel in the comparison. [Chiodo, Theorem 1.1.1 and §3.3, equations (52)–(53)](https://arxiv.org/pdf/math/0607324#page=25)

Each of the 108 parameters contributes two such nodes at each of its 2·3²² ramification points. The total correction to det Rπ_* is **−144·3²²**. On the original pencil the universal divisor on S×P¹ gives Rπ_*O_C(d)=χ(O_S(d))O−χ(O_S(d−2))O(−1) in K-theory. Hence the compactified determinant of the τ twist of the ordinary physical d-th hyperplane power has degree
```
χ(O_S(d−2)) N − 144·3²² .
```
For d=1 this gives 36·3²²−20=1,129,718,145,904. For d=4 it gives −45·3²²−11=−1,412,147,682,416. These are determinant-of-cohomology degrees on the general level pencil. The contact locus has codimension four in the five-dimensional level space, and its marked degree also depends on higher boundary intersections, the compact extension of the residue relations and the selected component.
