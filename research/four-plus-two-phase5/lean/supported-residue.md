# From divisor support to an actual rational current on the literal core

Let $H=H(a,b,c,d)$ be the catalogue core and $A=a^3c^3$. Over $\mathbf C$, the established geometric support theorem implies the existence of an ordinary rational solution to $v(H)-\tfrac32H\operatorname{div}v=A$. The implication applies directly to this fixed member. It places no restriction on the poles of the current apart from rationality.

## Geometric input

Set $T=(G=0)$, $G=u^3+H$, in $\mathbf P(1,1,1,1,2)$. The affine cone over $H$ has smooth projectivization by reduction modulo five: the reduction is $a^6+b^6+c^6+d^6$, and its partial derivatives have no projective common zero. The weighted hypersurface has no singularity either, because $G_u=3u^2$ and its weighted vertex lies outside $T$. The involution $\sigma$ negates both $c,d$ and fixes $u$.

The weighted threefold theorem puts $H^3(T,\mathbf Q)^-$ in coniveau one; its proof supplies a divisor supporting the anti-projector. This holds on every smooth invariant member by specialization of supported projectors. Alternatively, the fixed-member contact calculation identifies a proper curve correspondence with surjective cylinder map to $H^3(T)^-$, which is supported on its swept surface. The full scalar correspondence on a generic family gives the same conclusion on that family when its scalar is nonzero, since its image lies on the swept divisors. [Weighted threefold theorem, §§3B.1, 3B.5–6](https://github.com/hodge-computer/sextics/blob/2388542bd16eb4df849515418e1c2ab291e6dc6a/docs/03b-weighted-threefolds.md)

## The local residue, including its constant

Let $W$ be the weighted ambient space with the vertex removed. Its smooth divisor $T$ has a residue map $H^4(W\setminus T)\to H^3(T)$, with the Tate normalization immaterial to vanishing. Let $\Omega=\iota_E(da\wedge db\wedge dc\wedge dd\wedge du)$ with $E$ the weighted Euler field. Then $A\Omega/G^2$ descends to $W\setminus T$, since both numerator and denominator have weight twelve. The compact class $\alpha=\operatorname{Res}_T(A\Omega/G^2)$ is $\sigma$-odd, since $A$ is odd and $\Omega$ is invariant.

On $a=1$, set $x=b/a$, $y=c/a$, $z=d/a$, $w=u/a^2$, $h=H(1,x,y,z)$, $g=w^3+h$, and $\nu=dx\wedge dy\wedge dz$. Then $\Omega=\nu\wedge dw$. With $w\ne0$ put $\theta=y^3\nu/(3w^2)$. Direct differentiation on the ambient chart gives

$$\frac{y^3\nu\wedge dw}{g^2}=d\!\left(\frac{\theta}{g}\right)+\frac{dg}{g}\wedge\frac{2y^3\nu}{9w^5}.$$

Indeed $dg\wedge\theta=-y^3\nu\wedge dw$ and $d\theta=-2y^3dg\wedge\nu/(9w^5)$. The class $\alpha$ therefore restricts to $[2y^3\nu/(9w^5)]$ with the $d\log g$-first residue convention. In the Phase 4 coordinates $(x,y,z)=(r,t,st)$, $\nu=-t\,dr\wedge ds\wedge dt$; the chart numerator is correspondingly $-2t^4/9$. [Algebraic logarithmic residue convention, Stacks Project §50.15](https://stacks.math.columbia.edu/tag/0FMU)

## Algebraic primitive and descent of its character

Take a nonempty smooth affine open in $T$ outside the supporting divisor, $a=0$, and $w=0$. The restriction of $\alpha$ is zero. De Rham comparison and the affine calculation by regular differential forms therefore give an algebraic two-form $\eta$ with $d\eta=y^3\nu/w^5$, after multiplying by $9/2$. This is an existence statement in the function field of $T$. It uses neither the trace of marked points nor a choice of their individual residues. [Grothendieck, Theorems 1 and 1′](https://www.numdam.org/article/PMIHES_1966__29__95_0.pdf)

Write $K=\mathbf C(x,y,z)$ and $L=K(w)$ with $w^3=-h$. The nonconstant squarefree affine slice of the smooth branch hypersurface ensures $h$ is not a cube in $K$. Fix $\zeta^3=1$, $\zeta\ne1$, and $\rho(w)=\zeta w$. The form $y^3\nu/w^5$ has eigenvalue $\zeta$. Replacing $\eta$ in the function field by $\frac13\sum_{k=0}^2\zeta^{-k}(\rho^k)^*\eta$ preserves its derivative. Separability gives $\Omega^2_{L/\mathbf C}=L\otimes_K\Omega^2_{K/\mathbf C}$. Since $L_\zeta=w^{-2}K$, there is a unique $q=(q_x,q_y,q_z)\in K^3$ with $\eta=w^{-2}\iota_q\nu$ after projection. Thus

$$d\eta=w^{-5}\left[-h\operatorname{div}(q)+\tfrac23 q(h)\right]\nu,\qquad \tfrac23 q(h)-h\operatorname{div}(q)=y^3.$$

Here $q(w)=-q(h)/(3w^2)$, and $d\iota_q\nu=(\operatorname{div}q)\nu$. On the balanced four-variable field take $v_a=0$ and $(v_b,v_c,v_d)=(2a/3)q(b/a,c/a,d/a)$. Then $v(H)=(2/3)a^6q(h)$ and $\operatorname{div}_4(v)=(2/3)\operatorname{div}_3(q)$. The source is exactly $a^6y^3=a^3c^3$. The rational-field endpoint and homogeneous extraction already checked in Lean accept this equation.

The catalogue already records geometric certificate existence. This argument identifies the specific rational vector field existence theorem required by the new Lean endpoint. The formal proof of coniveau and its comparison with a rational primitive remains an additional Lean obligation. [Catalogue context](https://github.com/hodge-computer/sextics/blob/2388542bd16eb4df849515418e1c2ab291e6dc6a/data/four-plus-two-residual/why.md)

## A complete fixed-denominator-support search

Suppose $f(x,y,z)\ne0$ vanishes on the projection of a supporting divisor in the $a=1$ chart, which exists because the cubic cover is finite. Put $B=hf$. On $\operatorname{Spec}R$, $R=\mathbf C[x,y,z,1/B]$, the cover $R[w]/(w^3+h)$ is finite étale and $w$ is invertible. Its two-forms are obtained by base change. Applying the same affine argument and the character projector on the whole cover gives $q_i\in R$. Thus some $n\ge0$ and polynomials $N_i$ solve

$$\frac23 B\sum_i N_ih_i-h\left(B\sum_i \partial_iN_i-n\sum_iN_iB_i\right)=y^3B^{n+1}.$$

For each $n$ and each bound on the degrees of $N_i$, this is a finite linear system. Diagonal enumeration of both bounds halts with a solution for a valid $f$. When $B$ is rational, row reduction over $\mathbf Q$ suffices. The projected swept surface in the contact construction gives one possible $f$; multiplying it by other factors preserves the guarantee. Its actual implicit equation remains a computational input to this approach.

A solution also exists over $\mathbf Q(x,y,z)$ even before computing a supporting projection. Clear an arbitrary complex solution, fix its numerator and denominator degree bounds, and normalize a chosen nonzero denominator coefficient to one. The resulting finitely many coefficient equations have rational coefficients. The weak Nullstellensatz supplies an algebraic-number solution of the same system. Choose a finite Galois field $E/\mathbf Q$ containing its coefficients and replace its denominator by the product of its conjugates, which is rational. With that denominator fixed, the differential identity is a consistent finite *linear* system over $\mathbf Q$ for numerator coefficients of a sufficiently large finite degree. Gaussian elimination supplies a rational solution. This coefficient argument neither sums an incidence fiber nor assumes its degree. It makes enumeration of rational denominator polynomials a provably terminating alternative in principle.

For a more directed computation of $f$, start with the normalized 25 contact equations in 26 parameters, the two marked curve equations in $(x,y,z)$, and the chart inverses. Localize at a parameter Jacobian minor of size 25 that is generically nonzero on the chosen branch. On the smooth marked contact locus further localize at a nonzero two-coordinate Jacobian minor of the marked curve equations; one such minor is nonzero generically on the chosen family. Every component on that open has parameter dimension one and marked-curve dimension one. Its projection to three-space therefore has dimension at most two, and the closure contains the projection of the chosen swept surface. Any nonzero element of the resulting elimination ideal in $\mathbf Q[x,y,z]$ is a valid $f$ (enlarge by the sign-conjugate if the cylinder was separately projected). One can produce such an element by solving the empty generic fiber over $\mathbf Q(x,y,z)$ for a polynomial Bézout identity and clearing its coefficient denominators. This computation is still potentially large, but it asks only for one annihilating polynomial; it does not require decomposing the ratio fiber or recovering current moments.
