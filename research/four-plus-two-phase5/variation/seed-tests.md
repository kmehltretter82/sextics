# Pulling all 96 compact tests back to the marked contact

The two cyclic variations $S=b^6,d^6$ give the 48 degree-twelve classes needed alongside the 48 degree-six classes. Their curve representatives can be obtained by differentiating the canonical contact residue. This supplies a finite linear system for the local pairing with $a^3c^3$. A separate global interpolation condition governs which polar parts can be used in the current. [Exact spanning minors](README.md)

## Differentiating on a rational slice

Let $z_0$ be the lifted smooth contact above 11, $v$ its physical tangent, and $E$ the 70 pivot coordinate vectors of the differential of the presentation map $\Phi$. Use the two-dimensional slice
$$z=z_0+\tau v+E x,\qquad \Phi(z)=G+\lambda S.$$
A nonzero free coordinate of $v$ expresses $\tau$ rationally. The square matrix $J=(d\Phi)_z E$ is invertible at the seed. The field of this slice is finite separable over the two parameter coordinates on the selected component. The lifts $T=\partial_\tau$ and $D=\partial_\lambda$ commute, with
$$Dz=EJ^{-1}S,\qquad Tz=v-EJ^{-1}(d\Phi)_zv.$$

For $h=A_{\rm even}+B_{\rm odd}$, let $e_R$ be the normalized Artinian residue of $R\,Th$ for the five generators $(f_1,f_2,h,g_1,g_2)$. Differentiate the raw cubic first, before substituting for $u$. Differentiating the entire coefficient gives the test $(D e_R)_{\lambda=0}\,d\tau$ corresponding to $-2\operatorname{Res}(RS\Omega/G^3)$. On the smooth open this follows by naturality of the Gauss–Manin connection under the proper universal contact curve over the parameter open. Equivalently, transport a compact loop in that open and integrate over its cylinder of curve members. The derivative equals the period of the differentiated threefold residue. This identifies the class with the restriction of the compact correspondence pullback, so its meromorphic representative has zero residues at the omitted points. Changing the rational lift $D$ by $fT$ adds the exact form $d(fe_R)$, as follows by evaluating the Lie derivative on $T$. [Canonical forms on contact curves](../../../docs/03b-weighted-threefolds.md)

All 70 differentiated equations are checked coefficientwise. Two first-order embeddings over the residue field use shifts by $\tau^{33}$ and $\tau^{44}$, which have zero $\tau$ derivative modulo 11. The blocks agree through order 27 on all 60 raw inputs. A separate calculation evaluates $\lambda=0,11$ modulo 121 and recovers the first eleven derivative coefficients by taking their difference and dividing by 11. For $R=a^3c^3$ these are $(7,2,2,2,8,0,4,9,5,3,6)$ along $b^6$ and $(0,0,2,3,3,2,2,4,7,8,0)$ along $d^6$. The complete matrix $X$ of Taylor rows of orders 0 through 95 on the selected 96 representatives has determinant 5 modulo 11. Its unit determinant also certifies invertibility in characteristic zero at the lifted point. [Jet and lift verification](../quadric/canonical_variation.py)

## Matching the prescribed numerator locally

Let $H_p=\sum_{k=1}^{96}c_k\tau^{-k}$ be a proposed polar primitive of the parameter input. The pairing functional is $\operatorname{Res}_p(H_p\omega)$. Fix the common projective normalization by assigning $\operatorname{Res}_{\partial G}(\det\operatorname{Hess}G)=1250$. Holomorphic tests have zero pairing with $A$. A complementary test from $S,R$ has prescribed pairing $-2\operatorname{Res}_{\partial G}(A R S)$. Only the 32 products with one power of $u$ contribute, since the graded Jacobian socle has $u$-degree one.

Solving $X^{\mathsf T}c=y$ gives $c_1=4$, $c_{96}=8$ modulo 11. Hence order 96 is necessary for a primitive of order at most 96 in this local test, with this choice of representatives. The script provides all coefficients. Its solve uses no division by a polar index. In particular $c_{11},c_{22},\ldots,c_{88}$ reduce to $(6,9,1,8,2,8,3,7)$ even though their coefficients in $dH_p$ reduce to zero. A different common sign or projective normalization rescales the nonzero vector. [Principal coefficients and target pairings](../quadric/principal_solution.py)

With the patched conic/cubic representative $K_i$ the single seed contributes
$$-\tfrac43\operatorname{Res}_p(dH_p\,K_i)
=\tfrac43\sum_k k c_k[\tau^k]K_i.$$
Before division by the generic fiber degree, the three coordinates at $(r,s)=(4,6),(1,2),(2,3)$ are respectively $(3,3,6),(8,0,10),(6,8,5)$ modulo 11. The implementation checks the cubic kernel separately by its adjugate and the difference between representatives by a residue at infinity. [Six checked samples](../quadric/principal_source.py)

For a differential vanishing to the required order at the fixed complementary poles, one must also impose $\operatorname{Res}_p(H_p\,df)=0$ for all $f\in H^0(\Gamma,\mathcal O(B))$ for the chosen vanishing divisor $B$. For example, if $B$ contains the polar divisor of $\tau$, the function $f=\tau$ already gives $c_1=4$ modulo 11. The order-96 candidate fails that strengthened localization. Determining the minimal boundary conditions and their solutions remains open. [Global requirements for moving the trace to one point](../one-point-residue.md)
