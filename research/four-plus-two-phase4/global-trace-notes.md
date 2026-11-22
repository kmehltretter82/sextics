## Reduction to a residue on the parameter curve

Write the moving equations as
$$
q=U+V\alpha^2,\qquad h=a_0+b_1\alpha+a_2\alpha^2+b_3\alpha^3.
$$
The coefficients are functions of the normalized parameter curve $\Gamma$ and of the ratios $r,s$. Put
$$
E=Va_0-Ua_2,\qquad O=Vb_1-Ub_3,\qquad g=UO^2+VE^2=\operatorname{Res}_{\alpha}(q,h).
$$
On $VO\ne0$, a common root satisfies $\alpha=-E/O$. For a tangent derivation $D$ of $\Gamma$, the moving Jacobian satisfies $j=Dg/(VO)$ on the root. Indeed $Vh-(a_2+b_3\alpha)q=E+O\alpha$ and $q(-E/O)=g/O^2$. If $S\equiv S_e+\alpha S_o\pmod q$, then, for generic ratios on the selected component after removal of fixed base points, with $L=V(OS_e-ES_o)$ the unnormalized core trace of $\ell_R S/j$ is the sum over the ratio fiber of the residues of $\eta_R L/g$. The current convention multiplies by $4/3$, then divides the full sum by the degree. On the proper smooth $\Gamma$ this is the negative sum of residues away from the ratio fiber.

The projective tuple of the 49 coefficients of the bidegree-$(6,6)$ polynomial $g$ defines a morphism on the normalized projective curve after the common fixed divisor is removed. Call its pullback of $\mathcal O(1)$ the line bundle $M$. Linear combinations of the coefficients represent sections of $M$, including $g(\cdot;r,s)$ for generic ratios. For each such combination $J$, $\eta_R J/g$ is a global meromorphic differential with poles only over the generic ratio fiber. Its total residue there is zero. This gives a subspace of the kernel of the residue sum.

At a generic parameter the leading entries $V,b_3$ have no common projective root in $s$. Near the seed $V=\tau V_1+O(\tau^2)$ and $\operatorname{Res}_s(V_1,B)=8+6v\ne0$ modulo 23. Likewise $U,a_0$ at the seed are $ab,a^3+b^3$, without common projective root. Thus neither $\alpha=\infty$ nor $\alpha=0$ contributes a spurious component dominating the full ratio plane. Parameters at exceptional chart points are finite in number and, after removing the common fixed divisor, their ratios lie on proper closed sets. A common root on the remaining generic fiber gives an actual point on the contact, with $u$ determined by $f_1$. The generic contact is smooth and connected by the upstream branch construction. It and its sign conjugate are distinct, so the locus $q=A=B=0$ inside the universal contact cannot dominate the ratio plane. Thus on generic ratios with simple roots the one-variable residue formula applies directly. The degree of the generic resultant divisor is the finite ratio degree, counted with multiplicities.

An exact check on any future full trace comes from taking the section $J=\partial_r g$. At a simple fiber root, differentiating $g=0$ shows that the residue of $\eta_R(\partial_r g)/g$ is minus the $dr$ coefficient of the moving differential. The sum of those coefficients over the full fiber is zero, by the residue theorem for sections above. Over $r=0$, the three distinguished sheets alone give an unnormalized sum
$$
-3\ell_R(0)\,\frac{U^2}{U^3+V^3/B^2}.
$$
Indeed each sheet contributes $-\ell_R(0)/(U+V\alpha^2)$. For input $A$ this is generically nonzero. Contributions outside this factor must cancel it. This gives a full-fiber check independent of source matching.

## Moving all boundary poles to a fixed finite set

There are two representatives for the same moving residue. Let $S(\alpha)$ be one of the three source polynomials with the common $4\ell_R/3$ omitted. Set
$$
K_q=-\sum_{q=0}\operatorname{Res}_\alpha\frac{S\,d\alpha}{qh},
\qquad K_h=\sum_{h=0}\operatorname{Res}_\alpha\frac{S\,d\alpha}{qh}.
$$
For the first, $K_q=L/g$ as above. For the second, put $b=b_3$, work in the power basis of $h/b$, and let $Q$ denote multiplication by $q$ in that cubic algebra. Since $\det Q=g/b^2$, the second representative is
$$
K_h=\frac{b\,[\alpha^2]\big(\operatorname{adj}(Q)\,[S]_h\big)}{g}.
$$
Here $[S]_h$ is the coefficient vector of the remainder modulo $h/b$. At a simple common root one has $D(q|_{h=0})=j/h_\alpha$, giving the same parameter residue with the positive sign. The identity extends as one of rational functions. The difference $K_h-K_q$ is minus the residue at $\alpha=\infty$ by the residue theorem on the projective line. In particular it is regular at generic simple moving fiber roots. At the nodal seed for generic $r,s$ with $rB(s)\ne0$, $q=r$ and $b=B(s)$. Here $K_h=[\alpha^2]([S]_h)/(b r)$ is regular. The nonzero boundary residues of $K_q$ recorded below therefore depend on the chosen representative.

The coefficients $V$ and $b$ controlling these two charts are binary forms of degrees two and three in $s$ and depend only on the parameter point. Form their Sylvester resultant $T=\operatorname{Res}_s(V,b)$, a function on $\Gamma$. In the original chart it has order three at the seed and its leading coefficient reduces to $8+6v$. In the intrinsic chart the leading coefficient is the same because the chart-change scalars at the seed are one. Thus $T$ is nonzero. The Sylvester adjugate supplies polynomials $A_1,A_2$ in $s$ such that $A_1V+A_2b=T$. Work away from the fixed zeros and poles of $T$ and from poles of the parameter coordinates.

Choose $m,n\ge1$ so that $V^m gK_q$ and $b^n gK_h$ have regular coefficients there. The source numerators have degree at most 13 in $\alpha$. One may use $m=6,n=14$. For the cubic representative, the adjugate multiplication polynomial has degree at most two in $\alpha$ and denominator dividing $b^2$. Reducing its product with $S$ from degree at most 15 to at most two costs at most 13 more powers of $b$, and the leading factor $b$ recovers one. In the expansion of $(A_1V+A_2b)^{m+n-1}/T^{m+n-1}=1$, let $w$ be the sum of terms containing at least $m$ powers of $V$. Then $w$ is divisible by $V^m$, and $1-w$ is divisible by $b^n$. The differential
$$
\eta_R\big(w K_q+(1-w)K_h\big)
$$
has the desired fiber residues because the two representatives have the same simple principal part at each generic fiber point. The weights clear both moving chart denominators. On the stated parameter open its only possible poles lie over $g=0$. Remove the common fixed divisor of the coefficient sections of $g$ as well. Outside the generic ratio fiber, the poles on the proper normalization are now confined to a finite set of **fixed parameter points**, independent of $r,s$. Their principal parts and the degree normalization are the outstanding global inputs. The cubic representative can be used on its own near the seed, where $T$ vanishes, as long as $rB(s)\ne0$.

Over the field $k(s)$ there is also an existence argument for clearing the $T=0$ exceptions on the affine curve. A common zero of $V$ and $b$ there must lie above a fixed parameter with $T=0$. For transcendental $s$ this requires both binary forms to vanish identically. Restricting the norm equation to $X_0=X_1=0$ would then express the prescribed squarefree binary sextic as a cube of the remaining quadratic. This is impossible. Hence $V$ and $b$ generate the unit ideal on the affine parameter chart over $k(s)$. Computing coefficients for that global Bézout identity requires the selected curve's coordinate algebra.

The resultant sections have no common fixed base point at a finite intrinsic contact. If its conic and cubic had a common homogeneous factor $L$, the norm equation would give $H=Q_1^3+LG_1$ for a homogeneous $G_1$ of positive degree. The three forms $Q_1,L,G_1$ have a common point in projective three-space over the algebraic closure, forcing $H$ to be singular there. The fixed core is smooth, as follows from its smooth Fermat reduction at 5. Thus their common locus at a finite contact has dimension at most one and cannot dominate the two-dimensional ratio plane. A universal root at $\alpha=0$ is excluded by the unit $X_0X_1$ coefficient in $Q_2$, and a universal root at infinity by the binary restriction above. Combined with the unit-ideal argument, this confines the complementary poles of a suitably chosen representative to the compactification boundary of the affine normalized chart. The explicit Sylvester patch can retain auxiliary fixed poles at $T=0$.

## Size of the generic fiber and the line r = 0

The data underlying the coupled computations are the contact jets archived by `benediamond`, lifted at a specified place via unit recurrences. On $(r,s)\in\{0,\ldots,6\}^2$ and jet orders 0 through 48, the 49 resultant evaluations have determinant $15+5v\in\mathbf F_{23}[v]/(v^2-15)$. Its norm is 11. Evaluation is a tensor Vandermonde isomorphism on bidegree at most $(6,6)$. The resulting 49 coefficient functions are independent, including after lifting in characteristic zero. Consequently $h^0(M)\ge49$ and $\deg M\ge48$. On any chosen parameter normalization with genus at least 48, Riemann–Roch for nonspecial $M$ and Clifford for special $M$ strengthen the degree lower bound to 96.

At $r=0$, all seven section evaluations vanish through order two at the seed. Evaluation of orders 3 through 9 has determinant $14+v$. The seven sections stay independent after removing this and all other fixed base points. Their moving divisor therefore has degree at least six, disjoint from the fixed seed for generic $s$. On a parameter normalization of genus at least six the bound strengthens to twelve by the same Riemann–Roch and Clifford argument. At fixed $r=0$, a contact with a spurious solution at $\alpha=0$ for every $s$ gives a fixed base point of this linear system. These statements give lengths of moving divisors. Counting rational points over $\mathbf F_{23}$ is a separate calculation.

## Numerical coupled field checks and scope

The raw 90-coordinate eliminant differs from the intrinsic one by the common scalar $\delta^3 f_{1,u}^2$, with $\delta=f_{1,u}f_{2,ab}-f_{2,u}f_{1,ab}$. Elimination gives $q_{raw}=\delta q$ and $h_{raw}=f_{1,u}h+q\cdot \text{(linear form)}$ and the scalar is a unit at the seed. The check regenerates all 67 intrinsic jets and verifies the gauge on the whole 7 by 7 grid.

For the three coupled source entries of the supplied trace formula, use the above $L$ after reducing all powers of $\alpha$ using the conic. Clearing by $V^6$ across both source and section columns yields 52 independent columns at two samples. At $(r,s)=(1,2)$, the minor on orders 4 through 55 is $8+17v$; at $(4,6)$, orders 3 through 54 give $6+5v$. Thus independence of these representatives modulo the displayed section span persists at generic ratios in characteristic zero. Other residue cancellations and rational coboundaries remain to be determined.

For the prescribed input $R=A$, the reductions modulo 23 of the residues at the chosen parameter seed, including the factor $4/3$ and prior to normalization by the global fiber degree, are $(14+v,22+22v,0)$ at $(1,2)$ and $(16+8v,19+8v,11+22v)$ at $(4,6)$. These use the canonical input form through order two. Other poles must be included in the global residue sum, with the sign explained above. The seed boundary residue depends on the representative of the reduction. Its value alone does not determine a source class.

To reproduce, run `fiber_degree.py` with the path to the upstream `contact_span_coefficients.npz` and `test_residue_span.py` with the path to the unpacked `four-plus-two-current` attachment.
