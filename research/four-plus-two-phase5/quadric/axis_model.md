# A split-quadric model of the coupled norm-contact curve

Use $a,b,c,d=X_0,X_1,Y_0,Y_1$ and $\kappa=5/7$. Write the intrinsic norm equation as
$H-Q_1^3=A^2-B^2\bmod Q_2$. The coefficients of the quadrics are the ten parameters of the upstream chart, with the coefficient of $ab$ in $Q_2$ fixed at one. Both quadrics are sums of an $a,b$ form and a $c,d$ form. The cubic $A$ is even and $B$ odd under $(c,d)\mapsto(-c,-d)$.

Over the finite cover choosing an ordering of the two distinct roots of each binary block of $Q_2$, a smooth quadric has Segre coordinates $w,z$. On a finite-root chart put
$$
 Q_2=D(b-ma)(b-na)+e(d-uc)(d-vc),\qquad D=-1/(m+n).
$$
For $K=-e(v-u)^2/(D(n-m)^2)$ take
$$
 a=1+Kwz,\quad b=n+mKwz,\quad c=w+z,\quad d=vw+uz.
$$
This parametrization sends $\sigma(w,z)=(-w,-z)$ to the physical involution. The selected nodal seed requires a different factor chart, implemented in the checker. The same axis reconstruction applies on both charts.

The bidegree-$(6,6)$ section $f=(H-Q_1^3)|_{Q_2}$ is even under $\sigma$. The contact equation is equivalent to $f=h\cdot \sigma h$ for a section $h$ of bidegree $(3,3)$. Every such section is the restriction of a cubic by the Segre embedding. Its even and odd parts recover $A,B$ modulo $Q_2$.

On each of the two axes meeting at $(0,0)$, choose a cubic $u+vt+yt^2+xt^3$ whose product with its conjugate at $-t$ is $f_0+f_2t^2+f_4t^4+f_6t^6$. Use the same constant $u$ for the two axes. The possible axis cubics are given by
$$
 u^2=f_0,\quad x^2=-f_6,\quad
 (v^2+f_2)^2-8u^2xv-4u^2f_4=0,\quad
 y=(v^2+f_2)/(2u).
$$
These are finite algebraic choices. The derivative of the quartic with respect to $v$ is $-8u(ux-vy)$ at its root. Thus the cover is étale on $ux(ux-vy)\ne0$ for each axis, in characteristic different from two.

Write $h=\sum_{i=0}^3w^i h_i(z)$. The axis data determine $h_0(z)$ and all constants $h_i(0)$. For $i=1,2,3$ in order, subtract $[w^i](h_{<i}\sigma h_{<i})$ from row $i$ of $f$. Divide the difference by two. For odd $i$, compare coefficients $z,z^3,z^5$, and for even $i$ compare $z^2,z^4,z^6$. The matrices for the three unknown nonconstant coefficients of $h_i$, with $h_0=u+vz+yz^2+xz^3$, are
$$
 \begin{pmatrix}u&0&0\\y&-v&u\\0&-x&y\end{pmatrix},
 \qquad
 \begin{pmatrix}-v&u&0\\-x&y&-v\\0&0&-x\end{pmatrix}.
$$
Their determinants are $u(ux-vy)$ and $-x(ux-vy)$. All sixteen coefficients of the possible cubic now have rational expressions on the finite cover of the ten quadric coordinates. The surviving nine equations are the coefficients of $f-h\sigma h$ at $w^i z^j$ with $i=4,5,6$, $1\le j\le6$, and $i+j$ even. All other coefficients already agree. This gives equality of the contact schemes on this open, including nilpotents, by linear elimination.

For a physical point set $b/a=r,d/c=s,c/a=\alpha$. Let $Q_2=U(r)+V(s)\alpha^2$. On $V\ne0$, put $t=-U/V$. The fiber resultant on the quadric parameter space is
$$
 V^3\left(f_X(r)+\kappa rs t+2\kappa rs t^2+f_Y(s)t^3-(Q_{1,X}(r)+Q_{1,Y}(s)t)^3\right).
$$
It follows by evaluation of the homogeneous sextic and the identity $f=h\sigma h$. Adjoining a cube root of the bracket's first four terms allows elimination of the $a^2$ coefficient of $Q_1$ by a monic linear equation, leaving nine free quadric coordinates for the nine residual equations. For the original point fiber, lift $\alpha^2=t$ and keep the points with $h(\alpha)=0$. A covering trace must respect the choice of irreducible component, sign lifts, field degree and exceptional locus.

## Exact tests at the selected nodal branch

Pass the upstream intrinsic_branch_four_jets.npz to segre_test.py, axis_reduction.py and tangent_rank.py using Python with NumPy and SymPy. The latter imports intrinsic_norm_incidence.py from beside the archive.

At the selected mod-23 branch the $c,d$ block of $Q_2$ has a common factor $\tau$. Divide it out before splitting its two directions. The extension uses $v^2=15,\ \omega^2=14+13v$. The chart specializes to $b=0$ and parametrizes the smooth quadric for generic $\tau$. In $\mathbb F_{23^4}[\tau]/(\tau^3)$, both axis cubics have constant $u=1$, vanishing $v,y$ at the seed and nonzero $x$. All cubic coefficients recovered using only the axis seeds and the quadric jets agree through $\tau^2$ with the independent archive. All nine residuals vanish. Adding $\tau a^2$ to $Q_1$ produces a nonzero residual through the same order. The optional logarithmic differential test gives rank 35 for its 60-by-36 matrix at the seed and rank 36 on a smooth noncontact control.

In the intrinsic chart the differential is the remainder modulo $Q_2$ of
$$
 -3Q_1^2\,\delta Q_1-G_2\,\delta Q_2-2A\,\delta A+2B\,\delta B,\qquad
 G_2=(H-Q_1^3-A^2+B^2)/Q_2.
$$
It has rank 24 at the seed. Its sixteen cubic columns have rank sixteen there. The maximal minor omitting the $a^2$ coefficient of $Q_1$ has expansion $(5+19v)\tau\bmod\tau^2$, using the archived monomial order. This proves rank 25 along the generic formal branch and étaleness of the generic projection to this $Q_1$ coefficient. On the split smooth chart cubic reconstruction therefore leaves a rank-nine residual system. The characteristic-zero conclusions use the same good-place lifts as the upstream contact construction. The node itself has deficient rank; starting simple-root numerical continuation there would require a normalization or deflation first.

## Exhaustive split-chart specialization at 11

Compile split_fiber.cpp with c++ -std=c++17 -O3 -o fiber11 split_fiber.cpp. Run ./fiber11 4 6, ./fiber11 4 2, and ./fiber11 --census for the two fibers and the parameter census. The independent point check uses NumPy and SymPy; the lift and current calculations use NumPy.

The residue field is $\mathbb F_{11}$ and the coupled coefficient is $5/7=7$. Both binary parts of $Q_2$ are required to have distinct finite roots in this field. Order their representatives by $0\le m<n\le10$ and $0\le u<v\le10$. Require $m+n\ne0$ in the field and $e\ne0$. The search includes the contacts for which $h$ is defined over $\mathbb F_{11}$, its values at both pure-$X$ corners are nonzero, and the row reconstruction determinant $x(ux-vy)$ above is nonzero on the prescribed initial $z$ axis (the letters in this determinant refer to axis-cubic coefficients). The fiber runs additionally require $r\notin\{m,n\}$, $s\notin\{u,v\}$ and $-U/V$ a nonzero rational square. These conditions exclude chart boundaries, nonsplit quadrics and higher residue degrees. The count represents $h$ and $-h$ by the one with constant value $1,\ldots,5$. The sign conjugation $h\mapsto\sigma h$ is retained.

For each of the four coordinate ruling lines joining a pure-$X$ corner of ratio $i$ to a pure-$Y$ corner of ratio $j$, define $L_i=Q_{1,X}(1,i)$ and $M_j=Q_{1,Y}(1,j)$. At points $(1,i,t,jt)$ the four even coefficients of $H-Q_1^3$ are
$$
 f_X(i)-L_i^3,\quad \kappa ij-3L_i^2M_j,\quad
 2\kappa ij-3L_i M_j^2,\quad f_Y(j)-M_j^3 .
$$
A necessary condition is that each of the four sextics factors as $h_{\rm axis}(t)h_{\rm axis}(-t)$. The enumerator tabulates all cubic factors over the field, up to common sign. Four corner values determine the $X$ quadratic completely and leave one coefficient of the $Y$ quadratic. For a fixed moving fiber the remaining value $Q_{1,Y}(1,s)$ is determined by the cube condition because cubing permutes $\mathbb F_{11}$. For the parameter census this remaining value runs over the whole field, using any $s$ different from the two roots. The enumerator tests all factor pairs on the near axes and then all nine remaining equations. It checks the full 49 coefficients independently for every survivor.

| Sample | Quadrics | Possible $(Q_2,Q_1)$ | Axis-factor pairs tested | Contacts modulo common sign |
| --- | ---: | ---: | ---: | ---: |
| $r=4,s=6$ | 9,225 | 135,063,225 | 11,388,820 | 0 |
| $r=4,s=2$ | 9,225 | 135,063,225 | 11,558,400 | 2 |
| All parameters on specified chart | 27,500 | 4,428,902,500 | 373,592,120 | 2 |

The two contacts are mutual sign conjugates. They share $(m,n,u,v,e)=(2,6,1,7,7)$ and $Q_1=4a^2+8b^2+4c^2+3cd$. Their first cubic has coefficient matrix (rows indexed by $w$)
$$
 \begin{pmatrix}4&4&9&4\\6&0&6&9\\9&0&10&0\\7&3&9&6\end{pmatrix}.
$$
At $(r,s)=(4,2)$ the first cubic vanishes at $(w,z)=(2,10)$, with $\alpha=6$. Its conjugate has $\alpha=5$.

Run fiber_point.py for an independent verification. It evaluates the physical polynomials directly on a $7\times7$ grid. Its values determine all coefficients by interpolation because their degrees in each variable are at most six. The Jacobian ideals of the bidegree-$(3,3)$ cubic restriction on all four affine patches of the split quadric have Gröbner basis $[1]$ over the residue field. Its curve is therefore geometrically smooth. As an ample divisor it is connected, and adjunction gives genus four. Dual-number differentiation gives rank 25 for the contact equations, on ten quadric coordinates and sixteen cubic coefficients. Adding $w,z$ and the three equations $b-ra=d-sc=h=0$ gives rank 28. The reported 28-square minor is 8 modulo 11 on the first sheet. The selected 25 independent evaluations are an integral local basis of the 25 even coefficient equations: independence of their differentials forces the evaluation submatrix to be invertible. Hensel's lemma therefore gives a unique $\mathbb Q_{11}$ lift of this point on the chosen chart, with $r=4,s=2$, for the exact rational core. In fiber_lift.py finite differences modulo $11^2$ independently reproduce the minor; all 49 norm samples and three moving equations vanish through modulus $11^8$ under successive simple-root lifts.

Run fiber_current.py for the canonical single-sheet contribution. It lifts the displayed cubic to the physical 16-coefficient cubic and verifies $H-Q_1^3-A^2+B^2=Q_2G_2$ in all 44 monomials. It computes the one-dimensional kernel for the linearized contact equations, using a free tangent entry of one. In the intrinsic gauge $f_1=u+Q_1$, $g_1=u^2-uQ_1+Q_1^2$. Eliminating $u$ leaves the complete intersection $(Q_2,h,3Q_1^2,G_2)$. The degree-nine Macaulay matrix has a one-dimensional cokernel. Normalize its functional $\tau$ by $\tau(\det \operatorname{Jac})=96$ and compute $\ell=\tau(a^3c^3Dh)$. The two sheets respectively have $(\ell,j)=(8,10)$ and $(3,1)$ at $(4,2)$. Both give the single-sheet source triple
$$
 (9,9,6)\pmod{11},
$$
including the $4/3$ in the packet's convention. For the exact rational coupling both lifts give the same single-sheet triple $(13036912,13328335,93231)\bmod 11^7$. The checker recomputes the complete intersection, canonical residue, tangent and current over this residue ring using unit pivots. Evaluating the original flux density for each $t\ne\pm\alpha$ checks the current coefficient formula independently. These values provide a positive non-nodal test on the coupled specialization. The selected characteristic-zero component at the place above 23 has not been identified with either of these 11-adic sheets. This finite-field chart census computes neither a geometric degree nor a trace over the full fiber; the degree normalization, other sheets and rational correction remain open.

## Obstruction span of the new component

Run smooth_span.py with the directory containing the upstream contact_series_span.py. The script reuses the weighted polynomial bases, the fixed-Jacobian construction and the linear-algebra helpers in that file. It independently reconstructs the new physical seed and verifies the fixed equation and the cubic restriction. The weighted coefficient map has rank 70 of 70 at this smooth contact. Its first tangent comes from the physical 44-equation kernel used for the current. For successive coefficients of the 90-coordinate presentation, fix a set of 70 pivot columns of the differential, set higher free coefficients to zero and solve the coefficient equations. The script verifies all 70 equations through order 60.

The degree-six anti-invariant ideal matrix has rank 59 at this point. Normalize its one-dimensional left kernel by setting one fixed coordinate to one. Its series obeys a unit 59-by-59 recurrence. All 72 ideal-column annihilations and all 12 fixed-Jacobian-direction annihilations hold through order 60. The first 48 coefficient rows have rank 48 (pivot minor 9 modulo 11). The first 32 rows restricted to the $u^0$ sector have rank 32 (minor 1), and the first 16 restricted to $u^1$ have rank 16 (minor 7). The $u^2$ columns vanish. At the seed the one-dimensional raw degree-nine socle evaluated on $R\,Dh$ is five times the normalized obstruction row on all 60 weighted monomials. Its target entry is one, while the obstruction target entry is nine. Normalizing the Jacobian residue gives the value eight reported above.

The four physical homogeneous equations $(Q_2,h,3Q_1^2,G_2)$ generate all 286 monomials of degree ten, verifying the complete-intersection condition. The physical tangent is nontrivial in the normal space since the canonical residue on $a^3c^3$ is nonzero. As in the packet, the normal map is square at a smooth genus-four contact, so this tangent forces nonzero cokernel. Its invariant part vanishes by rank 70 of the invariant coefficient map. This gives the structural upper bound 59 for the anti-ideal rank alongside its unit 59-minor. The point lifts for the fixed characteristic-zero core, and the chosen 70-by-70 and 59-by-59 determinants are units at its place. The physical tangent lifts by the same linear algebra. Consequently the jet-rank lower bounds transfer to characteristic zero on the component through this smooth point. Its universal surface has dominant ratio projection by the invertible point-fiber differential above. This supplies a separate application of the packet's rank and dominance criteria, while its relationship to the component chosen at the place above 23 remains undetermined.

## A lower bound for the new projected swept surface

The plane $d=2c$ with $c=1$ contains three rational points of the central $\mathbb F_{11}$ contact, namely $(a,b)=(2,8),(4,7),(8,7)$. All three have invertible Jacobian in $a,b$. Run smooth_sweep.py UPSTREAM_CONTACT_SPAN_DIR 96 1664 with NumPy, SymPy and python-flint. The script computes 1,664 branch coefficients, independently multiplies all six weighted forms with FLINT, then lifts and verifies the three point pairs on the quadratic and cubic. On the 4,753 columns indexed by plane monomials $a^ib^j$ with $i+j\le96$, the resulting $3\cdot1,664$ rows have full rank over $\mathbb F_{11}$. The other three geometric points are unnecessary for this lower bound.

The full-rank minor persists on the characteristic-zero lift since every denominator is a unit. The swept image is a surface because the moving incidence dominates $(r,s)$, and it cannot be the plane $d=2c$. A defining homogeneous equation restricts nontrivially to this plane; dehomogenizing at $c=1$ would give a relation on the three lifted branches. Thus the projected surface on this component has degree at least 97. This bound addresses the defining equation of the projected surface. It gives no denominator bound for unrestricted rational coboundaries.
