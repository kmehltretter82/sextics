# A rank-three equation for the coupled core

The ratio chart of `four-plus-two-residual` has $w^3=-P$, $z=t^2$, $k=5/7$ and

$$
P=f+bz+2bz^2+gz^3,\quad b=krs,\quad
f=1+r^6+kr+2kr^5,\quad g=1+s^6+ks+3ks^5.
$$

The prescribed source is $z^2$ for $D_P(q)=\frac23q(P)-P\operatorname{div}q$. We reduce the fluxes whose finite fiber poles lie over $P=0$ or $t=0$, allowing arbitrary poles in the base field $K=\mathbf C(r,s)$. Under this hypothesis two base directions suffice to formulate a rank-three PDE. Both directions are necessary. Transverse moving base divisors and all members of $g=cf$ for $c\ne0$ have zero rational horizontal residue space. The outstanding choices include exceptional base supports and moving fiber poles such as the poles of the norm-contact current.

## Polynomial reduction and explicit reconstruction

Set $D_xh=\partial_xh-\frac23(P_x/P)h$. Then $D_P(q)=z^2$ means $\sum_xD_x q_x=-z^2/P$. The three operators commute. Symmetrizing for $t\mapsto-t$ gives $q_r,q_s$ even and $q_t$ odd. Over $K$, modulo the image of $D_t$ on odd rational functions with the specified poles, even coefficients have basis $1,z,z^2$. Geometrically multiply by $dt/w^2$.

Indeed, at each simple zero of $P$, $D_t$ removes poles of order $m+1\ge2$ using a primitive of order $m$ because $m+2/3\ne0$. These principal parts can be matched in $K[t]$ modulo the squarefree $P$. At $t=0$ the even pole orders are at least two and the connection is regular, so the same reduction applies there. For even polynomial degree $k\ge6$, $D_t(t^{k+1})$ has leading term $(k-3)t^k$ at infinity. The resulting remainder, with at most simple $P$ poles and polynomial degree at most four, has the form $U/P$ for $\deg_z U\le5$.

In the power basis $1,z,\ldots,z^5$, build a $6\times6$ matrix $M$ with first three columns $Pz^j$ and last three columns $z^j[(2j+1)P-\frac43zP_z]$, for $j=0,1,2$. The latter are $P D_t(tz^j)$. Its determinant is $64fg\Delta/27$ with $\Delta=\operatorname{disc}_zP$. Write $R$ for the first three rows of $M^{-1}$ and $[U]$ for the coefficient column. Thus $R[U]$ is the reduced class. Independence at arbitrary pole order follows directly. If $D_t h$ is a polynomial of degree at most four, $h$ has no finite poles and vanishes at all six roots of $P$. A nonzero odd $h$ then has degree $k\ge7$, giving leading term of degree $k-1\ge6$ with coefficient $k-4\ne0$.

For $x=r,s$ define $\Gamma_x=-\tfrac23 R[P_xz^j]_{j=0}^2$, $\nabla_x=\partial_x+\Gamma_x$, and $T=R[z^2]$. The exact equation is

$$\nabla_r B+\nabla_s C=-T,\qquad B,C\in K^3.$$

In fact, relative exact changes to $q_r,q_s$ pass to $q_t$ by commutation. Conversely write $B(z)=\sum B_jz^j$ and likewise $C(z)$. Form

$$
U=-z^2-P(\partial_r B(z)+\partial_s C(z))
   +\tfrac23(P_rB(z)+P_s C(z)).
$$

The PDE says that $M^{-1}[U]=(0,0,0,h_0,h_1,h_2)^T$. The flux $(B(z),C(z),t\sum h_jz^j)$ has exactly the prescribed source. Symbolic differentiation gives zero curvature and

$$\operatorname{tr}\Gamma=-\tfrac16(d\log f+7d\log g+d\log\Delta).$$

## Either single base direction fails at all degrees

Consider $\nabla_x V=-T$ over $\mathbf C(y)(x)$ for $(x,y)=(r,s)$ or $(s,r)$. The degree-18 polynomial $fg\Delta$ in $x$ is squarefree over the generic field, as certified at $y=1$. Both $T$ and $\Gamma_x$ have at most simple poles. At a simple zero, the determinant of $M$ has order one, so the residue has rank at most one. Its remaining eigenvalue, by the trace identity, is $-1/6$ at $f=0$ and $\Delta=0$, or $-7/6$ at $g=0$. A pole of positive order $n$ in $V$ would require the residue to have eigenvalue $n$ by comparing coefficients of order $-n-1$. At ordinary finite points a pole is likewise impossible. Hence $V$ is polynomial in $x$.

For $r=\infty$ use $t=r\,t'$, $w=r^2w'$. The limiting polynomial is $1+g(s)t'^6$. For $s=\infty$ use $t=t'/s$, $w=w'$, with limiting polynomial $f(r)+t'^6$. These fibers are smooth over the generic other coordinate and the normalized matrix $M$ is invertible. In the new frame the inhomogeneous form $(z^2/P)\,dt\,dx/w^2$ is regular and decays to order three in the local coordinate at infinity. Differentiating at fixed $t'$ in place of $t$ changes the covariant derivative only by a relatively exact term, by the derivative of a fiberwise change of coordinates. Thus the three coefficients must be regular in their new frames. The factors on the old fiber basis are $r^{2j-3}$ or $s^{-2j-1}$. The polynomial degree bounds for $V$ are $(3,1,-1)$ and $(1,3,5)$ respectively.

Multiplication by $fg\Delta$ now gives polynomial linear equations. At $y=1$ the augmented matrices have shapes $60\times7$ and $66\times13$. Full-column minors reduce to $643$ and $563$ modulo $1009$. This proves full augmented column rank over $\mathbf C(y)$, even if a hypothetical choice of rational coefficients had a pole at the test specialization. Thus neither single-direction primitive exists in this fiber-pole class. [Exact computations](gauss_manin.py)

## Moving base divisors

Suppose a curve with projective normalization $L$ in the $(r,s)$ base meets the finite loci $f=0$, $g=0$ and $\Delta=0$ at simple and mutually disjoint points. Suppose both ratios are nonconstant and, at infinity, any finite $r$ or $s$ has $f(r)\ne0$ or $g(s)\ne0$ respectively. The example $s=r+1$ satisfies these conditions. We show that $\nabla v=0$ on $L$ has only the zero rational section.

Let $R_\infty,S_\infty$ be the pole divisors of $r,s$, of degrees $m,n$. At a point at infinity, set $R_0=r$ if $r$ has a pole and $1$ otherwise, and similarly $S_0$. With $t=(R_0/S_0)t'$ and $w=R_0^2w'$, the leading fiber is a separable diagonal sextic. The normalized connection is regular. Transforming $t^{2j}dt/w^2$ gives the componentwise signed bound

$$\operatorname{div}(v_j)+(3-2j)R_\infty+(2j+1)S_\infty\ \ge 0.$$

At every finite point all $v_j$ are regular by the same residue test. At a simple $g$ zero the left kernel of $M$ extracts $[z^5]$. At a simple discriminant zero it is evaluation at the double root $z_*$. In each case the polar part of the first three rows of $M^{-1}$ is nonzero since the connection residue has nonzero trace. Regularity of $\nabla v=0$ therefore forces $v_2=0$ at each $g$ zero, and $v(z_*)=0$ at each discriminant zero. For the latter assertion, the local derivative of $P$ at the double root is nonzero by the simple discriminant hypothesis. There are $6n$ zeros of $g$ and the divisor bound for $v_2$ has degree $-m+5n$, whence $v_2\equiv0$.

At infinity the term $-27 f^2g^2$ in $\Delta$ has strictly larger pole order than its other terms, so $\Delta$ has $12(m+n)$ distinct zeros. Form the two rational functions obtained by homogenizing $P(-v_0/v_1)$ with $v_1^3$ and $P_z(-v_0/v_1)$ with $v_1^2$. They both vanish at these points, since $v_0+z_*v_1=0$. Their respective pole divisor bounds are $9R_\infty+9S_\infty$ and $6R_\infty+8S_\infty$. Hence both vanish identically. If $v_1\ne0$, this would give a multiple root on the generic smooth fiber. If $v_1=0$, $v_0$ must vanish at all the discriminant zeros, again exceeding its divisor bound. Thus $v=0$.

This test concerns residues of the base one-form $B\,ds-C\,dr$, whose covariant exterior derivative is the inhomogeneous target. At a smooth base divisor where the connection and target are regular, a closed leading polar part of order $h\ge2$ in a local normal coordinate can be removed by subtracting a covariant derivative of order $h-1$. The remaining logarithmic residue is a rational horizontal section along the divisor. Its vanishing excludes a nonzero residue there. This argument is local along the divisor and does not produce a global primitive.

## The special pencil $g=cf$

The previous transversality test does not cover simultaneous zeros of $f,g$. Yet the residue also vanishes on every irreducible component of $g=cf$ for each $c\in\mathbf C^\times$, including the reciprocal member $c=8$. On such a component the two ratios have equal poles with equal orders, with common divisor $R_\infty$. The $6\deg R_\infty$ zeros of $f$ are simple because both binary polynomials are squarefree. At those points $rs\ne0$.

Write $\lambda=b/f$ and $P=f(1+\lambda z+2\lambda z^2+cz^3)$. The connection becomes $-\frac23d\log f\ I+A(\lambda)d\lambda$. The discriminant of the normalized cubic is

$$D=4\lambda^4-(4c+32)\lambda^3+36c\lambda^2-27c^2,$$

with $\operatorname{disc}_\lambda D=-5038848 c^4(c-8)^4$. For $c\ne8$ all finite singularities of $A$ have residue eigenvalues $(-1/6,0,0)$. At $c=8$, $D=4(\lambda+2)(\lambda-6)^3$; the explicit matrix $A$ has only simple poles and the eigenvalues at $-2,6$ are $(-1/6,0,0)$ and $(-1/2,0,0)$. Pullback multiplies the residues by positive ramification indices, preserving regularity of rational horizontal sections. At the joint $f,g$ zeros the full pulled-back connection has residue

$$\begin{pmatrix}-1/6&0&0\\1/3&0&-5/12\\0&0&-7/6\end{pmatrix}.$$

Consequently $v_0,v_2$ both vanish at the $6\deg R_\infty$ joint zeros. All $v_j$ have infinity bound $4R_\infty$, so $v_0=v_2=0$ identically. Finally $A_{01}=(2\lambda^2-3c\lambda+9c)/D$ and $\lambda$ is nonconstant, having poles at the joint zeros. The first horizontal equation forces $v_1=0$. [Pencil identities and residues](gauss_manin.py)

## Finite tests with both base components

Allowing both $B,C$ with each of their six coefficients of the form $N/D_0$, $\deg_rN,\deg_sN\le12$, gives 1,014 unknown coefficients. For $D_0$ equal to $g-8f$, $\Delta$, $(g-8f)\Delta$, or $fg\Delta$, evaluation on 344 reproducibly chosen points mod 1009 gives 1,032 equations. In each case the augmented matrix has rank 1,015. The evaluated entries lift to rational numbers regular at 1009 and the nonzero maximal minor proves independence over $\mathbf C$. These bounds apply separately to each of the four denominator choices. [Reproduce with Python 3 and NumPy](base_search.py)
