# Exact parameter inputs and a one-point residue formula

Let $\Gamma$ be the proper smooth normalization of the chosen parameter curve. The patched source representatives in Phase 4 have only a fixed finite set of complementary parameter poles, with bounded orders over the generic ratio field. Given a particular regular closed parameter form $\eta$, choose a smooth point $p$ outside that set. Over characteristic zero one may replace $\eta$ by $\eta-d\phi$ so the complementary residues all vanish away from $p$. The global normalized **source** is then obtained from the local residue at $p$, with the sign, moving degree and $4/3$ factor of the original residue formula. This changes the representative by a rational coboundary which has an explicit formula without the moving fiber poles. Computing a global $\phi$ or an admissible principal part for the required class is still open. [Phase 4 residue formula and fixed pole set](../four-plus-two-phase4/global-trace-notes.md)

## The change of current

Write $w^3=-P(t)$, $\beta^3=-P(\alpha)$ and $u=(t-\alpha)/(t+\alpha)$. In the packet orientation the parameter input has primitive two-form $\mathcal B_\eta=(\beta^2/w^2)\,\eta\wedge d\log u$. For an exact input $d\phi$, subtract the exact two-form $d(\phi\,\beta^2w^{-2}d\log u)$ to obtain a primitive with the same differential,

$$
\mathcal B'_{\phi}=\frac{2\phi\beta^2}{3w^2}\,d\log\!\left(\frac{P(t)}{P(\alpha)}\right)\wedge d\log u.
$$

Set $B=(P(t)-P(\alpha))/(t^2-\alpha^2)$. This is polynomial in $t$ because $P$ is even. Expanding yields

$$
\mathcal B'_{\phi}=\frac{4\phi\beta^2}{3w^2P(t)}
\left[-2B\,dt\wedge d\alpha+
       (dB-B\,d\log P(\alpha))\wedge(\alpha\,dt-t\,d\alpha)\right].
$$

Thus $t=\alpha$ and $t=-\alpha$ have both disappeared from this difference of representatives on the generic chart. Trace is linear and commutes with exterior differentiation in the finite separable extension. This gives a rational correction for changing the parameter input by $d\phi$. It still requires traces of coefficients depending on $\phi$. The correction to the prescribed target has not been constructed. [Symbolic wedge check](exact_input.py)

## Which principal parts are allowed

Let $D=\sum_b M_b b$ be an effective divisor giving the required vanishing orders at the complementary poles, and let $E=\sum_b(M_b+1)b$. Locally choose primitives $F_b$ of $\eta$ at those points to the required order. For $N\deg(p)-\deg E>2g(\Gamma)-2$, Riemann–Roch makes the evaluation map from $H^0(\Gamma,\mathcal O(Np))$ to jets modulo $E$ surjective. Matching these primitives gives a global $\phi$ regular away from $p$ for which $\eta-d\phi$ vanishes to order at least $M_b$ at each $b$. For a rational point $p$ this gives the familiar bound $N\ge2g-1+\deg E$. After further base extension one may take all the fixed points rational.

There are constraints on independent choices of Laurent coefficients at $p$. For an arbitrary parameter differential vanishing at $D$ and with only pole at $p$, a proposed differential principal part $h$ is realizable exactly when $\operatorname{Res}_p(fh)=0$ for every $f\in H^0(\Gamma,\mathcal O(D))$. For the fixed cohomology class of $\eta$, a proposed function principal part $P_p$ of $\phi$ and the regular jets $F_b$ must instead satisfy

$$\operatorname{Res}_p(P_p\omega)+\sum_b\operatorname{Res}_b(F_b\omega)=0
\qquad(\omega\in H^0(\Gamma,\Omega^1(E))).$$

The integration constants in $F_b$ are free. These are the obstruction pairings from the principal-parts exact sequences and Serre duality. They explain precisely what a successful seed-jet computation must check globally. The traced flux itself has poles on the swept surface; removing complementary parameter residues in the source does not compute that flux.

At the smooth place above 11 the first 48 rows of the canonical differentials have determinant $6\bmod 11$. A global meromorphic function $\phi$ with a single pole at this point and principal part $\sum_k c_k\tau^{-k}$ satisfies $\sum_k c_k C_{k-1}(R)=0$ for every holomorphic canonical input, by the residue theorem applied to $\phi\,\eta_R$. The unit determinant lifts to characteristic zero. Consequently orders 1 through 48 are gaps in the pole orders of functions at this point, and $H^0(\Gamma,\mathcal O(48p))$ consists of constants. Any nonconstant one-point function used for the boundary interpolation here must have pole order at least 49. At larger orders the first 48 coefficients are already determined by the higher ones in these residue equations. [Canonical jet matrix and check](quadric/cartier_check.md)

If the local source kernel is $N_i(\tau;r,s)/g(\tau;r,s)$ and the seed image $g(0;r,s)=0$ is reduced, there is also a local obstruction to matching a regular target directly. Suppose $g_\tau(0)$ and $N_i(0)$ are generically nonzero along that image. For a principal part with leading coefficient $h_{-m}\ne0$, the highest-order pole in its seed residue is
$$
\frac{h_{-m}(-g_\tau(0))^{m-1}N_i(0)}{g(0)^m}
$$
up to the global residue sign. It comes from the Taylor coefficient of order $m-1$ of the kernel. Every lower principal-part coefficient has smaller pole order. Under these hypotheses, a rational chart correction is required to match the regular polynomial target.
