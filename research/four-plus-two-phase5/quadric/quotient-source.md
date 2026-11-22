# Source on the quotient of the coupled incidence

Fix generic ratios $(r,s)$, put $z=\alpha^2$, and abbreviate $p(z)=P(r,s,\alpha)$. Write $q(z)=Q_2(1,r,\alpha,s\alpha)=U+Vz$ and $c(z)=Q_1(1,r,\alpha,s\alpha)$. At the moving point put $Z=-U/V$, $L=c(Z)$ and $\widehat Q=Q_1/L$. On $VL\alpha\ne0$,

$$\widehat Q(1,r,\alpha,s\alpha)=1,\qquad L^3=p(Z).$$

Thus five coefficients of $Q_2$ and four free coefficients of $\widehat Q$ suffice. Form the bidegree-$(6,6)$ section $f=(H-p(Z)\widehat Q^3)|_{Q_2}$. Require that $f=h\,\sigma h$ for a section of bidegree $(3,3)$, where $\sigma$ negates $c,d$. The nine equations on the finite axis cover implement this factorization. At the moving point $f=0$ automatically. On the generic open where exactly one of $h,\sigma h$ vanishes at a given square root of $Z$, its choice selects the moving cubic. [Factor reconstruction](axis_model.md)

For a smooth integral cubic section $h$ distinct from $\sigma h$, unique factorization of the bihomogeneous section recovers its two projective factors up to order. Over the algebraic closure its compatible scalar choices and square roots give four pairs $(h,\alpha)$, permuted by $h\mapsto-h$ and $(h,\alpha)\mapsto(\sigma h,-\alpha)$. The three choices $Q_1\mapsto\zeta Q_1$, $\zeta^3=1$, give twelve presentations in the original chart over one generic quotient point, prior to selecting any component. This count uses the original normalized representatives for the cubics modulo $Q_2$.

For a polynomial input $R$ independent of $u$ and odd in $(c,d)$, the individual current and its source are invariant under all twelve choices. Choose an invariant tangent parameter locally, or use cancellation of its scale in the current. The parameter residue $\ell_R=\operatorname{Res}(R\,Dh)$ uses the four generators $(Q_2,h,3Q_1^2,G_2)$.

* The deck change $Q_1\mapsto\zeta Q_1$ multiplies this residue by $\zeta^{-2}$ and multiplies $\beta^2$ and every $S_i$ by $\zeta^2$.
* The generator change $h\mapsto-h$ changes both $Dh$ and the residue denominator sign, preserving $\ell_R$. It negates $j$ and all the moving minors, so $S_i/j$ is fixed.
* Under reflection, the change of physical variables has determinant one and sends $R$ to $-R$, hence $\ell_R$ changes sign. At the reflected moving root $j,E_r,E_s,E_Y$ change sign, and $E_Z$ and $S_i$ are preserved. Substitution in the flux and source formulas gives the same values.

Consequently the normalized trace of the source on an individual irreducible component may be computed on its image in the quotient, by trace transitivity. If its stabilizer is a proper subgroup, only the stabilizer degree is removed on that component. A trace over a union uses the actual field degrees. The component and its total fiber degree remain to be determined.

## An even moving equation

Set $k(z)=p(z)-c(z)^3$. At $q=h=0$ one has $dk=G_2\,dq+(\sigma h)\,dh$ in all four directions $D,r,s,\alpha$. Provided $(\sigma h)\ne0$ there, replacing the equation $h$ by $k$ multiplies the four moving minors by the same nonzero scalar. Their ratios in the source formula survive. Define, with evaluations at $z=Z$,

$$
J=(Dq)k_z-(Dk)q_z,\quad e_r=q_zk_r-k_zq_r,\quad e_s=q_zk_s-k_zq_s,\quad e_0=q_rk_s-q_sk_r .
$$

Let $\widehat c$ be the evaluation of $\widehat Q$ with independent moving variables. Subscripts on $q,p,k,\widehat c$ are moving-point derivatives with the parameter coefficients frozen. In particular $k_\nu=p_\nu-3L^3\widehat c_\nu$ for $\nu=r,s,z$ at this point. $D$ differentiates the parameter coefficients before substituting the point. Compute $Dk$ from $D(Q_1^3)$ before imposing the point equation.

For $i=1,2,3$ define

$$
T_i=Z\{3p_i(e_r\widehat c_s-e_s\widehat c_r+\widehat c_z e_0)
             -(e_rp_{i,s}-e_sp_{i,r})\}-i p_i e_0.
$$

The source in the packet's convention is

$$
{4\ell_R L^2\over3\alpha J}
(T_1+ZT_2+Z^2T_3,\ T_2+ZT_3,\ T_3).
$$

Indeed $j=2\alpha J$ and $S_i=2 L^2 T_i$ for the even equation. The scalar $\ell_R L^2/\alpha$ is invariant under all three transformations. It can be evaluated by reconstructing any lift. Global evaluation of these invariants and their field trace remains open. [Symbolic substitution](quotient_source.py)

At the smooth sample above 11, $(r,s)=(4,2)$ gives $Z=3$, $L=8$, $p(Z)=6$ and the normalized coefficients of $\widehat Q$ are $(6,1,6,10,0)$. Evaluation over $\mathbb F_{121}$ gives the same source $(9,9,6)$ in every one of the twelve presentations. [Exact orbit check](orbit_check.py)
