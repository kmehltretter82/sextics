# A finite cover by a surface times the Fermat elliptic curve

Let $h=H(1,x,y,z)$ and $T$ be the cyclic threefold $w^3=-h$. Put $A=y^3$ and $h_e=h+eA$. Work on the surface parametrizing $p=(x_0,y_0,z_0)$ and a direction $v=(1,q,r)$ such that the coefficients of $1,\tau,\tau^2$ in $h_e(p+\tau v)$ vanish. Use $b_1=x_0,b_2=y_0$ as parameters wherever the Jacobian in $z_0,q,r$ is invertible. The residual cubic and the volume coefficient are

$$C_e(\tau)=\tau^{-3}h_e(p+\tau v),\qquad J_2=q_{b_1}r_{b_2}-q_{b_2}r_{b_1}.$$

In the volume determinant, the constant term vanishes since $v$ is tangent to $h_e=0$. For the linear term choose one base tangent extending $v$ at the point. The only possible normal contribution is that of the derivative of the vector field $v$ in its own direction. Differentiating $dh_e(v)=0$ in that direction shows its normal component is $-d^2h_e(v,v)=0$. Thus

$$dX\wedge dY\wedge dZ=\tau^2 J_2\,db_1\wedge db_2\wedge d\tau,\qquad (X,Y,Z)=p+\tau v.$$

At $e=0$ and modulo $11$ take $p=(0,3,4)$ and $v=(1,8,4)$. The Jacobian determinant of the three contact equations, multiplied by seven prior to differentiation, is $2$. The residual cubic for the original $h$ is $2+\tau+7\tau^2+5\tau^3$. Its discriminant is $3$ and $J_2=6$. Hensel lifting at fixed $b_1,b_2$ proves there is a characteristic-zero component with nonzero volume determinant and separable residual cubic. Its line family dominates the three-dimensional physical space and the corresponding function-field extension is finite. [Arithmetic checker](asymptotic.py)

On this cover write $w=\tau u$, so $u^3=-C_e(\tau)$. After a finite coefficient extension there is a Möbius substitution $\tau=(a\xi+b)/(c\xi+d)$ and $u=\lambda V/(c\xi+d)$, with the fixed curve $E$ given by $V^3=-(1+\xi^3)$. Lifting a bijection of the three branch points and adjoining a scalar cube root gives the substitution. Modulo base differentials, $d\tau/u^2=\gamma\,d\xi/V^2$, with

$$\gamma=(ad-bc)/\lambda^2,\qquad \gamma^6=-27/\Delta,\qquad \Delta=\operatorname{disc}C_e.$$

The last equality follows from the binary cubic discriminant's scaling by the sixth power of the determinant and by the fourth power of the scalar. Thus $\Omega_e=dX\wedge dY\wedge dZ/w^2$ pulls back exactly to $J_2\gamma\, db_1\wedge db_2\wedge\omega_E$, where $\omega_E=d\xi/V^2$.

## The literal source and its exact correction

Differentiate at $e=0$ with $b_1,b_2,\xi,V$ held fixed, differentiating the algebraic contact and Möbius parameters implicitly. At fixed physical coordinates $\partial_e(w^{-2})=2A/(3w^5)$. The motion of the physical coordinates contributes $d\mathcal I$ by the Lie derivative formula, with

$$\mathcal I=w^{-2}(\dot X\,dY\wedge dZ-\dot Y\,dX\wedge dZ+\dot Z\,dX\wedge dY).$$

Every physical coordinate on the right is pulled back. In particular $(\dot X,\dot Y,\dot Z)=\dot p+\tau\dot v+\dot\tau v$. The Möbius derivative $\dot\tau$ at fixed $\xi$ is a polynomial of degree at most two in $\tau$ over the new base field. The exact identity is

$$\phi^*\!\left(\frac{A\,dX\wedge dY\wedge dZ}{w^5}\right)=\beta_A\wedge\omega_E-\frac32 d\mathcal I,\qquad
\beta_A=\frac32\gamma\!\left(\dot J_2-\frac{J_2\dot\Delta}{6\Delta}\right)db_1\wedge db_2.$$

The projected cubic action multiplies the elliptic differential by $\zeta_3$. Physical sign sends $e$ to $-e$ and leaves the volume form at $e=0$ invariant. One may lift it by preserving the ordered residual roots. The target, and consequently the displayed base coefficient with the elliptic form fixed, is sign-negative. The arithmetic checker differentiates the three contact equations in $b_1,b_2,e$, including the two mixed second derivatives. At the sample it gives $\dot J_2=7$, $\dot\Delta=5$, and the coefficient $\beta_A/(\gamma db_1\wedge db_2)=8\pmod {11}$. This proves the residual two-form is nonzero. The formula is meromorphic on the parameter surface and does not by itself determine the poles or its class on a smooth projective model.

## Base divisors and the slope test

If $\beta_A=d\eta$ in the base field, the displayed identity supplies the explicit form $\eta\wedge\omega_E-(3/2)\mathcal I$ before normalized trace. Away from the branch, the physical hyperplane at infinity and divisors over the base, this form is regular. Indeed all derivatives of the Möbius substitution have numerator of degree at most two in $\tau$ after division by the relative Jacobian; one may write the fixed $\omega_E$ as $\gamma^{-1}u^{-2}(d\tau-P_1(\tau)db_1-P_2(\tau)db_2)$ with $\deg P_i\le2$. Thus there are no additional finite moving poles on $E$. For a divisor of the surface parameter space, its image in the physical core is swept by a one-parameter family of the contact lines or is a single line. Over a base curve with separable residual cubic, the ruling is isotrivially elliptic; after a finite base extension its product with $E$ dominates, and the cubic-nontrivial part of the Albanese comes from $E$. A double-plus-simple residual cubic gives a geometrically irreducible rational ruling, hence trivial cubic character on the degree-one cohomology. At a residual cubic with a triple root the three rational components can be deck-permuted. If such lines are finite in the full line Grassmannian, their images in the threefold have codimension at least two. The smooth sextic surface itself has no positive-dimensional family of contained lines since it is of general type.

Under that finiteness condition, all permitted divisor-supported images in the target sector come from Tate twists of the $j=0$ elliptic curve. The project’s nonzero slope-two component of the prescribed class at inert primes excludes a base-field primitive. A slope argument merely using proper graphs of the rational map requires control of their exceptional divisors. The direct support argument identifies the degenerations that must be checked. A primitive with additional moving elliptic poles remains possible.

For the degree-one assertion, let $D$ be a surface swept over a parameter curve with smooth residual cubic. Pass to the finite cover $C'$ of the parameter curve which orders its three roots and supplies the scalar cube root. The cubic deck fixes $C'$ and acts on $E$ alone. There is an equivariant dominant generically finite rational map $C'\times E\dashrightarrow D$. Pullback embeds $H^1(\widetilde D,\mathbf Q)$ into $H^1(C',\mathbf Q)\oplus H^1(E,\mathbf Q)$: resolve the rational map, use the degree push-pull for the resulting morphism, and recall that blowups of smooth surfaces at points preserve $H^1$. Thus any cubic-nontrivial summand embeds into $H^1(E)$, even when $C'$ has large genus. For a double-plus-simple cubic replace $E$ by the rational normalization; the same proof gives only invariant $H^1$. Over a base divisor on a projective resolution, the underlying contact line has a unique limit because the Grassmannian is proper. If the base image is a point, the swept image lies over a single line. If it is a curve, the preceding generic-fiber classification applies. Finite trace preserves regular differential forms at divisors in characteristic zero by the local trace, so new horizontal divisorial poles are not introduced by descent.

An independent count is available for all components of the ordinary line congruence together. For a general physical point $Q$, its contact point $p$ satisfies $H(p)=D_QH(p)=D_Q^2H(p)=0$, of degrees $6,5,4$ in $p$. On the smooth surface the tangent directions satisfying the quadratic contact condition form a finite set except where the second fundamental form vanishes. The literal Hessian rank bound from the project's linear-support computation excludes such points. Hence the incidence of marked asymptotic lines has dimension two, and the general polar intersection is zero-dimensional by dominance. Bézout gives total multiplicity $120$. This count makes no irreducibility assertion about the selected component or the additional finite cover trivializing its cubic.

The [complete projective cube-line calculation](../lines/README.md) proves the finiteness assumption in the last paragraph on this exact sextic. The [class-specific slope calculation](../cartier/target-cartier.md) concerns the prescribed numerator.
