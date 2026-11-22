# Rational corrections on the four-plus-two chart

Use the coordinates and density sign of `benediamond`'s four-plus-two current. Write
$z=t^2$, $p=(f,\kappa rs,2\kappa rs,g)$ and $P=\sum_{i=0}^3p_i z^i$, where $\kappa=5/7$. A rational chart field has source
$$
\mathcal D_P(q)=\tfrac23 q(P)-P\operatorname{div}q .
$$
Normalize the protected source to $z^2$.

## An explicit correction test

Let $k=(k_0,k_1,k_2,k_3)\in\mathbb C(r,s)^4$ denote the coefficients of a prospective source. Form the matrix with rows
$$
M_i=(p_i,\tfrac{4i}{3}p_i,\tfrac23\partial_r p_i,\tfrac23\partial_s p_i).
$$
Its determinant is $27200/441$ at $(1,1)$. Therefore $x=(L,a,b,c)=M^{-1}k$ is uniquely defined in the ratio function field. Define
$$
J(k)=L+a+\partial_rb+\partial_sc .
$$
Then $k=\mathcal D_P(at,b,c)$ for a chart field with these three coefficient functions **if and only if** $J(k)=0$. Indeed, coefficient comparison gives $Mx=k$ and $L=-a-\partial_rb-\partial_sc$. In particular the inverse also gives the correction field. Finite separable trace commutes with both derivatives, so one may apply this test to the input or output of a normalized trace.

For $k=(0,0,1,0)$, $J(k)(1,1)=22379593/11895240$. For the 40 odd-odd monomial sextics transformed to chart sources $tS(1,r,t,st)$, this linear operator has rank exactly 36 over the constants. The script exhibits a $36\times36$ evaluation minor equal to 124 modulo 1009. Four independent elements of its kernel are $aH_c,bH_c,aH_d,bH_d$ in the original homogeneous coordinates $(a,b,c,d)$. Their ratio-chart fields have coefficient triples $(3/2,0,-3s/2),(3r/2,0,-3rs/2),(0,0,3/2),(0,0,3r/2)$. Thus the intersection of this rational chart family with the polynomial source space lies inside its eight polynomial-exact directions. The script checks an independent positive control with nonconstant poles $b=r/(1+s), c=s/(1+r)$.

For a field $q=\sum_n z^n(a_nt,b_n,c_n)$ of prescribed Laurent support, the corresponding condition is the existence of source vectors $u_n$ with $k(z)=\sum_n z^n u_n(z)$, $\deg u_n\le3$, and $J(u_n)+2n(M^{-1}u_n)_a=0$ for every $n$. This follows by differentiating $z^n t$ with respect to $t$. Allowing a second power introduces unknown functions through the overlap between the coefficient lists.

## Rational fields supported in one coordinate

There is also an exact check on the target for a field with just one nonzero component, with arbitrary rational dependence on all three variables. For each coordinate $x\in\{r,s,t\}$ the equation would be
$$
-P\,\partial_x b+\tfrac23 P_xb=t^4.
$$
View $P$ as a sextic over the rational function field in the other two coordinates. Its leading coefficient is nonzero and it is squarefree, as certified by specializing the other two coordinates to one. A pole of $b$ away from $P=0$ would give a pole of one higher order in $P\partial_x b$. At a simple root of $P$, a pole of order $m\ge1$ gives leading factor $(m+2/3)P_x$, which is nonzero. Thus $b$ is polynomial in $x$. A nonzero term of maximal degree $j$ gives leading coefficient proportional to $4-j$ in degree $j+5$. The target has degree at most four in $x$, which forces $j=4$.

Use $1,x,\ldots,x^4$ for the candidate, augment their source columns by $t^4$, and take the rows for degrees zero through five. At the same specialization the determinants for $r,s,t$ are respectively $14707736000/194481$, $22495616000/194481$, and $-20432613056/151263$. A nonzero specialized minor proves independence over the generic rational function field. This rules out each single-component solution of the prescribed equation. The calculation leaves fields with two or three nonzero components open.

## The fiber over the nodal contact

At the branch seed $r=0$ the quadric is $ab$, the cubic equation is $1+B\alpha^3=0$, and $\beta=-C\alpha^2$, where $C^3-B^2=g$. For tangent differentiation along a local parameter $\tau$, put $U+V\alpha^2=(\partial_\tau Q_2)(1,0,\alpha,s\alpha)$. For a regular differential $\eta=\ell(\tau)d\tau$, substitution in the packet gives three source moments proportional to $(\alpha^4,\alpha^2,1)$, with common coefficient
$$
w=\tfrac43\ell(0)\frac{\alpha^5 C(Cg'-3gC')}{U+V\alpha^2}.
$$
The normalized sum over the three sheets is their cubic algebra trace divided by three. One simple functional on the output is
$$
\left.(k_2/2-k_1)\right|_{r=0}=
-\tfrac43\ell(0)\frac{C(Cg'-3gC')}{B}
\frac{V(U+V/2)}{B^2 U^3+V^3}.
$$
For any base-field correction $(at,b,c)$ whose source extends regularly at generic $r=0$, this functional vanishes. In fact $ra=3(k_2/2-k_1)/(4\kappa s)$. Write its boundary value as $h$. The remaining rows of $M$ give $b,L$ regular there and $(rc)|_0=-6gh/g'$. Thus $(rJ)|_0=h-6(gh/g')'$. A nonzero rational solution of $h-6(gh/g')'=0$ would imply $(gh/g')^6$ is a nonzero constant multiple of $g$, contradicting squarefreeness. The latter follows, for instance, from the squarefree reduction $1+s^6$ modulo 5. The stored intrinsic tangent gives $U=8+21v$ and $V=18+5v+(17+12v)s+(20+6v)s^2$ modulo 23 with $v^2=15$. The displayed expression divided by $\ell(0)$ takes the values $11v,0,8+8v,6+10v$ at $s=1,2,3,4$.

A stronger calculation also allows all eight polynomial-exact directions. Apply $[r^0](rJ)$ and $[r^1](rJ)$ to the normalized three-sheet source. Only the two coefficients of $\eta=(u_0+u_1\tau+O(\tau^2))d\tau$ can contribute because $\tau=O(r)$ on these sheets. To see why only the first source jets enter, put $\bar a=ra$ and similarly for $L,b,c$. Then $\bar a$ depends directly on $k_2/2-k_1$, while $\bar b,\bar L$ are divisible by $r$ and their quotients depend regularly on the source. Now
$$
rJ=\bar L+\bar a+\partial_r\bar b-\bar b/r+\partial_s\bar c.
$$
Consequently an $O(r^2)$ change in the regular source changes $rJ$ by $O(r^2)$.

At $s=1,2,3,4$, each of the three seed roots is simple over $\mathbb F_{23^2}$ and the implicit Jacobian is a unit. The script solves $q=h=0$ by Newton iteration separately at all three roots in $\mathbb F_{23^2}[r,e]/(r^3,e^2)$ with $s=s_0+e$. It traces the original source formula, verifies every input against both defining equations, and computes both coefficients of $rJ$ using two distinct matrix inversions. The seven columns are $cH_a,dH_a,cH_b,dH_b$, the currents for $d\tau,\tau d\tau$, and the target. Taking both orders at $s=1,2,3$ and order zero at $s=4$ gives determinant $21+11v$ modulo 23, checked by the permutation formula. The remaining four exact directions are in the kernel of $J$. The nonzero minor lifts over the specified unit branch. Hence the target lies outside the sum of the three-sheet trace for any regular differential, all eight odd–odd polynomial-exact directions, and all base corrections of the displayed form. The full incidence trace can have additional sheets, and rational corrections can have additional dependence on $t$.

The input source is the original `four-plus-two-current` attachment by `benediamond` in [Discussion 120](https://github.com/hodge-computer/sextics/discussions/120#discussioncomment-19448881) and the Phase 4 handoff of `kmehltretter82`. Reproduce by running `correction.py`, `nodal.py`, and `local_jet.py /path/to/algebra/intrinsic_branch_four_jets.npz`. Regenerate the intrinsic jet file from the public attachment as described in the README.
