# Lines on which the catalogue core is a cube

Let $H$ be the four-variable homogeneous sextic of the four-plus-two instance, and work over the complex numbers. The scheme of pairs $(\ell,[Q])$ where $\ell$ is a projective line in $\mathbf P^3$, $Q$ is a nonzero binary quadratic up to scale and $H|_\ell$ is proportional to $Q^3$, is finite of length **2,628**. Each such line determines a unique projective $[Q]$. Some of the points have multiplicity three. The degree of the moving contact-curve fiber from Phase 4 and a rational certificate for the prescribed source remain open.

Put $k=5/7$, $f(r)=1+r^6+kr+2kr^5$ and $g(s)=1+s^6+ks+3ks^5$. Write $g_{\rm bin}(c,d)=c^6g(d/c)$. In the plane $b=ra$ the equation is

$$H=f(r)a^6+g_{\rm bin}(c,d)+kr a^2cd(a^2+2c^2).$$

## Contained lines

Write a basis for each line in reduced row echelon form as a $2\times4$ matrix. In each of the six pivot charts, substitution into $7H$ yields seven coefficient equations for the free entries. Modulo 11 their Gröbner basis is the unit ideal on each chart, including over the algebraic closure. The scheme of contained lines is closed in the projective Grassmannian over $\mathbf Z_{(11)}$, so its geometric generic fiber is empty. The checker generates the full Singular input. On all six cells its output was `J[1]=1`. [Restriction equations](no_lines.py)

## Intersection with one Plücker hyperplane

Let $L$ be $a=b=0$. A line meeting $L$ lies in one of the planes $b=ra$ or $a=0$. The line $L$ itself has squarefree sextic restriction. In the remaining cases a line is invariant under reflection of $c,d$ precisely when it joins one point on each coordinate line. In a fixed plane and the chart $a=1$, a noninvariant line can be written $d=sc+v$ or $c=v$, with $v\ne0$. At $a=0$ use $b=1$.

For the finite-slope case set $c=vX$, $d=v(sX+1)$. Division by $v^6$ gives

$$g_{\rm bin}(X,sX+1)+B X^3(sX+1)+C X(sX+1)+D,$$

where $B=2kr/v^2$, $C=kr/v^4$ and $D=f(r)/v^6$. In the plane $a=0$, $B=C=0$. If $g(s)=0$, squarefreeness makes the result degree five. Otherwise its quadratic cube root must give $g(s)(X^2+\alpha X+\beta)^3$. Put $q_j=g^{(j)}(s)/j!$. Matching powers in descending order gives

$$\begin{gathered}
\alpha=q_1/(3g),\qquad \beta_0=q_2/(3g)-\alpha^2,\qquad \beta=\beta_0+sB/(3g),\\
(1-2s\alpha)B=g\alpha^3+6g\alpha\beta_0-q_3,\\
C=3g\alpha\beta^2-q_5,\qquad D=g\beta^3-1,\\
3g(\beta^2+\alpha^2\beta)-q_4-sC=0.
\end{gathered}$$

The numerator of $1-2s\alpha$ is coprime to that of the right side of its equation for $B$. After discarding the factors of $g(1-2s\alpha)$, the final equation has degree 35 and is coprime to the numerators of $B$ and $C$. This excludes both planes where $B=C=0$. In the other planes we also require $r=B^2/(4kC)$ and $DB^3=8C^3 f(r)$. This last condition is coprime to the degree-35 polynomial. For the direction $c=v$, setting $d=vX$ leaves $g(X)$ plus a polynomial of degree at most one. Its two highest nonleading coefficients force $\alpha=k$ and $\beta=-k^2$. Its cubic coefficient would then be $-5k^3$, whereas it is zero in $H$. This exhausts the noninvariant cases. [Rational gcd calculations](schubert_lines.py)

For an invariant line $(a,b,c,d)=(X,rX,Y,sY)$ with finite nonzero $r,s$, its even restriction forces an even quadratic root. The two nonzero mixed coefficients give $s=6f(r)/(kr)$ and $g(s)=8f(r)$, with restricted sextic $f(r)(X^2+2Y^2)^3$. The eliminant in $r$ has 36 simple roots, all with $rf(r)\ne0$. Thus the two tangent directions preserving the coordinate blocks give zero tangent space to these equations. For the other two directions, varying $b$ by $Y$ and $d$ by $X$ gives columns $(f'(r),ks,2ks)$ and $(kr,2kr,g'(s))$ in the three odd sextic coefficients. The odd tangent of the quadratic cube locus is spanned by $(1,4,4)$. The determinant of these three columns is coprime to the eliminant. The 36 lines are therefore isolated and reduced in the full Grassmannian. [Mixed-line check](schubert_lines.py)

There are 24 remaining invariant lines. One ratio is zero or infinite and the other is a simple root of its binary sextic. The twelve zero endpoints are reduced. At the twelve infinite endpoints six of the seven coefficient equations formally solve for six chart variables, and the remaining equation has exact order three in the last variable. The unit test at order three works simultaneously over all the complementary roots. [Local equations and multiplicities](pure-lines.md)

The lines meeting $L$ form the Plücker hyperplane $p_{ab}=0$. Every positive-dimensional projective component meets this hyperplane, but all 60 points found there are isolated even in the full cube locus. The latter is closed by the projective incidence described below. Consequently the entire cube locus is finite.

## Total length

On $G=\operatorname{Gr}(2,4)$ take the tautological rank-two bundle $U$. Let $V=\operatorname{Sym}^2 U^*$ and $W=\operatorname{Sym}^6 U^*$. On the bundle $\pi:\mathbf P(V)\to G$ of one-dimensional subspaces, cubing the tautological quadratic gives a subbundle $\mathcal O(-3)\subset\pi^*W$. The section $H$ of the rank-six quotient cuts out the pairs in the statement. Because no line is contained in the surface and a nonzero cube has a unique projective root, its zero scheme is finite. On the smooth sixfold each local ideal is a complete intersection. Its total length is the degree of the top Chern class. [Stacks Project, Chern classes and sections](https://stacks.math.columbia.edu/tag/0FA8)

The projective bundle pushforward gives $\sum_{j=2}^6 3^j c_{6-j}(W)s_{j-2}(V)$, using $s(V)=c(V)^{-1}$. In the two Chern roots $x,y$ of $U^*$ this is

$$360(x^4+y^4)+2772(x^3y+xy^3)+5400x^2y^2.$$

Integration on $G$ extracts the coefficient of $x^3y^2$ after multiplication by $x-y$. This gives $5400-2772=2628$. Equivalently the top Schur class on $G$ is $s_{(2,2)}=x^2y^2$, and the coefficient extraction annihilates the Schur classes with first part greater than two. [Exact Chern calculation](chern.py)

From this directory run the following with SymPy and Singular installed.

```sh
python3 schubert_lines.py
python3 pure_lines.py
python3 chern.py
python3 no_lines.py /tmp/hodge-cube-lines.sing
Singular -q /tmp/hodge-cube-lines.sing
```
