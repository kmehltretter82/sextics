# A slope-two component of the prescribed class

For the weighted threefold $T:w^3+H=0$, the de Rham residue with numerator $A=a^3c^3$ has a crystalline Frobenius slope-two component at each of $p=11,17,23$. In particular this specific class cannot come entirely from Tate twists of elliptic curves with CM by $\mathbf Q(\zeta_3)$ via algebraic correspondences defined over number fields. Each calculation uses a 48-dimensional matrix over the prime field.

## Cartier on the projective complement

Use the ordinary sextic $Z:F=e^6+H(a,b,c,d)=0$ in $\mathbf P^4$. The map $Z\to T$ given by $w=e^2$ pulls the weighted class back to the residue of $2e a^3c^3\Omega/F^2$, with $\Omega$ the standard homogeneous projective volume form. All three primes are good for $Z$, since $6e^5$ forces $e=0$ at a singular point and the Jacobian of $H$ spans all monomials of degree 17 modulo each prime. [Smoothness verification](inert-traces.md)

Let $U=\mathbf P^4\setminus Z$. Its middle cohomology is $H^3_{\rm prim}(Z)(-1)$ under residue. In fact the map $H^2(Z)(-1)\to H^4(\mathbf P^4)$ in localization is surjective and $H^5(\mathbf P^4)=0$. On normalized projective coordinates $|F(x^p)-F(x)^p|\le p^{-1}$. For $p^{-1}<\eta<1$, the inverse image of $|F|>\eta$ under coordinatewise $p$th power is $|F|>\eta^{1/p}$, giving the finite trace on strict neighborhoods of the tube. On the $p$-adic tube of $U$ the finite map $\phi:[x_i]\mapsto[x_i^p]$ lifts Frobenius, since $F(x^p)\equiv F(x)^p\pmod p$. The pushforward on top forms is coefficient extraction $\mathcal C_p$ in logarithmic coordinates. For a Laurent monomial on the dense torus this follows by summing over the $p^4$ preimages, with the Jacobian $p^{-4}$ from four logarithmic differentials. The trace formula agrees on the entire projective open for forms regular there. Equivalently $\mathcal C_p(h(x^p)g(x))=h(x)\mathcal C_p(g(x))$, and $\mathcal C_p$ keeps only monomials whose exponents are multiples of $p$, dividing the exponents by $p$. On degree-four cohomology $\phi_*\phi^*=p^4$. Hence on the primitive residue classes the operator

$$V=\mathcal C_p/p=p^2 F_Z^{-1}$$

uses the crystalline Frobenius $F_Z$, with the extra factor $p$ in $F_U$ accounted for by the residue Tate twist. The calculation uses the projective rigid cohomology of the complement. Torus coordinates are a convenient way to evaluate its trace map. [Frobenius-equivariant residue, Proposition 2.4.1(b)](https://math.mit.edu/~roed/writings/research/picard.pdf) · [Cartier on rational forms](https://doi.org/10.1093/imrn/rnaa119)

Write $d\log x$ for the four-dimensional torus volume in the chart $x_5=1$. Interior monomials are the normalized forms

$$(k-1)!\,x^u F^{-k}\,d\log x,\qquad u_i\ge1,\quad \sum u_i=6k.$$

Their projective numerator is $x^{u-\mathbf1}$. Use the sector $u_c+u_d$ odd and $u_e\equiv2,4\pmod6$. It contains the pullback of $A$ with $u=(4,1,4,1,2)$, up to the unit 2. Taking the direct sum of both cubic characters makes it stable under $F_Z$ at a prime inert in $\mathbf Q(\zeta_3)$.

The $\mathbf Z_p$-span of the classes of the normalized forms is a lattice. Smoothness gives a regular sequence of the five degree-five partials modulo $p$. The graded Jacobian quotient lifts with a $\mathbf Z_p$-basis in every degree: the lifted regular sequence has the same Hilbert function, so representatives of a modular basis span by Nakayama and are independent over $\mathbf Q_p$. Every Jacobian multiple at pole order $k$ reduces with integral coefficients, since the Griffiths denominator $k-1$ cancels against $(k-1)!$. Larger pole orders reduce to finitely many degrees. The resulting primitive basis consists of complementary Jacobian monomials at each pole, by Griffiths reduction over $\mathbf Q_p$. Consequently a nonzero modular vector in the displayed basis detects a nonzero reduction in this lattice.

Put $s=\lceil k/p\rceil$ and $G=(F(x)^p-F(x^p))/p$. To apply $\mathcal C_p$ to $x^u/F^k$, multiply the numerator by $F^{ps-k}$ and expand $(F(x^p)+pG)^{-s}$. After normalization of the output at pole order $s+i$, the scalar of term $i$ is

$$\frac{(k-1)!}{(s-1)!}\frac{(-p)^i}{i!}.$$

The series converges $p$-adically, including after integral pole reduction. For $s\ge2$ the factorial ratio is a multiple of $p$. For $s=1$, the term $i=0$ vanishes identically: in this sector no strictly positive exponent vector of total degree six exists. Indeed even $u_e$ forces $u_e=2$ and the other four coordinates to be one, giving even $u_c+u_d$. All other terms contain a factor $p$. This proves $V$ preserves the lattice. At $k=2$ the terms $i\ge2$ vanish modulo $p$ after division by $p$. The $i=1$ term reduces to

$$V\left(\frac{x^u}{F^2}\right)\equiv -\frac1p\frac{\mathcal C_p(x^u F^{2p-2})}{F^2}\pmod p,$$

since $\mathcal C_p(x^u F^{p-2})=0$. In particular $V$ acts modulo $p$ within the pole-two part. This gives a matrix calculation with ordinary integer polynomial coefficients modulo $p^2$.

## Calculation on the literal numerator

At pole two, $u_e=2$ or $4$. The first block consists of 40 odd-odd core sextics; its eight Jacobian directions are $x_iH_j$ with $i,j$ from opposite variable pairs. The second consists of 16 odd-odd core quartics and has degree below the core Jacobian ideal. There is no lower pole in this sector. The script verifies that the eight relations are independent and preserved by the coefficient map. It then applies the resulting quotient operator 48 times to the coordinate vector of $a^3c^3e$.

| $p$ | $\operatorname{Tr}(V^2)\bmod p$ | $\operatorname{rank}(V^{48})$ | First coordinate of $V^{48}(A)$ |
|---:|---:|---:|---:|
| 11 | 8 | 30 | 8 |
| 17 | 16 | 32 | 7 |
| 23 | 21 | 32 | 19 |

The coordinate refers to the monomial and row-reduction orders in the script and omits the unit 2. The trace values independently agree with the direct counts over $\mathbf F_{p^2}$ divided by $p^2$, as predicted by Poincaré duality. The full lattice has a basis of pole orders at most four, so the same expansion shows the image of $V$ modulo $p$ lies at pole two. Run `python3 target_cartier.py 11 17 23` with NumPy. [Class-specific computation](target_cartier.py) · [Independent point counts](inert-traces.md) · [Independent multinomial and FLINT check](check_coeff.py) · [Check against canonical contact differentials](../quadric/cartier_check.md)

In a 48-dimensional vector space the image of $V^{48}$ is stable and $V$ is invertible on it. Since $V^{48}(A)\ne0$, its iterates never vanish modulo $p$. Thus the $p$-adic class has a $V$-unit-slope component (slope zero for $V$, two for $F_Z$). A class entirely in positive $V$ slopes would approach zero $p$-adically under iteration. Because the double cover is defined over the good integral model, the result applies to the original weighted class as well.

At these inert primes the CM elliptic curve with $j=0$ is supersingular. An algebraic correspondence from any power or twist of this elliptic curve has image of pure slope $3/2$ after the Tate shift. This remains true over any finite extension at $p$. For a correspondence defined over a number field, the functorial $p$-adic comparison identifies its de Rham image with that of the crystalline Galois representations of its good endpoints, even if a chosen representative of the correspondence has bad reduction. The slope-two component of the prescribed class rules out an expression of this type.
