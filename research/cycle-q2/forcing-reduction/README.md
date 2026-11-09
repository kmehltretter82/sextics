# A one-component forcing for the cyclic q=2 correction

This extends the exact connection presentation by `benediamond` in [PR #223](https://github.com/hodge-computer/sextics/pull/223#issuecomment-6712406587), building on `fqlx`'s unfinished [cycle draft #221](https://github.com/hodge-computer/sextics/pull/221). It supplies an equivalent correction equation and local connection checks. It does **not** construct a solution or remove the Lean admission. Prepared for `kmehltretter82` with OpenAI Codex assistance.

On the exact chart `X2=1`, set `(x,p,z,q)=(X0,X1,Y0,Y1)` and `y=Y2`. Retain all four independent base variables and put

$$
c=\frac{xz+p^3q}{2},\quad f=x^6+p^6+1+z^6+q^6+\frac{x^3pzq}{2},\quad a=x^2p^2z^2q^2,
\qquad F=y^6+cy+f,\quad g=6y^5+c.
$$

Use the previously checked paired direction `T=x∂x-z∂z`, `h=T(f)`, `h1=T(h)`, `h2=T(h1)`, and `lambda=3 h1^2-h h2`. Thus

$$h=6x^6-6z^6+x^3pzq,\qquad \sigma=\frac{a\lambda}{2h^4}.$$

The existing correction problem is `sum_j delta_j(t_j)=-sigma*y/g` in the algebra `K=C(x,p,z,q)[y]/(F)`, where `delta_j=partial_j-(F_j/g)partial_y`. The trace-zero basis is `e=(y,y²,y³,y⁴,y⁵+5c/6)`.

## Exact forcing reduction

The physical two-coordinate Jacobian satisfies

$$c_x f_z-c_z f_x=-h/2.$$

Consequently the rational base field

$$b_c=-\frac{2f_z}{h}\partial_x+\frac{2f_x}{h}\partial_z$$

has `b_c(c)=1` and `b_c(f)=0`, and its lift satisfies `delta_{b_c}(y)=-y/g`. Put `kappa=div(sigma*b_c)`. Direct product differentiation, with every coefficient and denominator derivative retained, gives the literal identity

$$\boxed{\sum_j\delta_j(\sigma (b_c)_j y)=\kappa y-\sigma y/g.}$$

This holds as a rational identity before reducing modulo `F`. Therefore translating each unrestricted trace-zero current by `t'_j=t_j-sigma*(b_c)_j*y` gives the equivalent equation

$$\boxed{\sum_j\delta_j(t'_j)=-\kappa y.}$$

Equivalently, for the unchanged rank-five connection matrices, the right side is now `-kappa*(1,0,0,0,0)^T`. There is no explicit discriminant denominator in this forcing presentation. The connection itself still has its discriminant poles, and no denominator bound for its solutions follows.

There is a compact exact expression `kappa=P/h^6`. With `M=a*lambda`,

$$P=h(-M_xf_z+M_zf_x)+5M(f_zh_x-f_xh_z).$$

The expanded `P` has 15 monomials. In particular, `h(1,1,1,1)=1` and `P(1,1,1,1)=42512272`, so the new residual has not vanished. The translation uses only powers of the same conical denominator `h`; these are legal because `h` is nonzero and independent of `y`, while `F` is monic of positive degree in `y`. It neither restricts the four-variable rational coefficients nor removes any trace-zero direction.

## Independent local connection checks

Write `Delta=46656*f^5-3125*c^6`. Independently reconstructing `B_c,B_f` from the quotient algebra verifies flatness and

$$\operatorname{tr}(B_t)=\frac{\partial_t\Delta}{2\Delta}.$$

At the generic discriminant point `c=-6r^5, f=5r^6`, `r!=0`, both transverse residues equal

$$R=v w^T,\qquad v=(1,r^{-1},r^{-2},r^{-3},r^{-4})^T,\qquad w=\frac1{30}(1,2r,3r^2,4r^3,5r^4)^T.$$

Hence `rank R=1`, `tr R=1/2`, and `R²=R/2`; `2R` is a projector. The residue of the coefficient vector of `y/g`, with transverse coordinate `f-5r^6`, is `-r*v/30`, in the image of `R`. These local identities do not prove global rational exactness.

## Reproduction

Run `python3 cycle_forcing_reduction.py` with SymPy installed. The standalone script reconstructs the universal connection independently, checks its flatness and residues, binds the actual physical coefficients, verifies the rational transgression exactly, and writes `cycle_forcing_reduction.json`. It uses exact rational arithmetic throughout and no search bounds. Checked with Python 3.14.6 and SymPy 1.14.0. No new Lean proof has been compiled, and no complete certificate is claimed.
