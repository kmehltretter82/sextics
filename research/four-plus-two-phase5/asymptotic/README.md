# Algebraic projections for the coupled four-plus-two instance

The two checks require Python 3 and SymPy. Run from this directory.

```sh
python asymptotic.py
python slanted_trace.py
```

* [Asymptotic cover and exact target identity](asymptotic.md). A surface times a fixed elliptic curve dominates the literal cyclic threefold after finite extension. The prescribed form becomes a computed base two-form times the elliptic differential, plus an explicit exact correction. The companion calculation verifies the finiteness of the cube lines used in the support argument.
* [Isolation of 24 pure endpoint lines](../lines/pure-lines.md). The full Grassmann chart gives 12 reduced points and 12 of local length three. These are the endpoint cases of the [global cube-line calculation](../lines/README.md).

There is also an elementary projection retaining the ratio $s$. Put $\lambda^6=-g(s)$, $p=r-\lambda t$ and $L=(6p+2\kappa)\lambda^5+2\kappa\lambda s$. The sextic $P(p+\lambda t,s,t)$ drops to degree five with leading coefficient $L$. On the generically degree-five map to $(p,s,w)$, Lagrange interpolation gives $\mathrm{Tr}(t^4/P_t|_p)=1/L$ (the denominator is the derivative in $t$ at fixed $p$). Its normalized pullback gives the exact source $(P_t+\lambda P_r)/(5L)$ with divergence-free flux $3(\lambda,0,1)/(10L)$. The residual source in the slanted coordinate has degree at most three over its fixed base after multiplying by $5L$, and its cubic coefficient is generically nonzero.

For the normalized six-direction trace define $A=6r+2\kappa$, $B=2\kappa s$, $D=6g(s)t$ and $G=g(s)$. Then

$$N=\operatorname{Norm}(D+A\lambda^5+B\lambda)=D^2(D^2+3ABG)^2+G(B^3+G^2 A^3)^2.$$

The rational flux is $(N_B,0,N_D)/(20N)$, where the partial derivatives are taken before the four displayed substitutions. Its divergence vanishes by $(\log N)_{AB}+G(\log N)_{DD}=0$. The residual still needs a correction. [Symbolic norm and derivative check](slanted_trace.py)
