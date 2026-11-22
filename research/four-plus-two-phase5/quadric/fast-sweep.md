# One smooth coupled branch gives projected surface degree at least 361

For the literal $\kappa=5/7$ core, use the smooth contact through the characteristic-11 point of the earlier [split-quadric calculation](https://github.com/kmehltretter82/sextics/blob/97fbaf87b5e238d50935c8474a3e1db2f212a2e0/research/four-plus-two-phase5/quadric/axis_model.md). On the physical plane $d=2c$ use $r=b/a$, $t=c/a$ and $Z=t^2$. At the selected marked point $(r,t)=(4,6)$. The following receipt uses **65,536 coefficients** (orders $0$ through $65,535$) of each series. With $x=r-4$, this sheet starts

$$t=6+9x+5x^2+8x^3+7x^4+O(x^5),\qquad Z=3+9x+9x^2+10x^3+O(x^5).$$

Over $\mathbf F_{11}$ the 65,341 polynomials $r^it^j$ of total degree $i+j\le360$ have linearly independent images in $\mathbf F_{11}[x]/x^{65536}$ on this **single branch**. For the same branch the following rectangular bounds use $Z$ as the second coordinate.

| Maximum $\deg_Z$ | Maximum $\deg_r$ with independent monomials |
|---:|---:|
| 1 | 32,767 |
| 2 | 21,844 |
| 4 | 13,106 |
| 8 | 7,280 |
| 16 | 3,854 |

These are finite-degree injectivity results; they do not establish the minimal polynomial of $t$ or $Z$. The first possible total-degree-361 truncated relations already follow by dimension count and require higher-order testing before use even as candidates.

## Coupled series generation and checks

Retain the existing 90-coefficient quadratic presentation $B(y,y)=H$ (70 equations), its smooth initial contact $y_0$, its physical tangent $y_1$, and the known invertible 70-column minor $J_0$. Freeze the other twenty presentation coordinates to their prescribed affine-linear functions. At order $k\ge2$ the nonzero unknown coefficients are the selected 70 entries of $y_k$, found by applying $-J_0^{-1}$ to $\sum_{i=1}^{k-1}B(y_i,y_{k-i})$.

The new recurrence accumulates this convolution online. For powers of two $m\ge16$, split the smaller positive index into $[m,2m)$ and the larger into consecutive intervals $[k-m,k)$ of length $m$, with $k\ge2m$. When $m$ divides $k$ both blocks are already known before $y_k$ is computed; add their ordered convolution to orders $k$ through $k+2m-2$. Use both orders for distinct blocks and one for equal blocks. Pairs with at least one index below 16 are added directly. Every ordered index pair occurs exactly once. The implementation uses FFTs on blocks and rounds to integers. **Every coefficient of every one of the 70 original quadratic equations is independently recomputed by exact FLINT polynomial arithmetic** at the requested precision. Together with the prescribed free entries and invertible constant Jacobian this certifies the branch even independently of numerical error estimates.

Use the physical point checks from the original construction: at $c=1,d=2$ there are three rational lifted points, initially $(a,b)=(2,8),(4,7),(8,7)$. Newton lifting checks both restricted quadric/cubic equations identically modulo $\tau^{65536}$ at each. In the first chart $r=4+7\tau+8\tau^2+\cdots$ and $Z=3+8\tau+7\tau^2+\cdots$. All three $r$-series have unit linear term. Reversion uses only public python-flint APIs, with $f(g)=\sum_{j=0}^{10}g^j f_j(g)^{11}$ for coefficient decimations $f_j$. Exact truncated composition verifies the inverse identity and recovers the original $Z(\tau)$ from $Z(r)$; independently composed $t$ satisfies $t^2=Z$ at full precision.

SHA-256 of the 65,536 individual coefficient bytes (first coefficient first, each in $0,\dots,10$) on the selected branch:

| Series | SHA-256 |
|---|---|
| $t(4+x)$ | `1a7993c218efa325f340ae806b00a78ec3f1330aee06828d671dd2a4d9289e6c` |
| $Z(4+x)$ | `c39c32e005ea4447cc4b7ae69b37162345cea47e4d92978ec671f8f10a16c075` |

## Structured rank proof

For $f=t$ or $Z$, fix $m$ and the column of series $F=(1,f,\ldots,f^{m-1})^\mathsf{T}$. At order $k$ maintain a polynomial basis of the module of row vectors $P$ with $PF=0\pmod{x^k}$. Start with the identity at order zero. Let the shifted degree of a row be $\max_j(\deg P_j+s_j)$ with $s_j=j$ for total degree and $s_j=0$ for rectangular degree. Break leading-position ties using the largest column index. More than one row can have the same degree. Among rows with nonzero discrepancy $[x^k]PF$, select the least shifted degree, then the lowest original row index. Subtract the scalar multiple of that pivot that kills the discrepancy in every other row; multiply the pivot by $x$.

These rows generate the next module: the constants of the polynomial coordinates of any approximant obey exactly the one discrepancy equation, and its new pivot coefficient is divisible by $x$. The leading positions stay distinct and equal to their original row indices, by the degree and index tie rule. Thus the basis has predictable shifted degree: for a row combination the maximum of multiplier degree plus row degree is attained without cancellation at the largest leading position of the terms with that maximum. The least row degree is therefore precisely the least degree of a nonzero approximant at the prescribed order. The checker needs only the residual series and these row degrees, avoiding the dense $65536$-row interpolation matrix.

For $t$ take $m=363$, $s_j=j$. At order 65,536 the minimal shifted degree is 361, maximum 362. This includes every physical polynomial in $r,t$ of degree at most 360, since $r=4+x$ is an invertible affine substitution. The C++ program has a scalar implementation and an optional ARM NEON implementation of the same bytewise modular subtraction. An independent dense FLINT audit checks the full nullity formula $\sum_i\max(0,D-d_i+1)$ on a pseudorandom and a rational test series, two shifts and all $D=0,\ldots,22$. The array of row degrees is computed entirely in exact modular arithmetic.

## Characteristic-zero meaning and scope

At the characteristic-zero 11-adic lift of the contact, the fixed 70-column minor, the two point-equation determinants, and the linear term of the chosen $r$-parameter are units. Hence solving the same formal coefficient equations in the valuation ring yields integral jets with the computed reduction. If a nonzero characteristic-zero plane polynomial of total degree at most 360 vanished on this lifted branch, its coefficients could be normalized to be integral with nonzero reduction (over a finite extension if necessary). The resulting nonzero coefficient vector would contradict the full-rank minor modulo 11. Equivalently, that minor remains nonzero over characteristic zero and every scalar extension. Since the moving surface dominates the $(r,s)$-base, it is not the plane $d=2c$. Its defining polynomial restricts nontrivially to this plane, and a degree-$D$ homogeneous polynomial dehomogenizes at $a=1$ to degree at most $D$ in $(b,c)=(r,t)$. **The characteristic-zero projected surface on this smooth component therefore has degree at least 361.** The earlier degree-97 lower bound used all three rational plane points. This conclusion concerns the defining equation of this particular swept surface. It gives no degree bound for an arbitrary rational primitive, no complete trace and no identification with the component chosen above 23.

## Reproduction

The generators reuse `smooth_span.py`, `smooth_sweep.py`, `fiber_current.py` and `fiber_point.py` from the linked `quadric` calculation, and the upstream directory containing `contact_series_span.py` and `contact_span_coefficients.npz`. The new programs require NumPy, SymPy and python-flint. The run used NumPy 2.3.5, python-flint 0.9.0 with FLINT 3.6.0, and a C++17 compiler. At order 65,536 the online convolution took about 30 seconds and exact checking of the 70 equations and three plane points another 2.3 seconds on the test Mac. The total-degree structured check took 57 seconds.

```sh
python3 fast_sweep.py UPSTREAM_CONTACT_SPAN_DIR --helpers . -n 65536 --out sweep.npz
python3 sweep_ratios.py sweep.npz ratios.npz --sheets 0 1 2
c++ -std=c++17 -O3 sweep_pade.cpp -o sweep_pade
python3 sweep_pade.py ratios.npz "$PWD/sweep_pade" --coordinate t -n 65536 --m 363 --weight 1 --audit
python3 sweep_pade.py ratios.npz "$PWD/sweep_pade" -n 65536 --m 2 3 5 9 17
```
