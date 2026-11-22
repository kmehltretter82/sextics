# The target direction in a varying contact family

For a contact correspondence defined over an open family of cyclic sextics, three coefficient variations force its character projection to act by a scalar on the 32-dimensional space containing the prescribed numerator. Thus the canonical input on that family can be the literal monomial $A=a^3c^3$. The scalar is nonzero when the correspondence pullback is injective and the projection uses its polarization adjoint. Determining the normalized trace as a rational function and producing its exact correction remain necessary for a certificate.

## The calculation at the catalogue polynomial

Put $D=\mathbf Q[a,b,c,d]/(H_a,H_b,H_c,H_d)$ for the literal four-variable core. For $G=u^3+H$, the weighted Jacobian ring is $D[u]/(u^2)$. Let the superscript minus denote odd parity under $(c,d)\mapsto(-c,-d)$. The sign-negative pole-two and pole-three graded spaces are $D_6^-\oplus uD_4^-$ and $D_{12}^-\oplus uD_{10}^-$, of respective dimensions $32+16$ and $16+32$. These dimensions follow from the equivariant complete-intersection Hilbert series. The core is smooth over the rationals by its Fermat reduction at five. The reductions used at 23 have the same Hilbert function. [Smooth reductions](../cartier/inert-traces.md)

Changing $H$ in an invariant sextic direction $S$ differentiates $R\Omega/G^2$ to $-2RS\Omega/G^3$. The multiplications by $b^6,d^6$ together give all 48 complementary tests. In core degree 12, 208 Jacobian columns and 16 selected $b^6$ products have a $224\times224$ integer determinant congruent to 19 modulo 23. In degree 10, 108 Jacobian columns and all 32 products of $b^6,d^6$ with the odd quartics have a $140\times140$ determinant congruent to 2. One cyclic deformation direction supplies at most 16 of those 32 quartic products. The receipt specifies the monomials and both integer determinants. [Generator](certify_span.py) · [Exact receipt](span.json)

For the target deck character put $L_S:D_6^-\longrightarrow D_{12}^-$ for multiplication by $S$. A flat endomorphism preserving the Hodge filtration induces maps $X,Y$ with $Y L_S=L_S X$. This follows by applying the Gauss–Manin connection and passing to the associated graded. For the cyclic family the deck projection also commutes with the connection. Work modulo 23 and stack $P=[L_{b^6};L_{d^6}]$. Its determinant is 20 in the computed quotient bases. Write $L_{a^6}P^{-1}=[C\;D']$. The first two intertwining equations give $PXP^{-1}=\operatorname{diag}(Y,Y)$ and the third makes $Y$ commute with $C,D'$.

The kernel of $C-I$ is one-dimensional. If $v$ is its specified generator, the determinant of $[v,D'v,\ldots,(D')^{15}v]$ is 10. Consequently $Y$ acts by the same scalar on a basis, forcing $X$ to be scalar as well. The intertwining equations are a linear system integral at 23 in the chosen quotient frames. The modular nullity bounds the characteristic-zero nullity from above and scalars always solve the equations. The same maximal-rank condition holds on a Zariski open in the three coefficient directions through the catalogue point. [Explicit eigenvector and finite-field check](certify_intertwiner.py)

The calculation applies on this open to any projected algebraic correspondence extending over those coefficient directions, including after a generically finite cover in characteristic zero. The generic total norm-contact incidence provides such a correspondence after compactifying the parameter curves and restricting the base. In the normalized chart, varying the pure-core coefficients surjects onto the 25 coordinates of the invariant sextic quotient by $Q_2$. Its universal parameter incidence over the chart coefficients therefore has affine-space fibers. The smooth contact of the fixed core deforms to nearby members, and its rank-48 jet minor remains a unit on an open. Taking the spreading families containing this contact, or their sum with conjugates, gives nonzero scalar action after projection. The argument identifies the input class on the generic family. Identifying the action on the individually selected special component still requires control of the proper family over the special fiber. [Contact spanning calculation](../quadric/finite_chart.md) · [Generic parameter family and compact correspondences](../../../docs/03b-weighted-threefolds.md)

## An independent check for the full weighted family

Allow $u$ times invariant quartics to vary as well. At five the literal Jacobian ring reduces to the monomial quotient by $(a^5,b^5,c^5,d^5,u^2)$. Three explicitly listed invariant sextics give maps $M_i$ between the two 48-dimensional graded pieces. Here $\det(M_0)=1$. The operator $M_0^{-1}M_1$ has a one-dimensional eigenspace for eigenvalue 2. The spans obtained by repeatedly applying $M_0^{-1}M_1$ and $M_0^{-1}M_2$ to it have dimensions $1,2,4,8,16,31,45,47,48$. The resulting basis determinant is 4. Any simultaneous intertwiner is scalar on this line and hence on the entire space. The specialization argument again gives the characteristic-zero result. [Weighted directions and orbit calculation](weighted_intertwiner.py)

## Reproduction

With Python and python-flint, run the following from this directory.

```sh
python3 certify_span.py --json /tmp/four-plus-two-span.json
cmp /tmp/four-plus-two-span.json span.json
python3 certify_intertwiner.py
python3 weighted_intertwiner.py
```
