# Homogeneous extraction for the four-plus-two core

Multiply the denominator and all four numerators by $x_0$. This preserves the rational current and makes the first nonzero homogeneous denominator component $h$ have positive degree $d$. Writing the rescaled denominator as $D=h+R$, the terms of $R$ begin in degree $d+1$.

Use $T=h^{d+2}$ and the finite sum $S=\sum_{j=0}^{d+1}(-R)^jh^{d+1-j}$. The identity $DS=T-(-R)^{d+2}$ makes the discrepancy between these two denominators start in sufficiently high degree. Homogeneous-component and derivative estimates then show that the degree-$2\deg T+6$ part of the cleared identity for $T$ and $SN$ is zero. Keeping the component of each $SN_i$ of degree $\deg T+1$ preserves precisely this identity. Since $h\ne0$ and $d>0$, the resulting witnesses satisfy every clause of the existing four-variable core certificate.

The finite geometric-sum construction comes from the Phase 4 packet. The proof adapts polynomial grading and differential lemmas from the Hodge Computer holdout work by `fqlx`. The new binding uses the mission's existing definition for its literal sextic and source. Local verification used Lean 4.33.1 and the published Mathlib revision.
