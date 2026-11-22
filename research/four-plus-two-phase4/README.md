# Coupled trace calculations for four-plus-two-residual

For the prescribed source $A=a^3c^3$, where $(a,b,c,d)=(X_0,X_1,Y_0,Y_1)$, the ratio chart requires the source triple $(0,0,1)$ in powers $(1,t^2,t^4)$. We computed the normalized trace on the three sheets through one nodal contact of the coupled core. The generic fiber on the selected parameter curve has degree at least 48, using the determinant below. The full trace, the constant combination of closed inputs, and a global rational correction to the prescribed source remain open.

These calculations extend the exact norm-contact family and current formulas of `benediamond`, the earlier core and suspension analyses by `kejace` and `DiamondJM`, and the Phase 2–3 work of `kmehltretter82` with Antigravity and ChatGPT. Codex prepared this continuation for `kmehltretter82`. [Original contact data and formulas](https://github.com/hodge-computer/sextics/discussions/120#discussioncomment-19448881) · [Phase 2–3 discussion](https://github.com/hodge-computer/sextics/discussions/120#discussioncomment-19791876)

## Binary control and coupled trace

For the uncoupled control $P=C^3+R^2+(D^3+S^2)t^6$, the normalized field trace is $(0,0,-4iJ_XJ_Y/9)$ with $J_X=2CR_r-3RC_r$ and $J_Y=2DS_s-3SD_s$. Tracing the three flux components and differentiating gives the same source. The two lower moments vanish in this control.

On the coupled nodal branch at $r=0$, write $C=C(1,s)$, $B=B(1,s)$, $g=C^3-B^2$, and $D_\tau Q_2=U+V\alpha^2$ at the moving point. For $\eta_R=\ell_R(\tau)d\tau$, the mean over all three roots of $1+B\alpha^3=0$ is

$$
\frac{4\ell_R(0) C(2CB'-3BC')}{3(B^2U^3+V^3)}
       (U^2,-UV,V^2).
$$

The formula holds on the open set where $B(B^2U^3+V^3)\ne0$, with the sign convention of the contact packet. The Lean theorems verify the numerator simplification and normalized companion-matrix traces. In a full trace normalized by a degree $d$, these three sheets carry weight $3/d$. [Lean proof](LocalTrace.lean) · [Symbolic checks](binary_and_resultant.py)

For $R=A$ we also computed the Taylor coefficients through $r^3$ at $s=1,2,3,4$ on these sheets, over $\mathbf F_{23}[v]/(v^2-15)$. At $(r,s)=(0,1)$ the source is $(13+15v,7+5v,8+5v)$. At $(0,2)$ all three entries vanish modulo 23 and their $r$ coefficients are $(13+15v,20+16v,9+v)$. The canonical input was obtained from the degree-nine complete-intersection residue, differentiating the cubic before eliminating the weighted variable. Its four coefficients are $5+5v,18+2v,20+2v,9+18v$ in ascending powers of $\tau$. The independently supplied obstruction rows agree through that order up to a single common unit for all 40 polynomial input monomials. [Deformed trace calculation](deformed_trace.py)

## Tests of source matching

A rational correction with chart components $(q_t,q_r,q_s)=(ta(r,s),b(r,s),c(r,s))$ is determined by a $4\times4$ rational matrix together with one differential condition. Even allowing any regular input differential on the chosen branch and all eight odd–odd polynomial-exact directions, the three-sheet contribution plus corrections of this form cannot give the protected source. A $7\times7$ minor reduces to $21+11v$ at the specified place. It lifts by the unit checks for the published branch. [Derivation and hypotheses](quotient.md) · [Independent matrix and jet checks](local_jet.py)

General rational corrections can match any finite jet at this line. For example, for a residual $r^m U(s,t)+O(r^{m+1})$, add $q_r=-r^{m+1}U/((m+1)P(0,s,t))$ with the other components zero. In characteristic zero its source has leading term $r^m U$. A global rational identity needs control away from the line. In fact each of the three choices of a single nonzero field component is excluded globally for the target, even allowing arbitrary rational dependence on all three variables. The valuation argument and three exact determinants are in the correction notes.

## The full fiber

Eliminating the moving point yields a bidegree-$(6,6)$ resultant on the normalized contact curve. Its 49 coefficients are independent. The determinant of their evaluations on a $7\times7$ grid of ratios against Taylor orders 0–48 reduces to $15+5v$. The associated moving divisor has degree at least 48. On $r=0$ at least six further degrees remain after removing the fixed contacts. [Fiber-degree argument](global-trace-notes.md) · [Determinant check](fiber_degree.py)

The full source sum can be calculated using residues of a meromorphic differential on the parameter curve. At $(r,s)=(4,6)$ the chosen seed contributes $(16+8v,19+8v,11+22v)$ modulo 23 for input $A$, with the current factor $4/3$ and before degree normalization. This is one boundary residue for the specified representative. A second representative, computed modulo the cubic moving equation, is regular at that seed for generic ratios. A Bézout identity combines the two representatives and confines all complementary poles to a finite set of parameter points independent of the ratios. Over the affine normalized chart, the unit-ideal argument further confines the required poles to its compactification boundary. Their principal parts remain to be computed. Two $52\times52$ minors also exclude absorption of the three reduced source numerators by the tested 49-dimensional space of sections. [Boundary calculation](test_residue_span.py) · [Cubic representative and Bézout patch](cubic_residue.py)

A direct Gröbner basis attempt at $(4,6)$ uses 19 auxiliary quartic coefficients to preserve the sparsity of the 25 norm equations and one inverse for the chart localization. This gives 47 equations in 47 variables over $\mathbf F_{23}$. The run of Singular 4.3.1 with `slimgb` reached its 120-second limit before returning a basis. Exact substitution with `elimpart` removed 16 variables. On the resulting 31-equation system, `slimgb` reported insufficient memory under a 5,500,000-KiB virtual-memory limit. On the same 31-variable system, `std` with a 12-GiB ceiling reached 150 seconds without a basis. This system contains all components in the stated chart. For the normalized full trace, one must isolate the selected component or account for its remaining boundary residues. [Singular input generator](global_fiber.py)

## Reproduction

Use Python 3 with NumPy and SymPy. The public source ZIP has SHA-256 `f73affdd3a4317b575d0417cf60675933ccca37385bc4ea59071b1586e92b48c`. The intrinsic four-jet file is generated by the upstream script. [Download the public contact packet](https://github.com/user-attachments/files/34648619/four-plus-two-current.zip)

Run the following from this directory with the extracted packet alongside it.

```sh
packet=./four-plus-two-current
python3 "$packet/algebra/intrinsic_norm_incidence.py"
python3 binary_and_resultant.py
python3 cubic_residue.py
python3 deformed_trace.py "$packet"
python3 fiber_degree.py "$packet/contact-span/contact_span_coefficients.npz"
python3 test_residue_span.py "$packet"
python3 correction.py
python3 nodal.py
python3 local_jet.py "$packet/algebra/intrinsic_branch_four_jets.npz"
python3 global_fiber.py "$packet/algebra/norm_system_circuit.json" > global-fiber.sing
```

From the research branch's `formal/lean` directory, run `lake env lean ../../research/four-plus-two-phase4/LocalTrace.lean`. Local and Linux VM checks used Lean 4.33.1 and Mathlib commit `0df444a360eaa60ab8c11dca51a86af692955474`. All four printed axiom lists are contained in `propext`, `Classical.choice` and `Quot.sound`.
