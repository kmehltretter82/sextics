# Complex continuation and certified source sums

At `(r,s)=(1/10,1)`, the original 27-equation circuit has 13 distinct high-precision numerical solutions recorded here. The first six start from one characteristic-zero nodal contact. Seven further solutions arise by continuing around closed paths in the complex `r`-plane. The largest absolute residual among the 13 refined solutions is `2.34e-71`; the minimum separation of their moving coordinates `alpha` is `0.00813`. All thirteen endpoints also have disjoint rigorous root enclosures. The orbit and the generic degree of the covering remain undetermined. An average over these solutions cannot yet be used as its normalized field trace.

## Root and residue certificates

Run the following in this directory with python-flint, NumPy and SymPy installed.

```sh
python3 certify_roots.py endpoints.json
python3 current.py endpoints.json --transport-test
```

Working with the 44 coefficients of the norm equation and the two moving equations gives a square polynomial system in 46 unknowns after introducing 19 quartic coefficients. At each center $z$, use a midpoint approximate inverse $R$ of its Jacobian. With complex sup radius $10^{-44}$, ball arithmetic certifies $\|R F(z)\|+10^{-44}\sup\|I-RJ\|<10^{-44}$ and $\sup\|I-RJ\|<1$. The complex moving coordinates separate all thirteen enclosures. These inequalities prove endpoint existence and uniqueness. The larger numerical loops used for discovery have empirical forward and reverse checks.

For each root we solve for the tangent by implicit differentiation. On the complete intersection $(Q_2,h,3Q_1^2,G_2)$ we choose 219 rows from the degree-nine ideal and append the polynomial Jacobian as the 220th row. An interval linear solve certifies invertibility. This also shows there is no common projective zero of the four forms, since at such a point Euler's formula would make the Jacobian vanish as well. The residue functional takes value 96 on the Jacobian. Applying it to each $R\,Dh$ gives a verified enclosure of the source from the current formula at each actual root.

For the unnormalized sum on these thirteen sheets at $(1/10,1)$, the columns for $bd^5,bcd^4,bc^2d^3$ have determinant approximately $-0.000768081286+0.000477629501i$, with zero excluded. Their unique combination giving $(0,0,1)$ there has coefficients approximately

$$(-2.28661001-0.76186008i,\quad -1.84441717+0.00887698i,\quad 0.48501918-3.08350003i).$$

To test constancy we enclose the same sheets along the segment to $r'=1/10+10^{-14}$, $s=1$. On a complex ball of radius $10^{-8}$ around each initial center, the contraction estimate holds uniformly for the entire parameter segment. A second radius-$10^{-44}$ enclosure around each terminal root lies in its original larger ball. With the original coefficients fixed, the real part of the first error at $r'$ lies between $1.45406506\cdot10^{-14}$ and $1.45406508\cdot10^{-14}$. Thus this three-input combination for this finite sum cannot give the prescribed source identically. The other sheets, other closed inputs and rational corrections are outside this test.

## Seed and continuation

The distinguished root `1<m<2` of the degree-40 binary shape polynomial gives the real contact `C³−B²=g`. The code selects the positive real roots `C(1,0)=m^(1/3)` and `B(1,0)=sqrt(m−1)`. At the central point the norm constraints are solved by `Q₁=C, Q₂=ab, A=a³+b³` and this binary `B`.

Adjoin the 19 invariant quartic coefficients of `T` and write all 44 coefficients of `H−Q₁³−A²+B²−Q₂T`. Their derivative in the 45 coefficients has numerical rank 43 at this contact. Project the quadratic part of the expansion onto the one-dimensional left cokernel and restrict it to the two-dimensional right kernel. Its two distinct tangent lines are ordered by the chart-coordinate ratio `δλ₂/δλ₀`, approximately `2.27706579` and `−1.87849130`.

At `s=1, r=0`, each tangent has three possible moving coordinates satisfying `1+B(1,1)alpha³=0`. Differentiate the two point equations to normalize the tangent using `dr=1` and compute `dalpha/dr`. Newton correction at `r=.002` gives six regular solutions. Continuing to `r=.1` gives records 0–5 below. The expanded 46-variable system uses the 26 original coefficients, the 19 invariant quartic coefficients and `alpha`. Eliminating the quartic through the exact packet's 19 pivot rows recovers the 27-equation system.

For the loops, the path begins at `.1`, travels via `.1+.06i` to `z+.06i` and then to `z=center+radius*exp(i*angle)`. It traverses 40 chords counterclockwise and returns along the same stem. Each accepted step bounds the correction from the tangent predictor by 8% of the predicted movement, halves rejected steps, and records residual, Jacobian condition and localization determinant. These diagnostics and numerical reversals test the paths. They do not certify the homotopy classes.

| Endpoint index | Seed index | Centre in the r plane | Radius | Endpoint alpha | Round-trip error |
|---:|---:|---|---:|---|---:|
| 6 | 3 | 0.4 | 0.38 | −0.256020430686 − 0.096313655189i | 4.42e−15 |
| 7 | 1 | 0.15 + 0.28i | 0.38 | 0.211129837442 − 0.407898897892i | 2.40e−15 |
| 8 | 4 | 0.25 + 0.20i | 0.32 | −0.041922820275 − 0.605996033117i | 2.67e−15 |

The stem angle here is `0.22`. Repeating each forward path with a step limit equal to 55% of the original reproduced all three endpoints within `2.19e−13` in the full coordinates. For the middle example the largest condition number was 34,672 and every step satisfied `|Δ|>0.0147`. Endpoints 9–12 arose from powers 2–5 of the middle loop. In the recorded double-precision run the sixth power reached the minimum step bound. An alternative based loop with angle `1.31` from the fifth iterate reaches `alpha≈−0.698420788842+0.192060197718i`; the reverse and finer forward discrepancies were below `2.3e−14`. This additional endpoint is separate from the fixed file of 13.

## Reproduction

Use Python 3 with NumPy, SymPy, mpmath and python-flint. Set `BASIS` to the characteristic-zero Singular output generated by the phase-5 binary shape calculation. Set `CIRCUIT` to the phase-4 packet's `algebra/norm_system_circuit.json`. Run these commands from this directory, with the phase-5 shared `verify_rur.py` in its parent directory.

```sh
python3 continue_incidence.py "$BASIS" --output continuation.json
python3 loops.py continuation.json 3 0.4 0.38 --output loop0.json
python3 loops.py continuation.json 1 0.15+0.28j 0.38 --output loop1.json
python3 loops.py continuation.json 4 0.25+0.20j 0.32 --output loop3.json
python3 cycle.py continuation.json 1 0.15+0.28j 0.38 5 cycle1.json
python3 verify_packet.py "$CIRCUIT" continuation.json reproduced.json loop0.json loop1.json loop3.json cycle1.json
python3 loops.py continuation.json 10 0.15+0.28j 0.38 --additional cycle1.json --angle 1.31 --step .01 --output alternate-stem.json
```

`endpoints.json` stores indices 0–12 as ordered above. Each record contains all 27 coordinates in the packet's order as pairs of decimal strings, along with `alpha`, `beta`, the determinant of the 27-by-27 Jacobian, the magnitude of the quartic pivot determinant, and the residuals. Refinement uses 72 decimal digits and cross-checks `β³+P(α)=0`. The written coordinates retain 57 significant digits.

The one-thread Google Cloud VM benchmark for the middle loop took 7.74 wall-clock seconds for the forward path, reversal and finer forward repetition together. On that host there were 523 accepted forward steps, one rejection, maximum residual `2.87e−14`, round-trip error `8.46e−15` and fine forward discrepancy `1.17e−14`. BLAS implementations can choose different adaptive step sequences.

To sample a bounded orbit on CPUs, start three single-threaded processes. The queue writes accepted endpoints, the graph of tested loop transitions, sample reverse checks, chart diagnostics, failures and the applied limits. The discovery order can vary with worker scheduling.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 orbit_queue.py continuation.json orbit.json loop0.json loop1.json loop3.json cycle1.json alternate-stem.json --workers 3 --tasks 240 --seconds 225
```

## Bounded VM orbit and interval verification

The three-worker run used 223.90 seconds, starting from the six radial solutions and eight subsequently discovered points. Of 240 submitted jobs, 153 paths completed, 35 reached the minimum Newton step, two hit the residual or chart guard, and 50 hit the time guard. Fourteen sample reversals returned to their starting points to within `2.40e−13`. The 153 completed paths gave 115 different endpoints including the starts. The raw graph records each tested generator, source, destination and error or diagnostics. [Orbit graph](orbit-vm.json.gz)

The shared Acb verifier encloses each of the 115 roots in a complex sup-norm ball of radius `1e−44`. The largest Newton displacement bound is `2.38e−57`, the largest Lipschitz bound is `7.13e−38`, and every quartic chart determinant has magnitude above `0.00117`. All 115 moving coordinates have pairwise disjoint enclosures. Complex conjugation gives at least 220 distinct values. The symmetry negating `alpha` and the eight coefficients of `B` fixes the point and norm constraints at `r=1/10, s=1`. Applying both symmetries gives at least **440 distinct moving coordinates** on the full specialized chart. Their membership in a single irreducible component remains open. [Certified coordinates](orbit-roots.json.gz) · [Interval summary](orbit-certificate.json)

The displayed coefficient system also admits `Q₁→ζQ₁` for each cube root of unity and simultaneous `(A,B)→(−A,−B)`. These two actions fix `alpha` and independently preserve the equations. Every certified box has a nonzero coefficient of `Q₁` and a nonzero coefficient of the cubic pair. The six coefficient tuples at each of the 440 distinct moving coordinates therefore give at least **2,640 coefficient lifts in the unquotiented specialized chart**. A component degree used in the normalized trace would require a separate argument.

To replay the interval checks on the recorded run, use the companion `certify_roots.py` in this directory.

```sh
python3 refine_orbit.py orbit-vm.json.gz verified-roots.json.gz verified-summary.json
```

Global completeness would require a degree or fibre argument for the relevant localized irreducible component and verification that its sheets have all been accounted for. The degree remains open.

The contact and trace circuit come from `benediamond` and the four-plus-two discussion. The code uses the exact binary shape coordinates in the enclosing phase-5 packet. [Contact family and chart](https://github.com/hodge-computer/sextics/discussions/120#discussioncomment-19448881)
