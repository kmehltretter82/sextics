import Mathlib.LinearAlgebra.Matrix.Notation
import Mathlib.LinearAlgebra.Matrix.Trace
import Mathlib.Tactic.FinCases
import Mathlib.Tactic.LinearCombination
import Mathlib.Tactic.Ring
import Mathlib.Tactic.FieldSimp

/-!
# The three sheets over the coupled nodal contact at r = 0

Writing x = α² in the cubic moving-point algebra gives x³ = k = 1/B².
The vector (x³,x²,x)/(u+v*x) multiplies a common seed residue.
The companion matrix below represents multiplication by x in the basis (1,x,x²).
-/

namespace FourPlusTwoSeedTrace

open Matrix

variable {K : Type*} [CommRing K]

def mulX (k : K) : Matrix (Fin 3) (Fin 3) K :=
  !![0, 0, k; 1, 0, 0; 0, 1, 0]

/-- Multiplication by the numerator of 1/(u + v*x). -/
def adjugateNumerator (k u v : K) : Matrix (Fin 3) (Fin 3) K :=
  !![u^2, k*v^2, -(k*u*v); -(u*v), u^2, k*v^2; v^2, -(u*v), u^2]

/-- The three regular-representation traces before division by 3 and by the norm. -/
theorem moment_traces (k u v : K) :
    let X := mulX k
    let adj := adjugateNumerator k u v
    (((X ^ 3 * adj).trace, (X ^ 2 * adj).trace, (X * adj).trace)) =
      (3*k*u^2, -(3*k*u*v), 3*k*v^2) := by
  simp [mulX, adjugateNumerator, Matrix.trace, Matrix.diag, Matrix.mul_apply,
    Fin.sum_univ_three, pow_succ]
  exact ⟨by ring, by ring, by ring⟩

/-- Clearing the inverse by the cubic norm u³ + k v³. -/
theorem clear_inverse (k u v : K) :
    (u • (1 : Matrix (Fin 3) (Fin 3) K) + v • mulX k) * adjugateNumerator k u v =
      (u^3+k*v^3) • (1 : Matrix (Fin 3) (Fin 3) K) := by
  ext i j
  fin_cases i <;> fin_cases j <;> simp [mulX, adjugateNumerator, Matrix.mul_apply, Fin.sum_univ_three] <;> ring

/-- At the central component q=r and h=1+α³ B the intrinsic top source
    S₃ simplifies to 3 α⁴ B C (2 C B' - 3 B C').  The Jacobian on a
    normalized branch is 3 α² B (u+v α²). -/
theorem top_seed_moment (a B C B' C' : K) (ha : a^3*B = -1) :
    let beta := -C*a^2
    let g := C^3-B^2
    let g' := 3*C^2*C'-2*B*B'
    let Er := -(3*a^2*B)
    let EZ := a^3*B'
    let EY := Er*(-C'*a^2) + (-(2*C*a))*EZ
    a*(3*g*beta*EY-beta^2*Er*g') - 6*g*beta^2*EZ =
      3*a^4*B*C*(2*C*B'-3*B*C') := by
  dsimp only
  linear_combination -(3*a^4*B*C*(2*C*B'-3*B*C')) * ha

section Field
variable {F : Type*} [Field F] [CharZero F]

/-- Normalized trace of the seed triple, with the common residue scalar removed.
    The companion matrix for x=α² has x³=1/B².  The inverse matrix has numerator
    adj and scalar denominator u³+k*v³ by clear_inverse. -/
theorem normalized_seed_trace (B u v : F) (hB : B ≠ 0) (hd : B^2*u^3+v^3 ≠ 0) :
    let k := 1/B^2
    let d := 3*(u^3+k*v^3)
    let X := mulX k
    let adj := adjugateNumerator k u v
    ((X^3*adj).trace / d, (X^2*adj).trace / d, (X*adj).trace / d) =
      (u^2/(B^2*u^3+v^3), -(u*v)/(B^2*u^3+v^3), v^2/(B^2*u^3+v^3)) := by
  let k := 1/B^2
  let d := 3*(u^3+k*v^3)
  let normalize (p : F × F × F) := (p.1/d, p.2.1/d, p.2.2/d)
  refine (congrArg normalize (moment_traces k u v)).trans ?_
  simp only [normalize, d, k, Prod.mk.injEq]
  exact ⟨by field_simp [hB, hd], by field_simp [hB, hd], by field_simp [hB, hd]⟩

end Field
end FourPlusTwoSeedTrace

#print axioms FourPlusTwoSeedTrace.top_seed_moment
#print axioms FourPlusTwoSeedTrace.clear_inverse
#print axioms FourPlusTwoSeedTrace.moment_traces

#print axioms FourPlusTwoSeedTrace.normalized_seed_trace
