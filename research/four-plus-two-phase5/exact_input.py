"""Check cancellation of the two moving poles for an exact parameter differential."""

import sympy as sp

t, alpha, b, beta, phi, Pa, Bt, Ba, Bm, Pt0, Pm0 = sp.symbols(
    "t alpha B beta phi P_alpha B_t B_alpha B_mu P0_alpha P0_mu"
)
divisor = t**2 - alpha**2
ddivisor = sp.Matrix([2 * t, -2 * alpha, 0])
direction = sp.Matrix([alpha, -t, 0])
dB = sp.Matrix([Bt, Ba, Bm])
dPa = sp.Matrix([0, Pt0, Pm0])
P = Pa + divisor * b
dP = dPa + b * ddivisor + divisor * dB
# Components of the wedge product in the auxiliary coordinates t, alpha, mu.
bracket = (dP - P * dPa / Pa).cross(direction)
explicit = -2 * b * sp.Matrix([0, 0, 1]) + (dB - b * dPa / Pa).cross(direction)
assert all(sp.cancel(v) == 0 for v in bracket - divisor * explicit)
for sign in (1, -1):
    assert all(sp.cancel(y.subs(t, sign * alpha)) == 0 for y in bracket)
print("Both moving poles cancel identically in the exact-input correction.")
