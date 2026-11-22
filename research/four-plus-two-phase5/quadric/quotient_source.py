"""Check the source after replacing the moving cubic by the even norm equation."""

import sympy as sp

alpha, scale, ell, J, er, es, ez, dcr, dcs, dcz, p, pr, ps, i = sp.symbols(
    "alpha L ell J er es ez c_r c_s c_z p p_r p_s i", nonzero=True
)
qr, qs, qz, dq, hr, hs, hz, dh, multiplier, conjugate = sp.symbols(
    "q_r q_s q_z Dq h_r h_s h_z Dh G hbar"
)
# Differentials at q=h=0 satisfy dk=G*dq+hbar*dh.
dotk, kr, ks, kz = [
    multiplier * a + conjugate * b for a, b in ((dq, dh), (qr, hr), (qs, hs), (qz, hz))
]
for left, right in [
    (dq * kz - dotk * qz, dq * hz - dh * qz),
    (qz * kr - kz * qr, qz * hr - hz * qr),
    (qz * ks - kz * qs, qz * hs - hz * qs),
    (qr * ks - qs * kr, qr * hs - qs * hr),
]:
    assert sp.expand(left - conjugate * right) == 0

# Here er,es,ez are the minors for k in coordinates (r,s,z).
beta = -scale
Er, Es, EZ = 2 * alpha * er, 2 * alpha * es, ez
EY = Er * (-scale * dcs) - Es * (-scale * dcr) + (-2 * alpha * scale * dcz) * EZ
S = (
    alpha * (3 * p * beta * EY - beta**2 * (Er * ps - Es * pr))
    - 2 * i * p * beta**2 * EZ
)
Z = alpha**2
T = Z * (3 * p * (er * dcs - es * dcr + dcz * ez) - (er * ps - es * pr)) - i * p * ez
assert sp.expand(S - 2 * scale**2 * T) == 0
old = ell * S / (2 * alpha * J)
reduced = (ell * scale**2 / alpha) * T / J
assert sp.cancel(old - reduced) == 0
print(
    "The four determinant substitutions and normalized quotient source identity hold."
)
