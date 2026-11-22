"""Prepare a three-variable system for the normalized binary decomposition."""

import sympy as sp
import argparse

p = argparse.ArgumentParser(description=__doc__)
p.add_argument("--prime", type=int, default=0)
args = p.parse_args()
m, v, w, x = sp.symbols("m v w x")
k = sp.Rational(5, 7)
g = 1 + k * x + 3 * k * x**5 + x**6
f = sp.Poly(m * (1 + v * x + w * x**2) ** 3 - g, x)
b = [sp.S(1)]
for i in range(1, 4):
    b.append(
        sp.factor((f.nth(i) / (m - 1) - sum(b[j] * b[i - j] for j in range(1, i))) / 2)
    )
eq = [
    sp.factor(
        (m - 1) * sum(b[j] * b[i - j] for j in range(max(0, i - 3), min(3, i) + 1))
        - f.nth(i)
    )
    for i in (4, 5, 6)
]
print(f"ring r = {args.prime},(iv,w,v,m),dp;")
print("option(redSB);")
print(
    "ideal i="
    + ",\n".join(str(sp.fraction(a)[0]).replace("**", "^") for a in eq)
    + ", iv*m*(m-1)*w-1;"
)
print('ideal j=std(i); print("dimension and length"); dim(j); vdim(j);')
print(
    "ring l="
    + str(args.prime)
    + ',(iv,w,v,m),lp; ideal jj=fglm(r,j); print("lex basis"); jj; quit;'
)
