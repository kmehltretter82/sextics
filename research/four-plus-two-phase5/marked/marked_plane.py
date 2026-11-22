"""Generate Singular checks for a complete plane of marked biquadrics."""

import argparse
import subprocess
import tempfile
from itertools import product
from pathlib import Path
from sympy import GF
from sympy.polys.rings import ring
from marked_axis import QMON, RES, seed, geometry
from marked_slice import add, neg, mul, power, value
from marked_axis_fiber import axis

p = 11


def inverse_matrix(rows):
    aug = [r[:] + [int(i == j) for j in range(3)] for i, r in enumerate(rows)]
    for i in range(3):
        k = next(j for j in range(i, 3) if aug[j][i])
        aug[i], aug[k] = aug[k], aug[i]
        u = power(aug[i][i], 119)
        aug[i] = [mul(u, v) for v in aug[i]]
        for j in range(3):
            if i != j:
                x = aug[j][i]
                aug[j] = [add(a, neg(mul(x, b))) for a, b in zip(aug[j], aug[i])]
    return [r[3:] for r in aug]


def field_string(t):
    a, b = t % p, t // p
    return str(a) if not b else f"({a}+{b}*ii)"


def seed_plane():
    b, Q, h = seed()
    F, P = geometry(b + [Q[m] for m in QMON], 4, 2, p)
    const = F[0, 0]
    R, w, z, l, m = ring("w,z,l,m", GF(p))
    old = sum(coef * w**i * z**j for (i, j), coef in Q.items())
    new = old + l * (w * z - w**2 * z**2) + m * (z**2 - w**2 * z**2)
    poly = (
        sum(coef * w**i * z**j for (i, j), coef in F.items()) + P * (old**3 - new**3)
    ) * pow(const, -1, p)
    polys = {}
    for (i, j, k, n), v in poly.terms():
        term = "*".join(
            (
                [str(int(v))]
                + (["l" if k == 1 else f"l^{k}"] if k else [])
                + (["m" if n == 1 else f"m^{n}"] if n else [])
            )
        )
        polys.setdefault((i, j), []).append(term)
    g = {k: t * pow(const, -1, p) % p for k, t in F.items()}
    f, _ = axis(p, g, True)
    T = [
        x for x in range(p * p) if value([int(f.nth(j)) % p for j in range(9)], x) == 0
    ]
    assert len(T) == 8
    return polys, g, T


def axis_coeffs(T, g):
    L = neg(g[6, 0])
    A = mul(6, add(add(mul(T, T), L), g[2, 0]))
    D = mul(2, mul(T, add(A, 1)))
    E = mul(
        add(add(mul(A, A), mul(mul(T, T), L)), add(mul(2, L), neg(g[4, 0]))),
        power(D, 119),
    )
    V = add(T, neg(E))
    Y = add(A, neg(mul(T, E)))
    assert mul(E, E) == L and add(mul(2, Y), neg(mul(V, V))) == g[2, 0]
    return V, Y, E


def script(T, polys, g):
    V, Y, E = axis_coeffs(T, g)
    ext = bool(T // p)
    ringdef = (
        "ring R=(11,ii),(l,m,v,e),dp; minpoly=ii^2+1;"
        if ext
        else "ring R=11,(l,m,v,e),dp;"
    )
    lines = [ringdef, "option(redSB);", "matrix G[7][7];", "matrix H[4][4];"]
    # First index in the two Singular matrices is z-degree; the second is w-degree.
    for (i, j), terms in polys.items():
        lines.append(f"G[{j + 1},{i + 1}]=" + "+".join(terms) + ";")
    for i, x in enumerate([1, V, Y, E]):
        lines.append(f"H[1,{i + 1}]={field_string(x)};")
    lines += ["H[2,1]=v;", "H[3,1]=(v^2+G[3,1])/2;", "H[4,1]=e;"]

    def norm(i, j):
        terms = []
        for k, a in product(range(4), repeat=2):
            if 0 <= i - k <= 3 and 0 <= j - a <= 3:
                terms.append(
                    ("-" if (i + j - k - a) % 2 else "+")
                    + f"H[{k + 1},{a + 1}]*H[{i - k + 1},{j - a + 1}]"
                )
        return "".join(terms).lstrip("+")

    def row(i, j):
        return f"({norm(i, j)})-G[{i + 1},{j + 1}]"

    eqs = [row(4, 0), row(6, 0)]
    for i in [1, 2, 3]:
        odd = i % 2
        indices = [1, 3, 5] if odd else [2, 4, 6]
        mat = (
            [[1, 0, 0], [Y, neg(V), 1], [0, neg(E), Y]]
            if odd
            else [[neg(V), 1, 0], [neg(E), Y, neg(V)], [0, 0, neg(E)]]
        )
        inv = inverse_matrix(mat)
        for j in indices:
            lines.append(f"poly k{i}{j}=-({row(i, j)})/2;")
        for idx, co in enumerate(inv):
            terms = [f"{field_string(c)}*k{i}{j}" for j, c in zip(indices, co) if c]
            lines.append(f"H[{i + 1},{idx + 2}]=" + "+".join(terms) + ";")
    eqs += [row(i, j) for i, j in RES]
    eqs += ["+".join(f"H[{i + 1},{j + 1}]" for i, j in product(range(4), repeat=2))]
    lines += ["ideal J=" + ",\n".join(eqs) + ";", CHECK.replace("NUMBER", str(T))]
    return "\n".join(lines) + "\n"


CHECK = """proc topweight(poly f) {
    int best=-1; int weightNow; intvec ex; poly topf=0;
    while(f!=0) {
       ex=leadexp(f); weightNow=ex[1]+2*ex[2]+ex[3]+3*ex[4];
       if(weightNow>best){topf=0; best=weightNow;}
       if(weightNow==best){topf=topf+lead(f);}
       f=f-lead(f);
    }
    return(topf);
}
ideal L=J[1],J[2],J[11],J[6];
ideal leading=topweight(L[1]),topweight(L[2]),topweight(L[3]),topweight(L[4]);
poly W=J[3]+2*J[4]+J[5]+7*J[7]+3*J[8]+2*J[9]+J[10];
ideal IW=W;
ideal IW2=W^2;
ideal coverBasis=std(L);
ideal C=std(J);
int dimCover=vdim(coverBasis);
int dimLeading=vdim(std(leading));
int dimW=vdim(std(L+IW));
int dimW2=vdim(std(L+IW2));
"Tw=NUMBER cover leading W W2 all";
dimCover,dimLeading,dimW,dimW2,vdim(C);
C;
quit;
"""


def run_all(dest, executable):
    polys, g, roots = seed_plane()
    for t in roots:
        source = dest / f"plane-{t}.sing"
        source.write_text(script(t, polys, g))
        if executable:
            result = subprocess.run(
                [executable, "-q", str(source)],
                capture_output=True,
                text=True,
                check=True,
            )
            lines = result.stdout.strip().splitlines()
            length = int(t == 6)
            expected = [
                f"Tw={t} cover leading W W2 all",
                f"60 60 {length} {length} {length}",
            ]
            expected += (
                ["C[1]=e+1", "C[2]=v+1", "C[3]=m", "C[4]=l"] if length else ["C[1]=1"]
            )
            if lines != expected or result.stderr.strip():
                raise RuntimeError(
                    f"Unexpected Singular output at axis {t}: {result.stdout} {result.stderr}"
                )
            print(lines[0], lines[1], flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--emit", type=Path, help="Write Singular inputs to this directory."
    )
    parser.add_argument("--singular", default="Singular")
    args = parser.parse_args()
    if args.emit:
        args.emit.mkdir(parents=True, exist_ok=True)
        run_all(args.emit, None)
    else:
        with tempfile.TemporaryDirectory() as path:
            run_all(Path(path), args.singular)
