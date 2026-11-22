"""Search for small bidegree algebraic relations using an exact approximant basis degree audit."""

from argparse import ArgumentParser
from pathlib import Path
import subprocess
from time import monotonic
from tempfile import TemporaryDirectory
import numpy as np
from flint import nmod_poly, nmod_mat


def run(executable, where, z, m, n, weight):
    f = nmod_poly(z[:n].tolist(), 11)
    p = nmod_poly([1], 11)
    a = np.zeros((m, n), dtype=np.uint8)
    for j in range(m):
        a[j, : len(p)] = list(p)
        p = p.mul_low(f, n)
    data = where / f"powers{m}.bin"
    a.tofile(data)
    r = subprocess.run(
        [str(executable), str(data), str(m), str(n), str(weight)],
        capture_output=True,
        text=True,
        check=True,
    )
    return [int(x) for x in r.stdout.strip().split(",")], a


def audit(executable, where):
    n = 87
    m = 6
    for case, z in [
        ("coupled", np.random.default_rng(8).integers(0, 11, n, dtype=np.uint8)),
        (
            "rational",
            np.array(
                list(nmod_poly([1, 1], 11).inverse_series_trunc(n)), dtype=np.uint8
            ),
        ),
    ]:
        for weight in [0, 2]:
            deg, rows = run(executable, where, z, m, n, weight)
            for d in range(23):
                cols = []
                for j in range(m):
                    for i in range(d - weight * j + 1):
                        c = np.zeros(n, dtype=np.uint8)
                        c[i:] = rows[j, : n - i]
                        cols.append(c)
                if not cols:
                    continue
                block = np.array(cols)
                rk = nmod_mat(*block.shape, block.flatten().tolist(), 11).rank()
                # The shifted row-reduced basis gives the Hilbert function of the kernel.
                assert len(cols) - rk == sum(max(0, d - e + 1) for e in deg)
    print(
        "dense FLINT rank agrees on two controls, two shifts, 23 degree bounds",
        flush=True,
    )


def main():
    ap = ArgumentParser(description=__doc__)
    ap.add_argument("data", type=Path)
    ap.add_argument("executable", type=Path)
    ap.add_argument("--sheet", type=int, default=0)
    ap.add_argument("--coordinate", choices=["Z", "t"], default="Z")
    ap.add_argument("-n", type=int, default=16384)
    ap.add_argument("--m", nargs="+", type=int, default=[2, 3, 5, 9, 17, 33, 65])
    ap.add_argument("--weight", type=int, default=0)
    ap.add_argument("--audit", action="store_true")
    args = ap.parse_args()
    z = np.load(args.data)[f"{args.coordinate}{args.sheet}"]
    assert args.n <= len(z)
    with TemporaryDirectory() as where:
        where = Path(where)
        if args.audit:
            audit(args.executable, where)
        for m in args.m:
            t = monotonic()
            deg, _ = run(args.executable, where, z, m, args.n, args.weight)
            assert sum(deg) == args.n + args.weight * m * (m - 1) // 2
            print(
                "sheet",
                args.sheet,
                "order",
                args.n,
                f"degree in {args.coordinate}",
                m - 1,
                "shift",
                args.weight,
                "least degree",
                min(deg),
                "max",
                max(deg),
                "seconds",
                round(monotonic() - t, 3),
                flush=True,
            )


if __name__ == "__main__":
    main()
