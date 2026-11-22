"""Track one monodromy cycle without treating it as the full fibre."""

import sys, json
from pathlib import Path
from loops import advance, route, vec, identify


def run(base, initial, c, rad, limit, out):
    d = json.load(open(base))
    s = d["s"]
    old = vec(d["branches"][initial]["end"])
    pts = [old]
    logs = []
    failure = None
    path = route(c, rad)
    for k in range(limit):
        try:
            w, stat = advance(old, path, s, 0.01)
        except ArithmeticError as e:
            failure = str(e)
            print("stop", failure)
            break
        i, dist = identify(w, pts)
        gap = min(abs(w[-1] - z[-1]) for z in pts)
        print(
            "power",
            k + 1,
            "alpha",
            w[-1],
            "distance",
            dist,
            "cond",
            stat["cond"],
            flush=True,
        )
        logs.append(
            dict(
                distance=dist,
                nearest=i,
                alpha_separation=gap,
                stats=stat,
                z=[[v.real, v.imag] for v in w],
            )
        )
        if dist < 1e-7:
            break
        pts.append(w)
        old = w
    Path(out).write_text(json.dumps(dict(reached=logs, stop=failure), indent=2) + "\n")


if __name__ == "__main__":
    run(
        sys.argv[1],
        int(sys.argv[2]),
        complex(sys.argv[3]),
        float(sys.argv[4]),
        int(sys.argv[5]),
        sys.argv[6],
    )
