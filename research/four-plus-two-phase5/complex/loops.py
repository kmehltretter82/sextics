"""Complex-parameter predictor-corrector around polygonal loops and their reversals."""

import json, argparse, time
from pathlib import Path
import numpy as np
from continue_incidence import rhs, jac, point


def vector(z, r, s):
    return np.linalg.solve(jac(z, r, s), -np.r_[np.zeros(44), point(z, r, s, dr=1)])


def advance(z, raw_path, s, maxstep=0.012, minstep=5e-8, deadline=None):
    z = z.copy()
    r = complex(raw_path[0])
    step = maxstep
    n = 0
    reject = 0
    stat = [0.0, 0.0, 0.0, 10.0, 0.0]
    for target in raw_path[1:]:
        while abs(target - r) > 3e-13:
            if deadline is not None and time.monotonic() > deadline:
                raise ArithmeticError("time limit")
            if n + reject > 30000:
                raise ArithmeticError("step budget")
            dr = (target - r) * min(1.0, step / abs(target - r))
            v = vector(z, r, s)
            guess = z + dr * v
            w = guess.copy()
            rr = r + dr
            conv = False
            for nit in range(7):
                e = rhs(w, rr, s)
                J = jac(w, rr, s)
                try:
                    corr = np.linalg.solve(J, -e)
                except np.linalg.LinAlgError:
                    break
                w += corr
                if np.max(abs(corr)) < 1e-12 * (1 + np.max(abs(w))):
                    conv = True
                    break
                if np.max(abs(corr)) > max(1, 5 * np.max(abs(dr * v))):
                    break
            err = np.linalg.norm(w - guess) / (1e-14 + np.linalg.norm(dr * v))
            if not conv or err > 0.08:
                step = abs(dr) * 0.5
                reject += 1
                if step < minstep:
                    raise ArithmeticError((r, step, float(err), nit))
                continue
            e1 = max(abs(rhs(w, rr, s)))
            step = min(maxstep, abs(dr) * (1.45 if err < 0.018 else 1))
            J = jac(w, rr, s)
            mult = w[5] * w[6]
            delta = abs((1 - mult) * (1 - 3 * mult) * (1 - 2 * mult) ** 3)
            stat = [
                max(stat[0], e1),
                max(stat[1], np.linalg.cond(J)),
                max(stat[2], err),
                min(stat[3], delta),
                max(stat[4], max(abs(w))),
            ]
            if e1 > 2e-9 or delta < 1e-6:
                raise ArithmeticError(("residual or chart pole", e1, delta, rr))
            z = w
            r = rr
            n += 1
    return z, dict(
        steps=n,
        rejects=reject,
        residual=stat[0],
        cond=stat[1],
        max_predictor_fraction=stat[2],
        min_chart_delta=stat[3],
        max_coordinate=stat[4],
    )


def route(center, rad, theta=0.22):
    v = rad * np.exp(1j * theta)
    first = center + v
    base = 0.1 + 0j
    return (
        [base, base + 0.06j, first + 0.06j, first]
        + [center + v * np.exp(2j * np.pi * i / 40) for i in range(1, 41)]
        + [first + 0.06j, base + 0.06j, base]
    )


def vec(x):
    return np.array([complex(*v) for v in x])


def identify(w, refs):
    distances = [np.max(abs(w - g)) for g in refs]
    return int(np.argmin(distances)), float(min(distances))


def go(doc, idx, center, rad, maxstep, save, angle, additional):
    refs = [vec(t["end"]) for t in doc["branches"]]
    if additional is not None:
        records = json.loads(Path(additional).read_text())
        for rec in records["reached"] if isinstance(records, dict) else records:
            refs.append(vec(rec["z"]))
    path = route(center, rad, angle)
    w, stats = advance(refs[idx], path, doc["s"], maxstep=maxstep)
    closest, dist = identify(w, refs)
    print(
        "loop",
        idx,
        center,
        rad,
        "alpha",
        w[-1],
        "match",
        closest,
        "distance",
        dist,
        "steps",
        stats,
        flush=True,
    )
    if idx == closest and dist < 1e-6:
        return
    rev, revstats = advance(w, path[::-1], doc["s"], maxstep=maxstep * 0.65)
    back = np.max(abs(rev - refs[idx]))
    fine, finestats = advance(refs[idx], path, doc["s"], maxstep=maxstep * 0.55)
    drift = np.max(abs(fine - w))
    print("reverse roundtrip", back, "independent forward", drift, flush=True)
    obj = dict(
        idx=idx,
        center=[center.real, center.imag],
        radius=rad,
        angle=angle,
        maxstep=maxstep,
        match=closest,
        distance=dist,
        stats=stats,
        reverse=revstats,
        refined=finestats,
        reverse_error=back,
        refined_error=drift,
        endpoint=[[x.real, x.imag] for x in w],
    )
    Path(save).write_text(json.dumps(obj, indent=2) + "\n")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("input")
    p.add_argument("idx", type=int)
    p.add_argument("center", type=complex)
    p.add_argument("radius", type=float)
    p.add_argument("--step", type=float, default=0.012)
    p.add_argument("--output", default="new-sheet.json")
    p.add_argument("--angle", type=float, default=0.22)
    p.add_argument(
        "--additional",
        type=Path,
        help="append saved cycle iterates as possible starting sheets",
    )
    a = p.parse_args()
    go(
        json.load(open(a.input)),
        a.idx,
        a.center,
        a.radius,
        a.step,
        a.output,
        a.angle,
        a.additional,
    )
