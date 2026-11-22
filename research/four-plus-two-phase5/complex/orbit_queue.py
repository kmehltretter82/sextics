"""Bounded, reproducible three-worker numerical exploration of the specialized incidence."""

import argparse
import concurrent.futures as cf
import json
import os
from pathlib import Path
import time
import numpy as np
from loops import advance, route, identify, vec

# Each path starts and ends at r=0.1 and uses a named based circle.
GENERATORS = [
    (0.15 + 0.28j, 0.38, 0.22),
    (0.4 + 0j, 0.38, 0.22),
    (0.25 + 0.20j, 0.32, 0.22),
    (0.15 + 0.28j, 0.38, 1.31),
    (-0.16 + 0.24j, 0.32, 0.8),
    (0.10 - 0.24j, 0.3, 0.3),
    (-0.28 - 0.10j, 0.34, 1.1),
    (0.23 - 0.31j, 0.26, 0.42),
]


def run_task(j):
    seq, idx, z, g, start, deadline = j
    path = route(*GENERATORS[g])
    # Reserve time for shutdown and partial-result collection.
    if time.monotonic() > deadline - 2:
        return seq, idx, g, None, {}, "time limit"
    try:
        end, info = advance(z, path, 1.0, 0.012, deadline=deadline)
        info["elapsed"] = time.monotonic() - start
        if seq % 11 == 0:
            back, binfo = advance(end, path[::-1], 1.0, 0.008, deadline=deadline)
            info["reverse_error"] = float(np.max(abs(back - z)))
            if info["reverse_error"] > 1e-6:
                raise ArithmeticError("failed reverse")
        return seq, idx, g, end, info, None
    except (ArithmeticError, np.linalg.LinAlgError) as e:
        return seq, idx, g, None, {}, str(e)


def run(initial, extras, out, maxtasks, duration, workers):
    start = time.monotonic()
    deadline = start + duration
    pts = [vec(rec["end"]) for rec in json.loads(Path(initial).read_text())["branches"]]
    for n in extras:
        rec = json.loads(Path(n).read_text())
        vals = (
            [rec["endpoint"]]
            if "endpoint" in rec
            else [x["z"] for x in rec.get("reached", [])]
            if isinstance(rec, dict)
            else [x["z"] for x in rec]
        )
        for zraw in vals:
            z = vec(zraw)
            if identify(z, pts)[1] > 1e-6:
                pts.append(z)
    n0 = len(pts)
    logs = []
    done = set()
    tasks = 0
    next_gen = [(i, g) for g in range(len(GENERATORS)) for i in range(n0)]
    with cf.ProcessPoolExecutor(workers) as pool:
        pending = {}
        while (tasks < maxtasks or pending) and time.monotonic() < deadline:
            while len(pending) < workers and tasks < maxtasks and next_gen:
                idx, g = next_gen.pop(0)
                if (idx, g) in done:
                    continue
                done.add((idx, g))
                seq = tasks
                tasks += 1
                f = pool.submit(run_task, (seq, idx, pts[idx], g, start, deadline))
                pending[f] = seq
            if not pending:
                break
            finished, _ = cf.wait(pending, timeout=0.5, return_when=cf.FIRST_COMPLETED)
            for f in finished:
                pending.pop(f)
                seq, idx, g, z, info, err = f.result()
                row = dict(
                    job=seq, source=idx, generator=g, diagnostics=info, error=err
                )
                if z is not None:
                    nearest, dist = identify(z, pts)
                    row["nearest_distance"] = dist
                    if dist > 1e-6:
                        nearest = len(pts)
                        pts.append(z)
                        # Front-load different generators on new sheets without losing the remaining BFS.
                        next_gen.append((nearest, (g + 1) % len(GENERATORS)))
                        next_gen.extend(
                            (nearest, h)
                            for h in range(len(GENERATORS))
                            if h != (g + 1) % len(GENERATORS)
                        )
                    row["target"] = nearest
                logs.append(row)
                if len(logs) % 15 == 0:
                    print(
                        len(logs),
                        "completed;",
                        len(pts),
                        "distinct numerical points",
                        flush=True,
                    )
            if time.monotonic() > deadline - 1:
                break
    out.write_text(
        json.dumps(
            dict(
                initial=n0,
                elapsed=time.monotonic() - start,
                submitted=tasks,
                generators=[[c.real, c.imag, r, t] for c, r, t in GENERATORS],
                paths=sorted(logs, key=lambda x: x["job"]),
                points=[[[v.real, v.imag] for v in z] for z in pts],
            ),
            separators=(",", ":"),
        )
        + "\n"
    )
    print(
        "finished",
        len(logs),
        "jobs",
        len(pts),
        "points",
        sum(x["error"] is not None for x in logs),
        "failures",
        sum("reverse_error" in x["diagnostics"] for x in logs),
        "reverse checks",
        flush=True,
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("initial")
    p.add_argument("output", type=Path)
    p.add_argument("extra", nargs="*")
    p.add_argument("--tasks", type=int, default=240)
    p.add_argument("--seconds", type=int, default=230)
    p.add_argument("--workers", type=int, default=3)
    a = p.parse_args()
    run(a.initial, a.extra, a.output, a.tasks, a.seconds, a.workers)
