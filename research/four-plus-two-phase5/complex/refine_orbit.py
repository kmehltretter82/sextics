"""Refine and certify a numerical orbit with the packet's shared Acb root checker."""

import argparse
import importlib.util
import gzip
import json
from pathlib import Path
from flint import acb, arb, ctx


def decimal_pair(z):
    return [z.real.mid().str(60, radius=False), z.imag.mid().str(60, radius=False)]


def disjoint_representatives(boxes):
    """Pick pairwise disjoint boxes. Overlapping projections could conceal more values."""
    result = []
    for i, box in enumerate(boxes):
        if all(any(not a.overlaps(b) for a, b in zip(box, boxes[j])) for j in result):
            result.append(i)
    return result


def symmetry(box, conjugate, negate):
    z = [v.conjugate() for v in box] if conjugate else list(box)
    if negate:
        # A is even and B is odd in (c,d). Negating the eight B coefficients and
        # alpha fixes both norm constraints and point constraints at fixed (r,s).
        z[18:26] = [-v for v in z[18:26]]
        z[26] = -z[26]
    return z


def certify(queue, checker, results, summary):
    spec = importlib.util.spec_from_file_location("phase5_root_checker", checker)
    cr = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(cr)
    ctx.prec = 320
    r, s = acb(1) / 10, acb(1)
    raw = queue.read_bytes()
    doc = json.loads(gzip.decompress(raw) if queue.suffix == ".gz" else raw)
    input_points = doc["points"]
    success, problems, boxes = [], [], []
    for index, point in enumerate(input_points):
        try:
            z = [acb(arb(str(a)), arb(str(b))).mid() for a, b in point]
            for _ in range(5):
                f, j = cr.system(z, r, s)
                step = j.solve(f)
                z = [(v - step[i, 0]).mid() for i, v in enumerate(z)]
            coordinates = z[:26] + z[-1:]
            record = dict(
                orbit_index=index,
                coordinates=[decimal_pair(v) for v in coordinates],
                alpha=decimal_pair(z[-1]),
            )
            # Match the public checker's center, including its exact quartic pivot recovery.
            box, eps, L = cr.root_box(record, r, s)
            p = box[5] * box[6]
            delta = (1 - p) * (1 - 3 * p) * (1 - 2 * p) ** 3
            assert delta.abs_lower() > 0
            record.update(
                newton_displacement_upper=float(eps),
                lipschitz_upper=float(L),
                localization_lower=float(delta.abs_lower()),
            )
            boxes.append(box[:26] + box[-1:])
            success.append(record)
        except (AssertionError, ArithmeticError, ValueError, ZeroDivisionError) as e:
            problems.append({"orbit_index": index, "error": repr(e)})
    transforms = []
    counts = {}
    for name, conjugate, negate in [
        ("identity", False, False),
        ("conjugates", True, False),
        ("sign", False, True),
        ("conjugate_sign", True, True),
    ]:
        transforms.extend(symmetry(b, conjugate, negate) for b in boxes)
        counts[name] = dict(
            cumulative_distinct_lifts_lower_bound=len(
                disjoint_representatives(transforms)
            ),
            cumulative_distinct_alpha_lower_bound=len(
                disjoint_representatives([[z[-1]] for z in transforms])
            ),
        )
    certificate = dict(
        numerical_points=len(input_points),
        certified=len(success),
        failed=problems,
        max_newton_displacement=max(v["newton_displacement_upper"] for v in success),
        max_lipschitz=max(v["lipschitz_upper"] for v in success),
        min_localization_lower=min(v["localization_lower"] for v in success),
        verified_complex_sup_radius="1e-44",
        symmetries=counts,
    )
    # Independent presentation symmetries fix alpha.  Scaling Q1 by a cube
    # root of unity leaves Q1^3 unchanged. Negating both A and B leaves
    # A^2-B^2 unchanged and rescales the moving cubic by -1.
    lambda_lower = min(max(v.abs_lower() for v in z[:5]) for z in boxes)
    cubic_lower = min(max(v.abs_lower() for v in z[10:26]) for z in boxes)
    assert lambda_lower > 0 and cubic_lower > 0
    certificate["unquotiented_fibre"] = dict(
        minimum_max_lambda_abs_lower=float(lambda_lower),
        minimum_max_cubic_coefficient_abs_lower=float(cubic_lower),
        distinct_lifts_lower_bound=6
        * counts["conjugate_sign"]["cumulative_distinct_alpha_lower_bound"],
        reason="independent Q1 times a cube root of unity and simultaneous A,B negation; each fixes alpha",
    )
    print(json.dumps(certificate, indent=2), flush=True)
    encoded = (json.dumps(success, separators=(",", ":")) + "\n").encode()
    results.write_bytes(
        gzip.compress(encoded, mtime=0) if results.suffix == ".gz" else encoded
    )
    summary.write_text(json.dumps(certificate, indent=2) + "\n")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("queue", type=Path)
    p.add_argument("endpoint_output", type=Path)
    p.add_argument("summary_output", type=Path)
    p.add_argument(
        "--checker", type=Path, default=Path(__file__).with_name("certify_roots.py")
    )
    a = p.parse_args()
    certify(a.queue, a.checker, a.endpoint_output, a.summary_output)
