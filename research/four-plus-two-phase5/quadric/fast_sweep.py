"""Dyadic online convolution for the quadratic marked-contact arc at p=11.
Usage: python fast_sweep.py UPSTREAM_CONTACT_SPAN_DIR --helpers PREVIOUS_SMOOTH_DIR -n 4096 --out sweep.npz
Requires numpy, python-flint, the existing smooth_span.py/smooth_sweep.py and their local imports.
"""

import argparse
import sys
from pathlib import Path
from time import monotonic
from functools import lru_cache
import numpy as np


def extend(seed, n, block):
    st = seed.upstream
    start, tangent = seed.seed()
    J = st.differential(start)
    _, cols = st.base_rref(J)
    inv = np.array(st.base_inverse(J[:, cols]), dtype=np.float64)
    signal = np.zeros((n, 90))
    signal[0] = start
    signal[1] = tangent
    accum = np.zeros((n, 70), dtype=np.int64)
    groups = []
    for row in range(70):
        sel = st.term_row == row
        groups.append((st.term_i[sel], st.term_j[sel], st.term_sign[sel]))

    @lru_cache(maxsize=None)
    def first_fft(m):
        return np.fft.rfft(signal[m : 2 * m].T, n=2 * m, axis=1)

    previous = []
    signal_float = signal
    for k in range(2, n):
        m = block
        while 2 * m <= k and k % m == 0:
            first = first_fft(m)
            same = 2 * m == k
            second = (
                first if same else np.fft.rfft(signal[k - m : k].T, n=2 * m, axis=1)
            )
            freq = np.zeros((70, m + 1), dtype=np.complex128)
            # Fourier convolutions across the time axis, grouped by physical monomial.
            for row, (ii, jj, ss) in enumerate(groups):
                z = first[ii] * second[jj]
                if not same:
                    z += second[ii] * first[jj]
                freq[row] = ss @ z
            values = np.fft.irfft(freq, n=2 * m, axis=1)
            integer = np.rint(values)
            assert np.max(np.abs(values - integer)) < 1e-4
            length = min(2 * m - 1, n - k)
            accum[k : k + length] += integer[:, :length].T.astype(np.int64)
            m *= 2
        if k < 2 * block:
            low = min(k - 1, block - 1)
            gram = signal_float[1 : low + 1].T @ signal_float[k - 1 : k - low - 1 : -1]
            other = min(k - block, block - 1)
            if other > 0:
                gram += (
                    signal_float[k - 1 : k - other - 1 : -1].T
                    @ signal_float[1 : other + 1]
                )
            recent = np.bincount(
                st.term_row,
                weights=st.term_sign * gram[st.term_i, st.term_j],
                minlength=70,
            )
        else:
            if isinstance(previous, list):
                previous = np.asarray(
                    [st.differential(signal[i].astype(int)) for i in range(1, block)],
                    dtype=np.float64,
                )
            recent = np.einsum(
                "ijk,ik->j", previous, signal_float[k - 1 : k - block : -1]
            )
        rhs = (recent + accum[k]) % 11
        signal[k, cols] = (-inv @ rhs) % 11
    out = signal.astype(np.int16)
    return out


def verify_contact(seed, arr, poly):
    st = seed.upstream
    n = len(arr)
    params = [poly(col.tolist(), 11) for col in arr.T]
    rows = [poly([], 11) for _ in range(70)]
    for i, j, row, sign in zip(st.term_i, st.term_j, st.term_row, st.term_sign):
        rows[row] += int(sign) * params[i].mul_low(params[j], n)
    c0 = seed.quadratic(arr[0], arr[0])
    assert all(f == int(c) for f, c in zip(rows, c0))
    return params


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("upstream")
    parser.add_argument("--helpers", required=True)
    parser.add_argument("-n", type=int, default=4096)
    parser.add_argument("--block", type=int, default=16)
    parser.add_argument("--out", type=Path, required=True)
    opts = parser.parse_args()
    if opts.n < 2 or opts.block < 2 or opts.block & (opts.block - 1):
        parser.error("n must be at least two and block must be a power of two >=2")
    sys.path.insert(0, opts.helpers)
    sys.argv[1] = opts.upstream
    import smooth_span as seed
    from smooth_sweep import restricted, points
    from flint import nmod_poly

    start = monotonic()
    arr = extend(seed, opts.n, opts.block)
    print(
        "generated",
        opts.n,
        "orders of 90 coordinates in",
        round(monotonic() - start, 2),
        "s",
        flush=True,
    )
    t = monotonic()
    params = verify_contact(seed, arr, nmod_poly)
    q, h = restricted(params, 2, opts.n)
    ps = list(points(q, h, opts.n))
    print(
        "all 70 contact equations and six plane point equations checked modulo x^"
        + str(opts.n),
        "in",
        round(monotonic() - t, 2),
        "s",
        flush=True,
    )
    outputs = {}
    for i, (a, b) in enumerate(ps):
        t = a.inverse_series_trunc(opts.n)
        r = t.mul_low(b, opts.n)
        z = t.mul_low(t, opts.n)
        for key, series in [("r", r), ("Z", z), ("a", a), ("b", b)]:
            buf = np.zeros(opts.n, dtype=np.uint8)
            buf[: len(series)] = list(series)
            outputs[f"{key}{i}"] = buf
        print("plane sheet", i, "first r,Z jets", list(r)[:3], list(z)[:3], flush=True)
    np.savez_compressed(opts.out, **outputs)


if __name__ == "__main__":
    main()
