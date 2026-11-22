"""Read the normalized algebra via the packet's shared verify_rur.py module."""

from pathlib import Path
import importlib.util

import numpy as np


def nodal_embedding(basis, verifier):
    from flint import arb, arb_poly, ctx

    spec = importlib.util.spec_from_file_location("phase5_verify_rur", Path(verifier))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    polynomial, v, w, _ = module.read(basis)
    with ctx.workprec(1024):
        embeddings = [
            x.real
            for x, _ in polynomial.complex_roots()
            if x.imag.is_zero() and 1 < x.real < 2
        ]
        (m,) = embeddings
        u = m.root(3)
        C = [u, u * arb_poly(v)(m), u * arb_poly(w)(m)]
        g = [1, arb(5) / 7, 0, 0, 0, arb(15) / 7, 1]
        remainder = arb_poly(C) ** 3 - arb_poly(g)
        B = [(m - 1).sqrt()]
        for k in range(1, 4):
            B.append(
                (remainder[k] - sum(B[j] * B[k - j] for j in range(1, k))) / (2 * B[0])
            )
        check = remainder - arb_poly(B) ** 2
        assert all(check[i].contains(0) for i in range(7))
        error = max(float(check[i].abs_upper()) for i in range(7))
        assert error < 1e-60
        return (
            float(m),
            np.array(list(map(float, C))),
            np.array(list(map(float, B))),
            error,
        )
