from itertools import product
import numpy as np
import sys


def indices(degree):
    return [
        (a, b, c, degree - a - b - c)
        for a in range(degree + 1)
        for b in range(degree - a + 1)
        for c in range(degree - a - b + 1)
    ]


def smooth(p):
    k = 5 * pow(7, -1, p) % p
    H = {
        (6, 0, 0, 0): 1,
        (0, 6, 0, 0): 1,
        (0, 0, 6, 0): 1,
        (0, 0, 0, 6): 1,
        (5, 1, 0, 0): k,
        (1, 5, 0, 0): 2 * k,
        (0, 0, 5, 1): k,
        (0, 0, 1, 5): 3 * k,
        (3, 1, 1, 1): k,
        (1, 1, 3, 1): 2 * k,
    }
    small = indices(12)
    large = {v: i for i, v in enumerate(indices(17))}
    a = np.zeros((4 * len(small), len(large)), dtype=np.int64)
    for j in range(4):
        for i, m in enumerate(small):
            for e, c in H.items():
                if e[j]:
                    key = tuple(m[n] + e[n] - int(n == j) for n in range(4))
                    a[len(small) * j + i, large[key]] = e[j] * c % p
    row = 0
    for col in range(a.shape[1]):
        avail = np.flatnonzero(a[row:, col])
        if not len(avail):
            continue
        pick = row + avail[0]
        a[[row, pick]] = a[[pick, row]]
        hot = np.flatnonzero(a[row + 1 :, col]) + row + 1
        pivot = a[row, col:] * pow(int(a[row, col]), -1, p) % p
        a[hot, col:] = (a[hot, col:] - a[hot, col, None] * pivot) % p
        row += 1
    print(
        p,
        "Jacobian matrix",
        a.shape,
        "rank",
        row,
        "smooth",
        row == a.shape[1],
        flush=True,
    )


for x in sys.argv[1:]:
    smooth(int(x))
