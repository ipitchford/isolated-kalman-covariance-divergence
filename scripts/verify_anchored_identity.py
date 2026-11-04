"""Exact noncommuting-block check of the anchored covariance transfer identity.

This checks factor order on a concrete positive definite rational example.
The general identity and period-by-period bound are proved in the manuscript.
Only the Python standard library is required.
"""
from fractions import Fraction as F
import json
import argparse
import hashlib
import os
import sys

CHECKS = {"passed": 0}


def check(condition, label):
    """Explicit check (works under python -O); exit with code 1 on failure."""
    if not condition:
        print(json.dumps({"status": "FAILED", "failed_check": label}))
        sys.exit(1)
    CHECKS["passed"] += 1


def output_dir():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", required=True, help="new directory (must not exist)")
    args = parser.parse_args()
    if os.path.exists(args.out):
        print(json.dumps({"error": "output directory exists; refusing to overwrite"}))
        sys.exit(2)
    return args.out


def emit(result, out):
    with open(__file__, "rb") as handle:
        result["script_sha256"] = hashlib.sha256(handle.read()).hexdigest()
    result["explicit_checks_passed"] = CHECKS["passed"]
    os.makedirs(out)
    text = json.dumps(result, indent=2) + "\n"
    with open(os.path.join(out, "verify_anchored_identity.json"), "w") as handle:
        handle.write(text)
    print(text, end="")


def mat(rows):
    return [[F(x) for x in row] for row in rows]


def det(a):
    if len(a) == 0:
        return F(1)
    if len(a) == 1:
        return a[0][0]
    return sum(((-1) ** j) * a[0][j] * det([
        row[:j] + row[j + 1:] for row in a[1:]
    ]) for j in range(len(a)))


def principal(a, ids):
    return [[a[i][j] for j in ids] for i in ids]


def is_spd(a):
    return all(det(principal(a, list(range(k)))) > 0 for k in range(1, len(a) + 1))


def string_matrix(p):
    return [[str(x) for x in row] for row in p]


def transpose(a):
    return [list(row) for row in zip(*a)]


def multiply(a, b):
    return [[sum(x * y for x, y in zip(row, col))
             for col in zip(*b)] for row in a]


def add(a, b):
    return [[x + y for x, y in zip(ar, br)] for ar, br in zip(a, b)]


def subtract(a, b):
    return [[x - y for x, y in zip(ar, br)] for ar, br in zip(a, b)]


def inverse(a):
    n = len(a)
    aug = [row[:] + [F(int(i == j)) for j in range(n)]
           for i, row in enumerate(a)]
    for k in range(n):
        pivot = next(i for i in range(k, n) if aug[i][k])
        aug[k], aug[pivot] = aug[pivot], aug[k]
        scale = aug[k][k]
        aug[k] = [x / scale for x in aug[k]]
        for i in range(n):
            if i != k:
                factor = aug[i][k]
                aug[i] = [x - factor * y for x, y in zip(aug[i], aug[k])]
    return [row[n:] for row in aug]


def main():
    out = output_dir()
    x = mat([[2, 1], [1, 3]])
    y = mat([[3, F(1, 2)], [F(1, 2), 2]])
    e = mat([[F(1, 2), F(1, 3)], [-F(1, 4), F(2, 5)]])
    r = mat([[1, F(1, 4)], [F(1, 4), F(3, 2)]])
    et = transpose(e)
    b = [x[i] + e[i] for i in range(2)] + [et[i] + y[i] for i in range(2)]
    check(is_spd(b) and is_spd(r), 'verify_anchored_identity:52')
    check(e != et, 'verify_anchored_identity:53')
    check(multiply(x, r) != multiply(r, x), 'verify_anchored_identity:54')
    s = subtract(subtract(add(add(x, y), r), e), et)
    check(is_spd(s), 'verify_anchored_identity:56')
    si = inverse(s)
    cjj = subtract(y, multiply(multiply(subtract(y, et), si), subtract(y, e)))
    left = subtract(add(x, r), cjj)
    right = multiply(multiply(subtract(add(x, r), e), si), subtract(add(x, r), et))
    check(left == right, 'verify_anchored_identity:61')
    check(is_spd(left) and is_spd(cjj), 'verify_anchored_identity:62')
    emit({
        "status": "all_exact_checks_passed",
        "prior": string_matrix(b),
        "measurement_noise": string_matrix(r),
        "prior_determinant": str(det(b)),
        "posterior_second_marginal": string_matrix(cjj),
        "gap_matrix": string_matrix(left),
        "gap_determinant": str(det(left)),
        "noncommuting_blocks": True,
        "identity": "X+R-Cjj=(X+R-E) S^(-1) (X+R-E^T)",
    }, out)


if __name__ == "__main__":
    try:
        main()
    except SystemExit:
        raise
    except Exception as exc:  # any unexpected failure becomes a JSON failure record
        print(json.dumps({"status": "FAILED", "error": type(exc).__name__, "detail": str(exc)}))
        sys.exit(1)
