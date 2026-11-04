"""Exact checks of the three-state pairwise-only IsoKF counterexample.

This is an independent recurrence calculation. The infinite-time conclusion
comes from the accompanying manuscript, not from the finite calculations here.
Only Python's standard library is required.
"""
from fractions import Fraction as F
from itertools import combinations
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
    with open(os.path.join(out, "verify_pairwise_divergence.json"), "w") as handle:
        handle.write(text)
    print(text, end="")


def pair_update(p, i, j, noise=F(1)):
    n = len(p)
    a, b, c = p[i][i], p[j][j], -p[i][j]
    check(a > 0 and b > 0 and 0 <= c < min(a, b), 'verify_pairwise_divergence:15')
    check(a * b - c * c > 0, 'verify_pairwise_divergence:16')
    u, v = a - c, b - c
    innovation = u + v + noise
    aa = a - u * u / innovation
    bb = b - v * v / innovation
    cc = c + u * v / innovation
    # Identities used in the proof, independently checked at each exact state.
    check(aa - cc == u * noise / innovation, 'verify_pairwise_divergence:23')
    check(bb - cc == v * noise / innovation, 'verify_pairwise_divergence:24')
    check(aa - a * b / (a + b) == (
        c * (2 * a * b - c * (a + b)) + a * a * noise
    ) / ((a + b) * innovation), 'verify_pairwise_divergence:25')
    check(aa >= a * b / (a + b) and bb >= a * b / (a + b), 'verify_pairwise_divergence:28')
    out = [row[:] for row in p]
    for k in range(n):
        if k not in (i, j):
            out[i][k] = out[k][i] = aa / a * p[i][k]
            out[j][k] = out[k][j] = bb / b * p[j][k]
    out[i][i], out[j][j] = aa, bb
    out[i][j] = out[j][i] = -cc
    check_cone(out)
    return out, innovation


def check_cone(p):
    for i in range(len(p)):
        check(p[i][i] > 0, 'verify_pairwise_divergence:42')
    for i, j in combinations(range(len(p)), 2):
        c = -p[i][j]
        check(0 <= c < min(p[i][i], p[j][j]), 'verify_pairwise_divergence:45')
        check(p[i][i] * p[j][j] - c * c > 0, 'verify_pairwise_divergence:46')


def det3(p):
    return (p[0][0] * (p[1][1] * p[2][2] - p[1][2] * p[2][1])
            - p[0][1] * (p[1][0] * p[2][2] - p[1][2] * p[2][0])
            + p[0][2] * (p[1][0] * p[2][1] - p[1][1] * p[2][0]))


def main():
    out = output_dir()
    p = [[F(int(i == j)) for j in range(3)] for i in range(3)]
    check(det3([[1, 1, 0], [0, 1, 1], [1, 0, 1]]) == 2, 'verify_pairwise_divergence:57')
    lower = F(1)
    snapshots = []
    for epoch in range(1, 5):
        p = [[9 * p[i][j] + F(int(i == j)) for j in range(3)]
             for i in range(3)]
        check_cone(p)
        innovations = []
        for i, j in [(0, 1), (1, 2), (2, 0)]:
            p, s = pair_update(p, i, j)
            check(s > 0, 'verify_pairwise_divergence:67')
            innovations.append(s)
        lower = (9 * lower + 1) / 5
        check(lower == F(5, 4) * F(9, 5) ** epoch - F(1, 4), 'verify_pairwise_divergence:70')
        check(all(p[i][i] >= lower for i in range(3)), 'verify_pairwise_divergence:71')
        snapshots.append({
            "epoch": epoch,
            "minimum_marginal_approx": float(min(p[i][i] for i in range(3))),
            "rigorous_lower_bound": str(lower),
            "global_determinant_sign": 1 if det3(p) > 0 else -1,
            "exact_numerator_max_bits": max(x.numerator.bit_length() for row in p for x in row),
        })
        if epoch == 1:
            first = {
                "after_first_epoch": [[str(x) for x in row] for row in p],
                "innovations": [str(s) for s in innovations],
            }

    # Fingerprint of the reciprocal-variance schedule comparison.
    matrix = [[2, 2, 1], [1, 1, 1], [2, 2, 1]]
    check(det3(matrix) == 0, 'verify_pairwise_divergence:87')
    check(sum(matrix[i][i] for i in range(3)) == 4, 'verify_pairwise_divergence:88')
    principal_minors_sum = sum(
        matrix[i][i] * matrix[j][j] - matrix[i][j] * matrix[j][i]
        for i, j in combinations(range(3), 2)
    )
    check(principal_minors_sum == -1, 'verify_pairwise_divergence:93')
    emit({
        "status": "all_exact_checks_passed",
        "observability_determinant": 2,
        "initial_covariance": "I_3",
        "state_transition": "3 I_3",
        "process_covariance": "I_3",
        "measurement_variances": [1, 1, 1],
        "comparison_matrix": matrix,
        "comparison_characteristic_polynomial": "lambda*(lambda^2-4*lambda-1)",
        **first,
        "finite_snapshots": snapshots,
        "proof_lower_bound": "min_i p_i(k) >= (5/4)*(9/5)^k - 1/4",
    }, out)


if __name__ == "__main__":
    try:
        main()
    except SystemExit:
        raise
    except Exception as exc:  # any unexpected failure becomes a JSON failure record
        print(json.dumps({"status": "FAILED", "error": type(exc).__name__, "detail": str(exc)}))
        sys.exit(1)
