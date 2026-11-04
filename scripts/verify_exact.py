"""Exact rational checks of the IsoKF well-posedness obstructions.

This independent implementation follows the final update equations (47)--(50)
of Jung, Luft and Weiss (2025), with scalar node states. It does not reproduce
software regularisation or rejection guards. Python standard library only.
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
    with open(os.path.join(out, "verify_exact.json"), "w") as handle:
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


def mv(a, x):
    return [sum(ai * xi for ai, xi in zip(row, x)) for row in a]


def dot(x, y):
    return sum(xi * yi for xi, yi in zip(x, y))


def principal(a, ids):
    return [[a[i][j] for j in ids] for i in ids]


def is_spd(a):
    return all(det(principal(a, list(range(k)))) > 0 for k in range(1, len(a) + 1))


def update(p, ids, h, r):
    """One scalar observation on selected participants, exact fractions."""
    ids = list(ids)
    h = [F(x) for x in h]
    r = F(r)
    b = principal(p, ids)
    v = mv(b, h)
    innovation = dot(h, v) + r
    if innovation == 0:
        raise ZeroDivisionError("The exact innovation variance is zero.")
    c = [[b[i][j] - v[i] * v[j] / innovation for j in range(len(ids))]
         for i in range(len(ids))]
    out = [row[:] for row in p]
    others = [i for i in range(len(p)) if i not in ids]
    for loc, i in enumerate(ids):
        factor = c[loc][loc] / b[loc][loc]
        for j in others:
            out[i][j] = factor * p[i][j]
            out[j][i] = out[i][j]
    for i, ii in enumerate(ids):
        for j, jj in enumerate(ids):
            out[ii][jj] = c[i][j]
    return out, innovation


def string_matrix(p):
    return [[str(x) for x in row] for row in p]


def main():
    out = output_dir()
    # A=2I, Q=I, P0=I, so the first predicted covariance is 5I.
    p = mat([[5, 0, 0], [0, 5, 0], [0, 0, 5]])
    p1, s1 = update(p, [0, 1], [1, -1], 1)
    check(is_spd(p1), 'verify_exact:75')
    p2, s2 = update(p1, [1, 2], [1, -1], 1)
    expected = mat([[F(30, 11), F(25, 16), 0],
                    [F(25, 16), F(15, 8), F(25, 16)],
                    [0, F(25, 16), F(205, 96)]])
    check(p2 == expected, 'verify_exact:80')
    for ids in combinations(range(3), 2):
        check(is_spd(principal(p2, ids)), 'verify_exact:82')
    check(det(p2) < 0, 'verify_exact:83')
    h = [F(2), F(-3), F(2)]
    v = mv(p2, h)
    check(v == [F(135, 176), F(5, 8), F(-5, 12)], 'verify_exact:86')
    check(dot(h, v) == -F(155, 132), 'verify_exact:87')
    check(det(mat([[1, -1, 0], [0, 1, -1], [2, -3, 2]])) == 1, 'verify_exact:88')
    singular_noise = F(155, 132)
    try:
        update(p2, [0, 1, 2], h, singular_noise)
    except ZeroDivisionError:
        pass
    else:
        raise AssertionError("Expected an exact singular innovation.")
    p3, s3 = update(p2, [0, 1, 2], h, F(5, 4))
    check(s3 == F(5, 66), 'verify_exact:97')
    check([p3[i][i] for i in range(3)] == [
        -F(7095, 1408), -F(105, 32), -F(5, 32)], 'verify_exact:98')
    thresholds = [v[i] ** 2 / p2[i][i] for i in range(3)]
    check(min(thresholds) == F(10, 123), 'verify_exact:101')

    # Sharp participant-family obstruction.
    g0 = mat([[1, 0, F(1, 2)], [0, 1, -F(1, 2)], [F(1, 2), -F(1, 2), 1]])
    check(is_spd(g0) and det(g0) == F(1, 2), 'verify_exact:105')
    g1, gs = update(g0, [0, 1], [1, -1], F(1, 4))
    expected_g1 = mat([[F(5, 9), F(4, 9), F(5, 18)],
                       [F(4, 9), F(5, 9), -F(5, 18)],
                       [F(5, 18), -F(5, 18), 1]])
    check(g1 == expected_g1, 'verify_exact:110')
    check(det(g1) == -F(7, 162), 'verify_exact:111')
    gh = [F(2), F(-2), F(-1)]
    check(dot(gh, mv(g1, gh)) == -F(1, 3), 'verify_exact:113')
    check(mv(g1, gh) == [-F(1, 18), F(1, 18), F(1, 9)], 'verify_exact:114')

    emit({
        "status": "all_exact_checks_passed",
        "first_innovations": [str(s1), str(s2)],
        "before_third": string_matrix(p2),
        "determinant_before_third": str(det(p2)),
        "third_measurement": [str(x) for x in h],
        "negative_quadratic_form": str(dot(h, v)),
        "singular_positive_noise": str(singular_noise),
        "robust_noise": "5/4",
        "robust_innovation": str(s3),
        "robust_negative_marginals": [str(p3[i][i]) for i in range(3)],
        "open_all_negative_noise_interval": [str(singular_noise), str(singular_noise + min(thresholds))],
        "participant_overlap_witness_after": string_matrix(g1),
        "participant_overlap_determinant": str(det(g1))
    }, out)


if __name__ == "__main__":
    try:
        main()
    except SystemExit:
        raise
    except Exception as exc:  # any unexpected failure becomes a JSON failure record
        print(json.dumps({"status": "FAILED", "error": type(exc).__name__, "detail": str(exc)}))
        sys.exit(1)
