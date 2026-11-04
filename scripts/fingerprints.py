"""Exact finite fingerprints for the documented formula-level novelty audit.

These matrices, their spectra, and the affine lower-bound sequence are not
claimed to be new objects. They identify the comparison used in the proof.
"""
from fractions import Fraction as F
from itertools import permutations
from hashlib import sha256
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
    with open(os.path.join(out, "fingerprints.json"), "w") as handle:
        handle.write(text)
    print(text, end="")


def identity(n):
    return [[F(int(i == j)) for j in range(n)] for i in range(n)]


def mm(a, b):
    return [[sum(x * y for x, y in zip(row, col)) for col in zip(*b)] for row in a]


def characteristic(a):
    """Faddeev--LeVerrier, with exact rational arithmetic."""
    n = len(a)
    b = identity(n)
    coefficients = [F(1)]
    for k in range(1, n + 1):
        ab = mm(a, b)
        coefficient = -sum(ab[i][i] for i in range(n)) / k
        coefficients.append(coefficient)
        b = [[ab[i][j] + (coefficient if i == j else 0)
              for j in range(n)] for i in range(n)]
    check(all(x == 0 for row in b for x in row), 'fingerprints:31')
    return coefficients


def rank(a):
    b = [row[:] for row in a]
    row = 0
    for col in range(len(b[0])):
        candidates = [i for i in range(row, len(b)) if b[i][col]]
        if not candidates:
            continue
        pivot = candidates[0]
        b[row], b[pivot] = b[pivot], b[row]
        denominator = b[row][col]
        b[row] = [x / denominator for x in b[row]]
        for i in range(len(b)):
            if i != row:
                factor = b[i][col]
                b[i] = [x - factor * y for x, y in zip(b[i], b[row])]
        row += 1
    return row


def encode(a):
    return json.dumps([[str(x) for x in row] for row in a], separators=(",", ":"))


def fingerprint(n, edges):
    a = identity(n)
    for i, j in edges:
        t = identity(n)
        t[i][j] = t[j][i] = 1
        a = mm(t, a)
    coefficients = characteristic(a)
    canonical = min(encode([[a[i][j] for j in order] for i in order])
                    for order in permutations(range(n)))
    return {
        "n": n,
        "ordered_edges_zero_based": edges,
        "comparison_matrix": [[int(x) for x in row] for row in a],
        "rank": rank(a),
        "trace": str(sum(a[i][i] for i in range(n))),
        "determinant": str((-1) ** n * coefficients[-1]),
        "characteristic_coefficients_descending": [str(x) for x in coefficients],
        "row_sums": [str(sum(row)) for row in a],
        "column_sums": [str(sum(row)) for row in zip(*a)],
        "sorted_entries": sorted(int(x) for row in a for x in row),
        "matrix_sha256": sha256(encode(a).encode()).hexdigest(),
        "transpose_sha256": sha256(encode(list(map(list, zip(*a)))).encode()).hexdigest(),
        "canonical_simultaneous_relabelling_sha256": sha256(canonical.encode()).hexdigest(),
    }


def main():
    outdir = output_dir()
    sequence = [F(5, 4) * F(9, 5) ** k - F(1, 4) for k in range(12)]
    scaled = [int(sequence[k] * 5 ** k) for k in range(len(sequence))]
    out = {
        "purpose": "Object identity and search fingerprints; not novelty claims for these standard ingredients.",
        "instances": [
            fingerprint(2, [(0, 1)]),
            fingerprint(3, [(0, 1), (1, 2), (2, 0)]),
            fingerprint(5, [(0, 1), (1, 2), (2, 3), (3, 4), (4, 0)]),
        ],
        "lower_bound_sequence": [str(x) for x in sequence],
        "lower_bound_recurrence": "b_(k+1)=(9*b_k+1)/5; b_0=1",
        "lower_bound_ordinary_generating_function": "(5-4*z)/(5-14*z+9*z^2)",
        "scaled_integer_sequence": scaled,
        "scaled_sequence_ordinary_generating_function": "(1-4*z)/(1-14*z+45*z^2)",
        "triangle_perron_root": "2+sqrt(5)",
        "triangle_perron_vector": ["(1+sqrt(5))/2", "1", "(1+sqrt(5))/2"],
    }
    check(out["instances"][1]["characteristic_coefficients_descending"] == ["1", "-4", "-1", "0"], 'fingerprints:102')
    emit(out, outdir)


if __name__ == "__main__":
    try:
        main()
    except SystemExit:
        raise
    except Exception as exc:  # any unexpected failure becomes a JSON failure record
        print(json.dumps({"status": "FAILED", "error": type(exc).__name__, "detail": str(exc)}))
        sys.exit(1)
