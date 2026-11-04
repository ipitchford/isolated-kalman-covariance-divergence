#!/bin/sh
# Re-run the interval certificates and numerical extension scripts.
# Usage: sh scripts/run_ext.sh RESULTS_DIR FLAG   (FLAG is "" or "-O"); run from the candidate root.
set -u
R="$1"; FL="${2:-}"
PY="uv run --with mpmath==1.3.0 python3 -I -B $FL"
mkdir -p "$R"
$PY scripts/cert_x1_tail.py --out "$R/x1_tail" --a2max 69/50 --step 1/500 --rad 0.006 --theta0 1e-6 > "$R/x1_tail.stdout.json"; echo "x1_tail exit=$?"
$PY scripts/cert_x1_entry.py --out "$R/x1_entry" --tail "$R/x1_tail/x1_tail.json" --a2 100000001/100000000 1000001/1000000 1001/1000 101/100 21/20 11/10 6/5 69/50 --c "10**12" --K 120 > "$R/x1_entry.stdout.json"; echo "x1_entry exit=$?"
$PY scripts/cert_x2_bounded.py --out "$R/x2" > "$R/x2.stdout.json"; echo "x2 exit=$?"
$PY scripts/x3_true_error.py --out "$R/x3" > "$R/x3.stdout.json"; echo "x3 exit=$?"
$PY scripts/x4_thresholds.py --out "$R/x4" > "$R/x4.stdout.json"; echo "x4 exit=$?"
