#!/usr/bin/env python3
"""X4 (NUMERICAL, double precision, not a proof): divergence thresholds in a2 for the triangle IsoKF
(A=aI, Q=R=I, P0=c I), standard rule (3) and the deletion variant (external cross-covariances set
to 0 after each update). Cone variables (p_i, deficits e_ij) with the exact cycle identity enforced.
Classification after at most --N epochs: 'div' if min_i p_i > 1e12, 'bdd' if the state has converged
(max relative change < 1e-13); otherwise 'div(slow)' if min p grew by >1% over the second half of the
run, else 'bdd?'. Brackets near a threshold are therefore resolution-limited. Bisection to --tol.
"""
import argparse, json, os, sys
PAIRS = [(0, 1), (1, 2), (2, 0)]

def epoch(p, e, a2, delete):
    p = [a2 * x + 1 for x in p]; e = {k: a2 * x + 1 for k, x in e.items()}
    for (i, j) in PAIRS:
        l = 3 - i - j
        u, v = e[(i, j)], e[(j, i)]; s = u + v + 1
        li, lj = u * u / s, v * v / s
        pi2, pj2 = p[i] - li, p[j] - lj
        if delete:
            e_li, e_lj, e_il, e_jl = p[l], p[l], pi2, pj2
        else:
            e_li = e[(l, i)] + (1 - e[(i, l)] / p[i]) * li; e_lj = e[(l, j)] + (1 - e[(j, l)] / p[j]) * lj
            e_il = e[(i, l)] * (pi2 / p[i]); e_jl = e[(j, l)] * (pj2 / p[j])
        e[(i, j)] = u / s; e[(j, i)] = v / s
        e[(l, i)], e[(l, j)], e[(i, l)], e[(j, l)] = e_li, e_lj, e_il, e_jl
        p[i], p[j] = pi2, pj2
        if not delete:
            e[(0, 2)] = e[(0, 1)] - e[(1, 0)] + e[(1, 2)] - e[(2, 1)] + e[(2, 0)]
    return p, e

def classify(a2, c, delete, N):
    p = [float(c)] * 3; e = {(i, j): float(c) for i in range(3) for j in range(3) if i != j}
    half = None
    for k in range(N):
        p0 = list(p); p, e = epoch(p, e, a2, delete)
        if min(p) > 1e12: return "div", k
        if max(abs(x - y) / y for x, y in zip(p, p0)) < 1e-13: return "bdd", k
        if k == N // 2: half = min(p)
    # undecided within N epochs: slow growth counts as divergence, otherwise 'bdd?'
    return ("div(slow)" if min(p) > 1.01 * half else "bdd?"), N

def bisect(c, delete, lo, hi, N, tol):
    if classify(lo, c, delete, N)[0] != "bdd" or classify(hi, c, delete, N)[0] == "bdd":
        return None
    while hi - lo > tol:
        m = (lo + hi) / 2
        if classify(m, c, delete, N)[0] in ("bdd", "bdd?"): lo = m
        else: hi = m
    return [lo, hi]

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--out", required=True)
    ap.add_argument("--N", type=int, default=20000); ap.add_argument("--tol", type=float, default=1e-4)
    a = ap.parse_args()
    if os.path.exists(a.out):
        print(json.dumps({"error": "output directory exists; refusing to overwrite"})); sys.exit(2)
    rows = []
    for c in [1, 3, 10, 100, 10**4]:
        rows.append({"rule": "IsoKF (3)", "P0": f"{c}*I", "a2_threshold_bracket": bisect(c, False, 1.0000001, 4.3, a.N, a.tol)})
    for c in [1, 10, 10**4]:
        rows.append({"rule": "deletion", "P0": f"{c}*I", "a2_threshold_bracket": bisect(c, True, 1.0000001, 12.0, a.N, a.tol)})
    spot = {f"a2={x},c={c}": classify(x, c, False, a.N)[0] for x in [1.0, 1.01, 1.05] for c in [1, 10**4, 10**8]}
    os.makedirs(a.out)
    out = {"script": "x4_thresholds.py", "label": "numerical (double precision; not a proof)",
           "N": a.N, "tol": a.tol, "thresholds": rows, "spot_checks": spot}
    with open(os.path.join(a.out, "x4_thresholds.json"), "w") as f: json.dump(out, f, indent=1)
    print(json.dumps(out))
    if any(r["a2_threshold_bracket"] is None for r in rows): sys.exit(1)

if __name__ == "__main__":
    try:
        main()
    except SystemExit:
        raise
    except Exception as exc:  # any unexpected failure becomes a JSON failure record
        print(json.dumps({"status": "FAILED", "error": type(exc).__name__, "detail": str(exc)}))
        sys.exit(1)
