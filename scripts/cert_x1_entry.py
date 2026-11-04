#!/usr/bin/env python3
"""X1 entry certificate: rigorous interval simulation (mpmath.iv, outward rounding) of the triangle
IsoKF recurrence in exact cone variables (p_i, e_ij) from priors P0 with diagonal entries in
[c, c(1+w)] and off-diagonal entries in [-d, 0] (i.e. c_ij in [0,d]), Q=R=I, at exact rational a2.
After K epochs it checks that the state lies in the certified tail box B_m (read from the JSON of
cert_x1_tail.py) for the a2-interval containing a2, and that min p >= pi_a = max(1/theta0,
(L_m-1)/(a2-1)). Exit 1 on failure. Formulas: see cert_x1_tail.py / EXTENSION_REPORT.md (X.2).
"""
import argparse, json, os, sys
from fractions import Fraction as Fr
from mpmath import iv, mp

PAIRS = [(0, 1), (1, 2), (2, 0)]
FREE = [(0, 1), (1, 0), (1, 2), (2, 1), (2, 0)]
def LO(x): return mp.make_mpf(x._mpi_[0])
def HI(x): return mp.make_mpf(x._mpi_[1])
def to_frac(x):
    s, man, ex, _ = x._mpf_
    v = Fr(man) * (Fr(2) ** ex)
    return -v if s else v

def frac_iv(f): return iv.mpf(f.numerator) / f.denominator

def meet(a, b):
    lo, hi = max(a.a, b.a), min(a.b, b.b)
    if lo > hi: raise SystemExit("empty intersection: enclosure inconsistent")
    return iv.mpf([lo, hi])

def epoch(p, e, a2):
    p = [a2 * x + 1 for x in p]; e = {k: a2 * x + 1 for k, x in e.items()}
    for (i, j) in PAIRS:
        l = 3 - i - j
        u, v = e[(i, j)], e[(j, i)]; s = u + v + 1
        if not s.a > 0: raise SystemExit("innovation not certified positive")
        li, lj = u * u / s, v * v / s
        e_li = e[(l, i)] + (1 - e[(i, l)] / p[i]) * li
        e_lj = e[(l, j)] + (1 - e[(j, l)] / p[j]) * lj
        e_il = e[(i, l)] * (1 - li / p[i]); e_jl = e[(j, l)] * (1 - lj / p[j])
        e[(i, j)] = u / s; e[(j, i)] = v / s
        e[(l, i)], e[(l, j)], e[(i, l)], e[(j, l)] = e_li, e_lj, e_il, e_jl
        p[i] = p[i] - li; p[j] = p[j] - lj
        # exact identity: e01-e10+e12-e21+e20-e02 = 0 and e_ij - e_ji = p_i - p_j
        e[(0, 2)] = meet(e[(0, 2)], e[(0, 1)] - e[(1, 0)] + e[(1, 2)] - e[(2, 1)] + e[(2, 0)])
        for (x, y) in [(0, 1), (1, 2), (2, 0)]:
            e[(x, y)] = meet(e[(x, y)], e[(y, x)] + p[x] - p[y])
    return p, e

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True); ap.add_argument("--tail", required=True)
    ap.add_argument("--a2", nargs="+", required=True); ap.add_argument("--c", default="10**12")
    ap.add_argument("--w", default="0"); ap.add_argument("--d", default="0")
    ap.add_argument("--K", type=int, default=80); ap.add_argument("--dps", type=int, default=60)
    a = ap.parse_args()
    if os.path.exists(a.out):
        print(json.dumps({"error": "output directory exists; refusing to overwrite"})); sys.exit(2)
    iv.dps = a.dps
    tail = json.load(open(a.tail))
    if not tail["summary"]["all_ok"]: print(json.dumps({"error": "tail certificate not ok"})); sys.exit(1)
    th0 = Fr(tail["summary"]["theta0"])
    c = Fr(eval(a.c, {"__builtins__": {}})); w = Fr(a.w); d = Fr(a.d)
    results = []; allok = True
    for s in a.a2:
        A2 = Fr(s); a2 = frac_iv(A2)
        m = [r for r in tail["intervals"] if Fr(r["a2_exact"][0]) <= A2 <= Fr(r["a2_exact"][1])]
        rec = {"a2": s, "c": str(c), "w": str(w), "d": str(d)}
        if not m or A2 <= 1:
            rec.update(ok=False, reason="a2 not covered by tail certificate"); results.append(rec); allok = False; continue
        r = m[0]
        sh = mp.mpf(10)**-12  # shrink printed (15-digit) box inward: rigorous containment in the certified box
        B = [iv.mpf([mp.mpf(lo) + sh, mp.mpf(hi) - sh]) for lo, hi in r["box"]]
        L = max(Fr(x) for x in r["loss_upper"])  # decimal strings -> exact fractions (upper bounds, rounded up below)
        L = L + Fr(1, 10**9)
        pi_a = max(1 / th0, (L - 1) / (A2 - 1))
        cc = frac_iv(c)
        p = [iv.mpf([LO(cc), HI(cc * (1 + frac_iv(w)))])] * 3
        cij = iv.mpf([0, HI(frac_iv(d))]) if d > 0 else iv.mpf(0)
        e = {k: p[k[0]] - cij for k in FREE + [(0, 2)]}
        hit = None
        try:
            for k in range(1, a.K + 1):
                p, e = epoch(p, e, a2)
                inB = all(e[key].a > b.a and e[key].b < b.b for key, b in zip(FREE, B))
                if inB and all(to_frac(LO(x)) > pi_a for x in p):  # exact, strict
                    hit = k; break
        except SystemExit as ex:
            rec.update(ok=False, reason=str(ex)); results.append(rec); allok = False; continue
        rec.update(ok=hit is not None, entry_epoch=hit, tail_interval=r["a2_exact"],
                   pi_a=float(pi_a), min_p_lower=mp.nstr(min(LO(x) for x in p), 12),
                   deficits=[[mp.nstr(LO(e[key]), 10), mp.nstr(HI(e[key]), 10)] for key in FREE])
        if hit is None: allok = False
        results.append(rec)
    os.makedirs(a.out)
    with open(os.path.join(a.out, "x1_entry.json"), "w") as f: json.dump(results, f, indent=1)
    print(json.dumps({"script": "cert_x1_entry.py", "all_ok": allok,
                      "cases": [(r["a2"], r.get("entry_epoch"), r["ok"]) for r in results]}))
    if not allok: sys.exit(1)

if __name__ == "__main__":
    try:
        main()
    except SystemExit:
        raise
    except Exception as exc:  # any unexpected failure becomes a JSON failure record
        print(json.dumps({"status": "FAILED", "error": type(exc).__name__, "detail": str(exc)}))
        sys.exit(1)
