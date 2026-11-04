#!/usr/bin/env python3
"""X2 certificate: at a2 = 21/20, Q=R=I, P0=I the triangle IsoKF covariance converges.
Variables x=(p0,p1,p2,c01,c12,c20), P_ii=p_i, P_ij=-c_ij. Exact recursions = manuscript (1)-(3), (9).
Phase A: outward-rounded interval simulation (mpmath.iv) of K epochs from P0=I.
Phase B: box B = xhat +- rad; mean-value form with interval forward-mode derivatives gives
  enclosures E_k of F^k(B), k=1..N, and an interval enclosure J of D(F^N) on B.  Checks:
  E_1 c B, E_N c int B, ||J||_inf < 1, every innovation and every p_i certified >0.
By the mean-value theorem F^N is then an inf-norm contraction of B into itself, so it has a unique
fixed point x* in B; F(x*) in E_1 c B is also fixed by F^N, hence F(x*)=x*; and the trajectory
(which enters B at epoch K) converges to x*. Prints JSON; exit 1 on any failed check.
"""
import argparse, json, os, sys
from fractions import Fraction as Fr
from mpmath import iv, mp

NV = 6
def LO(x): return mp.make_mpf(x._mpi_[0])
def HI(x): return mp.make_mpf(x._mpi_[1])

class D:
    __slots__ = ("v", "g")
    def __init__(s, v, g=None):
        s.v = v if isinstance(v, type(iv.mpf(0))) else iv.mpf(v)
        s.g = g if g is not None else [iv.mpf(0)] * NV
    def __add__(s, o):
        o = o if isinstance(o, D) else D(o); return D(s.v + o.v, [a + b for a, b in zip(s.g, o.g)])
    __radd__ = __add__
    def __sub__(s, o):
        o = o if isinstance(o, D) else D(o); return D(s.v - o.v, [a - b for a, b in zip(s.g, o.g)])
    def __rsub__(s, o): return D(o) - s
    def __mul__(s, o):
        o = o if isinstance(o, D) else D(o)
        return D(s.v * o.v, [s.v * b + o.v * a for a, b in zip(s.g, o.g)])
    __rmul__ = __mul__
    def __truediv__(s, o):
        o = o if isinstance(o, D) else D(o)
        if not o.v.a > 0: raise ZeroDivisionError("denominator not certified positive")
        q = s.v / o.v
        return D(q, [(a - q * b) / o.v for a, b in zip(s.g, o.g)])

PAIRS = [(0, 1), (1, 2), (2, 0)]
CI = {(0, 1): 3, (1, 0): 3, (1, 2): 4, (2, 1): 4, (2, 0): 5, (0, 2): 5}

def epoch(x, a2):
    x = [xi * a2 + (1 if n < 3 else 0) for n, xi in enumerate(x)]
    for (i, j) in PAIRS:
        l = 3 - i - j
        p, q, c = x[i], x[j], x[CI[(i, j)]]
        u, v = p - c, q - c; s = u + v + 1
        p2, q2, c2 = p - u * u / s, q - v * v / s, c + u * v / s
        cil, cjl = x[CI[(i, l)]] * (p2 / p), x[CI[(j, l)]] * (q2 / q)
        x[i], x[j], x[CI[(i, j)]], x[CI[(i, l)]], x[CI[(j, l)]] = p2, q2, c2, cil, cjl
    return x

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True); ap.add_argument("--a2", default="21/20")
    ap.add_argument("--K", type=int, default=150); ap.add_argument("--N", type=int, default=4)
    ap.add_argument("--rad", default="1e-6"); ap.add_argument("--dps", type=int, default=250)
    a = ap.parse_args()
    if os.path.exists(a.out):
        print(json.dumps({"error": "output directory exists; refusing to overwrite"})); sys.exit(2)
    iv.dps = a.dps; mp.dps = a.dps
    A2 = Fr(a.a2); a2 = iv.mpf(A2.numerator) / A2.denominator
    checks = {}
    # float-ish fixed point (mp, not rigorous) for the box centre
    xf = [mp.mpf(1)] * 3 + [mp.mpf(0)] * 3
    a2f = mp.mpf(A2.numerator) / A2.denominator
    for _ in range(400):
        y = [D(iv.mpf(t)) for t in xf]; y = epoch(y, iv.mpf(a2f)); xf = [mp.mpf(LO(t.v)) for t in y]
    xhat = [iv.mpf(mp.nstr(t, 25)) for t in xf]
    rad = iv.mpf(a.rad)
    B = [iv.mpf([LO(t - rad), HI(t + rad)]) for t in xhat]
    # Phase A
    x = [D(iv.mpf(1))] * 3 + [D(iv.mpf(0))] * 3
    entry = None
    for k in range(1, a.K + 1):
        x = epoch(x, a2)
        if all(t.v.a > b.a and t.v.b < b.b for t, b in zip(x, B)): entry = k; break
    checks["trajectory_enters_B"] = entry is not None
    # Phase B
    E = list(B); J = [[iv.mpf(1 if m == n else 0) for m in range(NV)] for n in range(NV)]
    Pk = [D(t) for t in xhat]; hull = list(B); E1 = None
    for k in range(1, a.N + 1):
        Pk = epoch(Pk, a2)
        sd = [D(E[n], [iv.mpf(1 if m == n else 0) for m in range(NV)]) for n in range(NV)]
        out = epoch(sd, a2)
        M = [[o.g[m] for m in range(NV)] for o in out]
        J = [[sum((M[n][t] * J[t][m] for t in range(NV)), iv.mpf(0)) for m in range(NV)] for n in range(NV)]
        newE = []
        for n in range(NV):
            enc = Pk[n].v + sum((J[n][m] * (B[m] - xhat[m]) for m in range(NV)), iv.mpf(0))
            lo, hi = max(enc.a, out[n].v.a), min(enc.b, out[n].v.b)
            if lo > hi: print(json.dumps({"error": "inconsistent enclosure"})); sys.exit(1)
            newE.append(iv.mpf([lo, hi]))
        E = newE
        if k == 1: E1 = list(E)
        hull = [iv.mpf([min(h.a, t.a), max(h.b, t.b)]) for h, t in zip(hull, E)]
    checks["E1_subset_B"] = all(t.a >= b.a and t.b <= b.b for t, b in zip(E1, B))
    checks["EN_in_interior_B"] = all(t.a > b.a and t.b < b.b for t, b in zip(E, B))
    # row sums of |J| in interval arithmetic; use upper endpoints
    jn = max(HI(sum((abs(J[n][m]) for m in range(NV)), iv.mpf(0))) for n in range(NV))
    checks["contraction_norm_lt_1"] = bool(jn < 1)
    checks["cone_positive"] = all(LO(b) > 0 for b in B[:3]) and all(LO(b) >= 0 for b in B[3:]) and \
        all(LO(B[n]) - HI(B[3 + t]) > 0 for n, t in [(0, 0), (1, 0), (1, 1), (2, 1), (2, 2), (0, 2)])
    ok = all(checks.values())
    names = ["p1", "p2", "p3", "c12", "c23", "c31"]
    res = {"script": "cert_x2_bounded.py", "a2": a.a2, "Q": 1, "R": 1, "P0": "I", "entry_epoch": entry,
           "N": a.N, "contraction_bound_inf_norm_of_DF^N": mp.nstr(jn, 8),
           "certificate_box_B_cert_exact_binary_endpoints": {nm: [repr(LO(b)), repr(HI(b))] for nm, b in zip(names, B)},
           "fixed_point_box_20_digits_for_reporting": {nm: [mp.nstr(LO(b), 20), mp.nstr(HI(b), 20)] for nm, b in zip(names, B)},
           "min_marginal_at_fixed_point_approx": mp.nstr(min(xf[:3]), 15),
           "checks": checks, "all_ok": ok}
    os.makedirs(a.out)
    with open(os.path.join(a.out, "x2_bounded.json"), "w") as f: json.dump(res, f, indent=1)
    print(json.dumps(res))
    if not ok: sys.exit(1)

if __name__ == "__main__":
    try:
        main()
    except SystemExit:
        raise
    except Exception as exc:  # any unexpected failure becomes a JSON failure record
        print(json.dumps({"status": "FAILED", "error": type(exc).__name__, "detail": str(exc)}))
        sys.exit(1)
