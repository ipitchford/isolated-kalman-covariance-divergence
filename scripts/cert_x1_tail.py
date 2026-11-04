#!/usr/bin/env python3
"""X1 tail certificate (computer-assisted, rigorous interval arithmetic via mpmath.iv).

Triangle IsoKF, A=aI, Q=R=I (equivalently Q=R=qI after scaling P by 1/q), cone variables:
p_i = P_ii, c_ij = -P_ij, deficits e_ij = p_i - c_ij (ordered pair: 'deficit of i towards j').
Exact recursions (derived in EXTENSION_REPORT.md, eq. (X.2)):
  prediction: p -> a2 p + 1, e -> a2 e + 1
  pair (i,j), l the third node: u=e_ij, v=e_ji, s=u+v+1, l_i=u^2/s, l_j=v^2/s
     e_ij'=u/s, e_ji'=v/s, e_li'=e_li+(1-th_i e_il) l_i, e_lj'=e_lj+(1-th_j e_jl) l_j,
     e_il'=(1-th_i l_i) e_il, e_jl'=(1-th_j l_j) e_jl, p_i'=p_i-l_i, p_j'=p_j-l_j, th=1/p.
The 6 deficits satisfy the exact identity e01-e10+e12-e21+e20-e02=0 (cycle of p-differences);
e02 is eliminated with it, leaving 5 coordinates x.
For each a2-subinterval I_m and th in [0,TH0] (each use independent) the program verifies with
outward-rounded intervals, via the mean-value form F^N(B) c F^N(xhat)+J(B)(B-xhat):
  (i) F^N(B_m) is contained in the interior of B_m,   (ii) all intermediate states k<=N lie in a
  certified hull H_m with all deficits >0,   (iii) an upper bound Lmax_m on the per-epoch loss
  sum_{updates of node i} u^2/s for every node.
Output: JSON with per-interval boxes, hulls, loss bounds. Exit code 1 on any failed check.
"""
import argparse, json, os, sys
from mpmath import iv, mp, mpf

iv.dps = 30
PAIRS = [(0, 1), (1, 2), (2, 0)]
FREE = [(0, 1), (1, 0), (1, 2), (2, 1), (2, 0)]
NX = 5
NV = 6  # 5 free deficits + the parameter a2 (mean-value form also in a2)

def LO(x): return mp.make_mpf(x._mpi_[0])
def HI(x): return mp.make_mpf(x._mpi_[1])
def S(x): return mp.nstr(x, 15)

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
        if not o.v.a > 0: raise ZeroDivisionError("denominator interval not positive")
        q = s.v / o.v
        return D(q, [(a - q * b) / o.v for a, b in zip(s.g, o.g)])

def fix(e):
    e[(0, 2)] = e[(0, 1)] - e[(1, 0)] + e[(1, 2)] - e[(2, 1)] + e[(2, 0)]

def epoch(e, a2, th):
    e = {k: x * a2 + 1 for k, x in e.items()}
    fix(e)
    loss = [D(0), D(0), D(0)]
    for (i, j) in PAIRS:
        l = 3 - i - j
        u, v = e[(i, j)], e[(j, i)]
        s = u + v + 1
        li, lj = u * u / s, v * v / s
        e_li = e[(l, i)] + (1 - e[(i, l)] * th) * li
        e_lj = e[(l, j)] + (1 - e[(j, l)] * th) * lj
        e_il = e[(i, l)] * (1 - li * th)
        e_jl = e[(j, l)] * (1 - lj * th)
        e[(i, j)] = u / s; e[(j, i)] = v / s
        e[(l, i)], e[(l, j)], e[(i, l)], e[(j, l)] = e_li, e_lj, e_il, e_jl
        fix(e)
        loss[i] = loss[i] + li; loss[j] = loss[j] + lj
    return e, loss

def float_fixed_point(a2):
    # plain float iteration of the theta=0 map (only used to choose box centres)
    e = {k: 1.0 for k in FREE + [(0, 2)]}
    for _ in range(400):
        e = {k: x * a2 + 1 for k, x in e.items()}
        for (i, j) in PAIRS:
            l = 3 - i - j; u, v = e[(i, j)], e[(j, i)]; s = u + v + 1
            li, lj = u * u / s, v * v / s
            e[(l, i)] += li; e[(l, j)] += lj
            e[(i, j)] = u / s; e[(j, i)] = v / s
            e[(0, 2)] = e[(0, 1)] - e[(1, 0)] + e[(1, 2)] - e[(2, 1)] + e[(2, 0)]
    return [e[k] for k in FREE]

def certify(lo, hi, th0, rad, N):
    # lo, hi are Fractions; enclose them rigorously
    elo = iv.mpf(lo.numerator) / lo.denominator; ehi = iv.mpf(hi.numerator) / hi.denominator
    a2 = iv.mpf([elo.a, ehi.b]); th = iv.mpf([0, th0])
    xh = float_fixed_point((float(lo) + float(hi)) / 2)
    xhat = [iv.mpf(mp.mpf(str(round(x, 12)))) for x in xh]
    B = [iv.mpf([x.a - rad, x.b + rad]) for x in xhat]
    # Mean-value (centred) propagation: E_k encloses F^k(B); J_k encloses D(F^k) on B
    # (chain rule: J_k = DF(E_{k-1}) J_{k-1}); P_k encloses F^k(xhat) for all parameters.
    def seeded(box, a2box):
        d = {k: D(box[n], [iv.mpf(1) if m == n else iv.mpf(0) for m in range(NV)]) for n, k in enumerate(FREE)}
        return d, D(a2box, [iv.mpf(0)] * NX + [iv.mpf(1)])
    ahat = iv.mpf((lo + hi) / 2) if False else (iv.mpf((lo + hi).numerator) / (lo + hi).denominator) / 2
    ahat = iv.mpf(HI(ahat))  # a point inside the interval (exact binary number)
    if not (ahat.a >= a2.a and ahat.b <= a2.b): raise SystemExit("bad centre")
    C = list(xhat) + [ahat]
    Bfull = list(B) + [a2]
    Pk = {k: D(xhat[n]) for n, k in enumerate(FREE)}
    E = list(B); J = [[iv.mpf(1) if m == n else iv.mpf(0) for m in range(NV)] for n in range(NV)]
    hull = list(B); Lmax = [mpf(0)] * 3; Lmin = [mpf(10)**9] * 3; okpos = True
    for k in range(1, N + 1):
        Pk, _ = epoch(Pk, D(ahat), th)
        sd, ad = seeded(E, a2)
        out, loss = epoch(sd, ad, th)
        Lmax = [max(Lmax[i], HI(loss[i].v)) for i in range(3)]
        Lmin = [min(Lmin[i], LO(loss[i].v)) for i in range(3)]
        if not (out[(0, 2)].v.a > 0): okpos = False
        M = [[out[key].g[m] for m in range(NV)] for key in FREE] + [[iv.mpf(0)] * NX + [iv.mpf(1)]]
        J = [[sum((M[n][t] * J[t][m] for t in range(NV)), iv.mpf(0)) for m in range(NV)] for n in range(NV)]
        newE = []
        for n, key in enumerate(FREE):
            enc = Pk[key].v + sum((J[n][m] * (Bfull[m] - C[m]) for m in range(NV)), iv.mpf(0))
            nv_ = out[key].v
            lo_, hi_ = max(enc.a, nv_.a), min(enc.b, nv_.b)
            if lo_ > hi_: raise SystemExit("inconsistent enclosure")
            newE.append(iv.mpf([lo_, hi_]))
            if not lo_ > 0: okpos = False
        E = newE
        hull = [iv.mpf([min(h.a, x.a), max(h.b, x.b)]) for h, x in zip(hull, E)]
    final = E
    inside = all(f.a > b.a and f.b < b.b for f, b in zip(final, B))
    return dict(a2=[str(lo), str(hi)], theta0=str(th0), N=N, radius=str(rad),
                box=[[S(LO(b)), S(HI(b))] for b in B],
                image=[[S(LO(f)), S(HI(f))] for f in final],
                hull=[[S(LO(h)), S(HI(h))] for h in hull],
                loss_upper=[S(x) for x in Lmax],
                loss_lower=[S(x) for x in Lmin],
                positive=okpos, inside=inside, ok=bool(inside and okpos))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--a2max", default="3/2")
    ap.add_argument("--step", default="1/100")
    ap.add_argument("--theta0", default="1e-6")
    ap.add_argument("--N", type=int, default=4)
    ap.add_argument("--rad", default="0.006")
    a = ap.parse_args()
    if os.path.exists(a.out):
        print(json.dumps({"error": "output directory exists; refusing to overwrite"})); sys.exit(2)
    os.makedirs(a.out)
    from fractions import Fraction as Fr
    amax, st = Fr(a.a2max), Fr(a.step)
    res = []; lo = Fr(1)
    # upper endpoint of an outward-rounded enclosure of the exact rational theta0
    th0_q = Fr(a.theta0)
    th0 = (iv.mpf(th0_q.numerator) / th0_q.denominator).b
    while lo < amax:
        hi = min(lo + st, amax)
        rad = mpf(a.rad)
        r = certify(lo, hi, th0, rad, a.N)
        r["a2_exact"] = [str(lo), str(hi)]
        res.append(r); lo = hi
    allok = all(r["ok"] for r in res)
    Lglob = max(float(x) for r in res for x in r["loss_upper"])
    Llow0 = min(float(x) for x in res[0]["loss_lower"]) if res else None
    out = {"script": "cert_x1_tail.py", "intervals": len(res), "all_ok": allok,
           "global_loss_upper": Lglob, "first_interval_loss_lower": Llow0, "theta0": a.theta0, "N": a.N,
           "failed": [r["a2_exact"] for r in res if not r["ok"]]}
    with open(os.path.join(a.out, "x1_tail.json"), "w") as f:
        json.dump({"summary": out, "intervals": res}, f, indent=1)
    print(json.dumps(out))
    if not allok: sys.exit(1)

if __name__ == "__main__":
    try:
        main()
    except SystemExit:
        raise
    except Exception as exc:  # any unexpected failure becomes a JSON failure record
        print(json.dumps({"status": "FAILED", "error": type(exc).__name__, "detail": str(exc)}))
        sys.exit(1)
