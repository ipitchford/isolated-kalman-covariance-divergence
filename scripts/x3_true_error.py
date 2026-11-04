#!/usr/bin/env python3
"""X3 (NUMERICAL, high precision, no Monte Carlo): true estimation-error covariance of the IsoKF
estimator versus the reported IsoKF covariance, triangle schedule, A=aI, Q=R=I.
Estimator model (stated assumption): at each scalar pair observation y=x_i+x_j+v, the participants
i,j update their estimates with the local Kalman gain K=(B h^T)/S computed from the reported local
block B (h=(1,1), S=hBh^T+R); nonparticipant estimates are unchanged (zero gain); prediction
xhat <- A xhat. The true error e=x-xhat then obeys the exact linear recursion
  prediction: Sigma <- A Sigma A^T + Q;   update: Sigma <- F Sigma F^T + K R K^T, F=I-K h (embedded),
with Sigma_0 = P_0 (prior covariance is the true one). Gains are deterministic, so Sigma is exact
(up to mpmath rounding at --dps digits). Reported P uses cone variables (p, deficits) to avoid
cancellation: K_i = (p_i - c_ij)/S, K_j = (p_j - c_ij)/S, S = u+v+R.
Also reports the centralised KF posterior for comparison and the spectral radius of the asymptotic
epoch error-transition matrix (numerical). Prints JSON; writes to --out (refuses to overwrite).
"""
import argparse, json, os, sys
from mpmath import mp, mpf, matrix, eig

PAIRS = [(0, 1), (1, 2), (2, 0)]

def run(a2, c0, K, dps):
    mp.dps = dps
    a2 = mpf(a2); a = mp.sqrt(a2)
    p = [mpf(c0)] * 3; e = {(i, j): mpf(c0) for i in range(3) for j in range(3) if i != j}
    Sig = matrix(3, 3); Pk = matrix(3, 3)
    for i in range(3): Sig[i, i] = mpf(c0); Pk[i, i] = mpf(c0)
    rows = []; Phi_last = None
    for k in range(1, K + 1):
        p = [a2 * x + 1 for x in p]; e = {kk: a2 * x + 1 for kk, x in e.items()}
        Sig = a2 * Sig + matrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        Pk = a2 * Pk + matrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        Phi = matrix([[a, 0, 0], [0, a, 0], [0, 0, a]])
        for (i, j) in PAIRS:
            l = 3 - i - j
            u, v = e[(i, j)], e[(j, i)]; s = u + v + 1
            Ki, Kj = u / s, v / s
            F = matrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
            F[i, i] -= Ki; F[i, j] -= Ki; F[j, i] -= Kj; F[j, j] -= Kj
            Kv = matrix(3, 1); Kv[i] = Ki; Kv[j] = Kj
            Sig = F * Sig * F.T + Kv * Kv.T
            Phi = F * Phi
            # centralised KF (sequential scalar updates = batch KF)
            h = matrix(1, 3); h[i] = 1; h[j] = 1
            Sk = (h * Pk * h.T)[0] + 1; G = Pk * h.T / Sk; Pk = Pk - G * h * Pk
            # reported IsoKF (cone variables)
            li, lj = u * u / s, v * v / s
            e_li = e[(l, i)] + (1 - e[(i, l)] / p[i]) * li; e_lj = e[(l, j)] + (1 - e[(j, l)] / p[j]) * lj
            e_il = e[(i, l)] * (1 - li / p[i]); e_jl = e[(j, l)] * (1 - lj / p[j])
            e[(i, j)] = u / s; e[(j, i)] = v / s
            e[(l, i)], e[(l, j)], e[(i, l)], e[(j, l)] = e_li, e_lj, e_il, e_jl
            p[i] -= li; p[j] -= lj
            # enforce the exact cycle identity e01-e10+e12-e21+e20-e02=0 (only O(deficit) terms;
            # removes the non-physical rounding mode that otherwise grows like a2^k)
            e[(0, 2)] = e[(0, 1)] - e[(1, 0)] + e[(1, 2)] - e[(2, 1)] + e[(2, 0)]
        Phi_last = Phi
        if k in (1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000) or k == K:
            rows.append({"epoch": k, "reported_diag": [mp.nstr(x, 8) for x in p],
                         "true_error_diag": [mp.nstr(Sig[i, i], 8) for i in range(3)],
                         "centralised_KF_diag": [mp.nstr(Pk[i, i], 8) for i in range(3)],
                         "true_error_min_eig": mp.nstr(min(mp.re(x) for x in eig(Sig)[0]), 6)})
    ev = eig(Phi_last)[0]
    rho = max(abs(x) for x in ev)
    return {"a2": str(a2), "P0": f"{c0}*I", "epochs": K, "rows": rows,
            "spectral_radius_last_epoch_error_transition": mp.nstr(rho, 10),
            "final_gains_note": "gains K_i=u/s computed from reported deficits"}

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--out", required=True)
    ap.add_argument("--dps", type=int, default=80)
    a = ap.parse_args()
    if os.path.exists(a.out):
        print(json.dumps({"error": "output directory exists; refusing to overwrite"})); sys.exit(2)
    cases = [("9", 1, 200), ("1.05", 1, 2000), ("1.05", 10**4, 2000), ("1.01", 10**8, 2000), ("2", 1, 500)]
    res = [run(*c, a.dps) for c in cases]
    os.makedirs(a.out)
    with open(os.path.join(a.out, "x3_true_error.json"), "w") as f: json.dump(res, f, indent=1)
    summ = [{"a2": r["a2"], "P0": r["P0"], "final": r["rows"][-1], "rho": r["spectral_radius_last_epoch_error_transition"]} for r in res]
    print(json.dumps({"script": "x3_true_error.py", "label": "numerical (high precision)", "summary": summ}))
    bad = [r for r in res if not r["rows"]]
    if bad: sys.exit(1)

if __name__ == "__main__":
    try:
        main()
    except SystemExit:
        raise
    except Exception as exc:  # any unexpected failure becomes a JSON failure record
        print(json.dumps({"status": "FAILED", "error": type(exc).__name__, "detail": str(exc)}))
        sys.exit(1)
