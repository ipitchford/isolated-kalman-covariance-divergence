"""Illustrate coexistence at a^2 = 21/20 and the true estimation error (Figure 2).

Triangle schedule of Theorem 1 (pair sums 12, 23, 31, unit noise), A = aI with
a^2 = 21/20, Q = I.  Three priors P0 = c I with c = 1, 1e4, 1e12.  For each prior
the script computes, with Python Decimal at --precision digits:
  * the reported IsoKF covariance, by the dense recurrence (1)-(3) of the paper;
  * the covariance of the true estimation error of the IsoKF estimator, whose
    participant gains K = B h^T / S come from the reported block and whose
    nonparticipant gains are zero (exact linear recursion, no Monte Carlo);
  * the centralised Kalman posterior.
This is an ILLUSTRATION.  The bounded case (c = 1) and the divergent case
(c = 1e12) are proved in the paper by interval certificates (scripts
cert_x2_bounded.py, cert_x1_tail.py, cert_x1_entry.py); the true-error curves are
numerical findings only.  Explicit checks; exit code 1 on failure; works under -O.
"""
from decimal import Decimal, localcontext
from pathlib import Path
import argparse
import csv
import json
import os
import sys
import tempfile

CHECKS = {"passed": 0}
PAIRS = ((0, 1), (1, 2), (2, 0))


def check(condition, label):
    if not condition:
        print(json.dumps({"status": "FAILED", "failed_check": label}))
        sys.exit(1)
    CHECKS["passed"] += 1


def eye(c):
    return [[Decimal(c) if i == j else Decimal(0) for j in range(3)] for i in range(3)]


def mm(a, b):
    return [[sum(a[i][k] * b[k][j] for k in range(3)) for j in range(3)] for i in range(3)]


def tr(a):
    return [[a[j][i] for j in range(3)] for i in range(3)]


def run(a2, c, epochs):
    p = eye(c)          # reported IsoKF array
    sig = eye(c)        # true error covariance of the IsoKF estimator
    kf = eye(c)         # centralised Kalman posterior
    rows = []
    for k in range(epochs + 1):
        if k:
            p = [[a2 * p[i][j] + (1 if i == j else 0) for j in range(3)] for i in range(3)]
            sig = [[a2 * sig[i][j] + (1 if i == j else 0) for j in range(3)] for i in range(3)]
            kf = [[a2 * kf[i][j] + (1 if i == j else 0) for j in range(3)] for i in range(3)]
            for i, j in PAIRS:
                l = 3 - i - j
                u = p[i][i] + p[i][j]
                v = p[j][j] + p[i][j]
                s = u + v + 1
                check(u > 0 and v > 0 and s > 0, "local block positive")
                gi, gj = u / s, v / s
                # true error: e <- (I - K h) e - K w
                f = eye(1)
                f[i][i] -= gi; f[i][j] -= gi; f[j][i] -= gj; f[j][j] -= gj
                kv = [Decimal(0)] * 3
                kv[i], kv[j] = gi, gj
                sig = mm(mm(f, sig), tr(f))
                sig = [[sig[r][q] + kv[r] * kv[q] for q in range(3)] for r in range(3)]
                # reported IsoKF update (2)-(3)
                pi_new = p[i][i] - u * u / s
                pj_new = p[j][j] - v * v / s
                pij_new = p[i][j] - u * v / s
                new = [row[:] for row in p]
                new[i][i], new[j][j] = pi_new, pj_new
                new[i][j] = new[j][i] = pij_new
                new[i][l] = new[l][i] = (pi_new / p[i][i]) * p[i][l]
                new[j][l] = new[l][j] = (pj_new / p[j][j]) * p[j][l]
                p = new
                # centralised KF
                ph = [kf[r][i] + kf[r][j] for r in range(3)]
                sk = ph[i] + ph[j] + 1
                kf = [[kf[r][q] - ph[r] * ph[q] / sk for q in range(3)] for r in range(3)]
        rows.append({
            "epoch": k,
            "reported_min": min(p[i][i] for i in range(3)),
            "reported_max": max(p[i][i] for i in range(3)),
            "true_error_max": max(sig[i][i] for i in range(3)),
            "central_max": max(kf[i][i] for i in range(3)),
        })
    return rows


def render(series, folder):
    os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "isokf_matplotlib_cache"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "pdf.fonttype": 42, "ps.fonttype": 42})
    colours = {"1": "#49756D", "10000": "#B77738", "1000000000000": "#824C71"}
    names = {"1": r"$P_0=I$", "10000": r"$P_0=10^4 I$", "1000000000000": r"$P_0=10^{12} I$"}
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.6, 4.0), layout="constrained")
    for c, rows in series.items():
        xx = [r["epoch"] for r in rows]
        ax1.semilogy(xx, [float(r["reported_max"]) for r in rows], color=colours[c],
                     linewidth=2, label=names[c])
        ax2.plot(xx, [float(r["true_error_max"]) for r in rows], color=colours[c],
                 linewidth=2, label=names[c])
    rows = series["1"]
    xx = [r["epoch"] for r in rows]
    ax1.semilogy(xx, [float(r["central_max"]) for r in rows], color="black", linestyle=":",
                 linewidth=1.6, label="Centralised KF ($P_0=I$)")
    ax2.plot(xx, [float(r["central_max"]) for r in rows], color="black", linestyle=":",
             linewidth=1.6, label="Centralised KF ($P_0=I$)")
    ax1.set_title("(a) Reported IsoKF covariance: largest marginal", fontsize=10)
    ax2.set_title("(b) Actual error variance of the IsoKF estimate: largest", fontsize=10)
    ax2.set_ylim(0, 1.2)
    for ax in (ax1, ax2):
        ax.set_xlabel(r"Epoch $k$")
        ax.grid(axis="y", color="#ddd8d2", linewidth=0.6)
        ax.set_axisbelow(True)
        ax.set_xlim(0, xx[-1])
    ax1.legend(loc="lower right", bbox_to_anchor=(1.0, 0.07), frameon=False, fontsize=8.5)
    ax2.legend(loc="upper right", frameon=False, fontsize=8.5)
    for ext in ("png", "pdf"):
        fig.savefig(folder / f"coexistence.{ext}", dpi=220, facecolor="white")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", required=True, help="new directory (must not exist)")
    ap.add_argument("--epochs", type=int, default=300)
    ap.add_argument("--precision", type=int, default=150)
    ap.add_argument("--no-figure", action="store_true")
    args = ap.parse_args()
    if os.path.exists(args.out):
        print(json.dumps({"error": "output directory exists; refusing to overwrite"}))
        sys.exit(2)
    series = {}
    with localcontext() as ctx:
        ctx.prec = args.precision
        a2 = Decimal(21) / Decimal(20)
        for c in (1, 10 ** 4, 10 ** 12):
            series[str(c)] = run(a2, c, args.epochs)
    last = {c: rows[-1] for c, rows in series.items()}
    # Consistency with the certified statements (Theorem 5.2 bounded box; Theorem 5.1 growth).
    check(Decimal("1.6897") < last["1"]["reported_max"] < Decimal("1.6898"), "P0=I bounded near certified fixed point")
    check(last["1000000000000"]["reported_min"] > Decimal(10) ** 12, "P0=1e12 I has grown")
    check(last["10000"]["reported_min"] > Decimal(10) ** 6, "P0=1e4 I has grown (numerical)")
    for c in series:
        check(last[c]["true_error_max"] < 1, f"true error bounded at end, c={c}")
    folder = Path(args.out)
    folder.mkdir(parents=True)
    with (folder / "coexistence.csv").open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["prior_c", "epoch", "reported_min", "reported_max", "true_error_max", "central_max"])
        for c, rows in series.items():
            for r in rows:
                w.writerow([c, r["epoch"]] + [format(r[k], ".12e") for k in
                                               ("reported_min", "reported_max", "true_error_max", "central_max")])
    summary = {"script": "plot_coexistence.py", "label": "illustration (Decimal arithmetic)",
               "a2": "21/20", "epochs": args.epochs, "precision": args.precision,
               "final": {c: {k: format(v, ".10e") for k, v in r.items() if k != "epoch"} for c, r in last.items()},
               "explicit_checks_passed": CHECKS["passed"]}
    (folder / "coexistence_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    if not args.no_figure:
        render(series, folder)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    try:
        main()
    except SystemExit:
        raise
    except Exception as exc:  # any unexpected failure becomes a JSON failure record
        print(json.dumps({"status": "FAILED", "error": type(exc).__name__, "detail": str(exc)}))
        sys.exit(1)
