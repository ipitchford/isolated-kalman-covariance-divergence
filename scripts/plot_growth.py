"""Reproduce the illustrative covariance trajectories in the research manuscript.

The proof of divergence is analytic. This program is an independent, high-
precision illustration, not a formal verification of the infinite-time claim.
It uses Python's standard-library Decimal for the recurrence and matplotlib
only for rendering. Run from any directory with Python 3.10 or later.
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


def check(condition, label):
    """Explicit check (works under python -O); exit with code 1 on failure."""
    if not condition:
        print(json.dumps({"status": "FAILED", "failed_check": label}))
        sys.exit(1)
    CHECKS["passed"] += 1


def identity(n):
    return [[Decimal(int(i == j)) for j in range(n)] for i in range(n)]


def predict(p):
    return [[9 * p[i][j] + int(i == j) for j in range(3)]
            for i in range(3)]


def isolated_pair(p, i, j):
    """One unit-noise pair-sum update, with marginal-ratio external correction."""
    out = [row[:] for row in p]
    pp, qq, cc = p[i][i], p[j][j], -p[i][j]
    uu, vv = pp - cc, qq - cc
    ss = uu + vv + 1
    check(uu > 0 and vv > 0 and ss > 0, 'plot_growth:33')
    pp_new = pp - uu * uu / ss
    qq_new = qq - vv * vv / ss
    cc_new = cc + uu * vv / ss
    out[i][i], out[j][j] = pp_new, qq_new
    out[i][j] = out[j][i] = -cc_new
    for k in range(3):
        if k not in (i, j):
            out[i][k] = out[k][i] = (pp_new / pp) * p[i][k]
            out[j][k] = out[k][j] = (qq_new / qq) * p[j][k]
    for u in range(3):
        check(out[u][u] > 0, 'plot_growth:44')
        for v in range(u + 1, 3):
            check(0 <= -out[u][v] < min(out[u][u], out[v][v]), 'plot_growth:46')
    return out


def central_pair(p, i, j):
    """The same observation, updating all states and all covariance entries."""
    ph = [p[k][i] + p[k][j] for k in range(3)]
    ss = ph[i] + ph[j] + 1
    check(ss > 0, 'plot_growth:54')
    return [[p[u][v] - ph[u] * ph[v] / ss for v in range(3)]
            for u in range(3)]


def trajectory(epochs, precision):
    rows = []
    with localcontext() as ctx:
        ctx.prec = precision
        isolated = identity(3)
        central = identity(3)
        lower = Decimal(1)
        for epoch in range(epochs + 1):
            if epoch:
                isolated, central = predict(isolated), predict(central)
                for i, j in ((0, 1), (1, 2), (2, 0)):
                    isolated = isolated_pair(isolated, i, j)
                    central = central_pair(central, i, j)
                lower = (9 * lower + 1) / 5
            minimum = min(isolated[i][i] for i in range(3))
            maximum = max(central[i][i] for i in range(3))
            check(minimum >= lower, 'plot_growth:75')
            if epoch:
                check(maximum < Decimal("0.75"), 'plot_growth:77')
            rows.append({
                "epoch": epoch,
                "isokf_min_marginal": str(minimum),
                "analytic_lower_bound": str(lower),
                "central_max_marginal": str(maximum),
                **{f"isokf_p{i+1}": str(isolated[i][i]) for i in range(3)},
            })
    return rows


def render(rows, folder):
    os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "isokf_matplotlib_cache"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator

    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 11,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })
    xx = [r["epoch"] for r in rows]
    fig, ax = plt.subplots(figsize=(8.6, 4.8), layout="constrained")
    ax.semilogy(xx, [float(r["isokf_min_marginal"]) for r in rows],
                color="#824C71", linewidth=2.3, marker="o", markersize=3,
                label="Smallest IsoKF marginal (reported)")
    ax.semilogy(xx, [float(r["analytic_lower_bound"]) for r in rows],
                color="#B77738", linewidth=2, linestyle="--",
                label=r"Proved lower bound: $\frac{5}{4}(\frac{9}{5})^k-\frac{1}{4}$")
    ax.semilogy(xx, [float(r["central_max_marginal"]) for r in rows],
                color="#49756D", linewidth=2.3,
                label="Largest centralised Kalman marginal")
    ax.set_xlabel(r"Completed observation epochs, $k$")
    ax.set_ylabel("Reported marginal variance (logarithmic scale)")
    ax.set_xlim(0, xx[-1])
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.grid(axis="y", color="#ddd8d2", linewidth=0.6)
    ax.set_axisbelow(True)
    ax.legend(loc="upper left", frameon=False, fontsize=10)
    for extension in ("png", "pdf"):
        fig.savefig(folder / f"covariance_growth.{extension}", dpi=220,
                    facecolor="white", metadata={"Title": "Observable system with divergent isolated Kalman covariances"})
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", required=True, help="new directory (must not exist)")
    parser.add_argument("--epochs", type=int, default=20)
    parser.add_argument("--precision", type=int, default=100)
    parser.add_argument("--no-figure", action="store_true", help="skip matplotlib rendering")
    args = parser.parse_args()
    if args.epochs < 1 or args.precision < 50:
        parser.error("Use at least one epoch and at least 50 digits of precision.")
    if os.path.exists(args.out):
        print(json.dumps({"error": "output directory exists; refusing to overwrite"}))
        sys.exit(2)
    rows = trajectory(args.epochs, args.precision)
    folder = Path(args.out)
    folder.mkdir(parents=True)
    with (folder / "covariance_growth.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    (folder / "figure_metadata.json").write_text(json.dumps({
        "arithmetic": "Python Decimal",
        "decimal_precision": args.precision,
        "epochs": args.epochs,
        "model": "A=3I, Q=I, P0=I; pair-sum observations 12,23,31 with R=1",
        "scope": "Illustration only; the analytic induction establishes divergence.",
        "explicit_checks_passed": CHECKS["passed"],
    }, indent=2) + "\n")
    if not args.no_figure:
        render(rows, folder)
    print(json.dumps({"status": "illustration_generated", "epochs": args.epochs,
                      "explicit_checks_passed": CHECKS["passed"],
                      "first_four": rows[1:5], "last": rows[-1]}, indent=2))


if __name__ == "__main__":
    try:
        main()
    except SystemExit:
        raise
    except Exception as exc:  # any unexpected failure becomes a JSON failure record
        print(json.dumps({"status": "FAILED", "error": type(exc).__name__, "detail": str(exc)}))
        sys.exit(1)
